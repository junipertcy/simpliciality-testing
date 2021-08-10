#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# simplicial_test -- a python module to realize simplicial degree-size sequences
#
# Copyright (C) 2020-2021 Tzu-Chi Yen <tzuchi.yen@colorado.edu>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

from . import validators
from .utils import *
from .custom_exceptions import NoMoreBalls, SimplicialSignal, GoToNextLevel
from copy import deepcopy

try:
    from sage.all import Combinations as combinations
except ImportError:
    from more_itertools import distinct_combinations as combinations
else:
    import sys

    sys.setrecursionlimit(100000)


class Test:
    r"""Base class for [Simpliciality]Test.

    Parameters
    ----------
    degree_list : ``iterable`` or :class:`numpy.ndarray`, required
        Sequence of vertex degree distribution.

    size_list : ``iterable`` or :class:`numpy.ndarray`, required
        Sequence of facet size distribution.

    blocked_sets : ``list`` of integer-valued ``tuple`` objects (optional, default ``None``)
        Facets that must be respected for non-inclusion.

    verbose : ``bool`` (optional, default ``False``)
        If ``True``, progress information will be shown.

    **kwargs :  keyword arguments
        Keyword arguments to be passed to base type constructor.

    Examples
    --------
    >>> from simplicial_test import *
    >>> degree_list, size_list = ([2, 1, 3, 2, 1], [3, 2, 2, 2])
    >>> st = Test(degree_list, size_list, verbose=False)
    >>> is_simplicial, facets = st.is_simplicial()
    >>> print(is_simplicial, facets)
    (True, ((0, 1, 3), (0, 2), (0, 4), (1, 2)))

    """

    def __init__(self, degree_list, size_list, verbose=False, **kwargs):
        super().__init__()
        self._level = 0
        self.s_depot = kwargs.pop("s_depot", SimplicialDepot(
            sorted(degree_list, reverse=True), sorted(size_list, reverse=True)
        ))
        trim_ones_first = kwargs.pop("trim_ones_first", True)
        self.facets_to_append = []
        if trim_ones_first:
            self.size_list, self.degree_list, num_ones = pair_one_by_one(size_list, degree_list)
            for _ in range(num_ones):
                self.facets_to_append += [(len(self.degree_list) + _,)]
        else:
            self.degree_list = self.s_depot.degree_list
            self.size_list = self.s_depot.size_list
        self.blocked_sets = []

        self.depth = kwargs.pop("depth", np.infty)
        self.width = kwargs.pop("width", 1e3)
        self.s_depot.cutoff = kwargs.pop("cutoff", 1e5)
        self.verbose = verbose
        self.current_fids = []
        self.non_simplicial_signal = False
        if len(kwargs) > 0:
            raise ValueError(f"unrecognized keyword arguments: {str(list(kwargs.keys()))}")

    @staticmethod
    def get_distinct_selection(size, degs, blocked_sets, level_map):
        r"""Compute the generator for candidate facets (each of size ``size``) from a set of nodes.
        These facets are ordered according to the branching heuristic as explained in the paper,
        see https://arxiv.org/abs/2106.00185.

        Parameters
        ----------
        size : ``int``
            Number of vertices to make the facet.

        degs : ``list`` of ``int``
            Sequence of wanting (or residual) vertex degree sequence.

        blocked_sets : ``list`` of integer-valued ``tuple`` objects (optional, default ``None``)
            Facets that must be respected for non-inclusion.

        level_map : ``dict``

        Returns
        -------
        facet : ``tuple``
            The candidate facet.

        Notes
        -----

        Note that _ = vtx degree at original view, vid = vtx index at original view,
        level_map[vid] =  vtx index at current view.

        """
        # blocked_sets = [_ for _ in blocked_sets if len(_) >= size]
        equiv2vid = defaultdict(list)
        for vid, _ in enumerate(degs):
            if level_map[vid] != -1:
                key = tuple(get_indices_of_k_in_blocked_sets(blocked_sets, level_map[vid]) + [_])
                equiv2vid[key] += [level_map[vid]]
                equiv2vid["pool"] += [key]
        equiv_class_pool = combinations(equiv2vid["pool"], size).__iter__()
        while True:
            facet = []
            tracker = defaultdict(int)
            for equiv_class in next(equiv_class_pool):
                facet += [equiv2vid[equiv_class][tracker[equiv_class]]]  # vids_same_equiv_class[tracker[equiv_class]]
                tracker[equiv_class] += 1
            yield tuple(facet)

    def sample_candidate_facet(self, size, valid_trials=None):
        r"""

        Parameters
        ----------
        size : ``int``
            Number of vertices to make the facet.

        valid_trials

        Returns
        -------

        """
        if not valid_trials:
            self.s_depot.valid_trials[self._level - 1] = self.get_distinct_selection(
                size, self.s_depot.prev_d[1], self.blocked_sets, self.s_depot.level_map[self._level - 1])
        counter = 0
        while True:
            if counter >= self.width or self.non_simplicial_signal:
                raise NoMoreBalls
            try:
                facet = next(self.s_depot.valid_trials[self._level - 1])  # candidate_facet
            except RuntimeError:
                raise NoMoreBalls
            if validators.a_issubset_any_b(facet, self.blocked_sets):
                continue
            counter += 1
            if not self.validate(facet):
                self.non_simplicial_signal = self.s_depot.add_to_time_counter(self._level - 1)

    def validate(self, facet) -> (bool, str):
        r"""This function must return True in order for the candidate facet to be considered.

        Parameters
        ----------
        facet : ``tuple``
            The candidate facet.

        Returns
        -------
        token : ``bool``
            A tuple of the form ``(ind, reason)``, where ``ind`` is the indicator for whether
            the candidate facet passes the validations and ``reason`` explains why they fail.

        Raises
        ------
        NoMoreBalls : :class:`simplicial_test.custom_exceptions.NoMoreBalls`
            if we depleted this level and wanted to go up (i.e., from `lv` to `lv - 1`).

        SimplicialSignal : :class:`simplicial_test.custom_exceptions.SimplicialSignal`
            Signal that indicates a simplicial instance has been found.

        Notes
        -----

        """
        residual_sizes = self.size_list
        if np.sum(residual_sizes) == 0:
            raise SimplicialSignal([facet])  # Last facet explored.

        safe, residual_degs, blocked_sets = validators.apply(self.degree_list, residual_sizes, self.blocked_sets, facet)
        if not safe:
            return False

        residual_degs, residual_sizes, collected_facets, exempt_vids = validators.reduce(
            residual_degs, residual_sizes, [facet] + blocked_sets, verbose=self.verbose
        )
        if np.sum(residual_sizes) == 0:
            raise SimplicialSignal([facet] + [exempt_vids] + collected_facets)

        filtered = filter_blocked_facets([facet] + blocked_sets, exempt_vids)
        blocked_sets = sort_facets(self.blocked_sets)

        sorted_wanting_degs = sorted(residual_degs, reverse=True)
        _ = (tuple(sorted_wanting_degs), tuple(blocked_sets))
        if _[1] in self.s_depot.explored[_[0]]:
            raise NoMoreBalls
        self.s_depot.explored[_[0]].add(_[1])

        self.s_depot.mappers[self._level] = get_seq2seq_mapping(residual_degs)  # (mapping2shrinked, mapping2enlarged)
        self.s_depot.compute_level_map(self._level, self.s_depot.mappers[self._level][0])
        self.s_depot.collects[self._level] = collected_facets
        self.s_depot.exempts[self._level] = exempt_vids
        self.s_depot.candidates[self._level] = [facet]

        blocked_sets = transform_facets(filtered, self.s_depot.mappers[self._level][0], to="l+1")
        raise GoToNextLevel(sorted_wanting_degs, residual_sizes, sort_facets(blocked_sets))

    def is_simplicial(self):
        if np.sum(self.size_list) == 0:
            return True, self.__mark(True, self._assemble_simplicial_facets([]))
        if not validators.preprocess(self.degree_list, self.size_list):
            self.s_depot.add_to_time_counter(0)
            return False, self.__mark(False, tuple())
        while True:
            self._level += 1
            self._verbose_logging()
            self._stage_state()
            s = self.size_list.pop(0)
            if self.non_simplicial_signal:
                return False, self.__mark(False, tuple())
            try:
                self.sample_candidate_facet(s, valid_trials=self.s_depot.valid_trials[self._level - 1])
            except SimplicialSignal as e:
                return True, self.__mark(True, self._assemble_simplicial_facets(e.message))
            except GoToNextLevel as e:
                self.degree_list, self.size_list, self.blocked_sets = e.message
                # degree_list, size_list, blocked_sets = e.message
                # if not validators.preprocess(degree_list, size_list):
                #     self._rollback(self._level - 1)
                # else:
                #     self.degree_list, self.size_list = degree_list, size_list
                #     self.blocked_sets = sort_facets(blocked_sets)
            except NoMoreBalls:
                if self._level == 1:
                    return False, self.__mark(False, tuple())
                if self._level > self.depth:
                    self._rollback(self.depth)
                else:
                    self._rollback(self._level - 1)

    def _clean_valid_trials(self):
        for ind, _ in enumerate(self.s_depot.size_list):
            if ind >= self._level - 1:
                self.s_depot.valid_trials[ind] = None

    def _assemble_simplicial_facets(self, seed_facets):
        if self._level != 1:
            for lv in np.arange(self._level - 1, 0, -1):
                facets = transform_facets(seed_facets, self.s_depot.mappers[lv][1], to="l-1")
                facets = [self.s_depot.exempts[lv] + list(s) for s in facets]
                seed_facets = self.s_depot.candidates[lv] + facets + self.s_depot.collects[lv]
        return tuple(sort_callback(seed_facets) + self.facets_to_append)

    def __mark(self, simplicial, facets):
        self.s_depot.simplicial = simplicial
        self.s_depot.facets = sort_facets(facets)
        del self.s_depot.valid_trials
        del self.s_depot.explored
        return self.s_depot.facets

    def _rollback(self, level):
        if level == 0:
            self.non_simplicial_signal = True
            return
        self._clean_valid_trials()
        self.non_simplicial_signal = self.s_depot.add_to_time_counter(self._level - 1)
        for _ in np.arange(self._level, level, -1):
            self.s_depot.candidates[_] = dict()
            self.s_depot.exempts[_] = dict()
            self.s_depot.collects[_] = dict()
        self._level = level
        self.size_list = self.s_depot.prev_s[self._level]
        self.degree_list = self.s_depot.prev_d[self._level]
        self.blocked_sets = self.s_depot.prev_b[self._level]
        self._level -= 1

    def _stage_state(self):
        self.s_depot.prev_s[self._level] = deepcopy(self.size_list)
        self.s_depot.prev_d[self._level] = deepcopy(self.degree_list)
        self.s_depot.prev_b[self._level] = deepcopy(self.blocked_sets)

    def _verbose_logging(self):
        if self.verbose:
            print(f"========== (l={self._level}) ==========\n"
                  f"Size list: {[_ for _ in self.size_list if _ > 0]}\n"
                  f"Degree list: {[_ for _ in self.degree_list if _ > 0]}\n"
                  f"blocked_sets = {self.blocked_sets}\n"
                  # f"currents = {self.s_depot.candidates}\n"
                  # f"exempts = {self.s_depot.exempts}\n"
                  # f"collects = {self.s_depot.collects}"
                  )
