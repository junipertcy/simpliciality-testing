#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# simplicial_test -- a python module to verify simplicial sequences
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
from .enumeration import sort_facets
from .utils import *
from .custom_exceptions import NoMoreBalls, SimplicialSignal, GoToNextLevel
from itertools import combinations
from copy import deepcopy


class Test(SimplexRegistrar):
    r"""Base class for SimplicialCheck.

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

        self.depth = kwargs.pop("depth", 1e2)
        self.width = kwargs.pop("width", 1e2)
        self.s_depot.cutoff = kwargs.pop("cutoff", np.infty)
        self.verbose = verbose
        self.summary = {}
        self.current_fids = []
        if len(kwargs) > 0:
            raise ValueError(f"unrecognized keyword arguments: {str(list(kwargs.keys()))}")

    def fids2facets(self):
        r"""Translate current selected facet ids to facets."""
        facets = []
        for _id in self.current_fids:
            facets += [self.id2name[_id]]
        return facets

    @staticmethod
    def get_distinct_selection(size, degs, blocked_sets, level_map):
        r"""Select a candidate facet (of size ``size``).

        Parameters
        ----------
        size : ``int``
            Number of vertices to make the facet.

        Returns
        -------
        facet : ``tuple``
            The candidate facet.

        Notes
        -----

        Note that _ = vtx degree at original view, vid = vtx index at original view,
        level_map[vid] =  vtx index at current view.

        """
        equiv2vid = defaultdict(list)
        blocked_sets = [_ for _ in blocked_sets if len(_) >= size]
        for vid, _ in enumerate(degs):
            if level_map[vid] != -1:
                key = tuple(get_indices_of_k_in_blocked_sets(blocked_sets, level_map[vid]) + [_])
                equiv2vid[key] += [level_map[vid]]
                equiv2vid["pool"] += [key]
        equiv_class_pool = combinations(equiv2vid["pool"], size)
        explored_set = set()
        while True:
            facet = []
            tracker = defaultdict(int)
            for equiv_class in next(equiv_class_pool):
                facet += [equiv2vid[equiv_class][tracker[equiv_class]]]  # vids_same_equiv_class[tracker[equiv_class]]
                tracker[equiv_class] += 1
            if tuple(facet) not in explored_set:
                explored_set.add(tuple(facet))
                yield tuple(facet)

    def sample_candidate_facet(self, size, valid_trials=None):
        if not valid_trials:
            self.s_depot.valid_trials[self._level] = self.get_distinct_selection(
                size, self.s_depot.prev_d[1], self.blocked_sets, self.s_depot.level_map[self._level - 1])
        while True:
            try:
                facet = next(self.s_depot.valid_trials[self._level])  # candidate_facet
            except RuntimeError:
                self.s_depot.add_to_time_counter(self._level)
                raise NoMoreBalls("nomore")
            except StopIteration:
                self.s_depot.add_to_time_counter(self._level)
                raise NoMoreBalls("nomore")
            if validators.validate_issubset_blocked_sets(facet, self.blocked_sets):
                continue
            ind, reason = self.validate(facet)
            if ind:
                picked_facet, picked_facet_id = self.register(facet)
                self.current_fids += [picked_facet_id]
                return picked_facet
            else:
                self.s_depot.add_to_time_counter(self._level)
                if reason == "explored":
                    raise NoMoreBalls

    def validate(self, facet) -> (bool, str):
        r"""This function must return True in order for the candidate facet to be considered.

        Parameters
        ----------
        facet : ``tuple``
            The candidate facet.

        Returns
        -------
       token : ``(bool, str)``
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
        token = validators.simple_validate(self.degree_list, self.size_list, facet)
        if type(token[0]) == bool:
            return token
        else:
            wanting_degs, sizes, current_facets = token
        blocked_sets = simplify_blocked_sets(self.blocked_sets)
        for blocked_set in blocked_sets:
            if set(np.nonzero(wanting_degs)[0]).issubset(set(blocked_set)):  # useful
                return False, "Rejected b/c the remaining facets are doomed to fail the no-inclusion constraint."
        if Counter(wanting_degs)[len(sizes)] != 0:
            _ = validators.validate_reduced_seq(wanting_degs, sizes, current_facets, blocked_sets, verbose=self.verbose)
            if _[0]:  # useful
                return False, "Rejected while reducing the sequences"
            wanting_degs, sizes, collected_facets, exempt_vids = _[1]
            if np.sum(wanting_degs) == np.sum(sizes) == 0:
                for collected_facet in collected_facets:
                    if set(exempt_vids).issubset(collected_facet):
                        raise SimplicialSignal(current_facets + collected_facets)
                raise SimplicialSignal(current_facets + [exempt_vids] + collected_facets)
        else:
            exempt_vids = []
            collected_facets = []

        filtered = filter_blocked_facets(current_facets + blocked_sets, exempt_vids)
        blocked_sets = sort_facets(self.blocked_sets)
        _ = (tuple(sorted(wanting_degs, reverse=True)), tuple(blocked_sets))
        if _ in self.s_depot.explored[self._level]:
            return False, "explored"
        self.s_depot.explored[self._level].add(_)

        self.s_depot.mappers[self._level] = get_seq2seq_mapping(wanting_degs)  # (mapping2shrinked, mapping2enlarged)
        self.s_depot.collects[self._level] = collected_facets
        self.s_depot.exempts[self._level] = exempt_vids
        self.s_depot.currents[self._level] = current_facets
        blocked_sets = transform_facets(filtered, self.s_depot.mappers[self._level][0], to="l+1")
        raise GoToNextLevel(wanting_degs, sizes, blocked_sets, self.s_depot.mappers[self._level][0])

    def is_simplicial(self):
        ind = None
        facets = None
        anew = True
        s = None
        while True:
            if anew:
                self._level += 1
                self._verbose_logging()
                if not validators.validate_data(self.degree_list, self.size_list):
                    self.s_depot.add_to_time_counter(self._level - 1)
                    if self._level != 1:
                        anew = True
                        self._rollback(self._level - 1)
                        continue  # todo: unsure
                    ind = 0
                    break
                self._stage_prev_state()
                s = self.size_list.pop(0)
            try:
                self.sample_candidate_facet(s, valid_trials=self.s_depot.valid_trials[self._level])
            except NonSimplicialSignal:
                ind = 0
                break
            except SimplicialSignal as e:
                msg = e.message
                if self._level != 1:
                    for lv in np.arange(self._level - 1, 0, -1):
                        facets = transform_facets(msg, self.s_depot.mappers[lv][1], to="l-1")
                        facets = [self.s_depot.exempts[lv] + list(s) for s in facets]
                        msg = self.s_depot.currents[lv] + facets + self.s_depot.collects[lv]
                facets = tuple(sort_callback(msg) + self.facets_to_append)
                ind = 1
                break
            except NoMoreBalls as e:
                if self._level == 1:
                    ind = 0
                    break
                if e.message == tuple(["nomore"]):
                    anew = True
                    self.s_depot.valid_trials[self._level] = None
                    self._rollback(self._level - 1)
                    continue
                else:
                    if self._level <= self.depth and self.s_depot.time[self._level] <= self.width:
                        anew = False
                        continue
                    anew = True
                    self.s_depot.valid_trials[self._level] = None
                    self._rollback(self._level - 1)
                    continue
            except GoToNextLevel as e:
                degs, sizes, blocked_sets, mapping2shrinked = e.message
                self.s_depot.compute_level_map(self._level, self.s_depot.mappers[self._level][0])
                self.degree_list = sorted(degs, reverse=True)
                self.size_list = sizes
                self.blocked_sets = blocked_sets
                anew = True
                continue
            else:
                if sum(self.size_list) == 0:
                    msg = self.fids2facets()
                    if self._level != 1:
                        for lv in np.arange(self._level, 0, -1):
                            facets = transform_facets(msg, self.s_depot.mappers[lv][1], to="l-1")
                            facets = [self.s_depot.exempts[lv - 1] + list(s) for s in facets]
                            msg = self.s_depot.currents[lv - 1] + facets + self.s_depot.collects[lv - 1]
                    facets = tuple(sort_callback(msg) + self.facets_to_append)
                    ind = 1
                    break
                self._level -= 1
                anew = True
                pass

        if ind == 0:
            return False, self.__mark(False, tuple())["facets"]
        elif ind == 1:
            return True, self.__mark(True, facets)["facets"]

    def __mark(self, simplicial, facets):
        self.summary["simplicial"] = simplicial
        self.summary["conv_time"] = sum(self.s_depot.time)
        self.summary["time"] = tuple(self.s_depot.time)
        self.summary["degs"] = tuple(self.s_depot.degree_list)
        self.summary["sizes"] = tuple(self.s_depot.size_list)
        self.summary["facets"] = facets
        return self.summary

    def _rollback(self, level):
        self.s_depot.add_to_time_counter(self._level)
        for _ in np.arange(self._level, level, -1):
            self.s_depot.currents[_] = dict()
            self.s_depot.exempts[_] = dict()
            self.s_depot.collects[_] = dict()
        self._level = level
        self.size_list = self.s_depot.prev_s[self._level]
        self.degree_list = self.s_depot.prev_d[self._level]
        self.blocked_sets = self.s_depot.prev_b[self._level]
        self._level -= 1

    def _stage_prev_state(self):
        self.s_depot.prev_s[self._level] = deepcopy(self.size_list)
        self.s_depot.prev_d[self._level] = deepcopy(self.degree_list)
        self.s_depot.prev_b[self._level] = deepcopy(self.blocked_sets)

    def _verbose_logging(self):
        if self.verbose:
            print(f"========== (l={self._level}) ==========\n"
                  f"Size list: {self.size_list}\n"
                  f"Degree list: {self.degree_list}\n"
                  f"blocked_sets = {self.blocked_sets}\n"
                  f"currents = {self.s_depot.currents}\n"
                  f"exempts = {self.s_depot.exempts}\n"
                  f"collects = {self.s_depot.collects}"
                  )
