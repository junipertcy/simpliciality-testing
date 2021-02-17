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
from .custom_exceptions import NoMoreBalls, SimplicialSignal
from itertools import combinations

import sys

sys.setrecursionlimit(10 ** 6)


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

    def __init__(self, degree_list, size_list, blocked_sets=None, verbose=False, **kwargs):
        super().__init__()
        if blocked_sets is None:
            blocked_sets = list(tuple())
        self.kwargs = kwargs.copy()
        self.depth = kwargs.pop("depth", 1e2)
        self.width = kwargs.pop("width", 1e2)
        trim_ones_first = kwargs.pop("trim_ones_first", True)
        self._level = kwargs.pop("level", 1)
        self.kwargs["level"] = self._level + 1
        self.facets_to_append = []
        if self._level == 1:
            if trim_ones_first:
                size_list, degree_list, num_ones = pair_one_by_one(size_list, degree_list)
                for _ in range(num_ones):
                    self.facets_to_append += [(len(degree_list) + _,)]
        self.m = len(size_list)
        self.n = len(degree_list)
        self.DEGREE_LIST = np.array(sorted(degree_list, reverse=True), dtype=np.int_)
        self.SIZE_LIST = np.array(sorted(size_list, reverse=True), dtype=np.int_)
        self._mutable_size_list = np.array(deepcopy(self.SIZE_LIST), dtype=np.int_)

        if blocked_sets is None:
            self.blocked_sets = []
        else:
            self.blocked_sets = simplify_blocked_sets(blocked_sets)

        self.s_depot = kwargs.pop("s_depot", SimplicialDepot(self.DEGREE_LIST, self.SIZE_LIST))
        self.s_depot.cutoff = kwargs.pop("cutoff", np.infty)

        self.verbose = verbose
        if self.verbose:
            print(f"========== (l={self._level}) ==========\n"
                  f"Size list: {self.SIZE_LIST.tolist()}\n"
                  f"Degree list: {self.DEGREE_LIST.tolist()}\n"
                  f"blocked_sets = {self.blocked_sets}"
                  # f"len::blocked_sets = {list(map(lambda x: len(x), self.blocked_sets))}\n"
                  )
        self.summary = {}
        self.current_fids = []
        self.callback_data = {}  # level2facets map
        if len(kwargs) > 0:
            raise ValueError(f"unrecognized keyword arguments: {str(list(kwargs.keys()))}")

    def identifier2facets(self):
        facets = []
        for _id in self.current_fids:
            facets += [self.id2name[_id]]
        return facets

    def get_distinct_selection(self, size):
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
        blocked_sets = [_ for _ in self.blocked_sets if len(_) >= size]
        level_map = self.s_depot.level_map[self._level - 1]
        for vid, _ in enumerate(self.s_depot.degree_list):
            if level_map[vid] != -1:
                key = tuple(get_indices_of_k_in_blocked_sets(blocked_sets, level_map[vid]) + [_])
                equiv2vid[key] += [level_map[vid]]
                equiv2vid["pool"] += [key]
        equiv_class_pool = combinations(equiv2vid["pool"], size)
        while True:
            facet = []
            tracker = defaultdict(int)
            for equiv_class in next(equiv_class_pool):
                facet += [equiv2vid[equiv_class][tracker[equiv_class]]]  # vids_same_equiv_class[tracker[equiv_class]]
                tracker[equiv_class] += 1
            yield tuple(facet)

    def sample_simplex_greedy(self, size):
        valid_trials = self.get_distinct_selection(size)
        while True:
            try:
                facet = next(valid_trials)  # candidate_facet
            except RuntimeError:
                self.s_depot.add_to_time_counter(self._level)
                raise NoMoreBalls(
                    "May not be solvable with the greedy algorithm OR the recursive_is_simplicial must return False.")
            # Note that self.blocked_sets are passed from a higher level,
            # whereas larger_simplices is anew for each level.
            if validators.validate_issubset_blocked_sets(facet, self.blocked_sets):
                continue
            ind, reason = self.validate(facet)
            if ind:
                picked_facet, picked_facet_id = self.register(facet)
                self.current_fids += [picked_facet_id]
                return picked_facet
            else:
                self.s_depot.add_to_time_counter(self._level)

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
        token = validators.simple_validate(self.DEGREE_LIST, self._mutable_size_list, facet)
        if type(token[0]) == bool:
            return token
        else:
            wanting_degs, sizes, current_facets = token
        blocked_sets = self.blocked_sets
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
                        self.callback_data[self._level] = current_facets + collected_facets
                        raise SimplicialSignal
                self.callback_data[self._level] = current_facets + [exempt_vids] + collected_facets
                raise SimplicialSignal
        else:
            exempt_vids = []
            collected_facets = []
        filtered = filter_blocked_facets(current_facets + blocked_sets, exempt_vids)
        mapping2shrinked, mapping2enlarged = get_seq2seq_mapping(wanting_degs)
        blocked_sets = shrink_facets(filtered, mapping2shrinked)

        try:
            deeper_facet_found, facets = self._is_simplicial(
                wanting_degs, sizes, mapping2shrinked, mapping2enlarged, blocked_sets=blocked_sets
            )
        except NoMoreBalls as e:
            lv, (self.s_depot.time, self.s_depot.explored) = e.message
            self.s_depot.add_to_time_counter(self._level)
            if int(lv) <= self.depth:
                if self.s_depot.time[self._level] > self.width:
                    raise NoMoreBalls
                # keep look for solutions at this level, until get_valid_trials emits a NoMoreBalls
                return False, "keep looking for balls!"
            else:
                # forget about this level
                raise NoMoreBalls

        if deeper_facet_found:
            facets = [exempt_vids + list(s) for s in facets]
            self.callback_data[self._level] = current_facets + facets + collected_facets
            raise SimplicialSignal
        else:
            self.s_depot.add_to_time_counter(self._level)
            return False, "No valid facet found at a higher level."

    def _is_simplicial(self, degs, sizes, mapping2shrinked, mapping2enlarged, blocked_sets=None) -> (bool, list):
        lv = self._level
        self.s_depot.compute_level_map(lv, mapping2shrinked)
        blocked_sets = sort_facets(blocked_sets)
        _ = (tuple(sorted(degs, reverse=True)), tuple(blocked_sets))
        if _ in self.s_depot.explored[lv]:
            raise NoMoreBalls((lv, (self.s_depot.time, self.s_depot.explored)))
        self.s_depot.explored[lv] += [_]
        self.kwargs["s_depot"] = self.s_depot
        st = Test(degs, sizes, blocked_sets=blocked_sets, verbose=self.verbose, **self.kwargs)
        try:
            deeper_facet_is_simplicial, cb_data = st.is_simplicial()  # cb_data contains facets from a deeper level
        except NoMoreBalls as e:
            raise NoMoreBalls((lv, e.message))

        if deeper_facet_is_simplicial:
            # transform the facets collected from a deeper level
            return True, get_mapped_seq(mapping2enlarged, cb_data[lv + 1])
        else:
            return False, list()

    def is_simplicial(self):
        lv = self._level
        if not validators.validate_data(self.DEGREE_LIST, self._mutable_size_list):
            self.s_depot.add_to_time_counter(self._level)
            return False, self.__mark(False, tuple())["facets"]
        else:
            if sum(self.DEGREE_LIST) == sum(self._mutable_size_list) == 0:
                self.callback_data[lv] = self.identifier2facets()
                if lv - 1 == 0:
                    return True, self.__mark(True, tuple(self.facets_to_append))["facets"]
                return True, self.callback_data

        while True:
            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self._mutable_size_list) == 0:
                self.callback_data[lv] = self.identifier2facets()
                if lv - 1 == 0:  # the very first level
                    facets = tuple(sort_callback(self.callback_data[lv]) + self.facets_to_append)
                    return True, self.__mark(True, facets)["facets"]
                return True, self.callback_data

            s = self._mutable_size_list[0]
            self._mutable_size_list = np.delete(self._mutable_size_list, 0)
            try:
                self.sample_simplex_greedy(s)
            except SimplicialSignal:
                if lv - 1 == 0:  # the very first level
                    facets = tuple(sort_callback(self.callback_data[lv]) + self.facets_to_append)
                    return True, self.__mark(True, facets)["facets"]
                return True, self.callback_data
            except NonSimplicialSignal:
                if lv - 1 == 0:  # the very first level
                    return False, self.__mark(False, tuple())["facets"]
                raise NonSimplicialSignal
            except NoMoreBalls:
                if lv - 1 == 0:
                    return False, self.__mark(False, tuple())["facets"]
                raise NoMoreBalls((self.s_depot.time, self.s_depot.explored))

    def __mark(self, simplicial, facets):
        self.summary["simplicial"] = simplicial
        self.summary["conv_time"] = sum(self.s_depot.time)
        self.summary["time"] = tuple(self.s_depot.time)

        self.summary["degs"] = tuple(self.DEGREE_LIST)
        self.summary["sizes"] = tuple(self.SIZE_LIST)
        self.summary["facets"] = facets
        return self.summary
