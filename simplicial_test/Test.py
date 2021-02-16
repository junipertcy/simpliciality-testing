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
from .custom_exceptions import *
from functools import partial
from itertools import combinations

import sys

sys.setrecursionlimit(10 ** 6)


class Test(SimplexRegistrar):
    r"""Base class for SimplicialCheck.

    Parameters
    ----------
    degree_list : ``iterable`` or :class:`numpy.ndarray`, required

    size_list : ``iterable`` or :class:`numpy.ndarray`, required

    """

    def __init__(self, degree_list, size_list, blocked_sets=None, level_map=None, verbose=False, **kwargs):

        super().__init__()
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

        self._sorted_s = np.array(deepcopy(self.SIZE_LIST), dtype=np.int_)
        self._sorted_d = np.array(deepcopy(self.DEGREE_LIST), dtype=np.int_)

        self.current_fids = list()

        if blocked_sets is None:
            self.blocked_sets = list()
        else:
            self.blocked_sets = simplify_blocked_sets(blocked_sets)
        self.verbose = verbose

        self.callback_data = dict()  # level: facets

        self.level_map = level_map
        if self.level_map is None:
            self.level_map = defaultdict(partial(np.ndarray, max(len(degree_list), len(size_list)), int))
            self.level_map[-1] = self.DEGREE_LIST
            self.level_map[self._level].fill(-1)
            self.level_map["time"].fill(0)
            self.level_map["explored"] = defaultdict(list)
            for ind, _ in enumerate(degree_list):
                self.level_map[0][ind] = ind

        self.summary = {}
        if verbose:
            print(f"========== (l={self._level}) ==========\n"
                  f"Size list: {self.SIZE_LIST.tolist()}\n"
                  f"Degree list: {self.DEGREE_LIST.tolist()}\n"
                  f"blocked_sets = {self.blocked_sets}"
                  # f"len::blocked_sets = {list(map(lambda x: len(x), self.blocked_sets))}\n"
                  )

    def identifier2facets(self):
        facets = []
        for _id in self.current_fids:
            facets += [self.id2name[_id]]
        return facets

    def get_distinct_selection(self, size):
        r"""Select N=size vertices as the candidate facet.

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
        level_map = self.level_map[self._level - 1]
        for vid, _ in enumerate(self.level_map[-1]):
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
                self.level_map["time"][self._level - 1] += 1
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
                return
            else:
                self.level_map["time"][self._level - 1] += 1

    def validate(self, facet) -> (bool, str):
        r"""This function must return True in order for the candidate facet to be considered.

        Parameters
        ----------
        facet : ``tuple``
            The candidate facet.

        Returns
        -------
        NoMoreBalls :
            if we depleted this level and wanted to go up (i.e., from `lv` to `lv - 1`).

        Notes
        -----

        """
        token = validators.simple_validate(self._sorted_d, self._sorted_s, facet)
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
                        raise WeCanStopSignal
                self.callback_data[self._level] = current_facets + [exempt_vids] + collected_facets
                raise WeCanStopSignal
        else:
            exempt_vids = []
            collected_facets = []
        filtered = filter_blocked_facets(current_facets + blocked_sets, exempt_vids)
        mapping2shrinked, mapping2enlarged = get_seq2seq_mapping(wanting_degs)
        blocked_sets = shrink_facets(filtered, mapping2shrinked)

        try:
            deeper_facet_found, facets = self._is_simplicial(
                wanting_degs, sizes, mapping2shrinked, mapping2enlarged,
                blocked_sets=blocked_sets,
                level_map=self.level_map,
                verbose=self.verbose
            )
        except NoMoreBalls as e:
            lv, (self.level_map["time"], self.level_map["explored"]) = e.message
            self.level_map["time"][self._level - 1] += 1
            if int(lv) <= self.depth:
                if self.level_map["time"][self._level - 1] > self.width:
                    raise NoMoreBalls
                # keep look for solutions at this level, until get_valid_trials emits a NoMoreBalls
                return False, "keep looking for balls!"
            else:
                # forget about this level
                raise NoMoreBalls
        if deeper_facet_found:
            facets = [exempt_vids + list(s) for s in facets]
            self.callback_data[self._level] = current_facets + facets + collected_facets
            raise WeCanStopSignal
        else:
            self.level_map["time"][self._level - 1] += 1
            return False, "No valid facet found at a higher level."

    def _is_simplicial(self, degs, sizes, mapping2shrinked, mapping2enlarged,
                       blocked_sets=None, level_map=None, verbose=False) -> (bool, list):
        lv = self._level
        self.level_map = compute_level_map(self.level_map, lv, mapping2shrinked)
        blocked_sets = sort_facets(blocked_sets)
        _ = (tuple(sorted(degs, reverse=True)), tuple(blocked_sets))
        if _ in self.level_map["explored"][lv]:
            raise NoMoreBalls((lv, (self.level_map["time"], self.level_map["explored"])))
        self.level_map["explored"][lv] += [_]
        st = Test(degs, sizes, blocked_sets=blocked_sets, level_map=level_map, verbose=verbose, **self.kwargs)
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
        if np.sum(self.level_map["time"]) > self.kwargs.get("cutoff", np.infty):
            self.level_map["time"][lv - 1] += 1
            return False, self.__mark(False, tuple())["facets"]
        if not validators.validate_data(self._sorted_d, self._sorted_s):
            self.level_map["time"][lv - 1] += 1
            return False, self.__mark(False, tuple())["facets"]
        else:
            if sum(self._sorted_d) == sum(self._sorted_s) == 0:
                self.callback_data[lv] = self.identifier2facets()
                if lv - 1 == 0:
                    return True, self.__mark(True, tuple(self.facets_to_append))["facets"]
                return True, self.callback_data
        while True:
            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self._sorted_s) == 0:
                self.callback_data[lv] = self.identifier2facets()
                if lv - 1 == 0:  # the very first level
                    facets = tuple(sort_callback(self.callback_data[lv]) + self.facets_to_append)
                    return True, self.__mark(True, facets)["facets"]
                return True, self.callback_data
            s = self._sorted_s[0]
            self._sorted_s = np.delete(self._sorted_s, 0)
            try:
                self.sample_simplex_greedy(s)
            except WeCanStopSignal:
                if lv - 1 == 0:  # the very first level
                    facets = tuple(sort_callback(self.callback_data[lv]) + self.facets_to_append)
                    return True, self.__mark(True, facets)["facets"]
                return True, self.callback_data
            except NoMoreBalls:
                if lv - 1 == 0:
                    return False, self.__mark(False, tuple())["facets"]
                raise NoMoreBalls((self.level_map["time"], self.level_map["explored"]))

    def __mark(self, simplicial, facets):
        self.summary["simplicial"] = simplicial
        self.summary["conv_time"] = sum(self.level_map["time"])
        self.summary["time"] = tuple(self.level_map["time"])

        self.summary["degs"] = tuple(self.DEGREE_LIST)
        self.summary["sizes"] = tuple(self.SIZE_LIST)
        self.summary["facets"] = facets
        return self.summary
