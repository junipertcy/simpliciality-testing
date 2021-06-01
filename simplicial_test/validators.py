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

from .utils import *


def preprocess(sorted_d, sorted_s):
    if np.max(sorted_d) > len(sorted_s):
        return False
    if len(sorted_s) > 1 and np.max(sorted_s) >= len(sorted_d):
        return False
    if np.sum(sorted_d) != np.sum(sorted_s):
        return False
    if Counter(sorted_s)[1] > Counter(sorted_d)[1]:  # TODO: perhaps we could only enforce this at the very first level.
        return False

    _sorted_s = []
    counter = 0
    for _ in sorted_s:
        if counter < max(sorted_d):
            _sorted_s += [_ - 1]
            counter += 1
        else:
            _sorted_s += [_]

    if Counter(_sorted_s)[1] > len([_ for _ in sorted_d if _ > 0]) - 1:
        return False

    return True


def get_wanting_slots(degs, facet):
    """
    Used in the greedy case only.

    Parameters
    ----------
    degs: wanting degrees
    facet: candidate facet

    Returns
    -------

    """
    n = len(degs)
    l_1 = [degs[_] - 1 if _ in set(facet) else degs[_] for _ in range(n)]
    l_2 = [0 if _ in set(facet) else degs[_] for _ in range(n)]
    return np.array(l_1, dtype=np.int_), np.array(l_2, dtype=np.int_)


def simple_validate(degs, sizes, facet):
    wanting_degs, non_shielding = get_wanting_slots(degs, facet)  # wanting_degs (non_shielding & shielding)
    if min(wanting_degs) < 0:
        return False, "Rejected b/c added in min(wanting_degs) < 0"
    elif len(sizes) == 1:
        for _ in [facet]:
            if set(np.nonzero(wanting_degs)[0]).issubset(set(_)):
                return False, "Second last facet explored. " \
                              "But the very last facet is doomed to fail the no-inclusion constraint."
    if np.any(wanting_degs > len(sizes)):
        return False, "Some degrees require more facets than the actual remaining number."
    if np.count_nonzero(wanting_degs) < sizes[0]:
        return False, "Rejected b/c np.count_nonzero(wanting_degs) < sizes[0]."
    elif np.count_nonzero(wanting_degs) == sizes[0]:
        if np.sum(wanting_degs) > np.count_nonzero(wanting_degs):
            return False, "Rejected b/c np.sum(wanting_degs) > np.count_nonzero(wanting_degs)."
    if np.min(non_shielding) >= 0:
        if validate_nonshielding(sizes, wanting_degs, non_shielding):
            return False, "Rejected while validating non_shielding vertices."
    return True, wanting_degs


def validate_nonshielding(curent_sizes, wanting_degs, non_shielding):
    r"""Verify that the sum of non-shielding slots not be less than the number of the remaining facets.

    Parameters
    ----------
    curent_sizes
    wanting_degs
    non_shielding

    Returns
    -------

    """
    shielding = wanting_degs - non_shielding  # only shielding
    remaining = len(curent_sizes)
    if np.sum(non_shielding) < remaining:
        return True
    elif np.sum(non_shielding) == remaining:
        _ = [dummy - 1 for dummy in curent_sizes]
        if len(_) - 1 <= 0:  # only the last to-be-chosen facet remains
            return False
        if _[0] > np.count_nonzero(shielding):
            return True
    return False  # safe!


def basic_validations_degs_and_sizes(degs, sizes):
    if Counter(degs)[len(sizes)] == np.min(sizes):
        return False
    if len(degs) == np.max(sizes):
        return False
    return True


def validate_reduced_seq(wanting_degs, sizes, current_facets, blocked_sets, verbose=False) -> (bool, tuple):
    if verbose:
        print(f"----- (BEGIN: reducing seq) -----\n"
              f"Wanting degrees: {[_ for _ in wanting_degs if _ > 0]}\n"
              f"Size list: {sizes}\n"
              f"Current facets: {current_facets}\n"
              f"blocked_sets = {blocked_sets}"
              )
    collected_facets = []
    exempt_vids = []
    while Counter(wanting_degs)[len(sizes)] != 0:
        if len(sizes) > 1:
            if not basic_validations_degs_and_sizes(degs=wanting_degs, sizes=sizes):
                return False, tuple([None] * 4)
        if Counter(sizes)[0] != len(sizes):
            must_be_filled_vids = np.where(wanting_degs == len(sizes))[0].tolist()
            exempt_vids += must_be_filled_vids
        else:
            break
        sizes = [_ - Counter(wanting_degs)[len(sizes)] for _ in sizes]
        wanting_degs[wanting_degs == len(sizes)] = 0
        if min(sizes) < 0:
            return False, tuple([None] * 4)
        if Counter(sizes)[1] > 0:
            if Counter(sizes)[1] > Counter(wanting_degs)[1]:
                return False, tuple([None] * 4)
            shielding_facets = get_shielding_facets_when_vids_filled(
                current_facets, blocked_sets, must_be_filled_vids, exempt_vids=exempt_vids
            )
            shielding_facets = [set(_).difference(set(must_be_filled_vids)) for _ in shielding_facets]
            if len(shielding_facets) == 0:  # You are free to do "remove_ones" (at a later stage, perhaps)
                nonshielding_vids = set(np.nonzero(wanting_degs)[0])
            else:
                nonshielding_vids = get_nonshielding_vids(len(wanting_degs), shielding_facets)
            all_one_sites = np.where(wanting_degs == 1)[0]
            nonshielding_vids.intersection_update(set(all_one_sites))
            if not nonshielding_vids:
                return False, tuple([None] * 4)
            wanting_degs, removed_sites = remove_ones(sizes, wanting_degs, choose_from=nonshielding_vids)
            [sizes.remove(1) for _ in range(Counter(sizes)[1])]
            collected_facets += [exempt_vids + [_] for _ in removed_sites]  # collected_facets immer from removed_sites

    if np.sum(wanting_degs) == np.sum(sizes) == 0:
        if validate_issubset_blocked_sets(exempt_vids, blocked_sets=blocked_sets):
            return False, tuple([None] * 4)
    if verbose:
        print(f"----- (↓ Returning these data ↓) -----\n"
              f"Wanting degrees: {[_ for _ in wanting_degs if _ > 0]}\n"
              f"Size list: {sizes}\n"
              f"Collected facets: {collected_facets}\n"
              f"Exempt vertex ids: {exempt_vids}\n"
              f"----- (END: reducing seq) -----"
              )
    return True, (wanting_degs, sizes, collected_facets, exempt_vids)


def validate_issubset_blocked_sets(candidate_facet, blocked_sets=None):
    if blocked_sets is not None:
        for facet in blocked_sets:
            if set(candidate_facet).issubset(set(facet)):
                return True
    return False


def get_shielding_facets_when_vids_filled(current_facets, blocked_sets, must_be_filled_vids, exempt_vids=None):
    """
    TODO: this function can be further simplified, along with the function::validate_reduced_seq
    The function works when one have "must_be_filled_vids" -- it goes by searching already existing facets,
    And find out the slots that must not be chosen in order to avoid clashes.

    Parameters
    ----------
    current_facets
    blocked_sets
    must_be_filled_vids
    exempt_vids


    Returns
    -------

    """
    if exempt_vids is None:
        exempt_vids = []
    shielding_facets = []
    # if a facet contains these must_be_filled_vids (or 'mbfv')
    mbfv = set(must_be_filled_vids).union(set(exempt_vids))
    for facet in current_facets:  # for all existing facets
        if mbfv.issubset(set(facet)):
            shielding_facets += [facet]
    for facet in blocked_sets:  # for all existing facets
        if mbfv.issubset(set(facet)):
            shielding_facets += [facet]
    return shielding_facets  # then we must avoid the slots in these shielding_facets


def get_nonshielding_vids(n, shielding_facets):
    nonshielding_vids = set()
    for facet in shielding_facets:
        non_shielding_part = set(range(n)).difference(set(facet))  # non_shielding_part vertex_id's
        if len(nonshielding_vids) == 0:
            nonshielding_vids = non_shielding_part
        else:
            nonshielding_vids.intersection_update(non_shielding_part)
    return nonshielding_vids
