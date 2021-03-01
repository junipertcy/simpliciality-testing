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

import numpy as np
from dataclasses import dataclass, field
from .custom_exceptions import NonSimplicialSignal
from collections import defaultdict, Counter
from functools import partial
from typing import List


@dataclass
class SimplicialDepot:
    degree_list: np.ndarray
    size_list: np.ndarray
    depths: List = field(repr=False, default_factory=list)
    conv_time: int = field(repr=False, default=0)
    cutoff: int = field(repr=False, init=False)

    def __post_init__(self):
        self.explored = defaultdict(list)
        self.time = np.zeros(max(len(self.degree_list), len(self.size_list)) + 1, np.int_)
        self.level_map = defaultdict(partial(np.ndarray, max(len(self.degree_list), len(self.size_list)), int))
        self.level_map[1].fill(-1)
        for ind, _ in enumerate(self.degree_list):
            self.level_map[0][ind] = ind

    def compute_level_map(self, level, mapping2shrinked):
        n = len(self.degree_list)
        if level > 1:
            for _ in range(n):
                vtx_current_view = self.level_map[level - 1][_]
                if vtx_current_view == -1:
                    self.level_map[level][_] = -1
                    continue
                try:
                    self.level_map[level][_] = mapping2shrinked[vtx_current_view]
                except KeyError:
                    self.level_map[level][_] = -1
        else:
            for _ in range(n):
                try:
                    self.level_map[1][_] = mapping2shrinked[_]
                except KeyError:
                    self.level_map[1][_] = -1

    def add_to_time_counter(self, level):
        self.time[level] += 1
        self.conv_time += 1
        self.depths += [level]
        if self.conv_time >= self.cutoff:
            raise NonSimplicialSignal


class SimplexRegistrar(object):
    def __init__(self):
        self.pointer = 0
        self.name2id = dict()
        self.id2name = dict()
        self.facet_size_per_id = np.array([], dtype=np.int_)
        self.logbook = dict()

    def register(self, name) -> (tuple, int):
        name = tuple(sorted(name, reverse=True))
        if name not in self.name2id:
            self.name2id[name] = self.pointer
            self.id2name[self.pointer] = name
            self.pointer += 1
            self.facet_size_per_id = np.append(self.facet_size_per_id, [len(name)])
        return name, self.name2id[name]

    def log_forbidden(self, name, reason) -> None:
        self.logbook[tuple(name)] = {
            "is_simplicial": False,
            "reason": reason
        }


def flatten(nested_list):
    f"""flattens out nested lists."""
    return [item for sublist in nested_list for item in sublist]


def simplify_blocked_sets(bsets):
    r"""No more inclusive blocked sets.

    Parameters
    ----------
    bsets : ``iterable`` of integer-valued ``tuple`` objects
        Facets that must be respected for non-inclusion, or "blocked sets (of facets)".

    Returns
    -------
    simplified_bsets : ``list`` of integer-valued ``tuple`` objects
        Simplified blocked sets, of which redundancy is removed.

    Examples
    --------
    >>> bsets = ((7, 44, 109, 273), (7, 273), (44, 273), (109, 273), (2,), (4,), (44,), (56,))
    >>> simplified_bsets = simplify_blocked_sets(bsets)
    >>> print(simplified_bsets)
    [(7, 44, 109, 273), (2,), (4,), (56,)]

    """
    simplified_bsets = []
    for _bsets in sorted(bsets, key=lambda _: -len(_)):
        if not tuple(_bsets) in simplified_bsets:
            if len(simplified_bsets) > 0:
                to_add = True
                for _ in simplified_bsets:
                    if set(_bsets).issubset(set(_)):
                        to_add = False
                        break
                if to_add:
                    simplified_bsets += [_bsets]
            else:
                simplified_bsets += [_bsets]
    return simplified_bsets


def get_indices_of_k_in_blocked_sets(blocked_sets, k):
    indices = []
    for idx, _ in enumerate(blocked_sets):
        if k in _:
            indices += [idx]
    return indices


def transform_facets(facets, mapping, to="l+1") -> list:
    if to == "l-1":
        return list(map(lambda x: tuple(map(lambda y: mapping[y], x)), facets))
    elif to == "l+1":
        _facets = []
        for facet in facets:
            transformed_facet = []
            for vtx in facet:
                try:
                    transformed_facet += [mapping[vtx]]
                except KeyError:
                    pass
            if len(transformed_facet) > 0:
                _facets += [transformed_facet]
        return _facets
    else:
        raise ValueError(f"data content to={to} not understood")


def remove_ones(sizes, wanting_degs, choose_from=None):
    removed_vtx_sites = np.array(list(choose_from), dtype=np.int_)[
                        :Counter(sizes)[1]]  # todo, to figure out: you cannot do [-Counter(s)[1]:]
    wanting_degs[removed_vtx_sites] = 0
    return wanting_degs, removed_vtx_sites.tolist()


def pair_one_by_one(size_list, degree_list) -> (list, list):
    """
    This function is used to pair up the ones in size/deg lists, since this is the only plausible situation.
    Note that it is only applied when level=0.

    Parameters
    ----------
    size_list
    degree_list

    Returns
    -------

    """
    size_list = list(size_list)
    degree_list = list(degree_list)
    _ = min(Counter(size_list)[1], Counter(degree_list)[1])
    for __ in range(_):
        size_list.remove(1)
        degree_list.remove(1)
    return size_list, degree_list, _


def sort_helper(st) -> list:
    d = defaultdict()
    for idx, _ in enumerate(st.DEGREE_LIST):
        d[idx] = _

    inv_map = defaultdict(list)

    for _ in d:
        inv_map[d[_]] += [_]
    joint_seq = st.compute_joint_seq_from_identifier(sorted_deg=False)
    if len(joint_seq[0]) == 0:
        raise ValueError("Invalid input instance. Perhaps execute SimplicialTest.is_simplicial() first?")
    old2new = dict()
    for idx, _ in enumerate(st.compute_joint_seq_from_identifier(sorted_deg=False)[1]):
        old2new[idx] = inv_map[_].pop(0)

    translated = []
    for facet in st.identifier2facets():
        translated += [sorted([old2new[_] for _ in facet])]
    return translated


def sort_callback(facets):
    t = []
    for facet in facets:
        t += [tuple(sorted(facet, reverse=False))]
    return t


def get_seq2seq_mapping(degs):
    mapping2enlarged = dict()
    nonzero_degs = np.count_nonzero(degs)
    argsort = np.argsort(degs)[::-1]
    for idx in np.arange(0, nonzero_degs, 1):
        mapping2enlarged[idx] = argsort[idx]
    mapping2shrinked = {v: k for k, v in mapping2enlarged.items()}
    return mapping2shrinked, mapping2enlarged


def filter_blocked_facets(blocked_facets, exempt_vids):
    r"""Supposedly every existing facet will block potential facets in the next level; however, this only applies if
    it contains `exempt_vids` because these vertices will be shared by next-level candidates.

    Parameters
    ----------
    blocked_facets
    exempt_vids

    Returns
    -------

    """
    filtered = []
    for facet in blocked_facets:
        if set(exempt_vids).issubset(facet):
            filtered += [facet]
    return filtered


def compute_joint_seq(facets) -> (list, list):
    if len(facets) == 0:
        raise ValueError("Empty input facets. Are you sure the input is correct?")
    flattened_facets = [item for sublist in facets for item in sublist]
    n = max(flattened_facets) - min(flattened_facets) + 1
    degs = np.zeros(n, dtype=np.int_)
    sizes = []

    for facet in facets:
        sizes += [len(facet)]
        for vertex_id in facet:
            degs[vertex_id] += 1

    return sorted(sizes, reverse=True), sorted(degs, reverse=True)


def if_facets_simplicial(facets) -> bool:
    for idx1, facet1 in enumerate(sorted(facets, key=len)):
        for idx2, facet2 in enumerate(sorted(facets, key=len)):
            if idx2 > idx1:
                if set(list(facet1)).issubset(set(list(facet2))):
                    return False
    return True


def accel_asc(n):
    """from: http://jeromekelleher.net/generating-integer-partitions.html"""
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def get_partition(n=1, sortby="asc"):
    gen = accel_asc(n)
    partitions = []
    while True:
        try:
            partitions += [sorted(next(gen), reverse=True)]
        except StopIteration:
            break
    if sortby == "asc":
        return partitions
    else:
        return sorted(partitions, key=lambda x: max(x) - len(x), reverse=True)


def read_hyperedge_list(path, delimiter=","):
    with open(path, 'r') as f:
        edge_list = set()
        for line in f:
            if not line.lstrip().startswith("%"):
                e = line.strip().split(delimiter)
                _edge_list = list()
                for _e in e:
                    _edge_list += [int(_e) - 1]
                edge_list.add(tuple(_edge_list))
        return list(edge_list)


def write_simplicial_list(el, path):
    with open(path, 'w') as f:
        for _el in el:
            for __el in _el:
                f.write(str(__el) + " ")
            f.write("\n")

        # return sorted(partitions, key=lambda x: [Counter(x)[1], max(x) - len(x)] + list(x), reverse=True)

# def get_slimmest_d(s):
#     """
#
#     Parameters
#     ----------
#     s: size list
#
#     Returns
#     -------
#
#     """
#     s = sorted(s)
#     pool = set()
#     tentative_ = []
#     for _s in s:
#         tentative = tentative_
#         if len(tentative) == 0:
#             idx = 0
#             for _ in range(_s):
#                 tentative += [idx]
#                 idx += 1
#             pool.add(tuple(tentative))
#             tentative_ = tentative
#             continue
#         tentative[-1] += 1
#         idx = tentative[-1]
#         for _ in range(_s - len(tentative)):
#             idx += 1
#             tentative += [idx]
#
#         pool.add(tuple(tentative))
#         tentative_ = tentative
#     return sorted(Counter(flatten(pool)).values(), reverse=True)

# def update_deg_seq(deg_seq, facet, value):
#     if value not in [+1, -1]:
#         raise NotImplementedError
#     for _ in facet:
#         deg_seq[_] += value


# def shrink_degs(degs, inv_map):
#     d = dict()
#     for idx, _degs in enumerate(degs):
#         try:
#             d[inv_map[idx]] = _degs
#         except KeyError:
#             pass
#     d_list = list()
#     for idx, _d in enumerate(d.values()):
#         d_list += [d[idx]]
#     return np.array(d_list, dtype=np.int_)
