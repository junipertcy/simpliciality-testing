from itertools import combinations
import numpy as np
from collections import defaultdict
import random
from copy import deepcopy


class SimplicialTest(object):
    """Base class for SimplicialCheck.

    Parameters
    ----------
    degree_list : ``iterable`` or :class:`numpy.ndarray`, required

    size_list : ``iterable`` or :class:`numpy.ndarray`, required

    """

    def __init__(self, degree_list, size_list):
        # random.seed(42)
        self.size_list = size_list
        self.degree_list = degree_list
        self.facet2id = defaultdict(int)
        self.id2facet = defaultdict(tuple)
        self.available_facets = defaultdict(dict)
        self.vertex_id = defaultdict(int)
        self.logbook = defaultdict(dict)

        self.m = len(size_list)
        self.n = len(degree_list)

        self.sorted_s = sorted(size_list, reverse=True)
        self._sorted_s = deepcopy(self.sorted_s)
        self.sorted_d = sorted(degree_list, reverse=True)
        self._sorted_d = deepcopy(self.sorted_d)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.init_all_simplices()
        self.init_vertex_id()
        self.symmetry_breaker = None
        self.identifier = None
        pass

    def init_all_simplices(self):
        simplices = dict()
        all_simplices = []
        for m_ in range(self.m):
            simplices[m_] = list(map(tuple, combinations(range(self.n), m_ + 1)))
            all_simplices += simplices[m_]
        for _id, simplex in enumerate(all_simplices):
            self.facet2id[simplex] = _id
        self.id2facet = {v: k for k, v in self.facet2id.items()}
        self.available_facets[tuple()] = simplices

    def init_vertex_id(self):

        for _id in range(self.n):
            self.vertex_id[_id] = _id

    def fill_up_deg_seq(self, facet):
        for _ in facet:
            self.deg_seq[_] += 1

    def resume_deg_seq(self, facet):
        for _ in facet:
            self.deg_seq[_] -= 1

    def deg_checks(self):
        # check 1
        if max(self.degree_list) < max(self.deg_seq):
            return False
        else:
            return True

    def remove_inclusive_facets(self, facet, old_identifier, new_identifier):
        m = len(facet)
        facets = self.available_facets[old_identifier]
        self.available_facets[new_identifier] = defaultdict(dict)
        for m_ in range(m):
            self.available_facets[new_identifier][m_] = [i for i in facets[m_] if i not in combinations(facet, m_ + 1)]

    def write_to_available_facets(self, identifier):
        self.available_facets[tuple(identifier)]["sorted_deg_seq"] = sorted(self.deg_seq, reverse=True)
        self.available_facets[tuple(identifier)]["sorted_size_seq"] = sorted(
            [len(self.id2facet[_id]) for _id in identifier],
            reverse=True)
        self.available_facets[tuple(identifier)]["chosen"] = list()
        for _id in identifier:
            self.available_facets[tuple(identifier)]["chosen"] += [self.id2facet[_id]]

    def write_to_logbook(self, identifier):
        if self.logbook[tuple(identifier)] != {"is_simplicial": False}:
            for key in self.available_facets[tuple(identifier)]:
                if key not in ["sorted_deg_seq", "sorted_size_seq", "chosen", "why-explore-more"]:
                    self.logbook[tuple(identifier)][key] = list()
                    for facet in self.available_facets[tuple(identifier)][key]:
                        self.logbook[tuple(identifier)][key] += [self.facet2id[facet]]

    def _reduce_pool(self, identifier, facets_pool):
        _pool = []
        for facet in facets_pool:
            if not self.logbook[tuple(identifier + [self.facet2id[facet]])] == {"is_simplicial": False}:
                _pool += [facet]
        return _pool

    def _break_symmetry(self, s):
        picked_facet = random.choice(self.available_facets[tuple([])][s - 1])
        picked_facet_id = self.facet2id[picked_facet]
        self.fill_up_deg_seq(picked_facet)
        identifier = [picked_facet_id]
        self.symmetry_breaker = picked_facet_id
        self.remove_inclusive_facets(picked_facet, tuple([]), tuple(identifier))
        self.write_to_available_facets(identifier)
        return identifier

    def is_simplicial(self):
        if max(self.degree_list) > self.m:
            print("1. This can never be simplicial.")
            return False
        identifier = self._break_symmetry(self.sorted_s.pop(0))

        while True:
            self.write_to_logbook(identifier)
            s = self.sorted_s.pop(0)
            # This facets_pool contains facets that we are supposed to draw from...
            facets_pool = deepcopy(self.available_facets[tuple(identifier)][s - 1])

            # but we do not need to explore cases that lead to "impossible" combinations
            facets_pool = self._reduce_pool(identifier, facets_pool)
            if len(facets_pool) == 0:
                if identifier == [self.symmetry_breaker]:
                    return False
                self.logbook[tuple(identifier)] = {"is_simplicial": False}  # marking impossible facet combinations
                self.sorted_s = [s] + self.sorted_s

                unwanted_facet = self.id2facet[identifier.pop(-1)]
                self.sorted_s = [len(unwanted_facet)] + self.sorted_s
                self.resume_deg_seq(unwanted_facet)
                continue
            picked_facet = random.choice(facets_pool)
            picked_facet_id = self.facet2id[picked_facet]
            self.fill_up_deg_seq(picked_facet)
            if self.deg_checks():
                old_identifier = deepcopy(identifier)
                identifier += [picked_facet_id]
                self.remove_inclusive_facets(picked_facet, tuple(old_identifier), tuple(identifier))
                self.write_to_available_facets(identifier)
            else:
                # Backtrack
                self.available_facets[tuple(identifier)][
                    "why-explore-more"] = f"deg_checks not passed if choosing {picked_facet_id}: {picked_facet}"
                self.logbook[tuple(identifier + [picked_facet_id])] = {"is_simplicial": False}  # marking impossible facet combinations
                self.resume_deg_seq(picked_facet)
                self.sorted_s.insert(0, len(self.id2facet[picked_facet_id]))
                continue

            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self.sorted_s) == 0:
                if self.available_facets[tuple(identifier)]["sorted_deg_seq"] == self._sorted_d and \
                        self.available_facets[tuple(identifier)]["sorted_size_seq"] == self._sorted_s:
                    self.identifier = identifier
                    return True
                else:
                    # Backtrack
                    self.available_facets[tuple(identifier)][
                        "why-explore-more"] = "depleted sorted_s, but simplicial complex not found."
                    self.logbook[tuple(identifier)] = {"is_simplicial": False}  # marking impossible facet combinations
                    last_facet = self.id2facet[identifier.pop(-1)]
                    self.sorted_s = [len(last_facet)] + self.sorted_s
                    self.resume_deg_seq(last_facet)
        return False
