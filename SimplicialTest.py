from itertools import combinations
import numpy as np
from collections import defaultdict
import random
from copy import deepcopy
from math import comb
import time
import sys


class SimplexRegistrar(object):
    def __init__(self):
        self.pointer = 0
        self.name2id = dict()
        self.id2name = dict()
        self.facet_size_per_id = np.array([], dtype=np.int_)

        self.logbook = dict()

        pass

    def register(self, name):
        name = tuple(sorted(name, reverse=True))
        if name not in self.name2id:
            self.name2id[name] = self.pointer
            self.id2name[self.pointer] = name
            self.pointer += 1
            self.facet_size_per_id = np.append(self.facet_size_per_id, [len(name)])
        return name, self.name2id[name]

    def log_forbidden(self, name, reason_id):
        self.logbook[tuple(name)] = {
            "is_simplicial": False,
            "reason": reason_id
        }


def prune_ones(size_list, degree_list):
    size_list = list(size_list)
    degree_list = list(degree_list)
    _1 = count_summary(size_list)[1]
    _2 = count_summary(degree_list)[1]
    if _1 > _2:
        print("Too many 0-simplices. We cannot satisfy the inclusive constraint.")
    else:
        for _ in range(_1):
            size_list.remove(1)
            degree_list.remove(1)
    return size_list, degree_list


def count_summary(seq):
    d = defaultdict(int)
    for s in seq:
        d[s] += 1
    return d


class SimplicialTest(SimplexRegistrar):
    """Base class for SimplicialCheck.

    Parameters
    ----------
    degree_list : ``iterable`` or :class:`numpy.ndarray`, required

    size_list : ``iterable`` or :class:`numpy.ndarray`, required

    """

    def __init__(self, degree_list, size_list):
        super().__init__()

        self.random_seed = time.thread_time_ns()
        random.seed(self.random_seed)

        self.size_list = size_list
        self.degree_list = np.array(degree_list, dtype=np.int_)
        self.logbook = defaultdict(dict)

        self.m = len(size_list)
        self.n = len(degree_list)

        self.sorted_s = sorted(size_list, reverse=True)
        self._sorted_s = np.array(deepcopy(self.sorted_s), dtype=np.int_)
        self._sorted_d = np.array(sorted(degree_list, reverse=True), dtype=np.int_)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.symmetry_breaker = None
        self.identifier = None

        self._backtrack_steps = 0
        self._counter = 0
        self._len_logbook = 0

    def update_deg_seq(self, facet, value):
        if value not in [+1, -1]:
            raise NotImplementedError
        for _ in facet:
            self.deg_seq[_] += value

    def checkpoint_1(self):
        return np.all(np.sort(self.degree_list) - np.sort(self.deg_seq) >= 0)

    def checkpoint_2(self, identifier):
        """Obsolete"""
        return np.any(self.get_available_slots(identifier) > self.m - len(identifier))

    def _break_symmetry(self, greedy=False):
        m = self.sorted_s.pop(0)
        if greedy:
            picked_facet, picked_facet_id = self.sample_simplex_greedy([], m)
        else:
            picked_facet, picked_facet_id = self.register(random.sample(range(self.n), k=m))
        self.update_deg_seq(picked_facet, +1)
        identifier = [picked_facet_id]
        self.symmetry_breaker = picked_facet_id
        return identifier

    def ensure_valid_draw(self, identifier, size):
        """
        This function ensures that our choice of facet is valid (not inclusive of any larger exisiting facet) and
        is potentially a good one (not explored).

        Parameters
        ----------
        identifier
        size

        Returns
        -------

        """

        pass

    def _backward_greedy(self, identifier, size):
        """Obsolete"""
        deg = self.compute_joint_seq_from_identifier(identifier, sorted_deg=False)[1]
        larger_selected_simplex_ids = self.get_selected_facet_ids(identifier, size)
        candidate_facet = []
        shift = self.n - 1
        while len(candidate_facet) < size:
            if shift < 0:
                raise NotImplementedError("Not solvable with greedy::backward strategy.")

            if len(candidate_facet) == size - 1:  # the last vertex
                vacancy_per_vertex = self._sorted_d - self.deg_seq
                sorted_vpv = vacancy_per_vertex[:(shift + 1)].argsort()
                shift_ = 1
                for _id in larger_selected_simplex_ids:
                    while set(candidate_facet + [sorted_vpv[-shift_]]).issubset(set(self.id2name[_id])):
                        shift_ += 1
                        if shift_ > len(sorted_vpv):
                            raise NotImplementedError("Not solvable with greedy::backward strategy.")
                        continue
                if len(identifier) == 0:
                    while np.sum([self._sorted_d[_] for _ in range(self.n) if
                                  _ not in set(candidate_facet + [sorted_vpv[-shift_]])]) < self.m - 1:
                        shift_ += 1
                        continue
                shift = sorted_vpv[-shift_]
            if deg[shift] + 1 <= self._sorted_d[shift]:
                candidate_facet += [shift]
            shift -= 1
        return candidate_facet

    def sample_simplex_greedy_special(self, identifier, size):
        """
        TODO: In case any other strategy needs to be used.
        Parameters
        ----------
        identifier
        size

        Returns
        -------

        """

        return

    def sample_simplex_greedy(self, identifier, size):

        deg = self.compute_joint_seq_from_identifier(identifier, sorted_deg=False)[1]
        larger_selected_simplex_ids = self.get_selected_facet_ids(identifier, size)

        if len(identifier) == 0:
            candidate_facet = [_ for _ in range(size)]
            pivot = 1
            # print(f"initial candidate_facet: {candidate_facet}")
            # print(f"initial pivot: {pivot}")
            while np.sum([self._sorted_d[_] for _ in range(self.n) if _ not in set(candidate_facet)]) < self.m - 1:
                while candidate_facet != [self.n - 1 - size + _ for _ in range(self.n - 1 - size)]:
                    last = candidate_facet.pop(-pivot)
                    if last < self.n - pivot:
                        last += 1
                    candidate_facet.insert(len(candidate_facet) - pivot + 1, last)
                    # print(f"candidate_facet now becomes: {candidate_facet}")
                    if last == self.n - pivot:
                        pivot += 1
                        # print(f"pivot now becomes: {pivot}")
                    break
            picked_facet, picked_facet_id = self.register(candidate_facet)
            return picked_facet, picked_facet_id

        candidate_facet = []
        shift = 0
        while len(candidate_facet) < size:
            if shift >= self.n:
                raise NotImplementedError("Not solvable with the greedy strategy.")
            if len(candidate_facet) == size - 1:  # the last vertex
                # This part may look overly complicated, but the main goal is,
                # rather than finding the vertex_id to be added to the `candidate_facet` while avoiding inclusion,
                # we want to pick the vertex_id that has a higher vacancy first.
                vacancy_per_vertex = self._sorted_d - self.deg_seq
                sorted_vpv = vacancy_per_vertex[shift:].argsort()  # previous indices do not count
                shift_ = 1

                for _id in larger_selected_simplex_ids:
                    while set(candidate_facet + [shift + sorted_vpv[-shift_]]).issubset(set(self.id2name[_id])):
                        shift_ += 1
                        if shift_ > len(sorted_vpv):
                            raise NotImplementedError("Not solvable with the greedy strategy.")
                        continue
                shift += sorted_vpv[-shift_]
            if deg[shift] + 1 <= self._sorted_d[shift]:
                candidate_facet += [shift]
            shift += 1

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    def get_available_slots(self, identifier):
        """
        Only used in greedy setting.
        Parameters
        ----------
        identifier

        Returns
        -------

        """
        degree = self.compute_joint_seq_from_identifier(identifier, sorted_deg=False)[1]
        return self._sorted_d - degree

    def get_selected_facet_ids(self, identifier, size):
        return [index for index, i in enumerate(self.facet_size_per_id) if (index in identifier and i >= size)]

    def sample_simplex(self, identifier, size, greedy=False):
        """

        Parameters
        ----------
        identifier
        size
        greedy

        Returns
        -------

        """
        if greedy:
            return self.sample_simplex_greedy(identifier, size)

        # Here, we have a good criterion! We may not need to explore further if...
        deg_sequence = np.array(self.compute_joint_seq_from_identifier(identifier)[1])
        deg_sequence_goal = self._sorted_d
        if np.any(deg_sequence_goal - deg_sequence < 0):
            self.log_forbidden(identifier, "NO need to explore further - 1")
            self._backtrack_steps = 1
            return list(), -1
        if len(np.nonzero(deg_sequence_goal - deg_sequence)[0]) < size:
            self.log_forbidden(identifier, "NO need to explore further - 2")
            self._backtrack_steps = 1
            return list(), -1

        larger_selected_simplex_ids = self.get_selected_facet_ids(identifier, size)

        set_of_vertices = set(range(self.n))
        picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
        qualified_draw = False
        _ = 0  # I think this is stupid.... but let's try this for now...
        while not qualified_draw:
            qualified_draw = True
            for _id in larger_selected_simplex_ids:

                if set(picked_facet).issubset(set(self.id2name[_id])):
                    self.log_forbidden(identifier + [picked_facet_id], 1)
                    picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
                    qualified_draw = False
                    _ += 1
                    break
            explored = [key for key in self.logbook.keys() if
                        key[:len(identifier)] == tuple(identifier) and len(key) == len(identifier) + 1]

            if len(explored) == comb(self.n, size):
                if len(identifier) == 2:
                    self._backtrack_steps = 1
                else:
                    self._backtrack_steps = 2
                return list(), -1
            if _ > 100:  # TODO
                self._backtrack_steps = 1
                return list(), -1

        return picked_facet, picked_facet_id

    def is_simplicial(self, greedy=False):
        if max(self.degree_list) > self.m:
            print("1. This can never be simplicial.")  # TODO.... why??
            return False
        if max(self.size_list) >= self.n:
            print("2. This can not be simplicial.")
            return False
        if np.sum(self.degree_list) != np.sum(self.size_list):
            print("Failing the Gale–Ryser criterion (1957), the sequence is not bigraphic.")
            return False
        # TODO: there is a second part of the GR criterion, which is not coded yet.

        identifier = self._break_symmetry(greedy=greedy)
        if len(self.sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d:  # TODO: start from identifier
                self.identifier = identifier
                return True
            else:
                return False

        while self._counter < 1e4:
            if len(self.logbook) == self._len_logbook:
                self._counter += 1
            s = self.sorted_s.pop(0)
            picked_facet, picked_facet_id = self.sample_simplex(identifier, s, greedy=greedy)
            if len(picked_facet) == 0:
                self.sorted_s = [s] + self.sorted_s
                self._pull_the_plug(identifier, self._backtrack_steps)
                continue

            self.update_deg_seq(picked_facet, +1)

            if self.checkpoint_1():
                identifier += [picked_facet_id]
            else:
                # Backtrack
                self.log_forbidden(tuple(list(identifier) + [picked_facet_id]), 2)
                self.update_deg_seq(picked_facet, -1)
                self.sorted_s.insert(0, len(self.id2name[picked_facet_id]))
                continue

            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self.sorted_s) == 0:
                if sorted(self.deg_seq, reverse=True) == self._sorted_d.tolist():  # TODO: start from identifier
                    self.identifier = identifier
                    return True
                else:
                    # Backtrack
                    self.log_forbidden(identifier, 3)
                    self._pull_the_plug(identifier, 1)
            if len(self.logbook) > self._len_logbook:
                self._counter = 0
                self._len_logbook = len(self.logbook)
        return False

    def _pull_the_plug(self, identifier, times=1):
        if not len(identifier) > times:
            raise ValueError("You cannot backtrack more times than the length of the identifier.")

        for _ in range(times):
            unwanted_facet = self.id2name[identifier.pop(-1)]
            self.sorted_s = [len(unwanted_facet)] + self.sorted_s
            self.update_deg_seq(unwanted_facet, -1)

    def count_explored_branches(self, identifier):
        s = set()
        for _ in self.logbook.keys():
            try:
                s.add(_[len(identifier) + 1])
            except IndexError:
                pass
        return len(s)

    def compute_joint_seq(self):
        degs = np.zeros(self.n, dtype=np.int_)
        sizes = []

        for _id in self.identifier:
            sizes += [len(self.id2name[_id])]
            for vertex_id in self.id2name[_id]:
                degs[vertex_id] += 1
        return sorted(sizes, reverse=True), sorted(degs, reverse=True)

    def compute_joint_seq_from_identifier(self, identifier, sorted_deg=True):
        degs = np.zeros(self.n, dtype=np.int_)
        sizes = []

        for _id in identifier:
            sizes += [len(self.id2name[_id])]
            for vertex_id in self.id2name[_id]:
                degs[vertex_id] += 1

        if sorted_deg:
            return sorted(sizes, reverse=True), sorted(degs, reverse=True)
        else:
            return sorted(sizes, reverse=True), degs

    def identifier2facets(self, identifier):
        facets = []
        for _id in identifier:
            facets += [self.id2name[_id]]
        return facets
