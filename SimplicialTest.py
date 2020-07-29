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
        self.facet_sizes = np.array([], dtype=np.int_)

        self.logbook = dict()
        pass

    def register(self, name):
        name = tuple(sorted(name, reverse=True))
        if name not in self.name2id:
            self.name2id[name] = self.pointer
            self.id2name[self.pointer] = name
            self.pointer += 1
            self.facet_sizes = np.append(self.facet_sizes, [len(name)])
        return name, self.name2id[name]

    def log_forbidden(self, name):
        self.logbook[tuple(name)] = False


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
        self._sorted_s = deepcopy(self.sorted_s)
        self.sorted_d = sorted(degree_list, reverse=True)
        self._sorted_d = deepcopy(self.sorted_d)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.symmetry_breaker = None
        self.identifier = None

        self._backtrack_steps = 0
        self._global_hit_wall = 0

        pass

    def update_deg_seq(self, facet, value):
        if value not in [+1, -1]:
            raise NotImplementedError
        for _ in facet:
            self.deg_seq[_] += value

    def checkpoint_1(self):
        return np.all(np.sort(self.degree_list) - np.sort(self.deg_seq) >= 0)

    def _break_symmetry(self):
        m = self.sorted_s.pop(0)
        picked_facet, picked_facet_id = self.register(random.sample(range(self.n), k=m))

        self.update_deg_seq(picked_facet, +1)
        identifier = [picked_facet_id]
        self.symmetry_breaker = 0  # used to check if the routine ends
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

    def sample_simplex(self, identifier, size):
        """

        Parameters
        ----------
        identifier
        size

        Returns
        -------

        """
        # TODO: Fix things here!
        explored_but_forbidden_ids = []
        explored_but_forbidden_ids = [k[len(identifier):][0] for k in self.logbook.keys() if
                                      len(k[len(identifier):]) == 1]
        larger_selected_simplex_ids = [index for index, i in enumerate(self.facet_sizes) if
                                       (index in identifier and i >= size) or index in explored_but_forbidden_ids]
        larger_selected_simplex_ids = [index for index, i in enumerate(self.facet_sizes) if (index in identifier and i >= size)]

        set_of_vertices = set(range(self.n))
        picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
        qualified_draw = False
        _ = 0  # I think this is stupid.... but let's try this for now...
        while not qualified_draw:
            qualified_draw = True
            for _id in larger_selected_simplex_ids:
                if set(picked_facet).issubset(set(self.id2name[_id])):
                    self.log_forbidden(identifier + [picked_facet_id])
                    aaa = tuple(list(identifier) + [picked_facet_id])
                    picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
                    qualified_draw = False
                    _ += 1
                    break
            if _ > 10:

                self._global_hit_wall += 1
                if self._global_hit_wall > 10:
                    self._global_hit_wall = 0
                    self._backtrack_steps = 2
                else:
                    self._backtrack_steps = 1
                return list(), -1

        [set_of_vertices.remove(i) for i in picked_facet]
        if len(set_of_vertices) < size:
            # This is an impossible path
            self.log_forbidden(identifier)

            self._backtrack_steps = 1
            return list(), -1

        pool = set()
        pool_size = comb(len(set_of_vertices), size)
        already_in_identifier = 0
        while picked_facet_id in identifier or tuple(identifier + [picked_facet_id]) in self.logbook:
            picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
            pool.add(tuple(identifier + [picked_facet_id]))
            if len(pool) == pool_size - already_in_identifier:
                self._backtrack_steps = 1
                return list(), -1

        # print(f"WE HAVE CHOSEN, IN SAMPLE_SIMPLEX, {picked_facet, picked_facet_id}")
        return picked_facet, picked_facet_id

    def is_simplicial(self):
        if max(self.degree_list) > self.m:
            print("1. This can never be simplicial.")  # TODO.... why??
            return False
        if np.sum(self.degree_list) != np.sum(self.size_list):
            print("2. This can never be simplicial.")
            return False

        identifier = self._break_symmetry()
        if len(self.sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d:  # TODO: start from identifier
                self.identifier = identifier
                return True
            else:
                return False

        while True:
            s = self.sorted_s.pop(0)
            picked_facet, picked_facet_id = self.sample_simplex(identifier, s)
            if len(picked_facet) == 0:
                if identifier == [self.symmetry_breaker]:
                    return False
                self.sorted_s = [s] + self.sorted_s
                unwanted_facet = self.id2name[identifier.pop(-1)]
                self.sorted_s = [len(unwanted_facet)] + self.sorted_s
                self.update_deg_seq(unwanted_facet, -1)
                if self._backtrack_steps == 2:
                    unwanted_facet = self.id2name[identifier.pop(-1)]
                    self.sorted_s = [len(unwanted_facet)] + self.sorted_s
                    self.update_deg_seq(unwanted_facet, -1)
                    unwanted_facet = self.id2name[identifier.pop(-1)]
                    self.sorted_s = [len(unwanted_facet)] + self.sorted_s
                    self.update_deg_seq(unwanted_facet, -1)
                continue

            self.update_deg_seq(picked_facet, +1)

            if self.checkpoint_1():
                identifier += [picked_facet_id]
            else:
                # Backtrack
                self.log_forbidden(tuple(list(identifier) + [picked_facet_id]))
                aaa = tuple(list(identifier) + [picked_facet_id])
                self.update_deg_seq(picked_facet, -1)
                self.sorted_s.insert(0, len(self.id2name[picked_facet_id]))
                continue

            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self.sorted_s) == 0:
                if sorted(self.deg_seq, reverse=True) == self._sorted_d:  # TODO: start from identifier
                    self.identifier = identifier
                    return True
                else:
                    # Backtrack
                    self.log_forbidden(identifier)
                    last_facet = self.id2name[identifier.pop(-1)]

                    self.sorted_s = [len(last_facet)] + self.sorted_s
                    self.update_deg_seq(last_facet, -1)
        return False

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
