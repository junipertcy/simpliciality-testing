import numpy as np
from collections import defaultdict
from itertools import combinations
import random
from copy import deepcopy
from math import comb
import time
from operator import itemgetter


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


def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]


def gen_joint_sequence(m, poisson_lambda, n_max=0):
    _size_list = [0]
    while min(_size_list) == 0:
        _size_list = sorted(np.random.poisson(lam=poisson_lambda, size=m), reverse=True)

    if n_max == 0:
        # print(n := np.random.randint(max(_size_list) + 1, np.sum(_size_list) + 1))
        n_max = np.sum(_size_list)
    if n_max <= max(_size_list):
        n_max = max(_size_list) + 1

    facets = []
    for facet_size in _size_list:
        candidate_facet = random.sample(range(0, n_max), k=facet_size)
        qualified_draw = False

        count = 0
        while not qualified_draw:
            count += 1
            if count > 1e5:
                print("Re-sample joint sequence at a higher n_max!")
                n_max += 1
                return gen_joint_sequence(m, poisson_lambda, n_max)
            qualified_draw = True
            for facet in facets:
                if set(candidate_facet).issubset(facet):
                    candidate_facet = random.sample(range(0, n_max), k=facet_size)
                    qualified_draw = False
                    break
        facets += [candidate_facet]

    _degree_list = sorted(list(count_summary(flatten(facets)).values()), reverse=True)
    return _size_list, _degree_list


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

    def sample_icebreaker(self, size):
        """
        Parameters
        ----------
        size

        Returns
        -------

        """
        candidate_facet = [_ for _ in range(size)]
        pivot = 1
        while np.sum([self._sorted_d[_] for _ in range(self.n) if _ not in set(candidate_facet)]) < self.m - 1:
            # actually... if it were a greater sign, it wouldn't work as well.
            while candidate_facet != [self.n - 1 - size + _ for _ in range(self.n - 1 - size)]:
                last = candidate_facet.pop(-pivot)
                if last < self.n - pivot:
                    last += 1
                candidate_facet.insert(len(candidate_facet) - pivot + 1, last)
                if last == self.n - pivot:
                    pivot += 1
                break
        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    @staticmethod
    def borrow_from_non_shielding(shielding, non_shielding, m):
        """

        Parameters
        ----------
        shielding
        non_shielding
        m

        Returns
        -------

        """
        for _ in range(m):
            shielding += [non_shielding.pop(0)]
        return shielding, non_shielding

    def prioritize(self, identifier):
        s = count_summary(self.id2name[identifier[-1]])
        shielding = []
        non_shielding = []
        for _ in range(self.n):
            if self._sorted_d[_] - self.deg_seq[_] > 0:
                if s[_] == 1:
                    shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                else:
                    non_shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]

        non_shielding = sorted(non_shielding, key=itemgetter(1), reverse=True)
        shielding = sorted(shielding, key=itemgetter(2), reverse=True)

        return shielding, non_shielding

    @staticmethod
    def get_valid_trails(shielding, non_shielding, size):
        n_shielding = len(shielding)
        n_non_shielding = len(non_shielding)
        n_ = n_shielding + n_non_shielding
        valid_trials = iter(_ for _ in combinations(range(n_), size) if not set(_).issubset(set(range(n_shielding))))
        return valid_trials

    def sample_simplex_greedy(self, identifier, size):
        if len(identifier) == 0:
            return self.sample_icebreaker(size)

        larger_selected_simplex_ids = self.get_selected_facet_ids(identifier, size)

        candidate_facet = []
        shielding, non_shielding = self.prioritize(identifier)
        options = shielding + non_shielding
        valid_trials = self.get_valid_trails(shielding, non_shielding, size)

        non_stop = True
        while non_stop:
            non_stop = False
            try:
                candidate_facet = [options[_][0] for _ in next(valid_trials)]
            except StopIteration:
                raise NotImplementedError("May not be solvable with the greedy algorithm.")
            for _id in larger_selected_simplex_ids:
                _ = self.m - len(identifier) - 1
                cond_1 = set(candidate_facet).issubset(set(self.id2name[_id]))
                cond_2 = np.any(self.get_remaining_slots(identifier, candidate_facet) > _)
                cond_3 = np.sum(self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)) < _
                if cond_1 or cond_2 or cond_3:
                    non_stop = True
                    break

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    def get_remaining_slots(self, identifier, facet, only_non_shielding=False):
        """
        Used in the greedy case only.

        Parameters
        ----------
        identifier: current identifier (not including the candidate facet)
        facet: candidate facet
        only_non_shielding: 

        Returns
        -------

        """
        remaining = self._sorted_d - self.compute_joint_seq_from_identifier(identifier, sorted_deg=False)[1]
        if only_non_shielding:
            return np.array([0 if _ in set(facet) else remaining[_] for _ in range(self.n)])
        else:
            return np.array([remaining[_] - 1 if _ in set(facet) else remaining[_] for _ in range(self.n)])

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
        if len(self._sorted_d) == len(self._sorted_s) == 0:
            return True
        if np.max(self.degree_list) > self.m:
            print("1. This can never be simplicial.")  # TODO.... why??
            return False
        if np.max(self.size_list) >= self.n:
            if len(self.size_list) == 1 and np.max(self.size_list) == self.n:
                return True
            else:
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
