import numpy as np
from collections import defaultdict, Counter
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


def trim_ones(size_list, degree_list):
    size_list = list(size_list)
    degree_list = list(degree_list)
    _1 = Counter(size_list)[1]
    _2 = Counter(degree_list)[1]
    if _1 > _2:
        print("Too many 0-simplices. We cannot satisfy the inclusive constraint.")
    else:
        for _ in range(_1):
            size_list.remove(1)
            degree_list.remove(1)
    return size_list, degree_list


def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]


def get_slimmest_d(s):
    """

    Parameters
    ----------
    s: size list

    Returns
    -------

    """
    s = sorted(s)
    pool = set()
    tentative_ = []
    for _s in s:
        tentative = tentative_
        if len(tentative) == 0:
            idx = 0
            for _ in range(_s):
                tentative += [idx]
                idx += 1
            pool.add(tuple(tentative))
            tentative_ = tentative
            continue
        tentative[-1] += 1
        idx = tentative[-1]
        for _ in range(_s - len(tentative)):
            idx += 1
            tentative += [idx]

        pool.add(tuple(tentative))
        tentative_ = tentative
    return sorted(Counter(flatten(pool)).values(), reverse=True)


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
                # print("Re-sample joint sequence at a higher n_max!")
                n_max += 1
                return gen_joint_sequence(m, poisson_lambda, n_max)
            qualified_draw = True
            for facet in facets:
                if set(candidate_facet).issubset(facet):
                    candidate_facet = random.sample(range(0, n_max), k=facet_size)
                    qualified_draw = False
                    break
        facets += [candidate_facet]

    _degree_list = sorted(list(Counter(flatten(facets)).values()), reverse=True)
    return _size_list, _degree_list


def sort_helper(st):
    d = defaultdict()
    for idx, _ in enumerate(st._sorted_d):
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

        self.logbook = defaultdict(dict)

        self.m = len(size_list)
        self.n = len(degree_list)

        self.DEGREE_LIST = sorted(degree_list, reverse=True)
        self.SIZE_LIST = sorted(size_list, reverse=True)

        self._sorted_s = np.array(deepcopy(self.SIZE_LIST), dtype=np.int_)
        self._sorted_d = np.array(deepcopy(self.DEGREE_LIST), dtype=np.int_)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.identifier = list()

        self._backtrack_steps = 0
        self._counter = 0
        self._len_logbook = 0

        self.matching = {
            "appending": 0,  # append n vertices to each facet at the end
            "isolated": 0,  # append n 0-simplices at the end
            "heuristic": 0
        }

    def update_deg_seq(self, facet, value):
        if value not in [+1, -1]:
            raise NotImplementedError
        for _ in facet:
            self.deg_seq[_] += value

    def checkpoint_1(self):
        return np.all(np.sort(self._sorted_d) - np.sort(self.deg_seq) >= 0)

    def _break_symmetry(self, greedy=False):
        m = self._sorted_s[0]
        self._sorted_s = np.delete(self._sorted_s, 0)
        if greedy:
            picked_facet, picked_facet_id = self.sample_simplex_greedy(m)
        else:
            picked_facet, picked_facet_id = self.register(random.sample(range(self.n), k=m))
        self.update_deg_seq(picked_facet, +1)
        self.identifier = [picked_facet_id]

    def get_icebreaker_combinations(self, size):
        to_try = defaultdict(list)
        n = len(self._sorted_d)
        _size = n - size
        for s in sorted([(combs, self._sorted_d[list(combs)], np.sum(self._sorted_d[list(combs)])) for combs in
                         combinations(range(n), _size)], key=lambda _1: _1[2], reverse=True):
            if len(to_try[tuple(s[1])]) == 0:
                to_try[tuple(s[1])] = s[0]
            else:
                continue
        return to_try

    def sample_icebreaker(self, size):
        """
        Parameters
        ----------
        size
        q: TODO: when q is minus a large number, reject it.

        Returns
        -------

        """
        to_try = self.get_icebreaker_combinations(size)
        __ = to_try.popitem()[1]
        candidate_facet = [_ for _ in range(self.n) if _ not in __]
        self.validate([], candidate_facet)

        while not self.validate([], candidate_facet):
            try:
                __ = to_try.popitem()[1]
            except KeyError:
                return [], None

            candidate_facet = [_ for _ in range(self.n) if _ not in __]
            self.validate([], candidate_facet)
        self.validate([], candidate_facet)

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    def prioritize(self, identifier):
        s = Counter(self.id2name[identifier[-1]])
        shielding = []
        non_shielding = []
        for _ in range(self.n):
            if self._sorted_d[_] - self.deg_seq[_] > 0:
                if s[_] == 1:
                    shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                else:
                    non_shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
        non_shielding = sorted(non_shielding, key=itemgetter(2), reverse=False)  # originally: 1, true.... exp: 2, false

        shielding = sorted(shielding, key=itemgetter(2), reverse=True)  # originally: 2, true

        return shielding, non_shielding

    @staticmethod
    def get_valid_trails(shielding, non_shielding, size):
        n_shielding = len(shielding)
        n_non_shielding = len(non_shielding)
        n_ = n_shielding + n_non_shielding
        valid_trials = []
        i = 0
        while size + i <= n_:
            valid_trials += list(
                _ for _ in combinations(range(size + i), size) if not set(_).issubset(set(range(n_shielding))))
            i += 1
        valid_trials = iter(valid_trials)
        return valid_trials

    def sample_simplex_greedy(self, size):
        identifier = self.identifier
        if len(identifier) == 0:
            return self.sample_icebreaker(size)

        larger_selected_simplex_ids = self.get_selected_facet_ids(size)

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
                if not self.validate(identifier, candidate_facet, _id=_id):
                    non_stop = True
                    break

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    @staticmethod
    def validate_cond_0(d, s):
        d = deepcopy(d)
        s = deepcopy(s)

        d = np.array(d, dtype=np.int_)
        s = np.array(s, dtype=np.int_)
        while Counter(d)[len(s)] != 0:
            s = np.sort(s)
            d = np.sort(d)
            if Counter(d)[len(s)] < np.min(s):
                s -= Counter(d)[len(s)]
                s = s[s != 0]
            elif len(s) > 1 and Counter(d)[len(s)] == np.min(s):
                return True

            d = d[d != len(s)]

            _1 = Counter(s)[1]
            _2 = Counter(d)[1]
            if _1 <= _2:
                for _ in range(_1):
                    s = s[1:]
                    d = d[1:]

            if len(s) > 1 and len(d) == np.max(s):
                return True
        return False

    def validate_issubset(self, candidate_facet, _id):
        """
        Verify that the selected facet not be a subset of any higher facet.

        Parameters
        ----------
        candidate_facet
        _id

        Returns
        -------

        """
        if _id is not None:
            return set(candidate_facet).issubset(set(self.id2name[_id]))
        return False

    def validate_nonshielding(self, remaining, both, non_shielding, shielding, candidate_facet=None):
        """
        Verify that the sum of non-shielding slots not be less than the number of the remaining facets.

        Parameters
        ----------
        remaining
        non_shielding
        shielding

        Returns
        -------

        """
        if np.sum(non_shielding) < remaining:
            return True
        elif np.sum(non_shielding) == remaining:
            _ = self._sorted_s - 1
            if len(_) - 1 <= 0:  # only the last to-be-chosen facet remains
                return False
            if np.count_nonzero(shielding) == _[0]:  # There must be at least 2 facets that remain to be chosen.
                if Counter(non_shielding)[1] == 0:
                    return True
            if len(self.identifier2facets()) == 0:

                if Counter(_)[1] > Counter(both)[1]:  # per unknown_cases[52]
                    return True
            # if Counter(_)[1] > Counter(both[never_filled_filter])[1]:  # per unknown_cases[52]
            #     return True
        else:
            return False  # safe!
        return False

    def validate_9(self, candidate_facet, remaining, both, identifier):
        """

        Parameters
        ----------
        remaining
        both
        identifier

        Returns
        -------

        """
        if len(np.where(both == remaining)[0]) > 0:
            must_fill = len(np.where(both == remaining)[0])
            both_ = deepcopy(both)
            both_[both_ == remaining] = 0
            remaining_ = self._sorted_s - must_fill
            if Counter(remaining_)[1] > Counter(both_)[1]:
                return True
            else:
                if len(remaining_[remaining_ != 1]) < np.max(both_[both_ != 1]):
                    return True

            for vertex_id in np.where(both == remaining)[0]:  # vertex_id's when we had to fill up all remaining slots
                for facet in self.identifier2facets(identifier) + [candidate_facet]:  # find existing facets that contain this vertex_id
                    if vertex_id in facet:
                        non_shielding_part = set(range(self.n)).difference(set(facet))  # non_shielding_part vertex_id's
                        if remaining > np.sum(both[np.array(list(non_shielding_part))]):  # remaining number of facets must be able to "hide" in those non_shielding slots
                            return True
                        if Counter(remaining_)[1] > Counter(both_)[1]:
                            return True
        return False

    def validate(self, identifier, candidate_facet, _id=None):
        """
        This function must return True in order for the candidate facet to be considered.
        In other words, if any condition appears to be True, we will not accept the candidate facet.

        cond_3: when the next-to-be-selected facet has the number of empty slots equal to that of the remaining slots,
                the sum of the remaining slots cannot be larger than that number.

        cond_4: when it's not the very last round of selection,
                the number of the remaining slots cannot be less than that of the next-to-be-selected facet.
                (otherwise, we simply cannot find such a facet)

        cond_5: the number of remaining slots, for each vertex, cannot be larger than the number of unselected facets.

        Parameters
        ----------
        identifier
        candidate_facet
        _id

        Returns
        -------

        """
        _ = self.m - len(identifier) - 1
        _1 = self.get_remaining_slots(identifier, candidate_facet)  # all (non_shielding & shielding)
        _2 = self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)
        _3 = _1 - _2  # only shielding
        if self.validate_9(candidate_facet, _, _1, identifier):
            return False
        if self.validate_cond_0(_1[_1 != 0], self._sorted_s):
            return False
        if self.validate_issubset(candidate_facet, _id):
            return False
        if self.validate_nonshielding(_, _1, _2, _3, candidate_facet=candidate_facet):
            return False
        if np.any(_1 > _):
            return False
        if np.sum(_1) > np.count_nonzero(_1) == self._sorted_s[0]:
            return False
        if len(self._sorted_s) > 0 and np.count_nonzero(_1) < self._sorted_s[0]:
            return False
        return True

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

    def get_selected_facet_ids(self, size):
        identifier = self.identifier
        return [index for index, i in enumerate(self.facet_size_per_id) if (index in identifier and i >= size)]

    def sample_simplex(self, size, greedy=False):
        """

        Parameters
        ----------
        size
        greedy

        Returns
        -------

        """
        identifier = self.identifier

        if greedy:
            return self.sample_simplex_greedy(size)

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

        larger_selected_simplex_ids = self.get_selected_facet_ids(size)

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

    def preprocess(self):
        idx = 0
        while self._sorted_d[idx] == len(self._sorted_s):
            self._sorted_s -= 1
            self.matching["appending"] += 1
            self._sorted_d = np.delete(self._sorted_d, 0)
            self.n -= 1
            idx += 1
            # TODO: if any 0 is present in the sizes list, it's a different story...

        _s = Counter(self._sorted_s)[1]
        _d = Counter(self._sorted_d)[1]
        if _s > _d:
            print("Too many 0-simplices. We cannot satisfy the inclusive constraint.")
            return False
        else:
            for _ in range(_s):
                self._sorted_d = np.delete(self._sorted_d, -1)
                self._sorted_s = np.delete(self._sorted_s, -1)
                self.m -= 1
                self.n -= 1
                self.matching["isolated"] += 1

        _s = Counter(self._sorted_s)[1]
        _d = Counter(self._sorted_d)[1]
        if _s == 0 and _d > 0:
            for _ in range(_d):
                copy_sorted_s = deepcopy(self._sorted_s)
                copy_sorted_d = deepcopy(self._sorted_d)
                self._sorted_d = np.delete(self._sorted_d, -1)
                self._sorted_d[0] -= 1
                self._sorted_d = np.array(sorted(self._sorted_d, reverse=True))
                self._sorted_s = np.delete(self._sorted_s, -1)

                self.m -= 1
                self.n -= 1
                if self._validate_data() is False:
                    self._sorted_s = copy_sorted_s
                    self._sorted_d = copy_sorted_d
                    self.m += 1
                    self.n += 1
                    break
                else:
                    m = self._sorted_s[0]
                    self._sorted_s = np.delete(self._sorted_s, 0)
                    if self.sample_icebreaker(m)[1] is None:
                        self._sorted_s = copy_sorted_s
                        self._sorted_d = copy_sorted_d
                        self.m += 1
                        self.n += 1
                        break
                    self._sorted_s = np.insert(self._sorted_s, 0, m)
                    self.matching["heuristic"] += 1

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        return True

    def postprocess(self):
        pass

    def is_simplicial(self, greedy=False, preprocess=False):
        if self._validate_data():
            return True
        elif self._validate_data() is False:
            return False

        if preprocess:
            _ = self.preprocess()
            if _ is False:
                return False

        self._break_symmetry(greedy=greedy)
        if len(self._sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d:  # TODO: start from identifier
                return True
            else:
                return False

        while self._counter < 1e4:
            if len(self.logbook) == self._len_logbook:
                self._counter += 1
            s = self._sorted_s[0]
            self._sorted_s = np.delete(self._sorted_s, 0)
            picked_facet, picked_facet_id = self.sample_simplex(s, greedy=greedy)
            if len(picked_facet) == 0:
                self._sorted_s = np.insert(self._sorted_s, 0, s)
                self._pull_the_plug(self._backtrack_steps)
                continue

            self.update_deg_seq(picked_facet, +1)

            if self.checkpoint_1():
                self.identifier += [picked_facet_id]
            else:
                # Backtrack
                self.log_forbidden(tuple(list(self.identifier) + [picked_facet_id]), 2)
                self.update_deg_seq(picked_facet, -1)
                self._sorted_s = np.insert(self._sorted_s, 0, len(self.id2name[picked_facet_id]))
                continue

            # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
            if len(self._sorted_s) == 0:
                if sorted(self.deg_seq, reverse=True) == self._sorted_d.tolist():  # TODO: start from identifier
                    return True
                else:
                    # Backtrack
                    self.log_forbidden(self.identifier, 3)
                    self._pull_the_plug(1)
            if len(self.logbook) > self._len_logbook:
                self._counter = 0
                self._len_logbook = len(self.logbook)
        return False

    def _pull_the_plug(self, times=1):
        if not len(self.identifier) > times:
            raise ValueError("You cannot backtrack more times than the length of the identifier.")

        for _ in range(times):
            unwanted_facet = self.id2name[self.identifier.pop(-1)]
            self._sorted_s = np.insert(self._sorted_s, 0, len(unwanted_facet))
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

    def compute_joint_seq_from_identifier(self, identifier=None, sorted_deg=True):
        if identifier is None:
            identifier = self.identifier
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

    def identifier2facets(self, identifier=None):
        if identifier is None:
            identifier = self.identifier
        facets = []
        for _id in identifier:
            facets += [self.id2name[_id]]
        return facets

    def _validate_data(self):
        if len(self._sorted_d) == len(self._sorted_s) == 0:
            return True
        if np.max(self._sorted_d) > self.m:
            # print("1. This can never be simplicial.")  # TODO.... why??
            return False

        if np.max(self._sorted_s) >= self.n:
            if len(self._sorted_s) == 1 and np.max(self._sorted_s) == self.n:
                return True
            else:
                # print("2. This can not be simplicial.")
                return False
        if np.sum(self._sorted_d) != np.sum(self._sorted_s):
            # print("Failing the Gale–Ryser criterion (1957), the sequence is not bigraphic.")
            return False
        # TODO: there is a second part of the GR criterion, which is not coded yet.
