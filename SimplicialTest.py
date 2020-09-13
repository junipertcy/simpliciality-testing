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


def get_nshielding(both, curent_sizes, current_facets, candidate_facet):
    remaining = len(curent_sizes)
    n = len(both)
    _nshielding = set()
    flag = True
    for vertex_id in np.where(both == remaining)[0]:  # vertex_id's when we had to fill up all remaining slots
        for facet in current_facets + [candidate_facet]:  # find existing facets that contain this vertex_id
            if vertex_id in facet:
                non_shielding_part = set(range(n)).difference(set(facet))  # non_shielding_part vertex_id's
                if flag:
                    _nshielding = non_shielding_part
                    flag = False
                else:
                    _nshielding.intersection_update(non_shielding_part)
                # remaining number of facets must be able to "hide" in those non_shielding slots, at least
                if remaining > np.sum(both[np.array(list(non_shielding_part), dtype=np.int_)]):
                    return True, (flag, set())
    return False, (flag, _nshielding)

def accel_asc(n, size, counter):
    """
    Modified from: http://jeromekelleher.net/tag/integer-partitions.html

    Parameters
    ----------
    n
    size
    counter

    Returns
    -------

    """
    a = [0 for _ in range(n + 1)]
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
            if len(a[:k + 2]) == size:
                _counter = Counter(a[:k + 2])
                indicator = True
                for _ in _counter:
                    if _counter[_] > counter[_]:
                        indicator = False
                        break
                if indicator:
                    yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        if len(a[:k + 1]) == size:
            _counter = Counter(a[:k + 1])

            indicator = True
            for _ in _counter:
                if _counter[_] > counter[_]:
                    indicator = False
                    break
            if indicator:
                yield a[:k + 1]


def remove_ones(s, both):
    both = [_ for _ in both]
    for _ in range(Counter(s)[1]):
        both.remove(1)
    both = np.array(both)
    return both


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
            "isolated": 0,  # append n 0-simplices at the end
            "heuristic": 0
        }

        # print(f"SimplicialTest initiated. degree_list={self.DEGREE_LIST}; size_list={self.SIZE_LIST}.")

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

    def get_nonselected(self, ind):
        nonselected = set()
        pool = defaultdict(list)
        for _ in ind:
            if len(pool[_]) == 0:
                pool[_] = list(np.where(self._sorted_d == _)[0])
            nonselected.add(pool[_].pop(-1))
        return nonselected

    def get_icebreaker_combinations(self, size):
        n = len(self._sorted_d)
        _size = n - size
        _min = np.sum(self._sorted_d[-_size:])
        _max = np.sum(self._sorted_d[:_size])
        k = _min
        gen = accel_asc(k, _size, Counter(self._sorted_d))
        while k <= _max:
            try:
                yield next(gen)
            except StopIteration:
                k += 1
                gen = accel_asc(k, _size, Counter(self._sorted_d))

    def sample_icebreaker(self, size, method="heuristic"):
        """
        Parameters
        ----------
        size
        method: `heuristic` or `strict`

        Returns
        -------

        """
        n = len(self._sorted_d)
        if method == "heuristic":
            gen = combinations(range(n), n - size)
            __ = [n - 1 - _ for _ in next(gen)]
        elif method == "strict":
            gen = self.get_icebreaker_combinations(size)
            __ = self.get_nonselected(next(gen))
        else:
            raise AttributeError("method can only be `heuristic` or `strict`.")

        candidate_facet = tuple([_ for _ in range(self.n) if _ not in __])
        while not self.validate([], candidate_facet):
            try:
                if method == "heuristic":
                    __ = [n - 1 - _ for _ in next(gen)]
                elif method == "strict":
                    __ = self.get_nonselected(next(gen))
            except StopIteration:
                return [], None

            candidate_facet = [_ for _ in range(self.n) if _ not in __]
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
    def get_valid_trials(shielding, non_shielding, size):
        n_shielding = len(shielding)
        n_non_shielding = len(non_shielding)
        n_ = n_shielding + n_non_shielding
        i = 0

        gen = combinations(range(size + i), size)
        while size + i <= n_:
            try:
                _comb = next(gen)
            except StopIteration:
                i += 1
                gen = combinations(range(size + i), size)
            else:
                if not set(_comb).issubset(set(range(n_shielding))):
                    yield _comb

    def sample_simplex_greedy(self, size):
        identifier = self.identifier
        if len(identifier) == 0:
            return self.sample_icebreaker(size)

        larger_selected_simplex_ids = self.get_selected_facet_ids(size)

        candidate_facet = []
        shielding, non_shielding = self.prioritize(identifier)
        options = shielding + non_shielding
        valid_trials = self.get_valid_trials(shielding, non_shielding, size)

        non_stop = True
        while non_stop:
            non_stop = False
            try:
                candidate_facet = tuple([options[_][0] for _ in next(valid_trials)])
            except StopIteration:
                raise NotImplementedError("May not be solvable with the greedy algorithm.")
            for _id in larger_selected_simplex_ids:
                if not self.validate(identifier, candidate_facet, _id=_id):
                    non_stop = True
                    break

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

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

    @staticmethod
    def validate_nonshielding(curent_sizes, non_shielding, shielding):
        """
        Verify that the sum of non-shielding slots not be less than the number of the remaining facets.

        Parameters
        ----------
        curent_sizes
        non_shielding
        shielding

        Returns
        -------

        """
        remaining = len(curent_sizes)
        if np.sum(non_shielding) < remaining:
            return True
        elif np.sum(non_shielding) == remaining:
            _ = curent_sizes - 1
            if len(_) - 1 <= 0:  # only the last to-be-chosen facet remains
                return False
            if np.count_nonzero(shielding) == _[0]:  # There must be at least 2 facets that remain to be chosen.
                if Counter(non_shielding)[1] == 0:
                    return True
        return False  # safe!

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
        if self.validate_issubset(candidate_facet, _id):
            return False

        _both = self.get_remaining_slots(identifier, candidate_facet)  # both (non_shielding & shielding)
        _sizes = self._sorted_s
        if np.sum(_both) > np.count_nonzero(_both) == _sizes[0]:
            return False
        if len(_sizes) > 0 and np.count_nonzero(_both) < _sizes[0]:
            return False
        if np.any(_both > len(_sizes)):
            return False
        _non_shielding = self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)
        _shielding = _both - _non_shielding  # only shielding
        if self.validate_nonshielding(_sizes, _non_shielding, _shielding):
            return False

        _degrees = self._sorted_d
        _current_facets = self.identifier2facets(identifier)
        if len(_sizes) > 0 and Counter(_both)[len(_sizes)] > 0:
            reduced_seq = self.is_reduced_seq(_both, _sizes, _degrees)
            # print(type(_sizes))
            _ = self.validate_reduced_seq(_both, _sizes, _current_facets, candidate_facet)
            if _ is True:
                return False
            elif _ is False:
                pass
            else:
                degs, sizes = _
                bool_ = self._recursive_is_simplicial(degs, sizes)
                if reduced_seq and bool_:
                    # print("///// This is definitely simplicial! /////")
                    pass
                return bool_

        return True

    @staticmethod
    def is_reduced_seq(both, sizes, degrees):
        return Counter(np.equal(both, degrees))[True] == Counter(both)[len(sizes)]

    @staticmethod
    def validate_reduced_seq(both, curent_sizes, current_facets, candidate_facet):
        if np.any(curent_sizes < 0):
            return True

        indicator, (flag, _nshielding) = get_nshielding(both, curent_sizes, current_facets, candidate_facet)
        if indicator:
            return True

        _ns_both = both[np.array(list(_nshielding), dtype=np.int_)]
        curent_sizes = np.array(curent_sizes, dtype=np.int_)

        firsttime = True
        while Counter(both)[len(curent_sizes)] != 0:
            if len(curent_sizes) > 1:
                if Counter(both)[len(curent_sizes)] == np.min(curent_sizes) or len(both) == np.max(curent_sizes):
                    return True

            _ = Counter(both)[len(curent_sizes)]

            curent_sizes -= Counter(both)[len(curent_sizes)]
            curent_sizes = curent_sizes[curent_sizes != 0]
            if len(curent_sizes) == 0:
                return False

            if Counter(curent_sizes)[1] > 0 and Counter(curent_sizes)[1] > Counter(_ns_both)[1]:
                if _ == 1 and flag is False and firsttime is True:
                    return True

            both = both[both != len(curent_sizes)]
            both = both[both != 0]

            if Counter(curent_sizes)[1] > Counter(both)[1]:
                return True
            else:
                both = remove_ones(curent_sizes, both)
                curent_sizes = curent_sizes[curent_sizes != 1]
            if firsttime:
                firsttime = False
        return both, curent_sizes

    @staticmethod
    def _recursive_is_simplicial(degs, sizes):
        try:
            st = SimplicialTest(degs, sizes)
            bool_ = st.is_simplicial(greedy=True, preprocess=False)
        except (KeyError, NotImplementedError):
            return False
        if bool_:
            return True
        else:
            return False

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

    def is_simplicial(self, greedy=False, preprocess=False):
        if self._validate_data():
            return True
        elif self._validate_data() is False:
            return False
        if not self.match_ones():
            return False

        if preprocess:
            _ = self.preprocess()
            if _ is False:
                return False
            if len(self._sorted_s) == len(self._sorted_d) == 0:
                return True

        self._break_symmetry(greedy=greedy)
        if len(self._sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d.tolist():  # TODO: start from identifier
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
        if len(self._sorted_d) > 0 and np.max(self._sorted_d) > self.m:
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

    def match_ones(self):
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
        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        return True

    def preprocess(self):
        _s = Counter(self._sorted_s)[1]
        _d = Counter(self._sorted_d)[1]
        if _s == 0 and _d > 0:
            for _ in range(_d):
                copy_sorted_s = deepcopy(self._sorted_s)
                copy_sorted_d = deepcopy(self._sorted_d)
                self._sorted_d = np.delete(self._sorted_d, -1)
                self._sorted_d[0] -= 1
                self._sorted_d = self._sorted_d[self._sorted_d != 0]
                self._sorted_d = np.array(sorted(self._sorted_d, reverse=True))
                self._sorted_s = np.delete(self._sorted_s, -1)

                self.m -= 1
                self.n -= 1
                if len(self._sorted_s) == 0:
                    break
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
