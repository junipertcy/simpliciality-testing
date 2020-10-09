import time
import validators

from utils import *
from math import comb
from operator import itemgetter
from itertools import combinations


class SimplicialTest(SimplexRegistrar):
    """Base class for SimplicialCheck.

    Parameters
    ----------
    degree_list : ``iterable`` or :class:`numpy.ndarray`, required

    size_list : ``iterable`` or :class:`numpy.ndarray`, required

    """
    def __init__(self, degree_list, size_list, level=0, blocked_sets=None, verbose=False):
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
        self._level = level + 1

        self.matching = {
            "isolated": 0,  # append n 0-simplices at the end
        }
        if blocked_sets is None:
            self.blocked_sets = list()
        else:
            self.blocked_sets = blocked_sets
        self.verbose = verbose
        self.reduced_data = None
        self.collected_facets = list()
        self.debug_counter = 0

        # print(f"SimplicialTest initiated. degree_list={self.DEGREE_LIST}; size_list={self.SIZE_LIST}.")

    def update_deg_seq(self, facet, value):
        if value not in [+1, -1]:
            raise NotImplementedError
        for _ in facet:
            self.deg_seq[_] += value

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
        candidate_facet = [_ for _ in range(self.n) if _ not in __]
        is_validate, reason = self.validate([], candidate_facet)
        while not is_validate:
            try:
                if method == "heuristic":
                    __ = [n - 1 - _ for _ in next(gen)]
                elif method == "strict":
                    __ = self.get_nonselected(next(gen))
            except StopIteration:
                return [], None

            candidate_facet = [_ for _ in range(self.n) if _ not in __]
            is_validate, reason = self.validate([], candidate_facet)

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    def prioritize(self, identifier):
        s = Counter(self.id2name[identifier[-1]])
        shielding = []
        non_shielding = []
        non_shielding_1 = []
        for _ in range(self.n):
            if self._sorted_d[_] - self.deg_seq[_] > 0:
                if s[_] == 1:
                    shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                else:
                    if self._sorted_d[_] - self.deg_seq[_] != 1:
                        non_shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                    else:
                        non_shielding_1 += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
        non_shielding = sorted(non_shielding, key=itemgetter(1), reverse=True)  # originally: 1, true.... exp: 2, false
        non_shielding = non_shielding_1 + non_shielding
        shielding = sorted(shielding, key=itemgetter(1), reverse=True)  # originally: 2, true
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
                ind, reason = self.validate(identifier, candidate_facet, _id=_id)
                if not ind:
                    non_stop = True
                    # print(f"[l={self._level}] okay, here we learnt that {candidate_facet} is rejected due to {reason}")
                    break

        picked_facet, picked_facet_id = self.register(candidate_facet)
        return picked_facet, picked_facet_id

    def validate(self, identifier, candidate_facet, _id=None) -> (bool, str):
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
        current_facets = self.identifier2facets(identifier) + [candidate_facet]

        if validators.validate_issubset(self.id2name, candidate_facet, _id, blocked_sets=self.blocked_sets):
            return False, "1"
        _sizes = self._sorted_s
        if len(_sizes) == 0:
            return True, "0"

        _both = self.get_remaining_slots(identifier, candidate_facet)  # both (non_shielding & shielding)
        if np.min(_both) < 0:
            return False, "13"

        if len(_sizes) == 1:
            for facet in current_facets:
                if set(np.nonzero(_both)[0]).issubset(set(facet)):
                    return False, "12"

        if np.any(_both > len(_sizes)):
            return False, "4"
        if np.sum(_both) > np.count_nonzero(_both) == _sizes[0]:
            return False, "2"

        _non_shielding = self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)
        if np.min(_non_shielding) < 0:
            return False, "14"
        _shielding = _both - _non_shielding  # only shielding
        # print(_sizes, _non_shielding, _shielding, current_facets)
        if validators.validate_nonshielding(_sizes, _non_shielding, _shielding):
            return False, "5"
        if np.count_nonzero(_both) < _sizes[0]:
            return False, "3"
        for blocked_set in self.blocked_sets:
            if set(np.nonzero(_both)[0]).issubset(set(blocked_set)):
                return False, "11"

        # _degrees = self._sorted_d
        if len(current_facets) > 1:
            if Counter(_both)[len(_sizes)] > 0:
                # reduced_seq = self.is_reduced_seq(_both, _sizes, _degrees)
                _ = self.validate_reduced_seq(_both, _sizes, current_facets)
                if _:
                    return False, "6"
                if self.reduced_data is not None:
                    # print(f"We have self.reduced_data={self.reduced_data}")
                    degs, sizes, shielding_facets = self.reduced_data
                    if np.sum(degs) == np.sum(sizes) == 0:
                        for facet in current_facets:
                            for collected_facet in self.collected_facets:
                                if set(collected_facet).issubset(set(facet)):
                                    # print(f"collected_facet = {collected_facet} is a subset of {facet}")
                                    return False, "10"
                        # print(f"self.collected_facets = {self.collected_facets}")
                        # print(f"accepting a proposal = {candidate_facet}")
                        return True, "0"
                    mapping = get_reduced_seq_mapping(degs)
                    inv_map = {v: k for k, v in mapping.items()}
                    _blocked_sets = transform_facets(shielding_facets, inv_map)
                    # print(f"(l={self._level}) recursive call, input (degs, sizes) = {(degs, sizes)} current_facets={current_facets}; ")
                    bool_, facets = self._recursive_is_simplicial(degs, sizes, blocked_sets=_blocked_sets, verbose=False)
                    # if bool_:
                    #     print(f"Oh, yeah, the recursive version is simplicial. facets = {facets}")  # TODO: use case 38 to complete this part
                    #     print(
                    #     f"also, we have self.collected_facets = {self.collected_facets}")  # TODO: use case 38 to complete this part
                    # if reduced_seq and bool_:
                    #     # print("📝 This is definitely simplicial! 📝")
                    #     pass
                    if not bool_:
                        return bool_, "7"
                    else:
                        return True, "0"
            else:
                self.debug_counter += 1
                # print(f"(l={self._level}) Checking.... current facets: {current_facets}")
                # if sorted(candidate_facet) == [0, 1, 2, 3, 4, 5, 12] or sorted(candidate_facet) == [0, 1, 2, 3, 4, 7, 12]:
                mapping = get_reduced_seq_mapping(_both)
                inv_map = {v: k for k, v in mapping.items()}
                # print(f"... \n candidate_facet={candidate_facet} \n _both= {_both}; _sizes = {_sizes}")

                # print(f"inv_map: {inv_map}")
                _blocked_sets = transform_facets(current_facets, inv_map)
                # print(f"current facets: {current_facets}")
                # print(f"(l={self._level}) Checking.... to recurse: {_both, _sizes, _blocked_sets}")
                bool_, facets = self._recursive_is_simplicial(_both, _sizes, blocked_sets=_blocked_sets, verbose=False)
                # print(f"(l={self._level}) bool_ = {bool_}; facets={facets} \n")
                # if self.debug_counter == 20000:
                #     exit(0)
                if not bool_:
                    return bool_, "333"

        return True, "0"

    def validate_reduced_seq(self, both, curent_sizes, current_facets) -> bool:
        self.reduced_data = None
        self.collected_facets = list()
        # print(f"START validate_reduced_seq:: (both, curent_sizes, current_facets) = {both, curent_sizes, current_facets}")
        curent_sizes = np.array(curent_sizes, dtype=np.int_)
        shielding_facets = []
        exempt_vids = []
        while Counter(both)[len(curent_sizes)] != 0:
            if len(curent_sizes) > 1:
                if not basic_validations_degs_and_sizes(degs=both, sizes=curent_sizes):
                    # print("validate_reduced_seq - 1")
                    return True

            must_be_filled_vids = np.where(both == len(curent_sizes))[0]

            # print(f"START get_nshielding:: (both, curent_sizes, current_facets) = {both, curent_sizes, current_facets}")
            shielding_facets = get_shielding_facets_when_vids_filled(current_facets, must_be_filled_vids, exempt_vids=exempt_vids)
            # print(f"Before (get_nonshielding_vids) shielding_facets = {shielding_facets}; {np.nonzero(both)}")
            nonshielding_vids = get_nonshielding_vids(shielding_facets, both)
            # print(f"After (get_nonshielding_vids) nonshielding_vids = {nonshielding_vids}")
            shielding_facets = [set(_).difference(set(must_be_filled_vids)) for _ in shielding_facets]

            curent_sizes -= Counter(both)[len(curent_sizes)]
            # print(f"both = {both}; must_be_filled_vids={must_be_filled_vids}")

            if Counter(curent_sizes)[0] == len(curent_sizes):
                # print(f"both = {both}; must_be_filled_vids={must_be_filled_vids}")
                self.collected_facets += [exempt_vids + must_be_filled_vids.tolist()]
                # print(f"Adding {exempt_vids + [must_be_filled_vids]} to collected_facets")

                # print(3, curent_sizes, "###### BREAK #########")
                both[both == len(curent_sizes)] = 0
                break
            both[both == len(curent_sizes)] = 0

            if Counter(curent_sizes)[1] > Counter(both)[1]:
                # print("validate_reduced_seq - 2")
                return True
            if Counter(curent_sizes)[1] > 0:
                try:
                    if len(shielding_facets) == 0:  # You are free to do "remove_ones"
                        nonshielding_vids = set(np.nonzero(both)[0])
                    # print(f"curent_sizes={curent_sizes}; both={both}; nonshielding_vids = {nonshielding_vids}")
                    both, removed_sites = remove_ones(curent_sizes, both, choose_from=nonshielding_vids)
                except NoSlotError:
                    # print("validate_reduced_seq - 3")
                    return True
                else:
                    curent_sizes = curent_sizes[curent_sizes != 1]
                    # print(f"removed_sites = {removed_sites}; must_be_filled_vids={must_be_filled_vids}; exempt_vids={exempt_vids}")
                    to_be_added_facets = [exempt_vids + must_be_filled_vids.tolist() + [s] for s in removed_sites]
                    # print(f"to_be_added_facets = {to_be_added_facets}")
                    self.collected_facets += to_be_added_facets
                    # print(f"Adding {[list(removed_sites) + list(must_be_filled_vids)]} to collected_facets")
            exempt_vids += must_be_filled_vids.tolist()
            # print(f"exempt_vids = {exempt_vids}")
        # print(
        #     f"END validate_reduced_seq:: (both, curent_sizes, shielding_facets) = {both, curent_sizes, shielding_facets} \n")
        self.reduced_data = (both, curent_sizes, shielding_facets)
        # print(6)
        return False

    def _recursive_is_simplicial(self, degs, sizes, blocked_sets=None, verbose=False) -> (bool, list):
        try:
            st = SimplicialTest(degs, sizes, level=self._level, blocked_sets=blocked_sets, verbose=verbose)
            bool_ = st.is_simplicial(greedy=True)
        except (KeyError, NotImplementedError):
            return False, list()
        if bool_:
            return True, get_reduced_seq(degs, st.identifier2facets())
        else:
            return False, list()

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

    def get_selected_facet_ids(self, size) -> list:
        identifier = self.identifier
        return [index for index, i in enumerate(self.facet_size_per_id) if (index in identifier and i >= size)]

    def sample_simplex(self, size, greedy=False) -> (list, int):
        """

        Parameters
        ----------
        size
        greedy

        Returns
        -------

        """
        if greedy:
            return self.sample_simplex_greedy(size)

        identifier = self.identifier
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

    def is_simplicial(self, greedy=False) -> bool:
        if self._validate_data():
            return True
        elif self._validate_data() is False:
            return False

        # match_ones
        ind, self._sorted_s, self._sorted_d, matching_isolated = match_ones(self._sorted_s, self._sorted_d)
        if not ind:
            return False
        else:
            self.m = len(self._sorted_s)
            self.n = len(self._sorted_d)
            self.deg_seq = np.zeros(self.n, dtype=np.int_)
            self.matching["isolated"] = matching_isolated

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
            if validators.checkpoint_1(self._sorted_d, self.deg_seq):
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

    def _pull_the_plug(self, times=1) -> None:
        if not len(self.identifier) > times:
            raise ValueError("You cannot backtrack more times than the length of the identifier.")

        for _ in range(times):
            unwanted_facet = self.id2name[self.identifier.pop(-1)]
            self._sorted_s = np.insert(self._sorted_s, 0, len(unwanted_facet))
            self.update_deg_seq(unwanted_facet, -1)

    def count_explored_branches(self, identifier) -> int:
        s = set()
        for _ in self.logbook.keys():
            try:
                s.add(_[len(identifier) + 1])
            except IndexError:
                pass
        return len(s)

    def compute_joint_seq(self) -> (list, list):
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
