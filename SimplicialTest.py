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

    def __init__(
            self,
            degree_list,
            size_list,
            level=0,
            blocked_sets=None,
            _degree_list=None,
            level_map=None,
            verbose=False):
        super().__init__()

        self.random_seed = time.thread_time_ns()
        random.seed(self.random_seed)

        self.logbook = defaultdict(dict)

        self.m = len(size_list)
        self.n = len(degree_list)

        self.DEGREE_LIST = np.array(sorted(degree_list, reverse=True), dtype=np.int_)
        self.SIZE_LIST = np.array(sorted(size_list, reverse=True), dtype=np.int_)

        self._sorted_s = np.array(deepcopy(self.SIZE_LIST), dtype=np.int_)
        self._sorted_d = np.array(deepcopy(self.DEGREE_LIST), dtype=np.int_)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.identifier = list()

        self._backtrack_steps = 0
        self._counter = 0
        self._len_logbook = 0
        self._level = level + 1

        if blocked_sets is None:
            self.blocked_sets = list()
        else:
            self.blocked_sets = blocked_sets
        # print(f"(l={self._level}) initiated self.blocked_sets = {self.blocked_sets}")
        self.verbose = verbose
        self.collected_facets = list()
        self.current_facets = list()
        self.exempt_vids = list()
        self.callback_data = dict()  # level: facets

        self.mapping2shrinked = dict()
        self.mapping2enlarged = dict()
        self.level_toggle = False
        self._DEGREE_LIST = _degree_list
        if self._DEGREE_LIST is None:
            self._DEGREE_LIST = self.DEGREE_LIST
        self.level_map = level_map
        if self.level_map is None:
            self.level_map = dict()
            self.level_map[self._level] = dict()
            for _ in range(self.n):
                self.level_map[self._level][_] = None


    def prioritize(self):
        # NOTE: deg_seq is all zero because of the recursive nature... (todo: this should be fixed.)
        # print(f"(l={self._level}) _sorted_d, deg_seq = {self._sorted_d, self.deg_seq}")
        # print(f"(l={self._level}) blocked_sets = {self.blocked_sets} \n")
        # if self._level > 1:
        #     level_map = self.level_map[self._level - 1]
        # else:
        #     return [], range(self.n)
        # if self.verbose:
        #     print(
        #         f"(l={self._level}) PRIORITIZE:: _DEGREE_LIST = {self._DEGREE_LIST}; "
        #         f"blocked_sets = {self.blocked_sets}; "
        #         f"level_map: {level_map}. "
        #         f"level_map -> values() = {list(level_map.values())}\n")
        # d = Counter()
        # for _ in map(lambda x: Counter(x), self.blocked_sets):
        #     d += _
        # if self.verbose:
        #     print(f"(l={self._level}) PRIORITIZE:: d = {d} \n")





        # # USE THIS
        if self._level == 1:
            if len(self.blocked_sets) == 0:
                s = defaultdict(int)
            else:
                s = Counter(self.blocked_sets[0])

            # print(s, self.id2name[identifier[-1]], identifier)
            shielding = []
            non_shielding = []
            for _ in range(self.n):
                if self._sorted_d[_] - self.deg_seq[_] > 0:
                    if s[_] == 1:
                        shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                    else:
                        non_shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
            non_shielding = sorted(non_shielding, key=itemgetter(1), reverse=True)  # originally: 1, true.... exp: 2, false
            shielding = sorted(shielding, key=itemgetter(1), reverse=True)  # originally: 2, true
            non_shielding = list(map(lambda x: x[0], non_shielding))
            shielding = list(map(lambda x: x[0], shielding))
        else:
        # # USE THIS



        # NOTE: deg_seq is all zero because of the recursive nature... (todo: this should be fixed.)
        # print(f"(l={self._level}) _sorted_d, deg_seq = {self._sorted_d, self.deg_seq}")
        # print(f"(l={self._level}) blocked_sets = {self.blocked_sets} \n")
        # if self._level > 1:
        #     level_map = self.level_map[self._level - 1]
        # else:
        #     return [], range(self.n)

            if len(self.blocked_sets) == 0:
                s = defaultdict(int)
            else:
                s = Counter(self.blocked_sets[0])

            # print(s, self.id2name[identifier[-1]], identifier)
            shielding = []
            non_shielding = []
            for _ in self.level_map[1].keys():
                if self.level_map[self._level - 1][_] is None:
                    continue
                if s[self.level_map[self._level - 1][_]] == 1:
                    shielding += [(self.level_map[self._level - 1][_], self._DEGREE_LIST[_])]
                else:
                    non_shielding += [(self.level_map[self._level - 1][_], self._DEGREE_LIST[_])]
            non_shielding = sorted(non_shielding, key=itemgetter(1), reverse=True)
            shielding = sorted(shielding, key=itemgetter(1), reverse=True)

            non_shielding = list(map(lambda x: x[0], non_shielding))
            shielding = list(map(lambda x: x[0], shielding))
            print(f"(l={self._level}) shielding, non_shielding = {shielding, non_shielding}")







        return shielding, non_shielding

    # @staticmethod
    # def get_valid_trials(shielding, non_shielding, size):
    #     n_shielding = len(shielding)
    #     n_non_shielding = len(non_shielding)
    #     n_ = n_shielding + n_non_shielding
    #     # print(f"n_shielding, n_non_shielding = {n_shielding, n_non_shielding}")
    #
    #     i = 0
    #
    #     gen = combinations(range(size + i), size)
    #     # print(f"comb of {range(size + i), size}")
    #     while size + i <= n_:
    #         if size + i < n_shielding:
    #             i += 1
    #             continue
    #         try:
    #             _comb = next(gen)
    #         except StopIteration:
    #             i += 1
    #             gen = combinations(range(size + i), size)
    #             # print(f"comb of {range(size + i), size}")
    #         else:
    #             if not set(_comb).issubset(set(range(n_shielding))):
    #                 # print(f"_comb = {_comb}")
    #                 yield _comb
    #             else:
    #                 # print(f"_comb = {_comb} is rejected")
    #                 pass

    def get_valid_trials(self, shielding, non_shielding, size):
        shielding_ids = list(map(lambda x: x, shielding))
        non_shielding_ids = list(map(lambda x: x, non_shielding))
        if self.verbose:
            print(
                f"(l={self._level}) size={size}; shielding_ids = {shielding_ids} \nnon_shielding_ids = {non_shielding_ids}")
        n_s = len(shielding)
        n_ns = len(non_shielding)

        d = Counter()
        for _ in map(lambda x: Counter(x), self.blocked_sets):
            d += _
        if self.verbose:
            print(f"(l={self._level}) PRIORITIZE:: d = {d} \n")

        # print(f"n_s, n_ns = {n_s, n_ns}")
        # if size > n_ns:
        #     raise NoSlotError
        # print(f"remaining = {remaining}")
        # gen = combinations(range(n), size)
        # while True:
        #     _comb = next(gen)
        #     # print(f"(l={self._level}) we have {_comb} in the combinations")
        #     yield _comb

        for _ in np.arange(1, n_ns + 1, +1):
            # print(f"size = {size}; _ = {_}; n_s, n_ns = {n_s, n_ns}")
            # print(f"size - _ = {size - _}")

            if size - _ > n_s or size - _ < 0:
                continue
            # if size - _ < 0:
            #     break
            gen_ns = combinations(non_shielding_ids, _)
            gen = combinations(shielding_ids, size - _)

            while True:
                try:
                    _comb_ns = next(gen_ns)
                    if self.verbose:
                        print(f"(l={self._level}) mixing {_comb_ns} in the combinations")
                    # print(np.array(remaining)[list(_comb_ns)] \n)
                except StopIteration:
                    break

                while True:
                    try:
                        _comb = next(gen)
                    except StopIteration:
                        gen = combinations(shielding_ids, size - _)
                        break
                    if self.verbose:
                        print(f"(l={self._level}) ---> ", _comb + _comb_ns)
                    yield _comb + _comb_ns
        # print(f"(l={self._level}) get_valid_trials::raise NoMoreBalls")
        raise NoMoreBalls

    def sample_simplex_greedy(self, size):
        identifier = self.identifier
        # if len(identifier) == 0:
        #     # print(f"(l={self._level}) compute ice breaker...")
        #     return self.sample_icebreaker(size)
        # larger_selected_simplex_ids = self.get_selected_facet_ids(size)
        # larger_simplices = [self.id2name[_] for _ in larger_selected_simplex_ids]

        # candidate_facet = []
        shielding, non_shielding = self.prioritize()

        # options = shielding + non_shielding
        # remaining = self._sorted_d - self.compute_joint_seq_from_identifier(identifier, sorted_deg=False)[1]
        # print(f"before running simple simplex greedy: {remaining}")
        valid_trials = self.get_valid_trials(shielding, non_shielding, size)
        # print(f"(l={self._level}) blocked sites are = {self.blocked_sets}")
        while True:
            try:
                # print(f"(l={self._level}) q = {q}")
                # candidate_facet = tuple([options[_][0] for _ in q])
                candidate_facet = tuple([_ for _ in next(valid_trials)])
            except RuntimeError:
                # print(f"(l={self._level}) stop iteration in sample_simplex_greedy called; ")
                raise NoMoreBalls(
                    "May not be solvable with the greedy algorithm OR the recursive_is_simplicial must return False.")

            # Note that self.blocked_sets are passed from a higher level,
            # whereas larger_simplices is anew for each level.
            if validators.validate_issubset_blocked_sets(candidate_facet, self.blocked_sets):
                # print(f"{candidate_facet} is a subset of some {self.blocked_sets}")
                continue

            # if validators.validate_issubset_larger_simplices(candidate_facet, larger_simplices):
            #     continue

            # if candidate_facet == (0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15):
            #     print(f"******* (l={self._level}) select w/ size {size}: candidate_facet = {candidate_facet} ")
            #     self.verbose = True

            if self.verbose:
                print(
                    f"(l={self._level}) new candidate_facet selected = {candidate_facet}, testing it, possibly to a deeper level.")
            ind, reason = self.validate(identifier, candidate_facet)
            if ind:
                print(f"(l={self._level}) select w/ size {size}: candidate_facet = {candidate_facet} ACCEPTED!")
                picked_facet, picked_facet_id = self.register(candidate_facet)
                # print(
                #     f"(l={self._level}) picked_facet: {picked_facet}, having already {self.identifier2facets(identifier)}; blocked sites are = {self.blocked_sets}")
                return picked_facet, picked_facet_id
            else:
                for key in self.level_map.keys():
                    if key >= self._level:
                        self.level_map[key] = {}
                # print(
                #         f"(l={self._level}) select w/ size {size}: candidate_facet = {candidate_facet} REJECTED due to {reason}!")
                pass

    def validate(self, identifier, candidate_facet) -> (bool, str):
        """
        This function must return True in order for the candidate facet to be considered.
        In other words, if any condition appears to be True, we will not accept the candidate facet.

        Parameters
        ----------
        identifier
        candidate_facet

        Returns
        -------

        """
        current_facets = self.identifier2facets(identifier) + [candidate_facet]
        # print(f"(l={self._level}) Testing {candidate_facet} when {self.identifier2facets(identifier)} existed & {self.blocked_sets} blocked")

        # if validators.validate_issubset(self.id2name, candidate_facet, _id):
        #     return False, "1"
        _sizes = self._sorted_s
        if len(_sizes) == 0:
            return True, "0"


        _both = self.get_remaining_slots(identifier, candidate_facet)  # both (non_shielding & shielding)
        if np.min(_both) < 0:
            return False, "13"
        #
        if len(_sizes) == 1:
            for facet in current_facets:
                if set(np.nonzero(_both)[0]).issubset(set(facet)):
                    return False, "12"
        #
        if np.any(_both > len(_sizes)):
            return False, "4"
        if np.sum(_both) > np.count_nonzero(_both) == _sizes[0]:
            return False, "2"

        _non_shielding = self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)
        _shielding = _both - _non_shielding  # only shielding

        if np.min(_non_shielding) < 0:
            return False, "14"
        #
        if validators.validate_nonshielding(_sizes, _non_shielding, _shielding):
            return False, "5"
        if np.count_nonzero(_both) < _sizes[0]:
            return False, "3"
        if validators.validate_interm_degs(self._sorted_d, candidate_facet):
            return False, "8i9"
        for blocked_set in self.blocked_sets:
            if set(np.nonzero(_both)[0]).issubset(set(blocked_set)):
                return False, "11"
        # print(
        #     f"(l={self._level}) ->> _both={_both}")
        # if len(current_facets) > 1:
        if Counter(_both)[len(_sizes)] > 0:
            _ = validators.validate_reduced_seq(_both, _sizes, current_facets)
            if _[0]:
                return False, "6"
            _both, _sizes, self.collected_facets, self.exempt_vids = _[1]
            # print(f"(l={self._level}) ->> _both={_both}; len(_both) = {len(_both)}")
            if np.sum(_both) == np.sum(_sizes) == 0:
                for facet in current_facets:
                    for collected_facet in self.collected_facets:
                        if set(collected_facet).issubset(set(facet)):
                            return False, "10"
                return True, "0"
        self.mapping2shrinked, self.mapping2enlarged = get_seq2seq_mapping(_both)
        filtered = filter_blocked_facets(current_facets + self.blocked_sets, self.exempt_vids)
        _blocked_sets = shrink_facets(filtered, self.mapping2shrinked)
        # print(f"(l={self._level}) self.mapping2shrinked = {self.mapping2shrinked}")
        # print(f"(l={self._level}) self.level_map = {self.level_map} \n")

        if candidate_facet == (0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12, 14, 15):
            self.verbose = True
            # print(
            #     f"self.callback_data={self.callback_data}; current_facets={current_facets}; len(_sizes)={len(_sizes)}")
        # print(s_degs_whole, s_degs_ns, '\n')
        # exit(0)
        #
        # print(
        #     f"(l={self._level}) "
        #     f"candidate_facet = {candidate_facet}"
        #     f"recursive input (degs, sizes) = {(_both, _sizes)} current_facets={current_facets}; blocked_sets={_blocked_sets} \n"
        # )
        try:
            self._compute_level_map()
            bool_, facets = self._recursive_is_simplicial(
                _both, _sizes,
                blocked_sets=_blocked_sets,
                level_map=self.level_map,
                _degree_list=self._DEGREE_LIST,
                verbose=self.verbose
            )
        except NoMoreBalls as e:
            if int(str(e)) <= 2:
                # keep look for solutions at this level, until get_valid_trials emits a NoMoreBalls
                return False, "keep looking for balls!"
            else:
                # forget about this level
                raise NoMoreBalls

            # return False, f"we caught an error in sample_simplex_greedy, at a deeper level = {self._level + 1}"
        # print(
        #     f"(l={self._level}) "
        #     f"After _recursive_is_simplicial, we have bool_, facets = {bool_, facets};"
        # )
        if bool_:
            facets = [self.exempt_vids + list(s) for s in facets]
            self.current_facets = current_facets
            # print(f"(l={self._level}) self.exempt_vids = {self.exempt_vids}")
            # print(f"(l={self._level}) ???? = {[self.exempt_vids + list(s) for s in facets] }")
            all_facets = self.current_facets + facets + self.collected_facets
            self.callback_data[self._level] = all_facets
            print(
                f"(l={self._level}) WeCanStopSignal fired::"
                f"id2name: {self.id2name} "
                f"accept {candidate_facet}, callback_data: {all_facets}; \n"
            )
            raise WeCanStopSignal
        else:
            # print(
            #     f"(l={self._level}) "
            #     f"candidate_facet = {candidate_facet} is not accepted. \n"
            # )
            return bool_, "333"
        # return True, "0"

    def _recursive_is_simplicial(self,
                                 degs,
                                 sizes,
                                 blocked_sets=None,
                                 level_map=None,
                                 _degree_list=None,
                                 verbose=False) -> (bool, list):

        try:
            st = SimplicialTest(
                degs,
                sizes,
                level=self._level,
                blocked_sets=blocked_sets,
                level_map=level_map,
                _degree_list=self._DEGREE_LIST,
                verbose=verbose
            )
            bool_, facets = st.is_simplicial(greedy=True)  # facets: facets from a deeper level
        except NoMoreBalls:
            raise NoMoreBalls(self._level)
            # return False, list()
        if bool_:
            # transform the facets collected from a deeper level
            return True, get_enlarged_seq(self.mapping2enlarged, facets[self._level + 1])
        else:
            return False, list()

    def _compute_level_map(self):
        if self._level > 1:
            self.level_map[self._level] = dict()
            n = len(self.level_map[1].keys())
            for _ in range(n):
                if self.level_map[self._level - 1][_] is None:
                    self.level_map[self._level][_] = None
                    continue
                try:
                    self.level_map[self._level][_] = self.mapping2shrinked[self.level_map[self._level - 1][_]]
                except KeyError:
                    self.level_map[self._level][_] = None
        else:
            n = self.n
            for _ in range(n):
                try:
                    self.level_map[1][_] = self.mapping2shrinked[_]
                except KeyError:
                    self.level_map[1][_] = None

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

    def is_simplicial(self, greedy=False) -> (bool, dict):
        if validators.validate_data(self._sorted_d, self._sorted_s):
            return True, dict()  # TODO, it should be facets instead of dict()
        elif validators.validate_data(self._sorted_d, self._sorted_s) is False:
            return False, dict()

        if len(self._sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d.tolist():  # TODO: start from identifier
                return True, dict()  # TODO, it should be facets instead of dict()
            else:
                return False, dict()

        while self._counter < 10:
            # if len(self.logbook) == self._len_logbook:
            #     self._counter += 1
            s = self._sorted_s[0]
            self._sorted_s = np.delete(self._sorted_s, 0)
            try:
                picked_facet, picked_facet_id = self.sample_simplex_greedy(s)
            except WeCanStopSignal:
                if self._level - 1 == 0:  # the very first level
                    return True, self.callback_data[self._level]
                # print(f"(l={self._level}) Sending back {self.callback_data}")
                return True, self.callback_data
            except NoMoreBalls:
                # print(f"(l={self._level}) is_simplicial::raising NoMoreBalls")
                raise NoMoreBalls
            else:
                self.identifier += [picked_facet_id]
                # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
                if len(self._sorted_s) == 0:
                    self.callback_data[self._level] = self.identifier2facets()
                    if self._level - 1 == 0:  # the very first level
                        return True, self.callback_data[self._level]
                    print(f"(l={self._level}) Sending back {self.callback_data}; picked_facet={picked_facet}")
                    return True, self.callback_data
        return False, dict()
    #
    # def _pull_the_plug(self, times=1) -> None:
    #     if not len(self.identifier) > times:
    #         raise ValueError("You cannot backtrack more times than the length of the identifier.")
    #
    #     for _ in range(times):
    #         unwanted_facet = self.id2name[self.identifier.pop(-1)]
    #         self._sorted_s = np.insert(self._sorted_s, 0, len(unwanted_facet))
    #         self.update_deg_seq(unwanted_facet, -1)

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

    # def count_explored_branches(self, identifier) -> int:
    #     s = set()
    #     for _ in self.logbook.keys():
    #         try:
    #             s.add(_[len(identifier) + 1])
    #         except IndexError:
    #             pass
    #     return len(s)
    #
    # def sample_simplex(self, size) -> (list, int):
    #     """
    #
    #     Parameters
    #     ----------
    #     size
    #     greedy
    #
    #     Returns
    #     -------
    #
    #     """
    #     identifier = self.identifier
    #     # Here, we have a good criterion! We may not need to explore further if...
    #     deg_sequence = np.array(self.compute_joint_seq_from_identifier(identifier)[1])
    #     deg_sequence_goal = self._sorted_d
    #     if np.any(deg_sequence_goal - deg_sequence < 0):
    #         self.log_forbidden(identifier, "NO need to explore further - 1")
    #         self._backtrack_steps = 1
    #         return list(), -1
    #     if len(np.nonzero(deg_sequence_goal - deg_sequence)[0]) < size:
    #         self.log_forbidden(identifier, "NO need to explore further - 2")
    #         self._backtrack_steps = 1
    #         return list(), -1
    #
    #     larger_selected_simplex_ids = self.get_selected_facet_ids(size)
    #
    #     set_of_vertices = set(range(self.n))
    #     picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
    #     qualified_draw = False
    #     _ = 0  # I think this is stupid.... but let's try this for now...
    #     while not qualified_draw:
    #         qualified_draw = True
    #         for _id in larger_selected_simplex_ids:
    #
    #             if set(picked_facet).issubset(set(self.id2name[_id])):
    #                 self.log_forbidden(identifier + [picked_facet_id], 1)
    #                 picked_facet, picked_facet_id = self.register(random.sample(list(set_of_vertices), k=size))
    #                 qualified_draw = False
    #                 _ += 1
    #                 break
    #         explored = [key for key in self.logbook.keys() if
    #                     key[:len(identifier)] == tuple(identifier) and len(key) == len(identifier) + 1]
    #
    #         if len(explored) == comb(self.n, size):
    #             if len(identifier) == 2:
    #                 self._backtrack_steps = 1
    #             else:
    #                 self._backtrack_steps = 2
    #             return list(), -1
    #         if _ > 100:  # TODO
    #             self._backtrack_steps = 1
    #             return list(), -1
    #
    #     return picked_facet, picked_facet_id
    #
    # def update_deg_seq(self, facet, value):
    #     if value not in [+1, -1]:
    #         raise NotImplementedError
    #     for _ in facet:
    #         self.deg_seq[_] += value
