import time
import random
import validators

from utils import *
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
        self._level = level + 1
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

        self._counter = 0

        if blocked_sets is None:
            self.blocked_sets = list()
        else:
            self.blocked_sets = blocked_sets
        self.verbose = verbose
        self.collected_facets = list()
        self.current_facets = list()
        self.exempt_vids = list()
        self.callback_data = dict()  # level: facets

        self.mapping2shrinked = dict()
        self.mapping2enlarged = dict()

        self._DEGREE_LIST = _degree_list
        if self._DEGREE_LIST is None:
            self._DEGREE_LIST = self.DEGREE_LIST

        self.level_map = level_map
        if self.level_map is None:
            self.level_map = dict()
            self.level_map[self._level] = dict()
            for _ in range(self.n):
                self.level_map[self._level][_] = None
        # print(f"(l={self._level}) deg, size, blocked_sets = {self.DEGREE_LIST, self.SIZE_LIST, self.blocked_sets} \n")
        # print(f"(l={self._level}) {list(map(lambda x: len(x), self.blocked_sets))}")

    def prioritize(self):
        if len(self.blocked_sets) == 0:
            s = defaultdict(int)
        else:
            s = Counter(self.blocked_sets[-1])  # This is critical to efficiency.
        shielding = []
        non_shielding = []
        if self._level == 1:
            for _ in range(self.n):
                if self._sorted_d[_] - self.deg_seq[_] > 0:
                    if s[_] == 1:
                        shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
                    else:
                        non_shielding += [(_, self._sorted_d[_], self._sorted_d[_] - self.deg_seq[_])]
        else:
            remaining_degs = self.DEGREE_LIST - self.compute_joint_seq_from_identifier(sorted_deg=False)[1]
            for _ in self.level_map[1].keys():
                if self.level_map[self._level - 1][_] is None:
                    continue
                if remaining_degs[self.level_map[self._level - 1][_]] == 0:
                    continue
                if s[self.level_map[self._level - 1][_]] == 1:
                    shielding += [(self.level_map[self._level - 1][_], self._DEGREE_LIST[_])]
                else:
                    non_shielding += [(self.level_map[self._level - 1][_], self._DEGREE_LIST[_])]
        non_shielding = sorted(non_shielding, key=itemgetter(1), reverse=True)
        shielding = sorted(shielding, key=itemgetter(1), reverse=True)
        non_shielding = list(map(lambda x: x[0], non_shielding))
        shielding = list(map(lambda x: x[0], shielding))

        return shielding, non_shielding

    @staticmethod
    def get_valid_trials(shielding, non_shielding, size):
        shielding_ids = list(map(lambda x: x, shielding))
        non_shielding_ids = list(map(lambda x: x, non_shielding))
        n_s = len(shielding)
        n_ns = len(non_shielding)

        for _ in np.arange(1, n_ns + 1, +1):
            if size - _ > n_s or size - _ < 0:
                continue
            gen_ns = combinations(non_shielding_ids, _)
            gen = combinations(shielding_ids, size - _)

            while True:
                try:
                    _comb_ns = next(gen_ns)
                except StopIteration:
                    break

                while True:
                    try:
                        _comb = next(gen)
                    except StopIteration:
                        gen = combinations(shielding_ids, size - _)
                        break
                    yield _comb + _comb_ns
        raise NoMoreBalls

    def sample_simplex_greedy(self, size):
        identifier = self.identifier
        shielding, non_shielding = self.prioritize()
        valid_trials = self.get_valid_trials(shielding, non_shielding, size)

        while True:
            try:
                candidate_facet = tuple([_ for _ in next(valid_trials)])
            except RuntimeError:
                raise NoMoreBalls(
                    "May not be solvable with the greedy algorithm OR the recursive_is_simplicial must return False.")
            # Note that self.blocked_sets are passed from a higher level,
            # whereas larger_simplices is anew for each level.
            if validators.validate_issubset_blocked_sets(candidate_facet, self.blocked_sets):
                continue
            if tuple(sorted(list(candidate_facet))) in self.logbook:
                continue
            ind, reason = self.validate(identifier, candidate_facet)
            if ind:
                picked_facet, picked_facet_id = self.register(candidate_facet)
                return picked_facet, picked_facet_id
            else:
                for key in self.level_map.keys():
                    if key >= self._level:
                        self.level_map[key] = {}

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
        _sizes = self._sorted_s
        if len(_sizes) == 0:
            return True, "0"
        _both = self.get_remaining_slots(identifier, candidate_facet)  # both (non_shielding & shielding)
        if np.min(_both) < 0:
            return False, "1"
        #
        if len(_sizes) == 1:
            for facet in current_facets:
                if set(np.nonzero(_both)[0]).issubset(set(facet)):
                    return False, "2"
        #
        if np.any(_both > len(_sizes)):
            return False, "3"
        if np.sum(_both) > np.count_nonzero(_both) == _sizes[0]:
            return False, "4"

        _non_shielding = self.get_remaining_slots(identifier, candidate_facet, only_non_shielding=True)
        _shielding = _both - _non_shielding  # only shielding

        if np.min(_non_shielding) < 0:
            return False, "5"
        if validators.validate_nonshielding(_sizes, _non_shielding, _shielding):
            return False, "6"
        if np.count_nonzero(_both) < _sizes[0]:
            return False, "7"

        for blocked_set in self.blocked_sets:
            if set(np.nonzero(_both)[0]).issubset(set(blocked_set)):
                return False, "8"

        if Counter(_both)[len(_sizes)] > 0:
            _ = validators.validate_reduced_seq(_both, _sizes, current_facets)
            if _[0]:
                return False, "9"
            _both, _sizes, self.collected_facets, self.exempt_vids = _[1]
            if np.sum(_both) == np.sum(_sizes) == 0:
                for facet in current_facets:
                    for collected_facet in self.collected_facets:
                        if set(collected_facet).issubset(set(facet)):
                            return False, "10"
                return True, "0"
        else:
            self.exempt_vids = []
            self.collected_facets = []
        self.mapping2shrinked, self.mapping2enlarged = get_seq2seq_mapping(_both)
        filtered = filter_blocked_facets(current_facets + self.blocked_sets, self.exempt_vids)
        _blocked_sets = shrink_facets(filtered, self.mapping2shrinked)

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
            if int(str(e)) <= 3:
                # keep look for solutions at this level, until get_valid_trials emits a NoMoreBalls
                return False, "keep looking for balls!"
            else:
                # forget about this level
                raise NoMoreBalls
        if bool_:
            facets = [self.exempt_vids + list(s) for s in facets]
            self.current_facets = current_facets
            all_facets = self.current_facets + facets + self.collected_facets
            self.callback_data[self._level] = all_facets
            # print(
            #     f"(l={self._level}) WeCanStopSignal fired::"
            #     f"id2name: {self.id2name} "
            #     f"accept {candidate_facet}, callback_data: {all_facets};"
            #     f"... within which facets={facets}, current_facets={self.current_facets}, collected_facets={self.collected_facets}; \n"
            # )
            raise WeCanStopSignal
        else:
            return bool_, "11"

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
            bool_, facets = st.is_simplicial()  # facets: facets from a deeper level
        except NoMoreBalls:
            raise NoMoreBalls(self._level)
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

    def is_simplicial(self) -> (bool, dict):
        if validators.validate_data(self._sorted_d, self._sorted_s):
            return True, dict()  # TODO, it should be facets instead of dict()
        elif validators.validate_data(self._sorted_d, self._sorted_s) is False:
            return False, dict()

        if len(self._sorted_s) == 0:
            if sorted(self.deg_seq, reverse=True) == self._sorted_d.tolist():  # TODO: start from identifier
                return True, dict()  # TODO, it should be facets instead of dict()
            else:
                return False, dict()

        while self._counter < 1e3:
            self._counter += 1
            s = self._sorted_s[0]
            self._sorted_s = np.delete(self._sorted_s, 0)
            try:
                picked_facet, picked_facet_id = self.sample_simplex_greedy(s)
            except WeCanStopSignal:
                if self._level - 1 == 0:  # the very first level
                    return True, self.callback_data[self._level]
                return True, self.callback_data
            except NoMoreBalls:
                raise NoMoreBalls
            else:
                self.update_deg_seq(picked_facet, +1)
                if validators.validate_interm_degs(self.DEGREE_LIST, self.deg_seq):
                    self.identifier += [picked_facet_id]
                else:
                    # Backtrack
                    self.log_forbidden(tuple(sorted(list(picked_facet))), "failed validate_interm_degs")
                    self.update_deg_seq(picked_facet, -1)
                    self._sorted_s = np.insert(self._sorted_s, 0, len(self.id2name[picked_facet_id]))
                    continue
                # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
                if len(self._sorted_s) == 0:
                    self.callback_data[self._level] = self.identifier2facets()
                    if self._level - 1 == 0:  # the very first level
                        return True, self.callback_data[self._level]
                    # print(f"(l={self._level}) Sending back {self.callback_data}")
                    return True, self.callback_data
        return False, dict()

    def update_deg_seq(self, facet, value):
        if value not in [+1, -1]:
            raise NotImplementedError
        for _ in facet:
            self.deg_seq[_] += value

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
