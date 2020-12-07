from simplicial_test import validators
from simplicial_test.utils import *
from operator import itemgetter
from itertools import combinations

from .mongo import DB

import sys

sys.setrecursionlimit(10 ** 6)


class Test(SimplexRegistrar):
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
            _dlist=None,
            _slist=None,
            level_map=None,
            verbose=False):
        super().__init__()
        self.logbook = defaultdict(dict)

        self.m = len(size_list)
        self.n = len(degree_list)

        self.DEGREE_LIST = np.array(sorted(degree_list, reverse=True), dtype=np.int_)
        self.SIZE_LIST = np.array(sorted(size_list, reverse=True), dtype=np.int_)

        self._sorted_s = np.array(deepcopy(self.SIZE_LIST), dtype=np.int_)
        self._sorted_d = np.array(deepcopy(self.DEGREE_LIST), dtype=np.int_)

        self.deg_seq = np.zeros(self.n, dtype=np.int_)
        self.identifier = list()

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

        self._DEGREE_LIST = _dlist
        if self._DEGREE_LIST is None:
            self._DEGREE_LIST = self.DEGREE_LIST
        self._SIZE_LIST = _slist
        if self._SIZE_LIST is None:
            self._SIZE_LIST = self.SIZE_LIST

        self.conn = DB("simplicial", "deadends", self._DEGREE_LIST.tolist(), self._SIZE_LIST.tolist())
        if level == 0:
            self.conn.init()

        self._level = level + 1
        self.level_map = level_map
        if self.level_map is None:
            self.level_map = dict()
            self.level_map[self._level] = dict()
            for _ in range(self.n):
                self.level_map[self._level][_] = None
        self.num_failures = 0
        if verbose:
            print(f"---------- (l={self._level}) ----------\n"
                  f"Size list: {self.SIZE_LIST.tolist()}\n"
                  f"Degree list: {self.DEGREE_LIST.tolist()}\n"
                  # f"blocked_sets = {self.blocked_sets}\n"
                  f"len::blocked_sets = {list(map(lambda x: len(x), self.blocked_sets))}\n"
                  )

    def prioritize(self):
        if len(self.blocked_sets) == 0:
            s = defaultdict(int)
        else:
            # self.blocked_sets.sort(key=lambda x: np.sum(x))
            # s = Counter(self.blocked_sets[-1])  # This is critical to efficiency.
            self.blocked_sets.sort(key=lambda x: len(x))
            s = Counter(self.blocked_sets[-1])  # This is critical to efficiency.
        shielding = []
        non_shielding = []
        if self._level == 1:
            for _ in range(self.n):
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

    def get_valid_trials(self, size):
        shielding, non_shielding = self.prioritize()
        for _ in np.arange(1, len(non_shielding) + 1, +1):
            if size - _ > len(shielding) or size - _ < 0:
                continue
            gen_ns = combinations(non_shielding, _)
            gen = combinations(shielding, size - _)
            while True:
                try:
                    _comb_ns = next(gen_ns)
                except StopIteration:
                    break

                while True:
                    try:
                        _comb = next(gen)
                    except StopIteration:
                        gen = combinations(shielding, size - _)
                        break
                    yield _comb + _comb_ns
        raise NoMoreBalls("No more usable proposals.")

    def sample_simplex_greedy(self, size):
        identifier = self.identifier

        valid_trials = self.get_valid_trials(size)
        # valid_trials = self.new_get_valid_trials(size)

        while True:
            try:
                candidate_facet = tuple([_ for _ in next(valid_trials)])
            except RuntimeError:
                raise NoMoreBalls(
                    "May not be solvable with the greedy algorithm OR the recursive_is_simplicial must return False.")
            except NoMoreBalls:
                raise NoMoreBalls("No more usable proposals.")
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
                self.num_failures += 1
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
            return True, "Last facet explored."
        _both = self.get_remaining_slots(identifier, candidate_facet)  # both (non_shielding & shielding)
        if np.min(_both) < 0:
            return False, "Invalid facet choice. Negative remaining degrees detected."

        if len(_sizes) == 1:
            for facet in current_facets:
                if set(np.nonzero(_both)[0]).issubset(set(facet)):
                    return False, "Second last facet explored. " \
                                  "But the very last facet is doomed to fail the no-inclusion constraint."

        if np.any(_both > len(_sizes)):
            return False, "Some degrees require more facets than the actual remaining number."
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
                return False, "The remaining facets are doomed to fail the no-inclusion constraint."

        if Counter(_both)[len(_sizes)] > 0:
            _ = validators.validate_reduced_seq(_both, _sizes, current_facets)
            if _[0]:
                return False, "9"
            _both, _sizes, self.collected_facets, self.exempt_vids = _[1]
            if np.sum(_both) == np.sum(_sizes) == 0:
                for facet in current_facets:
                    for collected_facet in self.collected_facets:
                        if set(collected_facet).issubset(set(facet)):
                            return False, "Last facet explored after must-do pairing. " \
                                          "But no-inclusion constraint failed."
                return True, "Last facet explored after must-do pairing."
        else:
            self.exempt_vids = []
            self.collected_facets = []
        self.mapping2shrinked, self.mapping2enlarged = get_seq2seq_mapping(_both)
        filtered = filter_blocked_facets(current_facets + self.blocked_sets, self.exempt_vids)
        _blocked_sets = shrink_facets(filtered, self.mapping2shrinked)
        # if self.verbose:
        #     print(
        #         f"(l={self._level}) Selected {candidate_facet} (next validating at l={self._level + 1})"
        #     )
        try:
            # try:
            #     inv_level_map = {v: k for k, v in self.level_map[self._level - 1].items()}
            #     translated = []
            #     for vtx in candidate_facet:
            #         translated += [inv_level_map[vtx]]
            #     print(f"(l={self._level}) okay... we selected (in orig coord){translated}; bsets={self.blocked_sets}")
            # except:
            #     print(f"(l={self._level}) okay... we selected {candidate_facet}; bsets={self.blocked_sets}")
            #     pass
            self._compute_level_map()
            deeper_facet_found, facets = self._recursive_is_simplicial(
                _both, _sizes,
                blocked_sets=_blocked_sets,
                level_map=self.level_map,
                _dlist=self._DEGREE_LIST,
                _slist=self._SIZE_LIST,
                verbose=self.verbose
            )
        except NoMoreBalls as e:

            if int(str(e)) <= 10:
                # print(f"(l={self._level}) num_failures={self.num_failures}")
                if self.num_failures > 20:
                    self.conn.add_to_failed_attempts()
                    # print(f"(l={self._level}) to many fails; skip this level")
                    raise NoMoreBalls
                # keep look for solutions at this level, until get_valid_trials emits a NoMoreBalls
                # print(f"(l={self._level}) keep looking")
                self.conn.add_to_failed_attempts()
                return False, "keep looking for balls!"
            else:
                # forget about this level
                self.conn.add_to_failed_attempts()
                # print(f"(l={self._level}) skip this level")
                raise NoMoreBalls
        if deeper_facet_found:
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
            return False, "No valid facet found at a higher level."

    def _recursive_is_simplicial(self,
                                 degs,
                                 sizes,
                                 blocked_sets=None,
                                 level_map=None,
                                 _dlist=None,
                                 _slist=None,
                                 verbose=False) -> (bool, list):
        try:
            st = Test(
                degs,
                sizes,
                level=self._level,
                blocked_sets=blocked_sets,
                level_map=level_map,
                _dlist=self._DEGREE_LIST,
                _slist=self._SIZE_LIST,
                verbose=verbose
            )
            # facets: facets from a deeper level
            deeper_facet_is_simplicial, facets = st.is_simplicial()
        except NoMoreBalls:
            raise NoMoreBalls(self._level)
        if deeper_facet_is_simplicial:
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

    def is_simplicial(self):
        if not validators.validate_data(self._sorted_d, self._sorted_s):
            self.conn.add_to_failed_attempts()
            return False, tuple()

        while True:
            s = self._sorted_s[0]
            self._sorted_s = np.delete(self._sorted_s, 0)
            try:
                picked_facet, picked_facet_id = self.sample_simplex_greedy(s)
            except WeCanStopSignal:
                if self._level - 1 == 0:  # the very first level
                    self.conn.mark_simplicial()
                    return True, sort_callback(self.callback_data[self._level])
                return True, self.callback_data
            except NoMoreBalls:
                raise NoMoreBalls(f"No more usable proposals at level={self._level}; roll back or say no.")
            else:
                update_deg_seq(self.deg_seq, picked_facet, +1)
                if validators.validate_interm_degs(self.DEGREE_LIST, self.deg_seq):
                    # print(f"(l={self._level}) validate_interm_degs check!")
                    self.identifier += [picked_facet_id]
                else:
                    # Backtrack
                    # print(f"(l={self._level}) validate_interm_degs failed!")
                    self.log_forbidden(tuple(sorted(list(picked_facet))), "failed validate_interm_degs")
                    update_deg_seq(self.deg_seq, picked_facet, -1)
                    self._sorted_s = np.insert(self._sorted_s, 0, len(self.id2name[picked_facet_id]))
                    continue
                # Here, assuming our algorithm is all good, we want to check if we indeed find the simplicial complex
                if len(self._sorted_s) == 0:
                    self.callback_data[self._level] = self.identifier2facets()
                    if self._level - 1 == 0:  # the very first level
                        self.conn.mark_simplicial()
                        return True, sort_callback(self.callback_data[self._level])
                    # print(f"(l={self._level}) Sending back {self.callback_data}")
                    return True, self.callback_data
        return False, tuple()

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
