from simplicial_test.utils import *
from simplicial_test.sample import get_hitting_sets
from copy import deepcopy
from itertools import combinations_with_replacement
from collections import defaultdict


class EnumRegistrar(object):
    def __init__(self):
        self.pointer = 0
        self.facets_pointer = 0

        self.facet2id = dict()
        self.id2facet = dict()

        self.facets2id = dict()
        self.id2facets = dict()

        self.states = dict()

        self.facet_size_per_id = np.array([], dtype=np.int_)
        self.logbook = dict()

        self.facet_book = dict()

    def register_facet(self, facet) -> (tuple, int):
        facet = tuple(sorted(facet, reverse=False))
        if facet not in self.facet2id:
            self.facet2id[facet] = self.pointer
            self.id2facet[self.pointer] = facet
            self.pointer += 1
        return facet, self.facet2id[facet]

    def register_facets(self, facets):
        sorted_facets = self.sort_facets(facets)
        if sorted_facets not in self.facets2id:
            self.facets2id[sorted_facets] = self.facets_pointer
            self.id2facets[self.facets_pointer] = sorted_facets
            self.facets_pointer += 1
        return sorted_facets, self.facets2id[sorted_facets]

    def register_state(self, facets):
        facets = self.sort_facets(facets)
        if self.facets2id[facets] not in self.states:
            state = self._compute_state(facets)
            vtx_symm_class = self.groupby_vtx_symm_class(self.dict2tuple(state[1]))
            self.states[self.facets2id[facets]] = {
                "n_facets": len(facets),
                "deg_per_vid": state[0],
                "vid_per_facet": state[1],
                "vtx_symm_class": vtx_symm_class
            }

    @staticmethod
    def groupby_vtx_symm_class(d_input):
        res = {}
        for i, v in d_input.items():
            res[v] = [i] if v not in res.keys() else res[v] + [i]
        return res

    @staticmethod
    def sort_facets(facets):
        sorted_facets = []
        for facet in facets:
            sorted_facets += [tuple(sorted(facet, reverse=False))]
        return tuple(sorted_facets)

    @staticmethod
    def get_created_vids(facets):
        return set(flatten(facets))

    def _compute_state(self, facets) -> tuple:
        deg_per_vid = defaultdict(int)
        vid_per_facet = defaultdict(list)
        for facet in facets:
            for vid in facet:
                deg_per_vid[vid] += 1
                vid_per_facet[vid] += [self.facet2id[tuple(facet)]]
        return deg_per_vid, vid_per_facet

    @staticmethod
    def dict2tuple(d_input):
        d_input_tupled = dict()
        for d in d_input:
            d_input_tupled[d] = tuple(d_input[d])
        return d_input_tupled

    def log_forbidden(self, name, reason_id) -> None:
        self.logbook[tuple(name)] = {
            "is_simplicial": False,
            "reason": reason_id
        }

    def update_incrementally(self, facet, facets):
        _facets = deepcopy(facets)
        self.register_facet(facet)
        _facets += [facet]
        self.register_facets(_facets)
        self.register_state(_facets)


class Enum(EnumRegistrar):
    def __init__(self, size_seq):
        super().__init__()
        self.size_seq = sorted(size_seq, reverse=True)
        self.created_vids = set()
        pass

    def compute(self):
        s = self.size_seq.pop(0)
        facets = []
        facet = []
        for _ in range(s):
            facet += [_]
        self.update_incrementally(facet, facets)

        facets_count = 1
        while len(self.size_seq) > 0:
            next_size = self.size_seq.pop(0)
            pool = self.get_facet_ids_per_n(facets_count)

            for facets_id in pool:
                facets = list(self.id2facets[facets_id])
                for ns_vtx in range(0, next_size + 1):
                    self.created_vids = self.get_created_vids(facets)
                    if ns_vtx == 0:  # must check with hitting set routine
                        self.fill_wo_creating_new_vertices(next_size, facets)
                    elif ns_vtx == next_size:  # only a single state could result
                        self.fill_with_only_ones(next_size, facets)
                    else:  # freedom in choosing slots in the shielded region, must consider identical vertices
                        self.fill_w_creating_new_vertices(next_size, facets, ns_vtx)
            facets_count += 1

    def fill_wo_creating_new_vertices(self, next_size, facets):
        # print(f"Dealing with ns_vtx = {ns_vtx} of {next_size}")
        hs_list = get_hitting_sets(facets, self.created_vids)
        vsc = self.states[self.facets2id[self.sort_facets(facets)]]["vtx_symm_class"]
        vsc_deepcopy = deepcopy(vsc)

        if len(hs_list) == 0:
            pass
        else:
            for _hs_list in hs_list:
                if len(_hs_list) <= next_size:
                    vsc = deepcopy(vsc_deepcopy)
                    _next_size = next_size - len(_hs_list)
                    leftover = list(self.created_vids.difference(set(_hs_list)))
                    _iter = combinations_with_replacement(vsc.keys(), _next_size)
                    for _iter_combn in _iter:
                        facet = deepcopy(_hs_list)
                        tracker = defaultdict(int)
                        for symm_class in _iter_combn:
                            try:
                                facet += [vsc[symm_class][tracker[symm_class]]]
                                tracker[symm_class] += 1
                            except IndexError:
                                break

                        if len(set(facet)) == next_size:
                            self.update_incrementally(facet, facets)

    def fill_w_creating_new_vertices(self, next_size, facets, ns_vtx):
        facet = []
        for _ in range(ns_vtx):
            facet += [_ + len(self.created_vids)]

        facet_deepcopy = deepcopy(facet)

        _next_size = next_size - ns_vtx
        vsc = self.states[self.facets2id[self.sort_facets(facets)]]["vtx_symm_class"]
        vsc_deepcopy = deepcopy(vsc)

        _iter = combinations_with_replacement(vsc.keys(), _next_size)
        for _iter_combn in _iter:
            facet = deepcopy(facet_deepcopy)
            vsc = deepcopy(vsc_deepcopy)
            tracker = defaultdict(int)
            for symm_class in _iter_combn:
                try:
                    facet += [vsc[symm_class][tracker[symm_class]]]
                    tracker[symm_class] += 1
                except IndexError:
                    break
            if len(facet) == next_size:
                self.update_incrementally(facet, facets)

    def fill_with_only_ones(self, next_size, facets):
        facet = []
        for _ in range(next_size):
            facet += [_ + len(self.created_vids)]
        self.update_incrementally(facet, facets)

    def get_facet_ids_per_n(self, n):
        pool = []
        for _ in self.states:
            if self.states[_]["n_facets"] == n:
                pool += [_]
        return pool
