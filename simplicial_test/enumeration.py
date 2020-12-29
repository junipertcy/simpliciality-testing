from simplicial_test.utils import *
from simplicial_test.sample import get_hitting_sets
from copy import deepcopy
from itertools import combinations_with_replacement, starmap, product
from collections import defaultdict


def sort_facets(facets):
    sorted_facets = []
    for facet in facets:
        sorted_facets += [tuple(sorted(facet, reverse=False))]
    sorted_facets = set(sorted_facets)
    return tuple(sorted(sorted_facets, key=lambda _: [-len(_)] + list(_)))


def get_relabeled_facets(facets):
    equiv = defaultdict(list)
    deg_list = compute_dpv(facets, is_sorted=False)

    for idx, _ in enumerate(deg_list):
        equiv[_] += [idx]
    index = 0
    old2new = dict()
    for _ in sorted(equiv.items(), key=lambda kv: kv[0], reverse=True):
        for _v in _[1]:
            if _v not in old2new:
                old2new[_v] = index
                index += 1
    new_facets = []
    for facet in facets:
        _facet = []
        for vtx in facet:
            _facet += [old2new[vtx]]
        new_facets += [tuple(sorted(_facet))]
    return sort_facets(new_facets)


class EnumRegistrar(object):
    def __init__(self):
        self.pointer = 0
        self.facets_pointer = 0
        self.facets_count = 0
        self.facet_count_per_facets = defaultdict(int)
        self.ns_vtx = 0
        self.dfs_locator = None

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
        sorted_facets = sort_facets(facets)
        if sorted_facets not in self.facets2id:
            self.facets2id[sorted_facets] = self.facets_pointer
            self.id2facets[self.facets_pointer] = sorted_facets
            self.facets_pointer += 1
        return sorted_facets, self.facets2id[sorted_facets]

    def register_state(self, facets):
        """
        dpv: `degree per vertex`
        
        vpf: `vertex id per facet`

        vsc: `vertex symmetry class`

        n: `number of facets`

        Parameters
        ----------
        facets

        Returns
        -------

        """
        facets = sort_facets(facets)
        dpv = compute_dpv(facets)
        if dpv not in self.states:
            vpf = self.compute_vpf(facets)
            vsc = groupby_vtx_symm_class(self.dict2tuple(vpf))
            self.states[dpv] = {
                "m": len(facets),
                "fid": self.facets2id[facets],
                "vpf": vpf,
                "vsc": vsc,
                "loc": (self.facets_count, self.facet_count_per_facets[self.facets_count], self.ns_vtx),
                "dfs": self.dfs_locator
            }
            self.facet_count_per_facets[self.facets_count] += 1

    @staticmethod
    def get_created_vids(facets):
        return set(flatten(facets))

    def compute_vpf(self, facets) -> dict:
        vpf = defaultdict(list)
        for facet in facets:
            for vid in facet:
                vpf[vid] += [self.facet2id[tuple(facet)]]
        return vpf

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

    @staticmethod
    def get_dfs_navigator(size_seq):
        m = len(size_seq)
        dummy1 = np.zeros([m], dtype=np.int_)
        dummy1[1] = 1
        dummy2 = [_ + 1 for _ in size_seq]
        dummy2[0] = 1
        s_map = starmap(range, zip(dummy1, dummy2))
        return product(*s_map)

    def compute_dfs(self):
        nav = self.get_dfs_navigator(self.size_seq)
        s = self.size_seq.pop(0)
        facets = []
        facet = []
        for _ in range(s):
            facet += [_]

        self.dfs_locator = tuple([0])
        self.update_incrementally(facet, facets)

        for _nav in nav:
            _nav = list(_nav)
            _nav.pop(0)
            idx = 0
            dfs_locator = [0]
            self.dfs_locator = tuple(dfs_locator)
            while len(_nav) > 0:
                ns_vtx = _nav.pop(0)
                next_size = self.size_seq[idx]
                self.ns_vtx = ns_vtx

                pool = self.get_facet_ids_per_dfs_locator(self.dfs_locator)
                dfs_locator += [ns_vtx]
                self.dfs_locator = tuple(dfs_locator)

                for facets_id in pool:
                    facets = list(self.id2facets[facets_id])
                    self.created_vids = self.get_created_vids(facets)
                    if ns_vtx == 0:
                        self.fill_wo_creating_new_vertices(next_size, facets)
                    elif ns_vtx == next_size:
                        self.fill_with_only_ones(next_size, facets)
                    else:
                        self.fill_w_creating_new_vertices(next_size, facets, ns_vtx)
                idx += 1

    def compute(self):
        s = self.size_seq.pop(0)
        facets = []
        facet = []
        for _ in range(s):
            facet += [_]
        self.update_incrementally(facet, facets)

        while len(self.size_seq) > 0:
            self.facets_count += 1
            next_size = self.size_seq.pop(0)
            pool = self.get_facet_ids_per_m(self.facets_count)

            for facets_id in pool:
                facets = list(self.id2facets[facets_id])
                for ns_vtx in range(0, next_size + 1):
                    self.ns_vtx = ns_vtx
                    self.created_vids = self.get_created_vids(facets)
                    if ns_vtx == 0:  # must check with hitting set routine
                        self.fill_wo_creating_new_vertices(next_size, facets)
                    elif ns_vtx == next_size:  # only a single state could result
                        self.fill_with_only_ones(next_size, facets)
                    else:  # freedom in choosing slots in the shielded region, must consider identical vertices
                        self.fill_w_creating_new_vertices(next_size, facets, ns_vtx)

    @staticmethod
    def get_hs_identifier(vsc, hs):
        vsc_inv = {tuple(v): k for k, v in vsc.items()}
        _id = []
        for _hs in hs:
            for k, v in vsc_inv.items():
                if _hs in k:
                    _id += [v[0]]
        return tuple(sorted(_id))

    def fill_wo_creating_new_vertices(self, next_size, facets):
        # print(f"Dealing with ns_vtx = {ns_vtx} of {next_size}")
        hs_list = get_hitting_sets(facets, self.created_vids)
        vsc = self.states[compute_dpv(facets)]["vsc"]
        vsc_deepcopy = deepcopy(vsc)

        if len(hs_list) == 0:
            pass
        else:
            hs_ids = dict()
            for _hs_list in hs_list:
                try:
                    hs_ids[self.get_hs_identifier(vsc, _hs_list)]
                except KeyError:
                    hs_ids[self.get_hs_identifier(vsc, _hs_list)] = _hs_list
                else:
                    if sum(_hs_list) < sum(hs_ids[self.get_hs_identifier(vsc, _hs_list)]):
                        hs_ids[self.get_hs_identifier(vsc, _hs_list)] = _hs_list
            for _hs in hs_ids.values():
                if len(_hs) <= next_size:
                    vsc = deepcopy(vsc_deepcopy)
                    _next_size = next_size - len(_hs)
                    _iter = combinations_with_replacement([_ for _ in vsc.keys() if _ not in _hs], _next_size)
                    for _iter_combn in _iter:
                        facet = deepcopy(_hs)
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
        vsc = self.states[compute_dpv(facets)]["vsc"]
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

    def get_facet_ids_per_m(self, m):
        pool = []
        for _ in self.states:
            if self.states[_]["m"] == m:
                pool += [_]
        return pool

    def get_dpv_ids_per_m(self, m):
        pool = []
        for _ in self.states:
            if self.states[_]["m"] == m:
                pool += [_]
        return pool

    def get_facet_ids_per_dfs_locator(self, dfs):
        pool = []
        for _ in self.states:
            if self.states[_]["dfs"] == dfs:
                pool += [self.states[_]["fid"]]
        return pool

    def get_placed_order_per_fid(self, _id):
        flen = len(self.id2facets[_id])
        codes = []
        for c in range(flen):
            fid = tuple(self.id2facets[_id][:c + 1])
            codes += [self.states[self.facets2id[fid]]["loc"][2]]
        return codes
