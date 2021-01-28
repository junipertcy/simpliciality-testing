import numpy as np
from collections import defaultdict, Counter
from copy import deepcopy

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable


def compute_level_map(level_map, level, mapping2shrinked):
    n = len(level_map[-1])
    if level > 1:
        for _ in range(n):
            vtx_current_view = level_map[level - 1][_]
            if vtx_current_view == -1:
                level_map[level][_] = -1
                continue
            try:
                level_map[level][_] = mapping2shrinked[vtx_current_view]
            except KeyError:
                level_map[level][_] = -1
    else:
        for _ in range(n):
            try:
                level_map[1][_] = mapping2shrinked[_]
            except KeyError:
                level_map[1][_] = -1
    return level_map


def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]


def simplify_blocked_sets(bsets):
    data = []
    for _bsets in bsets[::-1]:
        if not tuple(_bsets) in data:
            if len(data) > 0:
                to_add = True
                for _ in data:
                    if set(_bsets).issubset(set(_)):
                        to_add = False
                        break
                if to_add:
                    data += [_bsets]
            else:
                data += [_bsets]
    return data


def compute_dpv(facets, n=None, is_sorted=True):
    if n is None:
        dpv = defaultdict(int)
    else:
        dpv = np.zeros([n], dtype=np.int_)

    for facet in facets:
        for vid in facet:
            dpv[vid] += 1
    if n is not None:
        return dpv

    if is_sorted:
        return tuple(sorted(list(dpv.values()), reverse=True))
    else:
        _dpv = []
        for _ in range(len(dpv.keys())):
            _dpv += [dpv[_]]
        return tuple(_dpv)


def get_remaining_slots(degs, facets, facet):
    """
    Used in the greedy case only.

    Parameters
    ----------
    facet: candidate facet

    Returns
    -------

    """
    n = len(degs)
    remaining = degs - compute_dpv(facets, n=n, is_sorted=False)
    return np.array([remaining[_] - 1 if _ in set(facet) else remaining[_] for _ in range(n)], dtype=np.int_), \
           np.array([0 if _ in set(facet) else remaining[_] for _ in range(n)], dtype=np.int_)  # non_shielding


def groupby_vtx_symm_class(d_input):
    res = {}
    for i, v in d_input.items():
        res[v] = [i] if v not in res.keys() else res[v] + [i]
    return res


def get_indices_of_k_in_blocked_sets(blocked_sets, k):
    indices = []
    for idx, _ in enumerate(blocked_sets):
        if k in _:
            indices += [idx]
    return indices


class NoSlotError(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f'NoSlotError {self.message}'
        else:
            return 'NoSlotError has been raised'


class WeCanStopSignal(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f'WeCanStopSignal {self.message}'
        else:
            return 'WeCanStopSignal has been raised'


class NoMoreBalls(Exception):
    def __init__(self, *args):
        if args:
            self.message = args[0]
        else:
            self.message = None

    def __str__(self):
        if self.message:
            return f'{self.message}'
        else:
            return 'NoMoreBalls has been raised'


class SimplexRegistrar(object):
    def __init__(self):
        self.pointer = 0
        self.name2id = dict()
        self.id2name = dict()
        self.facet_size_per_id = np.array([], dtype=np.int_)
        self.logbook = dict()

    def register(self, name) -> (tuple, int):
        name = tuple(sorted(name, reverse=True))
        if name not in self.name2id:
            self.name2id[name] = self.pointer
            self.id2name[self.pointer] = name
            self.pointer += 1
            self.facet_size_per_id = np.append(self.facet_size_per_id, [len(name)])
        return name, self.name2id[name]

    def log_forbidden(self, name, reason_id) -> None:
        self.logbook[tuple(name)] = {
            "is_simplicial": False,
            "reason": reason_id
        }


def update_deg_seq(deg_seq, facet, value):
    if value not in [+1, -1]:
        raise NotImplementedError
    for _ in facet:
        deg_seq[_] += value


def shrink_facets(facets, inv_map):
    _facets = []
    for facet in facets:
        transformed_facet = []
        for _ in facet:
            try:
                transformed_facet += [inv_map[_]]
            except KeyError:
                pass
        if len(sorted(transformed_facet)) > 0:
            _facets += [sorted(transformed_facet)]
    _facets.sort(key=len)
    return _facets


def shrink_degs(degs, inv_map):
    d = dict()
    for idx, _degs in enumerate(degs):
        try:
            d[inv_map[idx]] = _degs
        except KeyError:
            pass
    d_list = list()
    for idx, _d in enumerate(d.values()):
        d_list += [d[idx]]
    return np.array(d_list, dtype=np.int_)


def get_shielding_facets_when_vids_filled(current_facets, blocked_sets, must_be_filled_vids, exempt_vids=None):
    """
    TODO: this function can be further simplified, along with the function::validate_reduced_seq
    The function works when one have "must_be_filled_vids" -- it goes by searching already existing facets,
    And find out the slots that must not be chosen in order to avoid clashes.

    Parameters
    ----------
    wanting_degs
    curent_sizes
    current_facets
    exempt_vids

    Returns
    -------

    """
    if exempt_vids is None:
        exempt_vids = []
    shielding_facets = []
    # if a facet contains these must_be_filled_vids (or 'mbfv')
    mbfv = set(must_be_filled_vids).union(set(exempt_vids))
    for facet in current_facets:  # for all existing facets
        if mbfv.issubset(set(facet)):
            shielding_facets += [facet]
    for facet in blocked_sets:  # for all existing facets
        if mbfv.issubset(set(facet)):
            shielding_facets += [facet]
    return shielding_facets  # then we must avoid the slots in these shielding_facets


def get_nonshielding_vids(n, shielding_facets):
    nonshielding_vids = set()
    for facet in shielding_facets:
        non_shielding_part = set(range(n)).difference(set(facet))  # non_shielding_part vertex_id's
        if len(nonshielding_vids) == 0:
            nonshielding_vids = non_shielding_part
        else:
            nonshielding_vids.intersection_update(non_shielding_part)
    return nonshielding_vids


def basic_validations_degs_and_sizes(degs, sizes):
    if Counter(degs)[len(sizes)] == np.min(sizes):
        return False
    if len(degs) == np.max(sizes):
        return False
    return True


def remove_ones(sizes, wanting_degs, choose_from=None):
    removed_vtx_sites = np.array(list(choose_from), dtype=np.int_)[:Counter(sizes)[1]]  # todo, to figure out: you cannot do [-Counter(s)[1]:]
    wanting_degs[removed_vtx_sites] = 0
    return wanting_degs, removed_vtx_sites.tolist()


def pair_one_by_one(size_list, degree_list) -> (list, list):
    """
    This function is used to pair up the ones in size/deg lists, since this is the only plausible situation.
    Note that it is only applied when level=0.

    Parameters
    ----------
    size_list
    degree_list

    Returns
    -------

    """
    size_list = list(size_list)
    degree_list = list(degree_list)
    _1 = Counter(size_list)[1]
    _2 = Counter(degree_list)[1]
    _ = min(_1, _2)
    for __ in range(_):
        size_list.remove(1)
        degree_list.remove(1)
    return size_list, degree_list, _


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


def sort_helper(st) -> list:
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


def get_mapped_seq(mapping2enlarged, facets) -> list:
    # _, mapping = get_seq2seq_mapping(degs)
    reduced_facets = list(map(lambda x: tuple(map(lambda y: mapping2enlarged[y], x)), facets))
    # print(f"The facets above have been `enlarged` (i.e., higher-level) to: {reduced_facets}")
    return reduced_facets


def sort_callback(facets):
    t = []
    for facet in facets:
        t += [tuple(sorted(facet, reverse=False))]
    return t


def get_seq2seq_mapping(degs):
    old2new = dict()
    for idx, _deg in enumerate(degs):
        old2new[idx] = _deg
    inv_old2new = defaultdict(list)
    for key in old2new.keys():
        if old2new[key] != 0:
            inv_old2new[old2new[key]] += [key]
    _idx = 0
    _inv_old2new = deepcopy(inv_old2new)

    keys = sorted(inv_old2new.keys(), reverse=True)
    for key in keys:
        for idx, _ in enumerate(inv_old2new[key]):
            inv_old2new[key][idx] = _idx
            _idx += 1

    mapping2enlarged = dict()
    for key in inv_old2new.keys():
        for idx, new_key in enumerate(inv_old2new[key]):
            mapping2enlarged[new_key] = _inv_old2new[key][idx]

    mapping2shrinked = {v: k for k, v in mapping2enlarged.items()}
    return mapping2shrinked, mapping2enlarged


def compute_joint_seq(facets) -> (list, list):
    # print(facets)
    flattened_facets = [item for sublist in facets for item in sublist]
    n = max(flattened_facets) - min(flattened_facets) + 1
    degs = np.zeros(n, dtype=np.int_)
    sizes = []

    for facet in facets:
        sizes += [len(facet)]
        for vertex_id in facet:
            degs[vertex_id] += 1

    return sorted(sizes, reverse=True), sorted(degs, reverse=True)


def if_facets_simplicial(facets) -> bool:
    for idx1, facet1 in enumerate(sorted(facets, key=len)):
        for idx2, facet2 in enumerate(sorted(facets, key=len)):
            if idx2 > idx1:
                if set(list(facet1)).issubset(set(list(facet2))):
                    return False
    return True


def filter_blocked_facets(blocked_facets, exempt_vids):
    """

    Supposedly every existing facet will block potential facets in the next level, however, this only applies if
    it contains `exempt_vids` because these vertices will be shared by next-level candidates.

    Parameters
    ----------
    blocked_facets
    exempt_vids

    Returns
    -------

    """
    # print(f"in filter_blocked_facets (blocked_facets, exempt_vids) = {(blocked_facets, exempt_vids)}")
    filtered = []
    for facet in blocked_facets:
        if set(exempt_vids).issubset(facet):
            filtered += [facet]
    return filtered


def compute_joint_seq_from_facets(facets):
    n = max(flatten(facets)) + 1
    degs = np.zeros(n, dtype=np.int_)
    sizes = np.zeros(len(facets), dtype=np.int_)

    for idx, facet in enumerate(facets):
        for vid in facet:
            degs[vid] += 1
        sizes[idx] = len(facet)

    return degs, sizes


def accel_asc(n):
    """from: http://jeromekelleher.net/generating-integer-partitions.html"""
    a = [0 for i in range(n + 1)]
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
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


def paint_block(facets, output=None, dpi=300, ):
    g = np.zeros([len(facets), max(flatten(facets)) - min(flatten(facets)) + 1])
    for idx, facet in enumerate(facets):
        for vid in facet:
            g[idx][vid] = 1

    plt.figure(figsize=(16, 12), dpi=300)

    plt.imshow(g, cmap='Greys', interpolation='nearest')
    plt.xlabel("Vertex index")
    plt.ylabel("Facet index")

    if output is not None:
        plt.savefig(output, dpi=dpi, transparent=True)


def paint_landscape(mat, output=None, dpi=300, ):
    max_sizes_ind = max_degs_ind = int(mat.size ** 0.5)
    fig, ax = plt.subplots(1, 1, figsize=(12, 12))  # setup the plot

    colors_undersea = plt.cm.tab20c(np.linspace(0, 0.1999, 256))
    colors_land = plt.cm.Wistia(np.linspace(0., 1., 256))
    all_colors = np.vstack((colors_undersea, colors_land))
    cmap = mpl.colors.LinearSegmentedColormap.from_list('noname', all_colors)

    divnorm = colors.TwoSlopeNorm(vmin=min(mat.flatten()) - 0.5, vcenter=0, vmax=max(mat.flatten()) + 0.5)

    ims = ax.imshow(mat, norm=divnorm, cmap=cmap, origin='lower', extent=[1, max_sizes_ind, 1, max_degs_ind],
                    rasterized=True)
    plt.xlabel("Size sequence")
    plt.ylabel("Degree sequence")

    # scaled colorbar that aligns with the frame
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    fig.colorbar(ims, cax=cax, label="Failed attempts", shrink=0.6)
    ax.tick_params(axis="y", direction="in", length=8)
    ax.tick_params(axis="x", direction="in", length=8)
    ax.xaxis.set_ticks_position("bottom")

    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)

    # ax.set_xticks(np.arange(1.5, max_sizes_ind + 1, 5))
    # ax.set_yticks(np.arange(1.5, max_degs_ind + 1, 5))
    # ax.set_xticklabels(np.arange(0, max_sizes_ind + 1, 5))
    # ax.set_yticklabels(np.arange(0, max_degs_ind + 1, 5))
    if output is not None:
        plt.savefig(output, dpi=dpi, transparent=True)


def get_partition(n=1):
    gen = accel_asc(n)
    partitions = []
    while True:
        try:
            partitions += [sorted(next(gen), reverse=True)]
        except StopIteration:
            break
    return sorted(partitions, key=lambda x: [Counter(x)[1], max(x) - len(x)] + list(x), reverse=True)
