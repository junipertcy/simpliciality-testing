import numpy as np
from collections import defaultdict, Counter
from copy import deepcopy
import random


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


def transform_facets(facets, inv_map):
    sets = []
    for facet in facets:
        transformed_facet = []
        for _ in facet:
            try:
                transformed_facet += [inv_map[_]]
            except KeyError:
                pass
        sets += [transformed_facet]
    return sets


def get_shielding_facets_when_vids_filled(current_facets, must_be_filled_vids, exempt_vids=list()) -> (bool, set, list):
    """
    TODO: this function can be further simplified, along with the function::validate_reduced_seq
    The function works when one have "must_be_filled_vids" -- it goes by searching already existing facets,
    And find out the slots that must not be chosen in order to avoid clashes.

    Parameters
    ----------
    both
    curent_sizes
    current_facets

    Returns
    -------

    """
    # print("-> must_be_filled_vids are ", must_be_filled_vids)
    shielding_facets = list()
    for facet in current_facets:  # for all existing facets
        if set(must_be_filled_vids).issubset(set(facet)) and set(exempt_vids).issubset(set(facet)):  # if a facet contains these must_be_filled_vids
            shielding_facets += [facet]
    return shielding_facets  # then we must avoid the slots in these shielding_facets


def get_nonshielding_vids(shielding_facets, both):
    nonshielding_vids = set()
    n = len(both)
    # print("-> shielding_facets = ", shielding_facets)
    for facet in shielding_facets:
        non_shielding_part = set(range(n)).difference(set(facet))  # non_shielding_part vertex_id's
        if len(nonshielding_vids) == 0:
            nonshielding_vids = non_shielding_part
        else:
            nonshielding_vids.intersection_update(non_shielding_part)
        # remaining number of facets must be able to "hide" in those non_shielding slots, at least
        # if remaining > np.sum(both[np.array(list(nonshielding_vids), dtype=np.int_)]):
        #     print(f"While checking if the remaining number of facets being able to hide in the nonshielding slots, ")
        #     print(f"We found that remaining={remaining} but nonshielding_vids={nonshielding_vids} (i.e., available sites to hide = {np.sum(both[np.array(list(nonshielding_vids), dtype=np.int_)])})")
        #     print("Therefore, we reject that the nodes being filled this way.")
        #     raise NoSlotError("The custom error in get_nshielding.")
    return nonshielding_vids


def basic_validations_degs_and_sizes(degs, sizes):
    if Counter(degs)[len(sizes)] == np.min(sizes):
        return False
    if len(degs) == np.max(sizes):
        return False
    return True


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


def remove_ones(s, both, choose_from=set()):
    # orig_both = deepcopy(both)
    # _both = deepcopy(both)
    # both = [_ for _ in both]
    # for _ in range(Counter(s)[1]):
    #     both.remove(1)
    # both = np.array(both)
    # print(type(choose_from), f"np.where(both == 1)[0] = {np.where(both == 1)[0]}")
    if len(choose_from) == 0:
        raise NoSlotError
    all_one_sites = np.where(both == 1)[0]
    if choose_from is not None and type(choose_from) is set:
        choose_from.intersection_update(set(all_one_sites))
        if not choose_from:
            raise NoSlotError
    # print(f"choose_from is {choose_from}")
    removed_sites = np.array(list(choose_from))[:Counter(s)[1]]  # todo, to figure out: you cannot do [-Counter(s)[1]:]
    both[removed_sites] = 0
    return both, removed_sites.tolist()


def trim_ones(size_list, degree_list) -> (list, list):
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


def match_ones(sorted_s, sorted_d):
    _s = Counter(sorted_s)[1]
    _d = Counter(sorted_d)[1]
    matching_isolated = 0
    if _s > _d:
        print("Too many 0-simplices. We cannot satisfy the inclusive constraint.")
        return False, list(), list(), 0
    else:
        for _ in range(_s):
            sorted_s = np.delete(sorted_s, -1)
            sorted_d = np.delete(sorted_d, -1)
            matching_isolated += 1
    return True, sorted_s, sorted_d, matching_isolated


def get_reduced_seq(degs, facets) -> list:
    mapping = get_reduced_seq_mapping(degs)
    return list(map(lambda x: tuple(map(lambda y: mapping[y], x)), facets))


def get_reduced_seq_mapping(degs):
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

    mapping = dict()
    for key in inv_old2new.keys():
        for idx, new_key in enumerate(inv_old2new[key]):
            mapping[new_key] = _inv_old2new[key][idx]
    return mapping
