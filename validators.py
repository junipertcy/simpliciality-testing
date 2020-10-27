import numpy as np
from collections import Counter
from utils import NoSlotError, basic_validations_degs_and_sizes, get_shielding_facets_when_vids_filled, \
    get_nonshielding_vids, remove_ones


def validate_data(sorted_d, sorted_s):
    n = len(sorted_d)
    m = len(sorted_s)
    if len(sorted_d) == len(sorted_s) == 0:
        return True
    if len(sorted_d) > 0 and np.max(sorted_d) > m:
        # print("1. This can never be simplicial.")  # TODO.... why??
        return False
    if np.max(sorted_s) >= n:
        if len(sorted_s) == 1 and np.max(sorted_s) == n:
            return True
        else:
            # print("2. This can not be simplicial.")
            return False
    if np.sum(sorted_d) != np.sum(sorted_s):
        # print("Failing the Gale–Ryser criterion (1957), the sequence is not bigraphic.")
        return False
    # TODO: there is a second part of the GR criterion, which is not coded yet.


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


def validate_reduced_seq(both, curent_sizes, current_facets) -> (bool, tuple):
    collected_facets = list()
    curent_sizes = np.array(curent_sizes, dtype=np.int_)
    shielding_facets = []
    exempt_vids = list()
    while Counter(both)[len(curent_sizes)] != 0:
        if len(curent_sizes) > 1:
            if not basic_validations_degs_and_sizes(degs=both, sizes=curent_sizes):
                return True, tuple()

        must_be_filled_vids = np.where(both == len(curent_sizes))[0]

        shielding_facets = get_shielding_facets_when_vids_filled(
            current_facets, must_be_filled_vids, exempt_vids=exempt_vids
        )
        nonshielding_vids = get_nonshielding_vids(shielding_facets, both)
        shielding_facets = [set(_).difference(set(must_be_filled_vids)) for _ in shielding_facets]

        curent_sizes -= Counter(both)[len(curent_sizes)]
        if Counter(curent_sizes)[0] == len(curent_sizes):
            collected_facets += [exempt_vids + must_be_filled_vids.tolist()]
            both[both == len(curent_sizes)] = 0
            break
        both[both == len(curent_sizes)] = 0

        if Counter(curent_sizes)[1] > Counter(both)[1]:
            return True, tuple()
        if Counter(curent_sizes)[1] > 0:
            try:
                if len(shielding_facets) == 0:  # You are free to do "remove_ones"
                    nonshielding_vids = set(np.nonzero(both)[0])
                both, removed_sites = remove_ones(curent_sizes, both, choose_from=nonshielding_vids)
            except NoSlotError:
                return True, tuple()
            else:
                curent_sizes = curent_sizes[curent_sizes != 1]
                to_be_added_facets = [exempt_vids + must_be_filled_vids.tolist() + [s] for s in removed_sites]
                collected_facets += to_be_added_facets
        exempt_vids += must_be_filled_vids.tolist()
    reduced_data = (both, curent_sizes, collected_facets, exempt_vids)
    return False, reduced_data


def validate_issubset(id2name, candidate_facet, _id):
    """
    Verify that the selected facet not be a subset of any higher facet.

    Parameters
    ----------
    candidate_facet
    _id
    blocked_sets

    Returns
    -------

    """
    if _id is not None:
        if set(candidate_facet).issubset(set(id2name[_id])):
            return True
    return False


def validate_issubset_larger_simplices(candidate_facet, larger_simplices=None):
    if larger_simplices is not None:
        for simplex in larger_simplices:
            if set(candidate_facet).issubset(set(simplex)):
                return True
    return False


def validate_issubset_blocked_sets(candidate_facet, blocked_sets=None):
    # if candidate_facet == (0, 1, 2, 3, 4, 5, 6):
    #     verbose = True
    # else:
    #     verbose = False
    if blocked_sets is not None:
        for facet in blocked_sets:
            # if verbose:
            # print(candidate_facet, facet, set(candidate_facet).issubset(set(facet)))
            if set(candidate_facet).issubset(set(facet)):
                return True
    return False


def check_reduced_facets(facets, current_facets):
    """

    Parameters
    ----------
    facets
    current_facets

    Returns
    -------

    """
    for facet in facets:
        for larger_facets in current_facets:
            if set(facet).issubset(set(larger_facets)):
                return False, "8"
    return True, "0"


def is_reduced_seq(both, sizes, degrees) -> bool:
    return Counter(np.equal(both, degrees))[True] == Counter(both)[len(sizes)]


def validate_interm_degs(degs, intermediate_degs):
    return np.all(degs - intermediate_degs >= 0)
