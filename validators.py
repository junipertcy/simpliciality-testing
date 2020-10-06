import numpy as np
from collections import Counter


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


def validate_issubset(id2name, candidate_facet, _id, blocked_sets=None):
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
    if blocked_sets is not None:
        for facet in blocked_sets:
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


def checkpoint_1(degs, intermediate_degs):
    return np.all(np.sort(degs) - np.sort(intermediate_degs) >= 0)
