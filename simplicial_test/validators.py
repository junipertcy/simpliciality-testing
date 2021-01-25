from simplicial_test.utils import *


def validate_data(sorted_d, sorted_s):
    # n = len(sorted_d)
    m = len(sorted_s)

    if len(sorted_d) > 0 and np.max(sorted_d) > m:
        return False
    if np.sum(sorted_d) != np.sum(sorted_s):
        return False
    if Counter(sorted_s)[1] > Counter(sorted_d)[1]:  # TODO: perhaps we could only enforce this at the very first level.
        return False

    # for ind_d in range(0, len(sorted_d) + 1):
    #     _ = 0
    #     for s in range(2 + ind_d, sorted_s[0] + 1):
    #         _ += Counter(sorted_s)[s]
    #     if _ < sorted_d[ind_d]:
    #         return False
    return True


def simple_validate(degs, sizes, facets, facet, enumeration=False):
    current_facets = facets + [facet]
    wanting_degs, non_shielding = get_remaining_slots(degs, facets, facet)  # wanting_degs (non_shielding & shielding)
    if enumeration:
        if min(wanting_degs) < 0:
            # print("1")
            return False, "added in Jan 14"
    if len(sizes) == 0:
        return True, "Last facet explored."
    elif len(sizes) == 1:
        for _ in current_facets:
            if set(np.nonzero(wanting_degs)[0]).issubset(set(_)):
                # print("2")
                return False, "Second last facet explored. " \
                              "But the very last facet is doomed to fail the no-inclusion constraint."
    if np.any(wanting_degs > len(sizes)):  # useful
        # print("3")
        return False, "Some degrees require more facets than the actual remaining number."

    if np.count_nonzero(wanting_degs) < sizes[0]:  # useful
        return False, "5"
    elif np.count_nonzero(wanting_degs) == sizes[0]:
        if np.sum(wanting_degs) > np.count_nonzero(wanting_degs):
            return False, "4"
    if np.min(non_shielding) >= 0:
        if validate_nonshielding(sizes, wanting_degs, non_shielding):  # useful
            return False, "6"

    return wanting_degs, sizes, current_facets


def validate_nonshielding(curent_sizes, wanting_degs, non_shielding):
    """
    Verify that the sum of non-shielding slots not be less than the number of the remaining facets.

    Parameters
    ----------
    curent_sizes
    non_shielding

    Returns
    -------

    """
    shielding = wanting_degs - non_shielding  # only shielding
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


def validate_reduced_seq(wanting_degs, curent_sizes, current_facets, blocked_sets) -> (bool, tuple):
    collected_facets = list()
    curent_sizes = np.array(curent_sizes, dtype=np.int_)
    exempt_vids = list()
    while Counter(wanting_degs)[len(curent_sizes)] != 0:
        if len(curent_sizes) > 1:
            if not basic_validations_degs_and_sizes(degs=wanting_degs, sizes=curent_sizes):
                return True, tuple()

        must_be_filled_vids = np.where(wanting_degs == len(curent_sizes))[0]
        shielding_facets = get_shielding_facets_when_vids_filled(
            current_facets, blocked_sets, must_be_filled_vids, exempt_vids=exempt_vids
        )
        nonshielding_vids = get_nonshielding_vids(shielding_facets, wanting_degs)
        shielding_facets = [set(_).difference(set(must_be_filled_vids)) for _ in shielding_facets]

        curent_sizes -= Counter(wanting_degs)[len(curent_sizes)]
        if Counter(curent_sizes)[0] == len(curent_sizes):
            collected_facets += [exempt_vids + must_be_filled_vids.tolist()]
            wanting_degs[wanting_degs == len(curent_sizes)] = 0
            break
        wanting_degs[wanting_degs == len(curent_sizes)] = 0

        if Counter(curent_sizes)[1] > Counter(wanting_degs)[1]:
            return True, tuple()
        if Counter(curent_sizes)[1] > 0:
            try:
                if len(shielding_facets) == 0:  # You are free to do "remove_ones"
                    nonshielding_vids = set(np.nonzero(wanting_degs)[0])
                wanting_degs, removed_sites = remove_ones(curent_sizes, wanting_degs, choose_from=nonshielding_vids)
            except NoSlotError:
                return True, tuple()
            else:
                curent_sizes = curent_sizes[curent_sizes != 1]
                to_be_added_facets = [exempt_vids + must_be_filled_vids.tolist() + [s] for s in removed_sites]
                collected_facets += to_be_added_facets
        exempt_vids += must_be_filled_vids.tolist()
    reduced_data = (wanting_degs, curent_sizes, collected_facets, exempt_vids)
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


def is_reduced_seq(wanting_degs, sizes, degrees) -> bool:
    return Counter(np.equal(wanting_degs, degrees))[True] == Counter(wanting_degs)[len(sizes)]


def validate_interm_degs(degs, intermediate_degs):
    return np.all(degs - intermediate_degs >= 0)
