from simplicial_test.utils import *


def validate_data(sorted_d, sorted_s):
    m = len(sorted_s)
    if len(sorted_d) > 0 and np.max(sorted_d) > m:
        return False
    if np.sum(sorted_d) != np.sum(sorted_s):
        return False
    if Counter(sorted_s)[1] > Counter(sorted_d)[1]:  # TODO: perhaps we could only enforce this at the very first level.
        return False
    return True


def simple_validate(degs, sizes, facets, facet, enumeration=False):
    current_facets = facets + [facet]
    wanting_degs, non_shielding = get_remaining_slots(degs, facets, facet)  # wanting_degs (non_shielding & shielding)
    if enumeration:
        if min(wanting_degs) < 0:
            return False, "added in Jan 14"
    if len(sizes) == 0:
        return True, "Last facet explored."
    if len(sizes) == 1:
        for _ in current_facets:
            if set(np.nonzero(wanting_degs)[0]).issubset(set(_)):
                return False, "Second last facet explored. " \
                              "But the very last facet is doomed to fail the no-inclusion constraint."
    if np.any(wanting_degs > len(sizes)):  # useful
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


def validate_reduced_seq(wanting_degs, sizes, current_facets, blocked_sets, verbose=False) -> (bool, tuple):
    if verbose:
        print(f"----- (BEGIN: reducing seq) -----\n"
              f"Wanting degrees: {wanting_degs}\n"
              f"Size list: {sizes}\n"
              f"Current facets: {current_facets}\n"
              f"blocked_sets = {blocked_sets}"
              )
    collected_facets = []
    exempt_vids = []
    sizes = np.array(sizes, dtype=np.int_)
    while Counter(wanting_degs)[len(sizes)] != 0:
        if len(sizes) > 1:
            if not basic_validations_degs_and_sizes(degs=wanting_degs, sizes=sizes):
                return True, tuple()
        if Counter(sizes)[0] != len(sizes):
            must_be_filled_vids = np.where(wanting_degs == len(sizes))[0].tolist()
            exempt_vids += must_be_filled_vids
        else:
            break

        sizes -= Counter(wanting_degs)[len(sizes)]
        wanting_degs[wanting_degs == len(sizes)] = 0
        if min(sizes) < 0:
            return True, tuple()
        if Counter(sizes)[1] > 0:
            if Counter(sizes)[1] > Counter(wanting_degs)[1]:
                return True, tuple()
            shielding_facets = get_shielding_facets_when_vids_filled(
                current_facets, blocked_sets, must_be_filled_vids, exempt_vids=exempt_vids
            )
            shielding_facets = [set(_).difference(set(must_be_filled_vids)) for _ in shielding_facets]
            if len(shielding_facets) == 0:  # You are free to do "remove_ones" (at a later stage, perhaps)
                nonshielding_vids = set(np.nonzero(wanting_degs)[0])
            else:
                nonshielding_vids = get_nonshielding_vids(len(wanting_degs), shielding_facets)
            all_one_sites = np.where(wanting_degs == 1)[0]
            nonshielding_vids.intersection_update(set(all_one_sites))
            if not nonshielding_vids:
                return True, tuple()
            wanting_degs, removed_sites = remove_ones(sizes, wanting_degs, choose_from=nonshielding_vids)
            sizes = sizes[sizes != 1]
            collected_facets += [exempt_vids + [_] for _ in removed_sites]  # collected_facets immer from removed_sites

    if np.sum(wanting_degs) == np.sum(sizes) == 0:
        if validate_issubset_blocked_sets(exempt_vids, blocked_sets=blocked_sets):
            return True, tuple()
        # for facet in collected_facets:
        #     if validate_issubset_blocked_sets(facet, blocked_sets=blocked_sets):
        #         return True, tuple()
    if verbose:
        print(f"----- (↓ Returning these data ↓) -----\n"
              f"Wanting degrees: {wanting_degs}\n"
              f"Size list: {sizes}\n"
              f"Collected facets: {collected_facets}\n"
              f"Exempt vertex ids: {exempt_vids}\n"
              f"----- (END: reducing seq) -----"
              )
    return False, (wanting_degs, sizes, collected_facets, exempt_vids)


def validate_issubset_blocked_sets(candidate_facet, blocked_sets=None):
    if blocked_sets is not None:
        for facet in blocked_sets:
            if set(candidate_facet).issubset(set(facet)):
                return True
    return False
