import random
import numpy as np
from collections import Counter


def flatten(nested_list):
    return [item for sublist in nested_list for item in sublist]


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
