#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# simplicial_test -- a python module to realize simplicial joint sequences
#
# Copyright (C) 2020-2021 Tzu-Chi Yen <tzuchi.yen@colorado.edu>
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3 of the License, or (at your option) any
# later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.

"""
The `cm_count` module computes the number of networks in the configuration model.


References
----------
.. [courtney-generalized-2016] Anita Liebenau and Nick Wormald, "Asymptotic enumeration
   of digraphs and bipartite graphs by degree sequence",
   :arxiv:`2006.15797`

"""

import numpy as np
from scipy.special import gammaln


def compute_entropy(degs, sizes):
    m = np.sum(degs)

    l = len(sizes)
    sigma_s = np.var(sizes)

    n = len(degs)
    sigma_t = np.var(degs)
    mu = m / (2 * l * n)

    term_1 = - lbinom(l * n, m)
    term_2 = 0.
    for _ in sizes:
        term_2 += lbinom(n, _)

    term_3 = 0.
    for _ in degs:
        term_3 += lbinom(l, _)

    term_4 = - 1 / 2 * (1 - sigma_s ** 2 / mu / (1 - mu) / n) * (1 - sigma_t ** 2 / mu / (1 - mu) / l)
    return term_1 + term_2 + term_3 + term_4


def lbinom(n, k):
    """Return log of binom(n, k)."""
    if type(n) in [float, int, np.int64, np.float64]:
        n = np.array([n])
        k = np.array([k])
    return (gammaln(np.array([float(x) for x in n + 1])) -
            gammaln(np.array([float(x) for x in n - k + 1])) -
            gammaln(np.array([float(x) for x in k + 1])))
