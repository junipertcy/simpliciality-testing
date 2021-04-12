How to Use Simplicial-test
==========================

Base on the heuristic, our first guess is having `[a, c, d]` as the candidate facet.
However, this configuration makes the non-shielding vertices being only `b` and `e`,
whose degrees sum to 2. This is implausible because each of the remaining facets has
to put at least one vertex in this non-shielding region, making at least one facet
doomed to be a subset of `[a, c, d]`.

.. list-table::
    Can deg_seq={3, 2, 2, 1, 1} and size_seq={3, 2, 2, 2} make
    a simplicial instance? Here's our first attempt, choosing [a, c, d] as the
    first facet. Note that stars (*) represent the non-shielding region.
    :widths: 5 10 5 5 5 5 5
    :stub-columns: 1
    :header-rows: 1

    * -
      -
      - 3
      - 2
      - 2
      - 1
      - 1
    * -
      -
      - c
      - a
      - d
      - b
      - e
    * - 3
      - [a, c, d]
      - 1
      - 1
      - 1
      -
      -
    * - 2
      - [?, ?]
      -
      -
      - *
      - *
      - *
    * - 2
      - [?, ?]
      -
      -
      - *
      - *
      - *
    * - 2
      - [?, ?]
      -
      -
      - *
      - *
      - *