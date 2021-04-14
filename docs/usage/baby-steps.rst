How to Use Simplicial-test
==========================

[young-construction-2017]_

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


.. tip::

   If you need some actual joint degree sequences to start with
   (and pretend that you do not know if they are realizable),
   Prof. Benson has `a nice collection`_ of hypergraph datasets.
   To use them, download a dataset and apply these utility functions in order.

   #. :func:`simplicial_test.utils.read_hyperedge_list`
   #. :func:`simplicial_test.utils.prune_included_facets`
   #. :func:`simplicial_test.utils.compute_joint_seq`
   #. Now you have the :code:`degree_list` and :code:`size_list` to test!

.. _`a nice collection`: https://www.cs.cornell.edu/~arb/data/



----

.. [young-construction-2017] Jean-Gabriel Young, Giovanni Petri, Francesco Vaccarino, and Alice Patania,
   "Construction of and efficient sampling from the simplicial configuration model", Phys.
   Rev. E 96, 032312 (2017), :doi:`10.1103/PhysRevE.96.032312`, :arxiv:`1705.10298`.