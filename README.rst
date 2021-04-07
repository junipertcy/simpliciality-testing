simpliciality-test
==================

.. image:: https://img.shields.io/badge/python-3.8-blue.svg?style=flat
   :target: https://github.com/junipertcy/simpliciality_test/blob/master/COPYING
   :alt: python
.. image:: https://img.shields.io/badge/license-LGPL-green.svg?style=flat
   :target: https://github.com/junipertcy/simpliciality_test/blob/master/COPYING
   :alt: license


A constructive greedy algorithm to check whether a joint sequence of degree/size distributions is simplicial.

Hi. Nice to meet you!
---------------------
The library requires `SageMath`_ and a few PyPI libraries:

..

    pip install -r requirements.txt

To warm up, try:

..

    python is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt

Or, you may try with a bunch of unit tests (see `test/`_):

..

   pytest

This is a deterministic, backtracking-based search algorithm for solving the *simplicial realization problem*.
If your input joint sequence is not realizable (as a simplicial complex), often we would need
to traverse the entire search tree, which would take a huge amount of time!
Happily, more than 90% of the input joint sequences lies in the *polynomial regime*,
which means that they can be solved easily.

For example, you can try the `crime network dataset`_ from this inspiring `Phys. Rev. E paper`_,
which has 551 nodes and 194 facets.

For example,

..

    python is_simplicial.py -k datasets/crime_degs.txt -s datasets/crime_sizes.txt

Boom! It's rather fast, isn't it?

Related experiments
-------------------
TODO.

The three figures below show the solution found by the greedy algorithm.

Those black regions denote 1's, whereas white blocks denote 0's.
Theoretically, at a deeper level, we need to avoid more "blocking facets," but as it turned out,
we barely need to avoid any blocking facets at later recursive levels, just because of the greedy construction.

.. figure:: figures/first_example.png
   :width: 400
   :alt: small first example

   Fig. 1. Small first example.

   We should note that vertices with higher degree will be paired with larger facets, unless such pairing would prohibit
   a solution. Note also that not any facet is in the subset of another (larger) facet.

The following two datasets come from the inspirational `Phys. Rev. E paper`_. Each of these requires less than 1 second
to compute!

.. figure:: figures/crime.png
   :width: 400
   :alt: map to buried treasure

   Fig. 2. Crimes (nodes, n=829) and suspects, victims, and witnesses (facets, f=378) network in St. Louis.

.. figure:: figures/pollinator.png
   :width: 400
   :alt: map to buried treasure

   Fig. 3. Flower-visiting insects (nodes, n=679) and plants (facets f=57) network in Kyoto.


.. list-table:: Empirical networks used in this work.
   :widths: 20 10 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dataset
     - # of edges, E
     - # of facets, \|s\|
     - # of nodes, \|d\|
     - (s_max, s_min)
     - (d_max, d_min)
     - conv. time
     - conv. time per facet
     - f(k)
   * - Macaque brain dataset
     - 628
     - n/a
     - n/a
     - n/a
     - n/a
     - n/a
     - n/a
     - n/a
   * - contact-primary-school
     - 20615
     - 8010
     - 242
     - (5, 2)
     - (174, 20)
     - 115934
     - 14.5
     - n/a
   * - contact-high-school
     - 11770
     - 4862
     - 327
     - (5, 2)
     - (90, 2)
     - 3
     - 0.00064
     - n/a
   * - senate-committees
     - 5165
     - 275
     - 282
     - (31, 4)
     - (59, 1)
     - 0
     - 0
     - n/a
   * - house-committees
     - 11660
     - 302
     - 1290
     - (82, 3)
     - (43, 1)
     - 0
     - 0
     - n/a
   * - mathoverflow-answers
     - 131406
     - 5296
     - 73851
     - (1784, 2)
     - (169, 1)
     - 0
     - 0
     - n/a
   * - senate-bills
     - 92876
     - 3599
     - 294
     - (99, 3)
     - (1147, 1)
     - too slow
     - too slow
     - n/a
   * - house-bills
     - 927075
     - 23267
     - 1494
     - (399, 2)
     - (3824, 1)
     - too slow
     - too slow
     - n/a
   * - walmart-trips
     - 447347
     - 63687
     - 88860
     - (25, 2)
     - (5412, 1)
     - too slow
     - too slow
     - n/a


Installation
------------
[final goal] The package can be installed with pip.

::

   pip install simplicial-test

MISC notes (to clean up later)
------------------------------
* The graphical `Erdős–Gallai theorem`_.
* The number of partitions of n (the partition numbers): OEIS:`A000041`_.

Acknowledgement
---------------


.. _`Erdős–Gallai theorem`: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem
.. _`crime network dataset`: https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt
.. _`Phys. Rev. E paper`: https://doi.org/10.1103/PhysRevE.96.032312
.. _`A000041`: https://oeis.org/A000041
.. _`Travis CI tests`: https://travis-ci.org/github/junipertcy/simpliciality_test
.. _`SageMath`: https://www.sagemath.org/index.html
.. _`test/`: test/
