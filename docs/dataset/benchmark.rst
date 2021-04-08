Benchmark
=========

The three figures below show the solution found by the greedy algorithm.

Those black regions denote 1's, whereas white blocks denote 0's.
Theoretically, at a deeper level, we need to avoid more "blocking facets," but as it turned out,
we barely need to avoid any blocking facets at later recursive levels, just because of the greedy construction.

.. figure:: ../assets/first_example.png
   :width: 400
   :alt: small first example

   Fig. 1. Small first example.

   We should note that vertices with higher degree will be paired with larger facets, unless such pairing would prohibit
   a solution. Note also that not any facet is in the subset of another (larger) facet.

The following two datasets come from the inspirational `Phys. Rev. E paper`_. Each of these requires less than 1 second
to compute!

.. figure:: ../assets/crime.png
   :width: 400
   :alt: crime dataset

   Fig. 2. Crimes (nodes, n=829) and suspects, victims, and witnesses (facets, f=378) network in St. Louis.

.. figure:: ../assets/pollinator.png
   :width: 400
   :alt: pollinator dataset

   Fig. 3. Flower-visiting insects (nodes, n=679) and plants (facets f=57) network in Kyoto.


.. _`crime network dataset`: https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt
.. _`Phys. Rev. E paper`: https://doi.org/10.1103/PhysRevE.96.032312

Related experiments
-------------------
TODO.


.. list-table:: Empirical networks used in this work.
   :widths: 20 10 10 10 10 10 10 10
   :header-rows: 1

   * - Dataset
     - # of edges, E
     - # of facets, \|s\|
     - # of nodes, \|d\|
     - (s_max, s_min)
     - (d_max, d_min)
     - conv. time
     - conv. time per facet
   * - Macaque brain dataset
     - 628
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
   * - contact-high-school
     - 11770
     - 4862
     - 327
     - (5, 2)
     - (90, 2)
     - 3
     - 0.00064
   * - senate-committees
     - 5165
     - 275
     - 282
     - (31, 4)
     - (59, 1)
     - 0
     - 0
   * - house-committees
     - 11660
     - 302
     - 1290
     - (82, 3)
     - (43, 1)
     - 0
     - 0
   * - mathoverflow-answers
     - 131406
     - 5296
     - 73851
     - (1784, 2)
     - (169, 1)
     - 0
     - 0
   * - senate-bills
     - 92876
     - 3599
     - 294
     - (99, 3)
     - (1147, 1)
     - too slow
     - too slow
   * - house-bills
     - 927075
     - 23267
     - 1494
     - (399, 2)
     - (3824, 1)
     - too slow
     - too slow
   * - walmart-trips
     - 447347
     - 63687
     - 88860
     - (25, 2)
     - (5412, 1)
     - too slow
     - too slow
