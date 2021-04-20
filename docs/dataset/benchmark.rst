Benchmark
=========
We suggested in the paper that more than half of the (realizable) joint degree sequences
falls in the "polynomial regime". How does this hold in real-world data? Well, the answer is likely affirmative.

.. note::
   We use :math:`\tau_{\text{conv}} = \tau_{\text{r}} + \tau_{\text{b}}` to denote the time that
   is necessary to determine if the input integer sequences can be realized as a simplicial complex.
   We deem a joint sequence polynomial if :math:`\tau_{\text{conv}} = 0`.

   In that case, no backtracking (:math:`\tau_{\text{b}} = 0`) and
   no rejection  (:math:`\tau_{\text{r}} = 0`) is necessary,
   the algorithm either rejects simpliciality immediately or finds an instance in linear time.

The results are measured with :code:`simplicial-test v1.0`,
for a machine with `Intel(R) Core(TM) i7-8700B CPU @ 3.20GHz`_ (6 cores & 12 threads) and 64GB 2667MHz DDR4 memory.


.. _`Intel(R) Core(TM) i7-8700B CPU @ 3.20GHz`: https://ark.intel.com/content/www/us/en/ark/products/134905/intel-core-i7-8700b-processor-12m-cache-up-to-4-60-ghz.html


**TABLE---Summary of the joint degree sequences derived from empirical datasets.**
Shown are the convergence time :math:`\tau_\text{conv}`,
wall time, number of edges :math:`E`, number of nodes (facets) :math:`n` (:math:`m`),
maximal/minimal facet size, and maximal/minimal node degree.

.. list-table::
   :widths: 20 10 10 10 10 10 10 10 10
   :align: left
   :header-rows: 1

   * - Dataset
     - :math:`\tau_\text{conv}`
     - Wall time
     - :math:`E`
     - :math:`m`
     - :math:`n`
     - :math:`(\textbf{s}_\text{max}, \textbf{s}_\text{min})`
     - :math:`(\textbf{d}_\text{max}, \textbf{d}_\text{min})`
     - Reference
   * - crime
     - **0**
     - 621 ms
     - 805
     - 194
     - 551
     - (25, 1)
     - (14, 1)
     -
   * - flower and pollinator
     - **0**
     - 295 ms
     - 1206
     - 91
     - 679
     - (189, 1)
     - (25, 1)
     -
   * - human diseasome
     - **0**
     - 2.1 s
     - 1416
     - 752
     - 1100
     - (10, 1)
     - (15, 1)
     -
   * - senate committees
     - **0**
     - 1.57 s
     - 5143
     - 275
     - 282
     - (31, 4)
     - (59, 1)
     - [stewart-congressional-2017]_
   * - house committees
     - **0**
     - 6.53 s
     - 11640
     - 302
     - 1290
     - (82, 3)
     - (43, 1)
     - [stewart-congressional-2017a]_
   * - contact high school
     - 3
     - 1min 59s
     - 11770
     - 4862
     - 327
     - (5, 2)
     - (90, 2)
     - [mastrandrea-contact-2015]_
   * - contact primary school
     - 61465
     - 7min 22s  [*]
     - 20615
     - 8010
     - 242
     - (5, 2)
     - (174, 20)
     - [stehlé-high-2011]_
   * - senate bills
     - DNF
     - ??
     - 92876
     - 3599
     - 294
     - (99, 3)
     - (1147, 1)
     - [chodrow-hypergraph-2021]_
   * - mathoverflow answers
     - **0**
     - 6h 19min  [*]
     - 131406
     - 5296
     - 73851
     - (1784, 2)
     - (169, 1)
     - [veldt-minimizing-2020]_
   * - walmart trips
     - DNF
     - ???
     - 447347
     - 63687
     - 88860
     - (25, 2)
     - (5412, 1)
     - [amburg-clustering-2020]_
   * - house bills
     - DNF
     - ???
     - 927075
     - 23267
     - 1494
     - (399, 2)
     - (3824, 1)
     - [chodrow-hypergraph-2021]_

[*] We also run on a `Raspberry Pi 4 (8GB)`_ for the "contact primary school" dataset, the wall time is 34min 30s.
However, for the "mathoverflow answers" dataset, MemoryError is raised (it has a memory use peak around 32GB).

DNF: Did Not Finish under the default cutoff threshold.

.. _`Raspberry Pi 4 (8GB)`: https://www.raspberrypi.org/products/raspberry-pi-4-model-b/specifications/

----

.. [chodrow-hypergraph-2021] Philip S. Chodrow, Nate Veldt, and Austin R. Benson,
   "Hypergraph clustering: from blockmodels to modularity." 2021.
   :arxiv:`2101.09611`.
.. [stehlé-high-2011] Juliette Stehlé, Nicolas Voirin, Alain Barrat, Ciro Cattuto, Lorenzo Isella, Jean-François Pinton,
   Marco Quaggiotto, Wouter Van den Broeck, Corinne Régis, Bruno Lina, and Philippe Vanhems,
   "High-Resolution Measurements of Face-to-Face Contact Patterns in a Primary School." PLOS ONE, 2011.
   :doi:`10.1371/journal.pone.0023176`, :arxiv:`1109.1015`.
.. [mastrandrea-contact-2015] Rossana Mastrandrea, Julie Fournet, and Alain Barrat,
   "Contact Patterns in a High School: A Comparison between Data Collected Using Wearable Sensors,
   Contact Diaries and Friendship Surveys." PLOS ONE, 2015.
   :doi:`10.1371/journal.pone.0136497`, :arxiv:`1506.03645`.
.. [stewart-congressional-2017] Charles Stewart III and Jonathan Woon.
   "Congressional Committee Assignments, 103rd to 114th Congresses, 1993--2017: Senate."
.. [stewart-congressional-2017a] Charles Stewart III and Jonathan Woon.
   "Congressional Committee Assignments, 103rd to 114th Congresses, 1993--2017: House."
.. [veldt-minimizing-2020] Nate Veldt, Austin R. Benson, and Jon Kleinberg.
   "Minimizing Localized Ratio Cut Objectives in Hypergraphs."
   Proceedings of the ACM SIGKDD International Conference on Knowledge Discovery and Data Mining (KDD), 2020.
   :doi:`10.1145/3394486.3403222`, :arxiv:`2002.09441`.
.. [amburg-clustering-2020] Ilya Amburg, Nate Veldt, and Austin R. Benson.
   "Clustering in graphs and hypergraphs with categorical edge labels."
   Proceedings of the Web Conference (WWW), 2020.
   :doi:`10.1145/3366423.3380152`, :arxiv:`1910.09943`.

