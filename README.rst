simpliciality-test
==================

.. image:: https://img.shields.io/badge/python-3.8-blue.svg?style=flat
   :target: https://github.com/junipertcy/simpliciality_test/blob/master/COPYING
   :alt: python
.. image:: https://img.shields.io/badge/license-LGPL-green.svg?style=flat
   :target: https://github.com/junipertcy/simpliciality_test/blob/master/COPYING
   :alt: license
.. figure:: docs/assets/simplicial-test-logo.png
   :width: 200
   :alt: logo

`Simplicial-test` implements a recursive algorithm to check whether a joint degree sequence is simplicial.

This is the software repository behind the paper:

* *Construction of simplicial complexes with prescribed joint degree sequences*, `Tzu-Chi Yen`_ (2021).

Read it on: [`arXiv`_].

* For full documentation, please visit `this site`_.
* For general Q&A, ideas, or other things, please visit `Discussions`_.
* For software-related bugs, issues, or suggestions, please use `Issues`_.

.. _`this site`: https://docs.netscied.tw/simplicial-test/index.html
.. _`Discussions`: https://github.com/junipertcy/simplicial-test/discussions
.. _`Issues`: https://github.com/junipertcy/simplicial-test/issues
.. _`Tzu-Chi Yen`: https://junipertcy.info/
.. _`arXiv`:

First steps
-----------
Before we ship the library to PyPI, the program can be downloaded from GitHub:

..

   git clone https://github.com/junipertcy/simplicial-test.git

The library requires `SageMath`_ and a few PyPI libraries:

..

    pip install -r requirements.txt

To warm up, try:

..

    python is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt

Or, you may try with a bunch of unit tests (for details, see `tests/`_):

..

   pytest

This is a deterministic, backtracking-based search algorithm for solving the *simplicial realization problem*.
If your input joint sequence is not realizable (as a simplicial complex), often we would need
to traverse the entire search tree, which would take a huge amount of time!
Happily, more than 90% of the input joint sequences lies in the *polynomial regime*,
which means that they can be solved easily.

For example, you can assemble a simplicial complex from the joint sequence of the `crime network dataset`_,
from this inspiring `Phys. Rev. E paper`_, which has 551 nodes and 194 facets.

..

    python is_simplicial.py -k datasets/crime_degs.txt -s datasets/crime_sizes.txt

Boom! It's rather fast, isn't it?



Related links
-------------
* The implementation of the `Simplicial Configuration Model`_.
* The graphical `Erdős–Gallai theorem`_.
* The partition numbers: `A000041`_.


Acknowledgement
---------------


.. _`Erdős–Gallai theorem`: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem
.. _`crime network dataset`: https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt
.. _`Phys. Rev. E paper`: https://doi.org/10.1103/PhysRevE.96.032312
.. _`A000041`: https://oeis.org/A000041
.. _`Travis CI tests`: https://travis-ci.org/github/junipertcy/simpliciality_test
.. _`SageMath`: https://www.sagemath.org/index.html
.. _`tests/`: tests/
.. _`Simplicial Configuration Model`: https://github.com/jg-you/scm

