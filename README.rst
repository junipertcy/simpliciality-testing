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

:code:`simplicial-test` implements a recursive algorithm to check whether a joint degree sequence is simplicial.

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
:code:`simplicial-test` is on PyPI. To start, hit this command on your shell,

..

    pip install simplicial-test

Please check the docs for `its usage`_ in the Python interactive console.

----

Here, we assume you have instead pulled the library from GitHub, via,

..

    git clone https://github.com/junipertcy/simplicial-test.git

And entered the :code:`simplicial-test` folder. Let's ask: Can :code:`d = (2, 3, 1, 1, 4, 2, 2, 1)`
and :code:`s = (3, 3, 2, 1, 4, 3)` be the joint degree sequence of some simplicial complex?
Let's do a simplicial test! Try:

..

    python utils/is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt

Look, the program gives an affirmative answer, with a realization in standard output.

As an alternative, you may try with a bunch of unit tests (for details, see `tests/`_):

..

    pytest

:code:`simplicial-test` implements a deterministic, backtracking-based search algorithm for solving
the *simplicial realization problem*. If your input joint sequence is not realizable
(as a simplicial complex), often we would need to traverse the entire search tree,
which would take a huge amount of time!

Happily, more than 90% of the input joint sequences lies in the *polynomial regime*,
which means that they can be solved easily.

For example, you can assemble a simplicial complex from the joint sequence of the `crime network dataset`_,
from this inspiring `Phys. Rev. E paper`_, which has 551 nodes and 194 facets.

..

    python utils/is_simplicial.py -k datasets/crime_degs.txt -s datasets/crime_sizes.txt

Boom! It's rather fast, isn't it?



Related links
-------------
* The implementation of the `Simplicial Configuration Model`_.
* The graphical `Erdős–Gallai theorem`_.
* The partition numbers: `A000041`_.
* To see how the algorithm scales, check out `this benchmark`_.


Acknowledgement
---------------
The simplicial-test library is inspired and supported by the following great humans,
Josh Grochow, Jean-Gabriel Young, and Alice Patania.

We want to thank Stefanie Molin (`@stefmolin`_) for their Custom-Colormaps,
Iacopo Iacopini (`@iaciac`_) for their py-draw-simplicial-complex,
whose libraries make pretty figures,
and Austin Benson (`@arbenson`_) for their hypergraph datasets that helped identify bugs in larger systems.


.. _`Erdős–Gallai theorem`: https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem
.. _`crime network dataset`: https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt
.. _`Phys. Rev. E paper`: https://doi.org/10.1103/PhysRevE.96.032312
.. _`A000041`: https://oeis.org/A000041
.. _`Travis CI tests`: https://travis-ci.org/github/junipertcy/simpliciality_test
.. _`SageMath`: https://www.sagemath.org/index.html
.. _`tests/`: tests/
.. _`Simplicial Configuration Model`: https://github.com/jg-you/scm
.. _`this benchmark`: https://docs.netscied.tw/simplicial-test/dataset/benchmark.html
.. _`its usage`:
.. _`@stefmolin`: https://github.com/stefmolin
.. _`@iaciac`: https://github.com/iaciac
.. _`@arbenson`: https://github.com/arbenson
