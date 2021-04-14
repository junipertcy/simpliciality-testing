<div align="left">
  <img src="docs/assets/simplicial-test-logo.png" alt="logo" width=400>
</div>

[![python](https://img.shields.io/badge/python-3.8-blue.svg?style=flat)](https://github.com/junipertcy/simplicial-test/blob/master/COPYING)
[![license](https://img.shields.io/badge/license-LGPL-green.svg?style=flat)](https://github.com/junipertcy/simplicial-test/blob/master/COPYING)


**Simplicial-test** implements a deterministic, backtracking-based algorithm to check whether a joint degree sequence is simplicial.

This is the software repository behind the paper:

* *Construction of simplicial complexes with prescribed joint degree sequences*, [Tzu-Chi Yen](https://junipertcy.info/) (2021).

Read it on: [arXiv]().

* For full documentation, please visit [this site](https://docs.netscied.tw/simplicial-test/index.html).
* For general Q&A, ideas, or other things, please visit [Discussions](https://github.com/junipertcy/simplicial-test/discussions).
* For software-related bugs, issues, or suggestions, please use [Issues](https://github.com/junipertcy/simplicial-test/issues).


First steps
-----------
`Simplicial-test` is on PyPI. To start, hit this command on the shell:

```sh
$ pip install simplicial-test
```

In your Python console, you may use `Simplicial-test` with the following:

```python
>>> from simplicial_test import Test
>>> degree_list = [3, 3, 2, 2, 1, 1, 1, 1]
>>> size_list = [4, 3, 2, 2, 2, 1]
>>> st = Test(degree_list, size_list)  # visit the docs for other options, like setting a cutoff to give up.
>>> is_simplicial, facets = st.is_simplicial()  # actual computation
>>> is_simplicial
True
>>> facets
((0, 1, 2, 4), (0, 1, 3), (0, 5), (1, 6), (2, 3), (7,))
```

Alternatively, you may install `Simplicial-test` from source:

```sh
$ git clone https://github.com/junipertcy/simplicial-test.git
$ cd simplicial-test
$ python setup.py install  # if you do not want to install it yet, skip this step.
```

In the `simplicial-test` folder, 
there's a useful command-line script that allows you to do the simplicial test (even when it's not installed). 
We ask the following: Can `d = (2, 3, 1, 1, 4, 2, 2, 1)`
and `s = (3, 3, 2, 1, 4, 3)` be the joint degree sequence of some simplicial complex? 
These integer sequences are hard-coded in `datasets/00_degs.txt` and `datasets/00_sizes.txt`. 

Let's do the simplicial test! Try:

```sh
$ python utils/is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt
```

Look, the program gives an affirmative answer, with a realization in standard output.

You may also try with a bunch of unit tests (for details, see [tests/](tests/)), by running:

```sh
$ pytest
```
    
Now you have sensed that `Simplicial-test` implements a search algorithm for solving
the *simplicial realization problem*. If your input joint sequence is not realizable
(as a simplicial complex), we may need to traverse the entire search tree,
which would take a huge amount of time!

Happily, more than 90% of the input joint sequences lies in the *polynomial regime* (check out the paper),
which means that they can be solved easily.

For example, you can assemble a simplicial complex from the joint sequence of the [crime network dataset](https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt),
from this inspiring [Phys. Rev. E paper](https://doi.org/10.1103/PhysRevE.96.032312), which has 551 nodes and 194 facets.

```sh
$ python utils/is_simplicial.py -k datasets/crime_degs.txt -s datasets/crime_sizes.txt
```    

Boom! It's rather fast, isn't it? 

Lastly, remember to check out [the docs](https://docs.netscied.tw/simplicial-test/index.html) for the full documentation!



Related links
-------------
* [The wonderful paper](https://doi.org/10.1103/PhysRevE.96.032312) on the Simplicial Configuration Model (SCM) and [its arXiv version](https://arxiv.org/abs/1705.10298).
* The implementation of the [SCM sampler](https://github.com/jg-you/scm).
* The graphical [Erdős–Gallai theorem](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem).
* The partition numbers: [A000041](https://oeis.org/A000041).
* SageMath's [`random_element_uniform()`](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/partition.html#sage.combinat.partition.Partitions_n.random_element_uniform) that returns a random partition of n with uniform probability.
* To see how `Simplicial-test` scales, check out [this benchmark](https://docs.netscied.tw/simplicial-test/dataset/benchmark.html).


Acknowledgement
---------------
The simplicial-test library is inspired and supported by the following great humans,
Josh Grochow, Jean-Gabriel Young, and Alice Patania.

We want to thank Stefanie Molin ([@stefmolin](https://github.com/stefmolin)) for their Custom-Colormaps,
Iacopo Iacopini ([@iaciac](https://github.com/iaciac)) for their py-draw-simplicial-complex,
whose libraries make pretty figures,
and Austin Benson ([@arbenson](https://github.com/arbenson)) for their hypergraph datasets that helped identify bugs in larger systems.
