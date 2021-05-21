<div align="left">
  <img src="https://github.com/junipertcy/simplicial-test/blob/main/docs/assets/simplicial-test-logo.png?raw=true" alt="logo" width=400>
</div>

[![python](https://img.shields.io/pypi/v/simplicial-test)](https://pypi.org/project/simplicial-test/)
[![license](https://img.shields.io/badge/license-LGPL-green.svg?style=flat)](https://github.com/junipertcy/simplicial-test/blob/master/COPYING)


**Simplicial-test** implements a deterministic, backtracking-based algorithm to check whether a joint degree sequence is simplicial.

This is the software repository behind the paper:

* *Construction of simplicial complexes with prescribed joint degree sequences*, [Tzu-Chi Yen](https://junipertcy.info/) (2021).

Read it on: [arXiv](). (link will be updated later in May)

* For full documentation, please visit [this site](https://docs.netscied.tw/simplicial-test/index.html).
* For general Q&A, ideas, or other things, please visit [Discussions](https://github.com/junipertcy/simplicial-test/discussions).
* For software-related bugs, issues, or suggestions, please use [Issues](https://github.com/junipertcy/simplicial-test/issues).


First steps
-----------
`Simplicial-test` is on PyPI. To start, hit this command on the shell:

```sh
$ pip install simplicial-test
```

Here's a typical *simplicial realization problem*: Can `d = (3, 3, 2, 2, 1, 1, 1, 1)`
and `s = (4, 3, 2, 2, 2, 1)` be the joint degree sequence of some simplicial complex? 

In your Python console, `simplicial-test` is invoked using:

```python
>>> from simplicial_test import Test, compute_joint_seq, if_facets_simplicial
>>> degree_list = [3, 3, 2, 2, 1, 1, 1, 1]
>>> size_list = [4, 3, 2, 2, 2, 1]
>>> st = Test(degree_list, size_list)  # visit the docs for other options, like setting a cutoff to give up.
>>> is_simplicial, facets = st.is_simplicial()  # actual computation
>>> assert is_simplicial is True
>>> assert if_facets_simplicial(facets) is True
>>> assert compute_joint_seq(facets) == (sorted(degree_list, reverse=True), sorted(size_list, reverse=True))
>>> facets
((0, 1, 2, 4), (0, 1, 3), (0, 5), (1, 6), (2, 3), (7,))
```

Command-line Utility and Unit Tests 
-----------------------------------
Alternatively, you may install `Simplicial-test` from the source:

```shell
$ git clone https://github.com/junipertcy/simplicial-test.git
$ cd simplicial-test
$ python setup.py install  # if you do not want to install it yet, skip this step.
```

In the `simplicial-test` folder, 
there's a useful script that allows you to do the simplicial test (even when it's not installed), 
by hard-coded integer sequences in [`datasets/00_degs.txt`](datasets/00_degs.txt) 
and [`datasets/00_sizes.txt`](datasets/00_sizes.txt) as the input.

```shell
$ python utils/is_simplicial.py --help

Usage: is_simplicial.py [OPTIONS]

Options:
  -k, --degree_seq_file FILENAME  Path to degree sequence file.
  -s, --size_seq_file FILENAME    Path to size sequence file.
  -c, --cutoff INTEGER            Cutoff (max. steps) to the search.
  -w, --width INTEGER             Search width (see the docs).
  -d, --depth INTEGER             Backtrack depth (see the docs).
  --verbose                       Turn on the verbose mode.
  --help                          Show this message and exit.


```

To run the simplicial test on the command line:
```sh
$ python utils/is_simplicial.py -k datasets/00_degs.txt -s datasets/00_sizes.txt

Yes, the joint degree sequence is simplicial. 
The complex is: ((0, 1, 2, 3), (0, 1, 4), (0, 1, 5), (0, 2, 4), (3, 6), (7,))

```

Look, the program gives an affirmative answer, with a realization in the standard output.

`Simplicial-test` implements a search algorithm for solving
the *simplicial realization problem*. If your input joint sequence is not realizable
(as a simplicial complex), we may need to traverse the entire search tree,
which would take a huge amount of time!

Happily, more than 90% of the input joint sequences lies in the *polynomial regime* (check out the paper),
which means that they can be solved easily.

For example, you can assemble a simplicial complex from the joint sequence
of the [crime network dataset](https://github.com/jg-you/scm/blob/master/datasets/crime_facet_list.txt),
from this inspiring [Phys. Rev. E paper](https://doi.org/10.1103/PhysRevE.96.032312), 
which has 551 nodes and 194 facets.

```sh
$ python utils/is_simplicial.py -k datasets/crime_degs.txt -s datasets/crime_sizes.txt
```    

Boom! It's rather fast, isn't it? 

Lastly, remember to check out [the docs](https://docs.netscied.tw/simplicial-test/index.html) for the full documentation!

### Development
`Simplicial-test` uses [poetry](https://python-poetry.org/) for packaging 
and [tox](https://tox.readthedocs.io/en/latest/) for testing. 
Once you have these two libraries installed, 
you can run `tox` to standardize testing on various Python versions (cf. [tox.ini](./tox.ini)) 

You may also test with your local environment (for details, see [tests/](tests/)), by running `pytest`.

Related links
-------------
* [The wonderful paper](https://doi.org/10.1103/PhysRevE.96.032312) on the Simplicial Configuration Model (SCM) and [its arXiv version](https://arxiv.org/abs/1705.10298).
* The implementation of the [SCM sampler](https://github.com/jg-you/scm).
* The [Erdős–Gallai theorem](https://en.wikipedia.org/wiki/Erd%C5%91s%E2%80%93Gallai_theorem) for graphical realization.
* The [Havel–Hakimi algorithm](https://en.wikipedia.org/wiki/Havel%E2%80%93Hakimi_algorithm) for _constructive_ graphical realization.
* The partition numbers: [A000041](https://oeis.org/A000041).
* SageMath's [`random_element_uniform()`](https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/partition.html#sage.combinat.partition.Partitions_n.random_element_uniform) that returns a random partition of n with uniform probability.
* To see how `Simplicial-test` scales, check out [this benchmark](https://docs.netscied.tw/simplicial-test/dataset/benchmark.html).


Acknowledgement
---------------
The simplicial-test library is inspired and supported by Josh Grochow, Jean-Gabriel Young, and Alice Patania.

We want to thank Stefanie Molin ([@stefmolin](https://github.com/stefmolin)) for their Custom-Colormaps,
Iacopo Iacopini ([@iaciac](https://github.com/iaciac)) for their py-draw-simplicial-complex,
whose libraries make pretty figures,
and Austin Benson ([@arbenson](https://github.com/arbenson)) for their hypergraph datasets that helped identify bugs in larger systems.
