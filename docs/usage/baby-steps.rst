How does the algorithm work?
============================

Here, we briefly explain how the algorithm finds a realization. For Python syntax, please see the `GitHub README`_.

.. _`GitHub README`: https://github.com/junipertcy/simplicial-test

🙋🙋🙋 A typical simplicial realization problem 🙋🙋🙋
Can the (input) joint degree sequence :math:`\mathbf{t} = (\mathbf{d}, \mathbf{s})`
make a simplicial instance? Here, we take :math:`\mathbf{d} = \lbrace 3, 3, 2, 2, 1, 1, 1, 1 \rbrace`
and :math:`\mathbf{s} = \lbrace 4, 3, 2, 2, 2, 1 \rbrace`.

**Step 1 (Setup)** First, we use the *nonincreasing* joint degree sequence
to make an "empty" incidence matrix as below.
Here, the header row represents the node labels and the first row and column
denote the node degree sequence and the facet size sequence, respectively. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+

---

**Step 2 (Preprocessing step)** We pair the 1’s in both sequences to make a partial output.
At this point, we may reject the input entirely, following *Rule 1* described in the paper. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 3 (Heuristic policy)** We fill the blanks with a branching heuristic
that prioritizes nodes with higher degrees. In other words, higher degree nodes
will be consumed earlier to make larger facets.
We call :math:`\hat{\sigma} = \lbrace a, b, c, d \rbrace` the *candidate facet*. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ? | ? | ? | ? |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 4 (Validation rules)** We check if :math:`\hat{\sigma} = \lbrace a, b, c, d \rbrace`
breaks any "obvious" conditions; for full set of rules, see *Rules 1-3* of the paper.
Here, we found that if we chose :math:`\lbrace a, b, c, d \rbrace`,
we left with 4 more "non-shielding sites", but there are still 5 more facets to fill,
meaning that we are doomed to violate the no-inclusion constraint later. This is essentially *Rule 2*.
Therefore, we reject :math:`\lbrace a, b, c, d \rbrace` and proceed with the next candidate in heuristic order,
namely, :math:`\lbrace a, b, c, e \rbrace`. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ? | ? | ? |   | ? |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 5 (Recursion)** We find that :math:`\lbrace a, b, c, e \rbrace` violates no validation rules.
So we move on to the next facet. This step is essentially the previous *Step 3*.
Note that we do not even choose :math:`\lbrace a, b, c \rbrace`
as a candidate because it is a subset of a former facet. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ✓ | ✓ | ✓ |   | ✓ |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 | ? | ? |   | ? |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 5 (Recursion)** We find that :math:`\lbrace a, b, d \rbrace` violates no validation rules. Moving forward... ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ✓ | ✓ | ✓ |   | ✓ |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 | ✓ | ✓ |   | ✓ |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 | ? |   |   |   |   | ? |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 5 (Recursion)** We find that :math:`\lbrace a, f \rbrace` violates no validation rules. Moving forward... ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ✓ | ✓ | ✓ |   | ✓ |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 | ✓ | ✓ |   | ✓ |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 | ✓ |   |   |   |   | ✓ |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   | ? |   |   |   |   | ? |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   |   |   |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+

---

**Step 6 (Stop)** We find that :math:`\lbrace b, g \rbrace` violates no validation rules. Moving forward...
We have no option but to choose :math:`\hat{\sigma} = \lbrace c, d \rbrace` for the last facet.
Since this makes :math:`| \mathbf{d} | = | \mathbf{s} | = 0`, we accept the facet and then stop. ⬇

+---+---+---+---+---+---+---+---+---+
|   | a | b | c | d | e | f | g | h |
+===+===+===+===+===+===+===+===+===+
|   | 3 | 3 | 2 | 2 | 1 | 1 | 1 | 1 |
+---+---+---+---+---+---+---+---+---+
| 4 | ✓ | ✓ | ✓ |   |✓  |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 3 | ✓ | ✓ |   | ✓ |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 | ✓ |   |   |   |   |✓  |   |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   | ✓ |   |   |   |   | ✓ |   |
+---+---+---+---+---+---+---+---+---+
| 2 |   |   | ✓ | ✓ |   |   |   |   |
+---+---+---+---+---+---+---+---+---+
| 1 |   |   |   |   |   |   |   | ✓ |
+---+---+---+---+---+---+---+---+---+


.. Attention::

   #. Because we use "not subsetted" proposal facets exclusively, once the *stop criterion* is reached, we are guaranteed to reach a simplicial realization.
   #. The example above has convergence time :math:`\tau_{\text{c}} = 1`, collected at *Step 4*. For the definition of :math:`\tau_{\text{c}}`, please check the next page.


.. Note::

   If you need some actual joint degree sequences to start with
   (and pretend that you do not know if they are realizable),
   Prof. Benson has `a nice collection`_ of hypergraph datasets.
   To use them, download a dataset and apply these utility functions in order.

   #. :func:`simplicial_test.utils.read_hyperedge_list`
   #. :func:`simplicial_test.utils.prune_included_facets`
   #. :func:`simplicial_test.utils.compute_joint_seq`
   #. Now you have the :code:`degree_list` and :code:`size_list` to test!


In the paper, we combine the algorithm with a MCMC sampler for the simplicial
configuration model [young-construction-2017]_ to facilitate efficient sampling
of simplicial ensembles from arbitrary joint degree distributions.

.. _`a nice collection`: https://www.cs.cornell.edu/~arb/data/



----

.. [young-construction-2017] Jean-Gabriel Young, Giovanni Petri, Francesco Vaccarino, and Alice Patania,
    "Construction of and efficient sampling from the simplicial configuration model", Phys.
    Rev. E 96, 032312 (2017), :doi:`10.1103/PhysRevE.96.032312`, :arxiv:`1705.10298`.