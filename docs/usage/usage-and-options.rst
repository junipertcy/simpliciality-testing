Basic Parameters
================================
:code:`Simplicial-test` implements a depth-bounded tree search algorithm, which enables a number
of techniques that may make it faster for certain input classes.

Here, we explain the :code:`cutoff`, :code:`width`, and :code:`depth` parameters,
which is used to stop the search, to try more candidate facets, and to backtrack, respectively.
There's another :code:`verbose` option, which helps debug.

.. toctree::
   :maxdepth: 1

   cutoff
   width
   depth
   verbose

Please refer to the following code block for the syntax.


::

    st = Test(degree_list, size_list, depth=1e5, width=1e5, cutoff=1e5)
    is_realizable, facets = st.is_simplicial()


Other Data Attributes
~~~~~~~~~~~~~~~~~~~~~
The main data holder is :class:`simplicial_test.utils.SimplicialDepot`. Its "representation" has 4 attributes,
which are :code:`degree_list`, :code:`size_list`, :code:`simplicial`, and :code:`facets`. The last two appears after
the computation, which makes it easy to carry around the "realized object".

Here's a selection of other interesting attributes.

- **conv_time**
   In the paper, we use :math:`\tau_{\text{c}} = \tau_{\text{r}} + \tau_{\text{b}}` to denote the time that
   is necessary to determine if the input integer sequences can be realized as a simplicial complex.
   Every time when the proposed facet is rejected or when the searcher backtracks,
   we increment by 1 the value of :math:`\tau_{\text{c}}`. So, at the end of computation,
   this variable gives an sense of "how hard" the input instance is.

   ::

      st = Test(degree_list, size_list)
      is_realizable, facets = st.is_simplicial()
      print(f"Convergence time is {st.s_depot.conv_time}")



- **levels_traj**
   As you may notice, :math:`\tau_{\text{c}}` is a function of the branching stage :math:`i`,
   meaning that the searcher may experience different amount of rejections at each stage.
   The :code:`levels_traj` records the number of rejections and backtracks at each iterative level.
   It gives you a chance to trace the "trajectory" of the searcher.

- **valid_trials**
   The :code:`valid_trials` is a vector,
   where each element contains a generator for candidate facet proposals at each level.
   It makes backtracking easier because, once we fail at a deeper level and are forced to backtrack,
   we can continue from some unexplored branch (i.e., the next item in the generator).
