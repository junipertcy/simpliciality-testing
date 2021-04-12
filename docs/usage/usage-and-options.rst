Basic Parameters
================================
:code:`Simplicial-test` implements a depth-bounded tree search algorithm, which enables a number
of techniques that may make it faster for certain input classes.

::

    from simplicial_test import Test


    degree_list = (3, 3, 2, 2, 1, 1, 1, 1)
    size_list = (4, 3, 2, 2, 2, 1)
    st = Test(degree_list, size_list, depth=1e5, width=1e5, cutoff=1e5)
    is_realizable, facets = st.is_simplicial()

.. toctree::
   :maxdepth: 1
   :caption: Basic Parameters

   cutoff
   width
   depth





Other data attributes
~~~~~~~~~~~~~~~~~~~~~
The main data holder is :class:`simplicial_test.utils.SimplicialDepot`.