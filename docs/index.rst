=====================================
The :mod:`Simplicial-test` User Guide
=====================================

We use :code:`simplicial-test` to study the realizability of simplicial complexes
with a given pair of integer sequences,
representing the node degree distribution and facet size distribution, respectively.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. toctree::
   :maxdepth: 1
   :caption: Quick start

   usage/installation
   usage/baby-steps
   usage/usage-and-options
   usage/enumerate-and-sample

.. toctree::
   :maxdepth: 1
   :caption: Dataset

   dataset/format
   dataset/benchmark

.. toctree::
   :maxdepth: 1
   :caption: Module documentation

   src/cores
   src/draw
   src/sample
   src/enumeration

Please contact Tzu-Chi Yen <tzuchi.yen@colorado.edu> if you have any questions.


Acknowledgements
----------------
The :code:`simplicial-test` library is inspired and supported by these great humans,
`Josh Grochow`_, `Jean-Gabriel Young`_, and `Alice Patania`_.

We want to thank
`Stefanie Molin`_ for their `Custom-Colormaps`_, `Iacopo Iacopini`_ for their
`py-draw-simplicial-complex`_, whose libraries make pretty figures,
and `Austin Benson`_ for their `hypergraph datasets`_ that helped identify bugs in larger systems.


.. _`Josh Grochow`: https://www.cs.colorado.edu/~jgrochow/index.html
.. _`Jean-Gabriel Young`: https://www.jgyoung.ca/
.. _`Alice Patania`: https://alpatania.github.io/
.. _`Stefanie Molin`: http://www.columbia.edu/~snm2121/
.. _`Custom-Colormaps`: https://github.com/stefmolin/Custom-Colormaps
.. _`Iacopo Iacopini`: https://iaciac.github.io/
.. _`py-draw-simplicial-complex`: https://github.com/iaciac/py-draw-simplicial-complex
.. _`Austin Benson`: https://www.cs.cornell.edu/~arb/
.. _`hypergraph datasets`: https://www.cs.cornell.edu/~arb/data/