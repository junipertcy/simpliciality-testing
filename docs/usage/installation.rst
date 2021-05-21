Installation
============

*simplicial-test* can be installed via PyPI.

.. code:: bash

    pip install simplicial-test


.. tip::

    By default, the library uses SageMath's `Combinations`_ to generate proposal
    facets, primarily because it correctly handles the cases when the "pool of
    elements" has duplicate items. However, it cannot be installed through
    :code:`pip` but has to be downloaded /slash/ installed from `their website`_.

    If installing Sage is not feasible, we offer :code:`more-itertools` as the
    fallback option, which has a similar :code:`distinct_combinations()` method.
    This method is, on the average, faster than SageMath's,
    especially in large system sizes (e.g., :math:`E \geq 10^4`). Yet, we've found
    that the generator is not as stable, and in some inputs, it may take a
    longer-than-usual time to output a proposal facet. Therefore, in smaller and most
    cases, SageMath is preferred.


.. _`their website`: https://www.sagemath.org/download.html
.. _`Combinations`: https://doc.sagemath.org/html/en/reference/combinat/sage/combinat/combination.html

