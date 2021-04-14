:code:`cutoff` (to stop)
~~~~~~~~~~~~~~~~~~~~~~~~
The algorithm will stop and conclude that the input sequence is not simplicial,
if :math:`\tau_{\text{conv.}} \geq \text{cutoff}`. So this algorithm has false negatives,
i.e., concluding that some input is not simplicial when it's *actually* simplicial.

By default, :math:`\text{cutoff} = 10^5`.

.. warning::
   On average, it takes longer to verify if a input is non-simplicial.
   You may wish to make :code:`cutoff` larger,
   just to inspect thoroughly. But be warned that, the entire search tree grows combinatorially,
   as you might spend a long time, just to learn that the input is not simplicial.

How to set :code:`cutoff`?
--------------------------
You *should* set a :code:`cutoff`, to roughly :code:`1e5` or larger (or you may wait forever for larger inputs).

There is only a small fraction of simplicial input sequence
that will take a long time to determine, whereas there are a lot of hard non-simplicial inputs.
Setting a moderate :code:`cutoff` will avoid false negatives.
