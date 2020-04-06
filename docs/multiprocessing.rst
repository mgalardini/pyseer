Multiprocessing
===============

``pyseer`` supports the use of multiple CPUs through the ``--cpu`` option. This
sends batches of processed variants to a core, which will fit the chosen model
on all variants in the batch.

The constant ``--block-size`` controls the number of variants sent to each
core. The higher this is set the more efficient the use of CPUs will be (up to
a limit, set by the time spent reading the variant input) at the expense of
a roughly linear increase in memory usage. The default is 1000, using which on
8 cores required around 1.5Gb of memory for a 1.4x speedup with the mixed model.
Increasing this to 30000 while using 4 cores gave a similar (1.5x) speedup, but needed 12Gb of memory.

Depending on your computing architecture, you may wish to split the input and
run separate jobs. This will be more efficient, but is less convenient. This
can be done using GNU ``split``::

   split -d -n l/8 fsm_kmers.txt fsm_out

This would split the input k-mers into 8 separate files.

Prediction
----------
The ``--wg enet`` mode also supports CPUs, but can be very memory-hungry (memory
use scales linearly with number of cores). For large datasets, if you are running
out of memory, you may wish to try with just a single core.