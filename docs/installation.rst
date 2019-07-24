Installation
============
The easiest way to install pyseer and its dependencies is through conda::

    conda install pyseer

If you need conda, download `miniconda <https://conda.io/miniconda.html>`__
and add the necessary channels::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

pyseer can also be installed through pip::

   python -m pip install pyseer

If you want multithreading make sure that you are using a version 3 python interpreter::

   python3 -m pip install pyseer

Prerequisites
-------------
These modules are installed through the pip command above,
but if you have cloned the repository you will need to install the depdencies
yourself.

We used the following versions, though higher should also work:

* python 3+ (3.5.3)
* ``numpy`` (1.13.3)
* ``scipy`` (1.0.0)
* ``pandas`` (0.21.0)
* ``scikit-learn`` (0.19.1)
* ``statsmodels`` (0.8.0)
* ``pysam`` (0.13)
* ``glmnet_py`` (commit ``946b65c``)
* ``matplotlib`` (2.1.0) -- for scree plots
* ``DendroPy`` (4.3.0) -- for phylogeny distances
* ``pybedtools`` (0.7.10) -- for annotating k-mers
* ``bedtools`` (2.27.0) -- for annotating k-mers
* ``bedops`` (2.4.9) -- for annotating k-mers

Test installation
-----------------
Run unit tests::

   pytest -v tests

Test functions and output::

   cd tests/ && bash run_test.sh && cd ../

Other software
--------------
To count k-mers, you may find `fsm-lite <https://github.com/nvalimak/fsm-lite>`__
or the original seer package useful. These can easily be installed with conda, set up
as above::

    conda install fsm-lite
    conda install seer

Packages to count `unitigs <https://doi.org/10.1371/journal.pgen.1007758>`__ and
help with their interpretation are also available on bioconda::

   conda install unitig-counter unitig-caller
