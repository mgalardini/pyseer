.. pyseer documentation master file, created by
   sphinx-quickstart on Thu Feb  1 13:55:57 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyseer documentation
==================================
``pyseer`` is a python reimplementation of `seer <https://github.com/johnlees/seer>`_, which was written in C++.
``pyseer`` uses linear models with fixed or mixed effects to estimate the
effect of genetic variation in a bacterial population on a phenotype of
interest, while accounting for potentially very strong confounding population
structure. This allows for genome-wide association studies (GWAS) to be performed in
clonal organisms such as bacteria and viruses.

.. image:: pyseer_logo.png
   :alt: pyseer - python version of seer
   :align: right

The original version of ``seer`` used sequence elements (k-mers) to represent
variation across the pan-genome. ``pyseer`` also allows variants stored in VCF
files (e.g. SNPs and INDELs mapped against a reference genome) or Rtab files
(e.g. from `roary <https://sanger-pathogens.github.io/Roary/>`_ or
`piggy <https://github.com/harry-thorpe/piggy>`_ to be used too). There are also a greater range of association models
available, and tools to help with processing the output.

Testing shows that results (p-values) should be the same as the original
``seer``, with a runtime that is roughly twice as long as the optimised C++
code.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   installation.rst
   options.rst
   usage.rst
   tutorial.rst
   multiprocessing.rst
   api.rst

Index:
------

* :ref:`genindex`
* :ref:`search`
