pyseer
======

[SEER](https://github.com/johnlees/seer), reimplemented in python by
[Marco Galardini](https://github.com/mgalardini) and [John Lees](https://github.com/johnlees)

    pyseer --phenotypes phenotypes.tsv --kmers kmers.gz --distances structure.tsv --min-af 0.01 --max-af 0.99 --cpu 15 --filter-pvalue 1E-8

![Run tests](https://github.com/mgalardini/pyseer/workflows/Run%20tests/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/pyseer/badge/?version=master)](http://pyseer.readthedocs.io/)
[![Anaconda package](https://anaconda.org/bioconda/pyseer/badges/version.svg)](https://anaconda.org/bioconda/pyseer)

Motivation
----------

Kmers-based GWAS analysis is particularly well suited for bacterial samples,
given their high genetic variability. This approach has been
implemented by [Lees, Vehkala et al.](https://www.nature.com/articles/ncomms12797),
in the form of the [SEER](https://github.com/johnlees/seer) software.

The reimplementation presented here should be consistent with the
current version of the C++ seer (though we do not guarantee this for all
possible cases).

In this version, as well as all the original features, many new features (input types,
association models and output parsing) have been implemented. See the
[documentation](http://pyseer.readthedocs.io/) and
[paper](https://doi.org/10.1093/bioinformatics/bty539) for full details.

Citations
--------

Unitigs and elastic net preprint: `Lees, John A., Tien Mai, T., et al.` [Improved inference and prediction of bacterial genotype-phenotype associations using pangenome-spanning regressions.](https://www.biorxiv.org/content/10.1101/852426v1) `bioRxiv 852426 (2019) doi: 10.1101/852426`

pyseer and LMM implementation paper: `Lees, John A., Galardini, M., et al.` [pyseer: a comprehensive tool for microbial
pangenome-wide association studies.](https://academic.oup.com/bioinformatics/article/34/24/4310/5047751) `Bioinformatics 34:4310–4312 (2018). doi: 10.1093/bioinformatics/bty539`

Original SEER implementation paper: `Lees, John A., et al.` [Sequence element enrichment analysis to determine
the genetic basis of bacterial phenotypes.](https://www.nature.com/articles/ncomms12797) `Nature communications 7:12797 (2016). doi: 10.1038/ncomms12797`

Documentation
--------------------

Full documentation is available at [readthedocs](http://pyseer.readthedocs.io/).

You can also build the docs locally (requires `sphinx`) by typing:

    cd docs/
    make html

Prerequisites
-------------

Between parenthesis the versions the script was tested against:

* `python` 3+ (3.6.6)
* `numpy` (1.15.2)
* `scipy` (1.1.0)
* `pandas` (0.23.4)
* `scikit-learn` (0.20.0)
* `statsmodels` (0.9.0)
* `pysam` (0.15.1)
* `glmnet_python` (commit `946b65c`)
* `DendroPy` (4.4.0)

If you would like to use the `scree_plot_pyseer` script you will also need to have
`matplotlib` installed.
If you would like to use the scripts to map and annotate kmers, you will also need
`bwa`, `bedtools`,
`bedops` and `pybedtools` installed.

Installation
------------

The easiest way to install pyseer and its dependencies is through `conda`::

    conda install pyseer

If you need `conda`, download [miniconda](https://conda.io/miniconda.html)
and add the necessary channels::

    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

`pyseer` can also be installed through `pip`; download this repository
(or one of the [releases](https://github.com/mgalardini/pyseer/releases)), `cd` into it, followed by:

   `python -m pip install .`

If you want multithreading make sure that you are using a version 3 python interpreter:

   `python3 -m pip install .`

**If you want the next pre-release**, just clone/download this repository and run:

    python pyseer-runner.py

Copyright
---------

Copyright 2017 EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
