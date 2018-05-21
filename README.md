pyseer
======

[SEER](https://github.com/johnlees/seer), reimplemented in python by
[Marco Galardini](https://github.com/mgalardini) and [John Lees](https://github.com/johnlees)

    pyseer --phenotypes phenotypes.tsv --kmers kmers.gz --distances structure.tsv --min-af 0.01 --max-af 0.99 --cpu 15 --filter-pvalue 1E-8

[![Build Status](https://travis-ci.org/mgalardini/pyseer.svg?branch=master)](https://travis-ci.org/mgalardini/pyseer)
[![Documentation Status](https://readthedocs.org/projects/pyseer/badge/?version=master)](http://pyseer.readthedocs.io/)
[![PyPI version](https://badge.fury.io/py/pyseer.svg)](https://badge.fury.io/py/pyseer)
[![Anaconda package](https://anaconda.org/bioconda/pyseer/badges/version.svg)](https://anaconda.org/bioconda/pyseer/badges/version.svg)

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
[preprint](https://www.biorxiv.org/content/early/2018/02/15/266312) for full details.

Citation
--------

``Lees, John A., Galardini, M., et al. "pyseer: a comprehensive tool for microbial
pangenome-wide association studies". bioRxiv 266312 (2018).``

``doi: 10.1101/266312``

``Lees, John A., et al. "Sequence element enrichment analysis to determine
the genetic basis of bacterial phenotypes." Nature communications 7 (2016): 12797.``

``doi: 10.1038/ncomms12797``

Documentation
--------------------

Full documentation is available at [readthedocs](http://pyseer.readthedocs.io/).
You can also build the docs locally (requires `sphinx`) by typing:

    cd docs/
    make html

Prerequisites
-------------

Between parenthesis the versions the script was tested against:

* `python` 3+ (3.5.3)
* `numpy` (1.13.3)
* `scipy` (1.0.0)
* `pandas` (0.21.0)
* `scikit-learn` (0.19.1)
* `statsmodels` (0.8.0)
* `pysam` (0.13)
* `dendropy` (4.3.0)

If you would like to use the `scree_plot_pyseer` script you will also need to have
`matplotlib` installed.
If you would like to use the scripts to map and annotate kmers, you will also need
`bwa`, `bedtools`,
`bedops` and `pybedtools` installed.

Installation
------------

The easiest way to install `pyseer` is through `pip`:

    python -m pip install pyseer

If you want multithreading make sure that you are using a version 3 python interpreter:

    python3 -m pip install pyseer

You can also use the [conda](https://conda.io/docs/) package manager, using the [bioconda](https://bioconda.github.io/) channel:

    conda install -c bioconda pyseer

**For the impatient**, just clone/download this repository and run:

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
