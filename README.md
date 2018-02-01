pyseer
======

[SEER](https://github.com/johnlees/seer), reimplemented in python by
[Marco Galardini](https://github.com/mgalardini) and [John Lees](https://github.com/johnlees)

    python pyseer-runner.py --phenotypes phenotypes.tsv --kmers kmers.gz --distances structure.tsv --min-af 0.01 --max-af 0.99 --cpu 15 --filter-pvalue 1E-8

[![Build Status](https://travis-ci.org/mgalardini/pyseer.svg?branch=master)](https://travis-ci.org/mgalardini/pyseer)
[![PyPI version](https://badge.fury.io/py/pyseer.svg)](https://badge.fury.io/py/pyseer)
[![Documentation Status](https://readthedocs.org/projects/pyseer/badge/?version=latest)](http://pyseer.readthedocs.io/en/latest/?badge=version2)

Motivation
----------

Kmers-based GWAS analysis is particularly well suited for bacterial samples,
given their high genetic variability. This approach has been for the first
time formally implemented by [Lees, Vehkala et al.](https://www.nature.com/articles/ncomms12797),
in the form of the [SEER](https://github.com/johnlees/seer) software.

The reimplementation presented here should be consistent with the
current version of the C++ seer, though we do not guarantee this for all
possible cases.

In this version, as well as all the original features, many new features (input types,
association models andoutput parsing) have been implemented. See the
[documentation](http://pyseer.readthedocs.io/en/latest/) for full details.

Citation
--------

``Lees, John A., et al. "Sequence element enrichment analysis to determine
the genetic basis of bacterial phenotypes." Nature communications 7 (2016): 12797.``

``doi: 10.1038/ncomms12797``

Prerequisites
-------------

Between parenthesis the versions the script was tested against:

* python 3+ (3.5.3)
* `numpy` (1.13.3)
* `scipy` (1.0.0)
* `pandas` (0.21.0)
* `scikit-learn` (0.19.1)
* `statsmodels` (0.8.0)
* `pysam` (0.13)

If you would like to use the `scree_plot` script you will also need to have `matplotlib` installed

Installation
------------

The easiest way to install `pyseer` is through `pip`:

    pip install pyseer

If you want multithreading:

    pip3 install pyseer

**For the impatient**, just clone/download this repository and run:

    python pyseer-runner.py

Testing
-------

While waiting on proper unit tests to be implemented, you can check that the script doesn't crash
by running:

    cd test/ && bash run_test.sh && cd ../

Documentation
--------------------

Full documentation is available at [readthedocs](http://pyseer.readthedocs.io/en/latest/).

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
