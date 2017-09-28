pyseer
======

[SEER](https://github.com/johnlees/seer), reimplemented in python

Motivation
----------

Kmers-based GWAS analysis is particularly well suited for bacterial samples,
given their high genetic variability. This approach has been for the first
time formally implemented by [Lees, Vehkala et al.](https://www.nature.com/articles/ncomms12797),
in the form of the [SEER](https://github.com/johnlees/seer) software. While
wrapping my head around SEER's algorithm I figured I could gain a better
understanding if I tried to reimplement it in a language I'm
very familiar with (python). This also allowed me to get around a
few possible bugs related to low sample sizes.

The reimplementation presented here is solely based on my interpretation
of the methods section of the SEER paper, further constrained to what
could be easily (~1-2 days) reimplemented in python. Based on (very) few
tests I've run the results are going to be exactly the same to the
C++ implementation. However, **no guarantee whatsoever is given that the
script works in the same way as the original**, especially for corner cases
and across millions of kmers.

Citation
--------

As this script is a reimplementation of previously published work, if you end
up using it you should cite the original paper. They did the hard work and
they need to be recognized for it.

``Lees, John A., et al. "Sequence element enrichment analysis to determine
the genetic basis of bacterial phenotypes." Nature communications 7 (2016): 12797.``

``doi: 10.1038/ncomms12797``

Prerequisites
-------------

Between parenthesis the versions the script was tested against:

* python 3+ (3.5.3)
* numpy (1.13.1)
* scipy (0.19.1)
* pandas (0.20.3)
* statsmodels (0.8.0)

Installation
------------

Since everything fits in a script (`pyseer`), you can just place it in your `$PATH`,
make it executable (`chmod 755 pyseer`) and type:

    pyseer -h

**Note:** this assumes that your default python's interpreter is at least version 3.
If that is not the case then edit the first line to:

    #!/usr/bin/env python3

Population structure
--------------------

`pyseer` accepts as input a tab-delimited square matrix of distances between samples, with
the first row and column listing sample names. Such a square matrix can be easily obtained
by piping the `mash dist` command into the provided `square_mash` script:

    mash dist samples.msh samples.msh | python square_mash > mash.tsv

**Note:** `square_mash` extracts the sample names as the string after the last `/` character
and up to the first full stop (`.`).

Notes
-----

SEER's features present in this script:

* binary and continuos phenotypes
* population structure correction
* MAF filtering
* kmers prefiltering
* multi-threading

Absent features:

* automatic determination of binary/continuous phenotypes
* binary phenotypes: Firth regression upon failure of bfgs and Newton-Raphson
* user-defined covariates
* betas for covariates are not present in output
* reasons for failures are not reported

Additional features:

* Multidimensional scaling from squared [mash](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) matrix
* Uncompressed kmers input file is acceptable too
* List of samples without the kmer are also outputed

Copyright
---------

Copyright (C) <2017> EMBL-European Bioinformatics Institute

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the   
GNU General Public License for more details.

Neither the institution name nor the name pyseer
can be used to endorse or promote products derived from
this software without prior written permission.
For written permission, please contact <marco@ebi.ac.uk>.

Products derived from this software may not be called pyseer
nor may pyseer appear in their names without prior written
permission of the developers. You should have received a copy
of the GNU General Public License along with this program.
If not, see <http://www.gnu.org/licenses/>.

