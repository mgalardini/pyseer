Usage
=====
Quick start::

   pyseer --phenotypes phenotypes.tsv --kmers kmers.gz --distances structure.tsv --min-af 0.01 --max-af 0.99 --cpu 15 --filter-pvalue 1E-8 > pyseer.assoc

Will run the original ``seer`` model on given phenotypes and k-mers, using
MDS scaling of the pairwise distances provided to correct for population
structure. This will paralellize the analysis over 15 cores.

See the :doc:`best_practices` page for guidance on which options to use.

.. contents::
   :local:

Input
-----

``pyseer`` will automatically take the intersection of samples found in the
phenotype file and the population structure file. Only variation within these
samples will be considered. Information on this is printed to STDERR.

Phenotype and covariates
^^^^^^^^^^^^^^^^^^^^^^^^
The phenotype file is required to be supplied using the ``--phenotypes``
option. The format is tab-delimited, with the sample name in the first
column, and the phenotype in the last column. A header is required as the first
row::

   samples continuous      binary
   sample_1        1       0
   sample_2        2       1
   sample_3        3       1
   sample_4        4       1
   sample_5        5       1
   sample_6        6       1
   sample_7        7       0

The default column to use as the phenotype is the last column, but you can
provide an explicit value with ``--phenotype-column``.
Missing phenotypes can be supplied as 'NA'. If all values are 0 or 1 a binary
phenotype is assumed (only relevant for the fixed effect model), otherwise a continuous phenotype is used.
Use ``--continuous`` to force this behaviour.

.. warning:: Using numbers as the sample names has been reported to cause
   problems in some modes and versions of pyseer. While we have tried to fix
   this issue, if you run into trouble try chaning your sample names into a 
   string (e.g. by adding an underscore at the end of every name).

Covariate files (``--covariates``) must be tab-delimited with a header row, and the first column
must contain the sample names::

   samples      time       cluster
   sample_1        1       cluster1
   sample_2        2       cluster2
   sample_3        3       cluster0
   sample_4        4       cluster1
   sample_5        5       cluster2
   sample_6        6       cluster0
   sample_7        7       cluster1

Choose which covariates to use with ``--use-covariates``. Provide space
separated column numbers to use. The default is that the covariates are labels,
but for a quantitative covariate add 'q' after the column number. For the above
example ``--use-covariates 2q 3`` would be the correct argument.

k-mers
^^^^^^
Variable length k-mers counted by `fsm-lite <https://github.com/nvalimak/fsm-lite>`_
or `dsm-framework <https://github.com/HIITMetagenomics/dsm-framework>`_ are input with the
``--kmers`` option. This file is assumed to be gzipped, use the
``--uncompressed`` option if they aren't. If you wish to use `dsk <https://github.com/GATB/dsk>`_ to count
k-mers you will need to use ``combineKmers`` from the original ``seer``
installation to convert them to the correct input format.

If needed, both fsm-lite and seer can be installed through conda. See :doc:`installation` for
details.

.. note:: For common variation k-mers or unitigs should probably be your variant of choice.
   ``seer`` was mainly designed to work with k-mers, due to their ability to
   test variation across the pan-genome without the need to call variants
   against multiple references, or deal with the complexities of constructing
   accurate COGs for the whole population. We have included these input formats
   for convenience and flexibility.

   We would recommend the use of SNPs and genes *in addition* to k-mers, or for
   a quick first pass analysis.

unitigs
^^^^^^^
Unitigs are nodes in a compressed de Bruijn graph, and remove some of the redundancy present
in k-mer counting, as well as presenting fewer tests (and advantage both computationally and
statistically) and being easier to interpret thanks to their length and context provided by
the variation graph.

Unitigs can both be counted, and called consistently in new populations, using the
`unitig-caller <https://github.com/johnlees/unitig-caller>`__ package.

An older version of the package, giving the same results, is available as 
`unitig-counter <https://github.com/johnlees/unitig-counter>`__ (see documentation
in the ``README.md``). 

Usage is then identical to k-mers, providing input with the ``--kmers`` options, and ``--uncompressed`` if necessary.

.. note:: Both packages can be installed thorough conda, see :doc:`installation` for
   details.

SNPs and INDELs
^^^^^^^^^^^^^^^
Short variation (SNPs and INDELs) can be read from a VCF file using the ``PySAM`` module. Simply use
the ``--vcf`` option to read in your file.

If you have multiple VCF files (e.g. one per sample) you can combine them with
``bcftools``::

   bcftools merge -m none -0 -O z *.vcf.gz > merged.vcf.gz

Sample names are taken from the header row. Only one ``ALT`` variant per row is supported,
if you have multiple alternative variants use::

   bcftools norm -m - <in.vcf> > out.vcf

to split them into multiple rows otherwise they will be skipped. If ``FILTER``
fields are present only those with 'PASS' will be processed.

.. note::
   The ``GT`` field is used to determine variant presence/absence.
   '0' or '.' is absence, anything else is presence.

Genes and intergenic regions, or any other variant type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
COG or intergenic region variation is represented as an .Rtab file by `roary <https://sanger-pathogens.github.io/Roary/>`_ and
`piggy <https://github.com/harry-thorpe/piggy>`_::

   Gene	sample_1	sample_2
   COG1	1	1
   COG2	1	0

These can be used directly with ``--pres``, and this format can be used flexibly
to represent variants from other sources.

Rare variants
^^^^^^^^^^^^^

``pyseer`` supports burden testing of rare variants. Variants at low frequency
which are associated with the phenotype cannot be detected by a standard
regression model. A burden test groups sets of rare variants with the same
predicted biological effect, and then treats these sets like common variants.

.. note:: Group variants only with the same predicted functional effect.
   A good start would be all loss of function mutations (frameshift or stop
   gained/nonsense) within a gene. This can be expanded to operons or pathways,
   and to variants predicted as damaging (missense) or all variants. Burden
   tests assume all variants in a group have the same direction of effect, and
   will lose power if this assumption is broken.

To run a burden test, available under any of the association models below,
requires a VCF file of SNPs and INDELs. First predict the function of mutations
(using `VEP <https://useast.ensembl.org/info/docs/tools/vep/index.html>`_
or `bcftools csq <http://www.htslib.org/doc/bcftools.html#csq>`_) and filter the
VCF file appropriately on variant frequency and predicted effect::

   bcftools view -Q 0.01 -i 'CSQ[*] ~ "stop_gained" snps_indels.vcf.gz | CSQ[*] ~ "frameshift_variant"' | bgzip -c > low_freq_vars.vcf.gz

Then run ``pyseer`` providing a list of regions to group variants by to the
``--burden`` option and the filtered VCF file with ``--vcf``.
These regions are one per line, with their name and the bcftools style region co-ordinates::

   CDS1    FM211187:3910-3951
   CDS2    FM211187:4006-4057

Multiple regions can be specified for a single burden test, by separating each
region using a comma::

   pathway1    FM211187:4006-4057,FM211187:5673-5777

.. warning:: The same frequency filters as for common variants still apply.
   Only groups within the threshold will be tested. To ensure only rare
   variants enter the sets, you will need to pre-filter the VCF file with
   bcftools as shown above.


Filtering
^^^^^^^^^
Filtering on allele frequency is necessary, unless the input has already been
filtered. We would recommend only including variants with a minor allele count
of at least five. Use ``--min-af`` and ``--max-af`` to achieve this. The
default is to test variants with a MAF > 1%.

If computational resources are limited, you can use the unadjusted p-value as
a pre-filter ``--filter-pvalue``. :math:`10^{-5}` is a reasonable value, or
three orders of magnitude below your final significance threshold. If you just
want to plot the significant results, or save space in the output you can also
print just those passing a final threshold with ``--lrt-pvalue``.

.. warning:: We would recommend not filtering on p-value if possible.
   It is possible that variants not significant before correction may be
   significant afterwards, and taking a final threshold will prevent a Q-Q plot
   from being used to test for inflation of p-values.

Population structure
--------------------

To adjust for population structure, the fixed effects (:ref:`fixed_effects`) model
needs a matrix with distances between all pairs of samples in the analysis::

   	sample_1	sample_2	sample_3
   sample_1	0	0.0115761	0.0119383
   sample_2	0.0115761       0.0     0.0101878
   sample_3	0.0119383       0.0101878       0.0

This file is included with ``--distances``. The default is to perform classical MDS on
this matrix and retain 10 dimensions. The type of MDS performed can be changed
with the ``--mds`` option to metric or non-metric if desired. Once the MDS has run once,
the ``--save-m`` argument can be used to save the result to file. Subsequent runs can
then be provided with this decomposition directly using ``load-m`` rather than recomputing the MDS.

An alternative to using a distance matrix in the fixed effects analysis is to provide clusters of samples with the same genetic
background (e.g. from BAPS) as a categorical covariate with the
``--use-covariates`` option. In this case you should also add the ``--no-distances`` options
to allow running without one of the matrices below, which would define these
covariates twice.

The mixed effects model (:ref:`mixed_model`) needs a matrix with
covariances/similarities included with ``--similarities`` between all pairs of samples in the analysis::

   	sample_1	sample_2	sample_3
   sample_1	0.319	0.004	0.153
   sample_2	0.004	0.004	0.004
   sample_3	0.153	0.004	0.288

This is known as the kinship matrix :math:`K`.
Analagously to the MDS runs, the decomposition can be save with ``--save-lmm``
and loaded with ``--load-lmm`` in subsequent analysis rather than processing the
similarity matrix again.

Both types of matrix are necessarily symmetric. The entries along the diagonal of a pairwise distance
matrix are zeros. The matrices can be generated in three ways.

mash
^^^^
`mash <http://mash.readthedocs.io/en/latest/>`_ can be used to rapidly estimate distance between samples.
First of all create a sketch of all your samples (assuming assembled contigs in fasta
files)::

   mash sketch -s 10000 -o samples *.fa

Calculate the pairwise distances and create a distance matrix::

   mash dist samples.msh samples.msh | square_mash > mash.tsv

These distances can only be used with the fixed effects model.

Phylogeny based
^^^^^^^^^^^^^^^
If you have a high quality phylogeny (removing recombination, using a more
accurate model of evolution) using this to calculate pairwise distances may be more accurate than mash.
For the fixed effects model you can extract the
patristic distances between all samples. Using a newick file::

   python scripts/phylogeny_distance.py core_genome.tree > phylogeny_distances.tsv

For use with :ref:`mixed_model` add the ``--calc-C`` or ``--lmm`` option (which are equivalent).
This calculates the similarities based on the shared branch length between each pair's MRCA and
the root (as PDDIST)::

   python scripts/phylogeny_distance.py --lmm core_genome.tree > phylogeny_similarity.tsv

If you want to ignore branch lengths (not usually recommended) use the
``--topology`` option. Other tree formats supported by `dendropy <https://pypi.python.org/pypi/DendroPy>`_
can be used by specifying ``--format``.

Genotype matrix
^^^^^^^^^^^^^^^
For a mixed model association the FaST-LMM default is to use the genotype
matrix (design matrix) of variant presence absence to calculate the kinship
matrix :math:`K = GG^T`. To use this method for the ``--similarity`` option use
the similarity script with any valid pyseer input variant type::

   similarity_pyseer --vcf core_gene_snps.vcf sample_list.txt > genotype_kinship.tsv

Where ``sample_list.txt`` is a file containing sample names to keep, one on
each line.

.. warning:: Choose the input to this command carefully.
   Using too few variants or those which don't represent vertical evolution may
   be inaccurate (e.g. the roary gene presence/absence list). Choosing too many
   will be prohibitive in terms of memory use and runtime (e.g. all k-mers).
   A VCF of SNPs from the core genome is a good tradeoff in many cases.

No population structure correction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
You can run the fixed effects model without a population structure correction.
As this is generally not recommended you need to add the ``--no-distances``
option to allow the analysis to run.

Situations where this may be desirable are when you are using population
structure(/lineage) as the phenotype i.e. looking for k-mers which define
lineages, or if you are correcting for population structure manually using
covariates such as cluster IDs.

Association models
------------------

Symbols used:

=====================  =======
Symbol                 Meaning
=====================  =======
:math:`y`              A vector containing the phenotype for each sample.
:math:`W`              A design matrix containing the covariates, and the MDS components if SEER's model is used.
:math:`a`              Fixed effects for the covariates.
:math:`X`              A design matrix (/vector) containing the variant presence/absence.
:math:`b`              Fixed effects for the variant (also known as beta/effect size).
:math:`K`              The kinship matrix of relations between all pairs of samples.
:math:`G`              The genotype matrix of all variant presence/absence.
:math:`u`              Random effects for each row of the kinship matrix.
=====================  =======

.. _fixed_effects:

Fixed effects (SEER)
^^^^^^^^^^^^^^^^^^^^

If provided with a valid phenotype and variant file this is the default
analysis run by ``pyseer``. In summary, a generalized linear model is run on each
k-mer (variant), amounting to multiple linear regression for continuous
phenotypes and logistic regression for binary phenotypes. Firth regression is
used in the latter case when large effect sizes are predicted.
For details see the `original publication <https://www.nature.com/articles/ncomms12797>`_.

.. math::
   y \sim Wa + Xb

The most important adjustment to this analysis is choosing the number of MDS
components with the ``--max-dimensions`` argument. Once you have your
``--distances`` matrix, draw a scree plot::

   scree_plot_pyseer mash.tsv

This will show the variance explained (the eigenvalues of each MDS component)
for the first 30 dimensions (increased using ``--max-dimensions`` to
``scree_plot_pyseer``). You can pick a value at the 'knee' of this plot, or
choose to include much of the total variation. Consider choosing around the
first 30 components.

.. _mixed_model:

Mixed model (FaST-LMM)
^^^^^^^^^^^^^^^^^^^^^^
A linear mixed model (LMM) of fixed and random effects can be fitted by
adding the ``--lmm`` option, as well as either ``--similarities`` or
``--load-lmm`` from a previous analysis.

.. math::
   y \sim Wa + Xb + Ku

We use `FaST-LMM's <http://dx.doi.org/10.1038/nmeth.1681>`_ likelihood calculation
to compute this model in linear time for each variant. The phenotype is always
treated as continuous, which in the case of case/control data may cause some
loss of power.

The main advantage of this model is that all relationships are implicitly
included and selection of the number of components to retain is not necessary.
In comparison to the fixed effect model this has shown to better control inflation of
p-values (https://elifesciences.org/articles/26255).

In addition this model will output the narrow sense heritability :math:`h^2`, which is the
proportion of variance in phenotype explained by the genetic variation when
maximizing the log-likelihood:

.. math::
   LL(\sigma^2_E, \sigma^2_G, \beta) = \log N (y | X\beta; \sigma^2_GK + \sigma^2_EI) \\
   h^2 = \frac{\sigma^2_G}{\sigma^2_G + \sigma^2_E}

This assumes effect sizes are normally distributed, with a variance proportional
to the total genetic variance (the GCTA model). See
`this paper <http://dx.doi.org/10.1093/molbev/msx328>`_ for more information on
the heritability of pathogen traits.

.. warning:: pyseer will print the :math:`h^2` estimate to STDERR, but it will
   only be valid under the assumptions of the model used. You may wish to
   compare estimates from other software, and particular care should be taken
   with binary phenotypes.

Whole genome models (elastic net)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
All variants can be included at once with the ``--wg`` mode. Currently only the elastic
net is implemented, but more models will be included in future.

An elastic net can be fitted to all the variants at once by providing the ``--wg enet``
option, using the `glmnet <https://web.stanford.edu/~hastie/glmnet/glmnet_alpha.html>`__
package to solve the following problem:

.. math::
     \min_{b_0, b}\frac{1}{N} \sum_{i=1}^N w_i l(y_i, b_0+ b^T x_i)^2+\lambda \left[ (1-\alpha)||b||_2^2/2 + \alpha||b||_1\right]

with the link function :math:`w_i l()` set by the phenotype error distribution.

In this mode, all the variants are read into an object in memory, a correlation-based
filter is applied, the model is fitted, then those variants with non-zero :math:`b`
are printed in the output. The model is fit by ten-fold cross-validation to pick the
:math:`\lambda` which gives the lowest deviance when compared to the true phenotypes. Higher
:math:`\lambda` leads to smaller fitted :math:`b` values.
These values, along with the corresponding best :math:`R^2` will be written to ``STDERR``.
Setting :math:`\alpha` closer to one will remove more variants from the model by giving
them zero beta.

.. tip:: Population structure can be included using ``--sequence-reweighting`` and
      ``--lineage-clusters``. Use of the latter will also use these clusters to give
      a more representative cross-validation accuracy. See :doc:`predict` for more details.

Cross-validation uses ``--cpu`` threads, which is recommended for better performance.

.. warning:: As all variants are stored in memory, and potentially copied, very large
    variant files will cause this method to run out of RAM. We therefore do not recommend
    running on k-mers, but to use unitigs instead. SNPs and genes work fine.

By default, the top 75% of variants correlated with the phenotype are included in the fit.
Variants will include the unadjusted single-variate p-values, if distances have been provided
with either ``--distances`` or ``--load-m`` the adjusted p-values will also be present.

=====================  =======
Option                 Use
=====================  =======
``--save-vars``        Save the object representing all objects to disk. Useful for reruns, or using multiple phenotypes.
``--load-vars``        Load the variants saved to disk, the most time-consuming step.
``--save-model``       Save the fitted model so that one can perform :ref:`enet-predict` on samples with unobserved phenotypes.
``--alpha``            Sets the mixing between ridge regression (0) and lasso regression (1) in the above formula. Default is 0.0069 (closer to ridge regression)
``--n-folds``          Number of folds in cross validation (samples removed to test prediction accuracy). Default is 10.
``--cor-filter``       Set the correlation filter to discard the variants with low correlation to the phenotype. Default is 0.25 (keeping the top 75% variants correlated with phenotype).
=====================  =======

.. note:: When using ``--load-vars`` you still need to provide the original variant file with
    ``--vcf``, ``--kmers`` or ``--pres`` as this is read again to output the selected variants. pyseer will
    test that the checksums of this files is identical to that used with ``--save-vars``, and will
    warn if any difference is detected.

.. _enet-predict:

Prediction with the elastic net
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
If ``--wg`` was used with ``--save-model`` this fit can be used to attempt to predict the
phenotype of new samples without a phenotype label::

    enet_predict --vcf new_snps.vcf.gz old_snps.lasso_model.pkl samples.list > lasso.predictions.txt

Provide the samples you wish to predict the phenotype of in ``samples.list`` along with
comparable variants and covariates to that which were used in the original model. If any
variant or covariate is not found in the new input this will be noted on ``STDERR`` and the
mean values (the originally observed allele frequency) will be used instead. Use
``--ignore-missing`` to turn this off.

See :doc:`predict` for more examples.

.. information:: :math:`\beta` in the output is coded in terms of the minor variant, so
      when making predictions observation vectors will need to be coded in the same manner
      as the reference panel. Using `enet_predict` does this automatically.

Lineage effects (bugwas)
^^^^^^^^^^^^^^^^^^^^^^^^
`Earle et al <https://www.nature.com/articles/nmicrobiol201641>`_ introduced
the distinction between 'lineage' and 'locus' effects. Also see `this review <https://figshare.com/articles/The_background_of_bacterial_GWAS/5550037>`_.
The p-values output by ``pyseer`` are aimed at finding 'locus' effects. To find
lineage effects Earle et al proposed ordering variants by those associated with
both the phenotype and a lineage highly associated with a phenotype. They
performed this by decomposing the random effects to find the principal
component each variant was most associated with, and then order variants by
those principal components most associated with the phenotype.

To perform a similar analysis in ``pyseer``, add the ``--lineage`` option. This
first checks the lineages most associated with the phenotype:

.. math::
   y \sim Wa

writing the results to ``--lineage_file``, ordered by the most associated
lineage. For each variant, after the main regression the lineage the variant
belongs to is chosen by the most significant when regressing the variant
presence/absence on the lineages:

.. math::
   X \sim Wa

To pick lineage effects, those variants assigned to a lineage highly associated
with the phenotype in the ``--lineage_file`` and with a significant p-value
should be chosen. A Manhattan plot, with the x-axis order defined by the
lineage column in the output, can be created.

The default is to use the MDS components to define lineage effects, but you
can supply custom lineage definitions such as BAPS clusters with the
``--lineage-clusters`` options::

   sample_1        BAPS_3
   sample_2        BAPS_16
   sample_3        BAPS_27
   sample_4        BAPS_3

.. note:: One of these clusters will be removed to ensure the regressions are of full rank.
   Therefore there is one cluster variants will never be assigned to. This
   is chosen as the cluster least associated with the phenotype.

Output
------

``pyseer`` writes output to STDOUT, which you can redirect with a pipe ``>``.
The format is tab separated, one line per variant tested and passing filtering,
with the first line as a header. Add ``--print-samples`` to print the k-samples
and nk-samples fields.

Fields for a fixed effect analysis:

=====================  =======
Field                  Meaning
=====================  =======
variant                sequence of k-mer or ID of variant from VCF or Rtab.
af                     allele frequency.  The proportion of samples the variant is present in.
filter-pvalue          association of the variant with the phenotype, unadjusted for population structure.
lrt-pvalue             the p-value of association, adjusted for population structure. This corresponds to the LRT p-value of ``seer``.
beta                   the effect size/slope of the variant. For a binary phenotype, exponentiate to obtain the odds-ratio.
beta-std-err           the standard error of the fit on beta.
intercept              the intercept of the regression.
PCX                    the slope each fixed effect (covariate and MDS component).
k-samples (optional)   the samples the variant is present in (comma separated).
nk-samples (optional)  the samples the variant is not present in (comma separated).
lineage (optional)     the lineage the variant is most associated with.
notes                  notes about the fit.
=====================  =======


Fields for a mixed model analysis:

======================  =======
Field                   Meaning
======================  =======
variant                 sequence of k-mer or ID of variant from VCF or Rtab.
af                      allele frequency.  The proportion of samples the variant is present in.
filter-pvalue           association of the variant with the phenotype, unadjusted for population structure.
lrt-pvalue              the p-value from the mixed model association, as given by FaST-LMM.
beta                    the effect size/slope of the variant. For a binary phenotype, exponentiate to obtain the odds-ratio.
beta-std-err            the standard error of the fit on beta.
variant_h2              the variance in phenotype explained by the variant. The :math:`h^2` for this variant alone.
k-samples (optional)    the samples the variant is present in
nk-samples (optional)   the samples the variant is not present in
lineage (optional)      the lineage the variant is most associated with.
notes                   notes about the fit.
======================  =======


Notes field
^^^^^^^^^^^^

Possible 'notes' are:

===========================  =======
Note                         Meaning
===========================  =======
af-filter                    Variant failed set allele frequency filters ``--min-af`` or ``--max-af``.
pre-filtering-failed         Variant failed ``filter-pvalue`` filter .
lrt-filtering-failed         Variant failed ``lrt-pvalue`` filter.
bad-chisq                    :math:`\chi^2` test was invalid, suggesting either a very high effect size or low allele frequency. Firth regression used.
high-bse                     SE of fit was >3, which may imply a high effect size. Firth regression used.
perfectly-separable-data     Variant presence and phenotype exactly correlate, so regression cannot be fitted.
firth-fail                   Firth regression failed (did not converge after 1000 iterations).
matrix-inversion-error       A pseudo-inverse could not be taken, preventing model from being fitted. This likely implies nearly separable data.
===========================  =======

Number of unique patterns
^^^^^^^^^^^^^^^^^^^^^^^^^
One way to pick the threshold for significance is to use a Bonferroni
correction with the number of unique variant patterns as the number of multiple
tests. When running ``pyseer`` add the ``--output-patterns`` option to write
a file with hashes of the patterns.

Then run the ``count_patterns.py`` script on this output::

   python scripts/count_patterns.py --alpha 0.05 --cores 4 --memory 1000 --temp /tmp patterns.txt

This will return the number of unique patterns and the significance threshold.
``--alpha`` is the unadjusted  significance threshold to use. The other options interface
to GNU ``sort`` to speed up the calculation, and control the amount of data
stored in main memory/where to store on disk.

Processing k-mer output
-----------------------

See the :doc:`tutorial` for full concrete examples.

Mapping to references (phandango)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

K-mers can be mapped to reference genomes using the provided script and a fasta
file of the reference::

   phandango pyseer_kmers.assoc reference_1.fa reference_1.plot

These ``.plot`` files can be dragged and dropped into `phandango <http://jameshadfield.github.io/phandango/#/>`_
along with a reference annotation file (the ``.gff`` file corresponding to the
fasta reference file). Phandango will display the length of the k-mer as well
as its position. The y-axis is :math:`-\mathrm{log}_{10}(p)`.

.. warning:: If all the k-mers are plotted performance will be slow.
   It is computationally challenging to render tens of millions of k-mers with
   a real time interface, so we recommend filtering out those with a p-value
   below a threshold value for interactive performance.

.. _annotate-kmers:

Annotating k-mers
^^^^^^^^^^^^^^^^^

K-mers can also be annotated with the gene they are in, or nearby. This
requires a list of annotations. Trusted references are used first, and
allow a close match of k-mer (using ``bwa mem``). Draft annotations, ideally
those the k-mers were counted from, are used second, and require an exact match
of the k-mer (using ``bwa fastmap``).

K-mers will be iteratively mapped to references in the order provided, either until all the
references are used, or all k-mers have been mapped::

   annotate_hits_pyseer pyseer_kmers.assoc references.txt kmer_annotation.txt

The ``references.txt`` file contains the sequence, annotation and type of the
references to be used::

   D39.fa	D39.gff	ref
   TIGR4.fa	TIGR4.gff	ref
   sample1.fa	sample1.gff	draft
   sample2.fa	sample2.gff	draft

To map all of the k-mers, and ensure good quality annotation where possible, provide
a few trusted references as the first lines in this file. You can then list all of the assemblies
used as input after this, designated as draft.

For each k-mer, each match will be returned in the format 'contig:pos;gene_down;gene_in;gene_up'
i.e. the closest downstream gene, the gene the k-mer is in (if it is), the closest
upstream gene. The gene name will be chosen if in the GFF, otherwise the gene
ID will be used.

.. note:: This analysis uses bedtools to find overlapping and nearby genes.
   A working installation of bedtools is therefore required. The construction
   of each query is slow, so only significant k-mers should be annotated in
   this manner.

To summarise these annotations over all significant k-mers, use the
``summarise_annotations.py`` script::

   python scripts/summarise_annotations.py kmer_annotation.txt

For each gene name, the number of overlapping significant k-mers, maximum p-value, average
MAF and average effect size will be reported. This is ideal input for plotting with
`ggplot2 <http://ggplot2.tidyverse.org/reference/>`_.

Processing unitig output
------------------------

As unitigs are sequence elements of variable length, identical steps can be taken as for k-mers,
as described above.

Additionally, ``cdbg-ops`` provided by installing ``unitig-counter`` can be used to extend
short unitigs leftwards and rightwards by following the neightbouring nodes in the de Bruijn graph.
This can help map sequences which on their own are difficult to align in a specific manner.

Create a file ``unitigs.txt`` with the unitigs to extend (probably your significantly associated hits) and run::

   cdbg-ops extend --graph output/graph --unitigs unitigs.txt > extended.txt

The output ``extended.txt`` will contain possible extensions, comma separated, with lines corresponding
to unitigs in the input. See the help for more options.
