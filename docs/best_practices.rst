Best practices
==============

Doing a genome-wide association study?
--------------------------------------
If you want to evaluate the effect of individual genomic variants on a phenotype of interest,
while accounting for possible confounders, you will want to run a genome-wide association study.
This will give a p-value for every genomic variant, comparing the alternative hypothesis that
the variant does have an effect on phenotype (has an effect size :math:`\beta > 0`) with the 
null hypothesis that the variant has no effect.

For this mode you will need at least three input files:

* Genetic variants (``--kmers``, ``--vcf`` or ``--pres``).
* A phenotype (``--phenotypes``).
* A representation of the population structure (``--distances`` or ``--similarity``).

For a starting point, have a look at :doc:`tutorial`.

Current 'best-practice' GWAS recommendations:

* Use the ``--lmm`` mode.
* Use a phylogeny to generate the ``--similarity`` matrix.
* Use `unitigs <https://github.com/johnlees/unitig-caller>`__ as the input,
  provided with the ``--kmers`` option. End-to-end analysis is identical to k-mers.
* If you have covariates, provide them with `--covariates` and `--use-covariates`.

Once this works, you may also wish to also add the following extra analyses:

* A burden test with ``--vcf`` and ``--burden``.
* Tests of other forms of variation (genes, structural variants from 
  `panaroo <https://gtonkinhill.github.io/panaroo/#/>`__.
* Extract lineage effects with ``--lineage``.

Trying to predict a phenotype from genetics?
--------------------------------------------
If you want to predict a phenotype in new samples where it is unmeasured, or look at the 
power of genetic variants to predict a phenotype, you'll want to use a whole-genome model.

You will need:

* Genetic variants (``--kmers``, ``--vcf`` or ``--pres``).
* A phenotype (``--phenotypes``). 

A good starting place is to read :doc:`predict`.

Current 'best-practice' prediction recommendations:

* Use ``--wg enet --save-vars`` and ``--wg enet --load-vars`` to save time in future runs.
* Use unitigs, if you can.
* For large variant sets, use a small number of ``--cpu`` to keep memory use manageable.
* Divide the population into strains with `PopPUNK <https://www.poppunk.net>`__ and use
  these definitions with ``--lineage-clusters`` and ``--sequence-reweighting``.
* Turn the correlation filter off with ``--cor-filter 0``.

Trying to calculate heritability?
---------------------------------
If you want an estimate of what proportion of the phenotype variance can be explained
by genomic variation, known as the heritability :math:`h^2`, you can use either of the
above modes to do this.

With ``--lmm`` an estimate for :math:`h^2` will be printed to stderr, based on the GCTA
model (all variants affect the phenotype, with normally distributed effect sizes).

With ``--wg enet`` and estimate for :math:`h^2` will also be printed to stderr, based on
the average prediction accuracy :math:`R^2` in held-out samples during cross-validation.

For a comparison of these approaches, see:

:emphasis:`Lees, John A., Mai, T. T., et al. Improved inference and prediction of bacterial genotype-phenotype associations 
using interpretable pangenome-spanning regressions. (2020)`

Preprint: `<https://doi.org/10.1101/852426>`__
