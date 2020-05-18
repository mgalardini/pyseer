Prediction tutorial
===================

.. |nbsp| unicode:: 0xA0
   :trim:

This page describes how to fit whole genome models with ``--wg``, and how they can be used
to predict the phenotype for new samples. This tutorial starts from the same dataset as the
:doc:`tutorial`, which is described at the top of that page.

.. note:: Presently only the elastic net is implemented, which is the method used in this
        tutorial. Future methods will include random forests and best linear unbiased
        predictors (BLUPs).

.. contents::
   :local:

Fitting a whole-genome model
------------------------------------------------
The first step to performing prediction is to a train a model on genetic data with a known
phenotype. The trained models in pyseer can also be used for association purposes, as the individual
variants associated with the phenotype are reported.

Here we will try and find SNPs which can predict penicillin resistance in *S. pneumoniae*. It
would also be possible to use unitigs by changing ``--vcf`` to ``--kmers``. We will use the same
data as from the GWAS :doc:`tutorial` -- instructions on how to download this data can be found
at the top of that page. 

Variants are loaded, the model is fitted and saved. 
This can all be done in a single step::

    pyseer --vcf snps.vcf.gz --phenotypes resistances.pheno --wg enet \
    --save-vars output/ma_snps --save-model penicillin.lasso --cpu 4 --alpha 1 > selected.txt

    Read 603 phenotypes
    Detected binary phenotype
    Reading all variants
    198248variants [04:46, 691.24variants/s]
    Saved enet variants as output/ma_snps.pkl
    Applying correlation filtering
    100%|████████████████████████████████████████████████████| 89703/89703 [00:51<00:00, 1742.25variants/s]
    Fitting elastic net to top 67277 variants
    [status]	Parallel glmnet cv with 4 cores
    Best penalty (lambda) from cross-validation: 2.09E-02
    Best model deviance from cross-validation: 0.405 ± 4.57E-02
    Best R^2 from cross-validation: 0.822
    Finding and printing selected variants
    Saved enet model as penicillin.lasso.pkl
    198248 loaded variants
    130971 filtered variants
    67277 tested variants
    35 printed variants

.. warning:: You may see warnings about variants with no observations. In this case the VCF has
        many missing calls causing this, which can be ignored. In other settings this often
        points to a mismatch between sample labels in the variant and phenotype files.

``selected.txt`` now contains the selected variants in a GWAS-like format::

    variant	af	filter-pvalue	lrt-pvalue	beta	notes
    26_31771_C_T	3.48E-02	8.06E-02		1.10E-01
    26_292628_G_A	3.78E-01	6.12E-94		9.92E-02
    26_292653_T_C	6.63E-02	2.95E-16		9.56E-01	bad-chisq

To calculate an adjusted p-value you can add ``--distances`` as one would do for
GWAS with the SEER fixed effects model, or create a new variant file with just the
selected variants, then run pyseer again.

Differences of this approach from the univariate GWAS approach covered in :doc:`tutorial`:

* ``--wg enet`` fits an elastic net to all variants with ``--n-folds`` cross-validation (default 10-fold).
* ``--save-vars`` saves the variants loaded by ``--vcf`` in an efficient sparse matrix format,
  which can be quickly loaded for new model fitting.
* ``--save-model`` saves the fitted model so it can be used for prediction.
* ``--cpu`` uses four cores efficiently during cross-valdation.

``--alpha`` controls the mixing between ridge regression and lasso regression. Above we have used a
value of 1, which is lasso regression, selecting just a few variants. We can use a value closer to ridge
regression if desired, which will select more variants with smaller effect sizes::

    pyseer --vcf snps.vcf.gz --phenotypes resistances.pheno --wg enet \
    --load-vars output/ma_snps --save-model penicillin.001 --alpha 0.01 > selected.txt

    Read 603 phenotypes
    Detected binary phenotype
    Reading all variants
    Analysing 603 samples found in both phenotype and loaded npy
    Applying correlation filtering
    100%|████████████████████████████████████████████████████| 89703/89703 [01:03<00:00, 1421.87variants/s]
    Fitting elastic net to top 67275 variants
    Best penalty (lambda) from cross-validation: 8.26E-01
    Best model deviance from cross-validation: 0.402 ± 4.45E-02
    Best R^2 from cross-validation: 0.815
    Finding and printing selected variants
    Saved enet model as penicillin.001.pkl
    198248 loaded variants
    130973 filtered variants
    67275 tested variants
    3523 printed variants

We can load the variants saved previously which saves a lot of time. The variant file is needed
to print the selected variants at the end -- this is checked to ensure it is the same as the one
originally provided.

Loading the variants can also be used when just a subset of ``--phenotypes`` is provided, which
is useful for training-test validation.

Accounting for population structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
As the model includes all genetic variants at once, covariance between them from population
structure can implicitly be included already. However, it is possible to include an explicit
correction for population structure which may improve prediction accuracy in new populations.

This correction is based on providing discrete definitions of lineages/strains. Prepare a file
``lineages.txt`` with the following format::

    7001_3#17	0
    6999_7#9	0
    7622_5#50	0
    6999_1#2	0
    7622_4#1	0
    ...
    7622_2#40	59
    7622_3#86	60
    7622_5#61	61

.. important:: Rare lineages must be represented correctly, i.e. in their own cluster rather
        than being grouped in a 'bin'. One method we recommend to do this is
        `PopPUNK <https://poppunk.readthedocs.io/en/latest/>`__. Connecting samples together
        which are below a certain distance threshold will also work.

If you need to convert from PopPUNK output to this format with the tutorial, you can use the 
following code::

    import csv, re
    reader = csv.DictReader(open("poppunk/poppunk_clusters.csv"))
    writer = csv.DictWriter(open("lineages.txt", "w"), delimiter=' ', fieldnames=reader.fieldnames)
    for row in reader:
        row['Taxon'] = re.match(r'.*/(.*)\.contigs_velvet\.fa', row['Taxon']).group(1)
        writer.writerow(row)

Now add this to the analysis::

    pyseer --vcf snps.vcf.gz --phenotypes resistances.pheno --wg enet \
    --load-vars output/ma_snps --lineage-clusters lineages.txt --sequence-reweighting

    Read 603 phenotypes
    Detected binary phenotype
    Reading all variants
    Analysing 603 samples found in both phenotype and loaded npy
    Applying correlation filtering
    100%|████████████████████████████████████████████████████| 89703/89703 [00:59<00:00, 1513.70variants/s]
    Fitting elastic net to top 67275 variants
    Fitting elastic net to top 67275 variants
    Best penalty (lambda) from cross-validation: 1.17E+00
    Best model deviance from cross-validation: 0.572 ± 8.76E-02
    Best R^2 from cross-validation: 0.815
    Predictions within each lineage
    Lineage	Size	R2	TP	TN	FP	FN
    0	96	0.820	35	57	4	0
    1	55	0.182	2	48	0	5
    ...
    8	18	-0.200	0	15	0	3
    9	18	1.000	0	18	0	0
    Finding and printing selected variants
    198248 loaded variants
    130973 filtered variants
    67275 tested variants
    4357 printed variants

Adding ``--lineage-clusters`` has two effects. Cross-validation will be performed by leaving one strain
out. This will usually take longer as there are more strains than folds, but may help reduce the number
of lineage effects included. Also, training predition accuracy for each lineage will be reported,
making it easier to see whether there are some parts of the data where the model is performing better.
For binary phenotypes :math:`R^2` can be difficult to interpret, so true/false positives/negatives are
also reported.

Adding ``--sequence-reweighting`` has one further effect. Within each lineage, the weight :math:`w_i`
given to each sample in the loss function

.. math::
     \min_{b_0, b}\frac{1}{N} \sum_{i=1}^N w_i l(y_i, b_0+ b^T x_i)^2+\lambda \left[ (1-\alpha)||b||_2^2/2 + \alpha||b||_1\right]

is set by

.. math::
    \frac{1}{u_i} = \sum_{j=1}^N [j \in C(i)] \\
    w_i = u_i \cdot \frac{N}{\sum_{j=1}^N u_j}

where :math:`C(x)` is the lineage cluster of :math:`x`.

This sets the weights as being inversely proportional to the size of the cluster, and rescales all
weights to sum to :math:`N`. Without this option :math:`w_i = 1 \; \forall \; i`.

Using the model to predict phenotype in new samples
---------------------------------------------------
The elastic net models can be used to predict phenotypes in new samples. We will first
split the samples into training and test sets::

    head -500 resistances.pheno > train.pheno
    cat <(head -1 resistances.pheno) <(tail -104 resistances.pheno) > test.pheno
    cut -f 1 test.pheno | sed '1d' > test.samples

.. warning:: This is a random split of the samples, unlikely to be equivalent to different
        sample collections made up of different proportions of strains. Accuracy is likely
        overestimated, but within strain accuracies can be useful.

We will use lasso regression as fewer variants are selected, so if they were uncalled
in the test set this should be less of a problem (but is still an important concern). Fit
a model to the training set::

    pyseer --vcf snps.vcf.gz --phenotypes train.pheno --wg enet \
    --load-vars output/ma_snps --alpha 1 --save-model test_lasso --cpu 4 \
    --lineage-clusters lineages.txt --sequence-reweighting

    Read 499 phenotypes
    Detected binary phenotype
    Reading all variants
    Analysing 499 samples found in both phenotype and loaded npy
    Applying correlation filtering
    100%|████████████████████████████████████████████████████| 89703/89703 [00:56<00:00, 1597.01variants/s]
    Fitting elastic net to top 67277 variants
    [status]	Parallel glmnet cv with 4 cores
    Best penalty (lambda) from cross-validation: 3.38E-02
    Best model deviance from cross-validation: 0.605 ± 1.01E-01
    Best R^2 from cross-validation: 0.788
    Predictions within each lineage
    Lineage	Size	R2	TP	TN	FP	FN
    0	74	0.753	24	46	4	0
    1	41	0.219	2	35	0	4
    10	12	1.000	0	12	0	0
    11	9	1.000	8	1	0	0
    12	8	1.000	8	0	0	0
    13	11	1.000	11	0	0	0
    14	9	1.000	3	6	0	0
    15	9	1.000	0	9	0	0
    16	10	1.000	0	10	0	0
    17	7	-0.167	0	6	0	1
    18	6	1.000	0	6	0	0
    19	5	-0.250	0	4	0	1
    2	35	1.000	0	35	0	0
    20	3	1.000	3	0	0	0
    21	6	-0.200	0	5	0	1
    22	7	1.000	0	7	0	0
    23	6	1.000	0	6	0	0
    24	7	-0.167	0	6	0	1
    25	7	1.000	0	7	0	0
    26	6	-0.200	0	5	0	1
    27	5	1.000	0	5	0	0
    28	5	1.000	2	3	0	0
    29	5	1.000	5	0	0	0
    3	36	-0.059	34	0	2	0
    30	3	1.000	0	3	0	0
    31	4	-0.333	0	3	0	1
    32	4	1.000	0	4	0	0
    33	3	-0.500	0	2	0	1
    34	3	1.000	3	0	0	0
    35	3	1.000	0	3	0	0
    36	3	1.000	0	3	0	0
    37	3	1.000	3	0	0	0
    38	1	1.000	0	1	0	0
    39	1	1.000	0	1	0	0
    4	24	-0.043	23	0	1	0
    40	2	1.000	2	0	0	0
    41	2	1.000	0	2	0	0
    42	1	1.000	0	1	0	0
    43	2	1.000	0	2	0	0
    44	1	1.000	0	1	0	0
    45	1	1.000	0	1	0	0
    46	2	1.000	0	2	0	0
    47	1	1.000	1	0	0	0
    48	1	1.000	0	1	0	0
    49	1	1.000	0	1	0	0
    5	24	1.000	24	0	0	0
    50	1	1.000	0	1	0	0
    51	1	1.000	1	0	0	0
    52	1	1.000	0	1	0	0
    53	1	1.000	1	0	0	0
    54	1	1.000	1	0	0	0
    55	1	1.000	1	0	0	0
    56	1	1.000	1	0	0	0
    57	1	1.000	0	1	0	0
    58	1	1.000	0	1	0	0
    59	1	1.000	0	1	0	0
    6	18	-0.059	0	17	0	1
    7	18	-0.200	12	3	0	3
    8	18	-0.200	0	15	0	3
    9	16	1.000	0	16	0	0
    Finding and printing selected variants
    Saved enet model as test_lasso.pkl
    198248 loaded variants
    130971 filtered variants
    67277 tested variants
    32 printed variants

The prediction accuracy is pretty similar across lineages, which is good. As the
test set is a similar makeup of lineages hopefully prediction accuracy will be similar.

``enet_predict`` is used to make the predictions::

    enet_predict --vcf snps.vcf.gz --lineage-clusters lineages.txt --true-values test.pheno \
    test_lasso.pkl test.samples > test_predictions.txt

    Reading variants from input
    198248variants [00:11, 17657.99variants/s]
    Overall prediction accuracy
    R2: 0.8668373879641486
    tn: 69
    fp: 2
    fn: 1
    tp: 32
    Predictions within each lineage
    Lineage	Size	R2	TP	TN	FP	FN
    0	22	1.000	11	11	0	0
    1	14	-0.077	0	13	0	1
    10	3	1.000	0	3	0	0
    11	5	1.000	2	3	0	0
    12	5	-0.250	4	0	1	0
    13	1	1.000	1	0	0	0
    14	2	1.000	1	1	0	0
    15	2	1.000	0	2	0	0
    17	2	1.000	0	2	0	0
    18	2	1.000	0	2	0	0
    19	4	1.000	0	4	0	0
    2	11	1.000	0	11	0	0
    20	4	-0.333	3	0	1	0
    21	1	1.000	0	1	0	0
    23	1	1.000	0	1	0	0
    26	1	1.000	0	1	0	0
    27	1	1.000	0	1	0	0
    3	8	1.000	8	0	0	0
    30	1	1.000	0	1	0	0
    33	1	1.000	0	1	0	0
    39	1	1.000	0	1	0	0
    4	1	1.000	1	0	0	0
    42	1	1.000	0	1	0	0
    44	1	1.000	0	1	0	0
    45	1	1.000	0	1	0	0
    5	1	1.000	1	0	0	0
    6	2	1.000	0	2	0	0
    60	1	1.000	0	1	0	0
    61	1	1.000	0	1	0	0
    7	1	1.000	0	1	0	0
    9	2	1.000	0	2	0	0

The required options are a variant file, in this case the same ``--vcf`` contains
calls for the test samples, but this could be a new file, as long as the variant labels
match (non-trivial!). ``test_lasso.pkl`` is the saved model and ``test.samples`` are
the names of samples appearing in the variants file to produce predictions for.

Here, providing ``--true-values`` is needed to give the prediction accuracies. Providing
``--lineage-clusters`` in addition gives the per lineage prediction accuracy. For the reasons
noted above, the test accuracy is pretty similar to the training set.

The predictions are in ``test_predictions.txt``::

    Sample  Prediction      Link    Probability
    7622_3#79       1.0     1.1723387708055686      0.7635674993396665
    7622_3#80       1.0     2.828167499490956       0.9441790988402875
    7622_3#81       1.0     2.2308130622857987      0.9029826106893201
    7622_3#82       0.0     -0.7572524088945985     0.3192430949937001

For a binary phenotype:

* a 0/1 prediction at ``--threshold`` on the probability.
* Link is the value of the linear sum of the model betas, before entering the logit link function.
* Probability is a continuous prediction (after taking logit).

Generating consistent unitig calls
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
It has been mentioned many times above that it is necessary that variant calls match between
the inputs of the training and test data. This was ensured above as all variants were called
together and merged into a single file. Generally this may not be possible, especially if testing
prediction accuracy in a new cohort. If a variant in the model is missing its mean slope value will be
used for all samples, which may significantly reduce accuracy.

One way around this issue is to use unitigs. However, sequences which are unitigs in the DBG of
one population may not be unitigs in the DBG of a different sample set, even if they are present.
So simply running ``unitig-counter`` on both training and test datasets will result in many missing calls.

You should instead use `unitig-caller <https://github.com/johnlees/unitig-caller>`__ to make variant calls in
the test population using the same unitigs definitions as in the training population. Full usage and details
are given in the ``README.md``, but briefly::

    gzip -d -c unitigs.txt.gz | cut -f 1 > queries.txt
    unitig-caller --mode simple --strains strain_list.txt --unitigs queries.txt --output calls.txt

Will write a file of sequence elements for the samples in ``strain_list.txt`` to ``calls.txt``, which
is guaranteed to overlap with the original training set calls, and can therefore be used with ``enet_predict``.
