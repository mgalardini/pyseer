# Copyright 2017 Marco Galardini and John Lees
# Functions here are based on FaST-LMM Copyright 2014 Microsoft Corporation

'''LMM interface implementations'''

import os
import sys
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import math
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

import pyseer.classes as var_obj
from limix.heritability import estimate as estimate_h2
from limix.qtl import qtl_test_lmm, qtl_test_glmm

from .model import pre_filtering
from .model import fit_lineage_effect


def initialise_lmm(p, cov, K_in, continuous=False, lineage_samples=None):
    """Initialises LMM using the similarity matrix
    see _internal_single in fastlmm.association.single_snp

    Args:
        p (pandas.Series)
            Phenotypes vector (n, 1)
        cov (pandas.DataFrame)
            Covariance matrix (n, m)
        K_in (str)
            Similarity matrix filename
        continuous (bool)
            If the phenotype is continuous (otherwise binary)

            [default = False]
        lineage_samples (list or None)
            Sample names used for lineage (must match K_in)

    Returns:
        p (pandas.Series)
            Phenotype vector with the samples present in the similarity matrix
        cov (numpy.array)
            Covariates for use in LMM
        K (numpy.array)
            Kinship matrix for use in LMM
        h2 (float)
            Estimate of narrow-sense heritability
    """
    # read and normalise K
    K = pd.read_table(K_in,
                        index_col=0)
    K.index = K.index.astype(str)
    sys.stderr.write("Similarity matrix has dimension " + str(K.shape) + "\n")

    # If using lineages, check compatible with LMM
    if lineage_samples is not None and set(K.index) != set(lineage_samples):
        sys.stderr.write("Lineage file and similarity matrix contain different sets"
                            " of samples\n")
        sys.exit(1)

    intersecting_samples = p.index.intersection(K.index)
    sys.stderr.write("Analysing " + str(len(intersecting_samples)) + " samples"
                        " found in both phenotype and similarity matrix\n")
    p = p.loc[intersecting_samples]
    y = np.reshape(p.values, (-1, 1))
    K = K.loc[p.index, p.index]
    if cov.shape[0] == p.shape[0]:
        cov = cov.loc[intersecting_samples]
        covar = np.c_[cov.values, np.ones((p.shape[0], 1))]
    else:
        covar = np.ones((p.shape[0], 1))

    factor = float(len(p)) / np.diag(K.values).sum()
    if abs(factor-1.0) > 1e-15:
        K *= factor

    if continuous:
        model_type = 'normal'
    else:
        model_type = 'bernoulli'
    h2 = estimate_h2(pheno=y, lik=model_type, K=K, covs=covar, verbose=False)

    return(p, covar, K, h2)


def fit_lmm(covar, K, variants, variant_mat, lineage_effects,
            lineage_clusters, covariates, continuous,
            filter_pvalue, lrt_pvalue):
    """Fits LMM and returns LMM tuples for printing

    Args:
        covar (numpy.array)
            Covariates from initialise_lmm
        K (numpy.array)
            Kinship matrix from initialise_lmm
        variants (iterable)
            Tuples with variant object, phenotype vector and variant vector
            (pyseer.classes.LMM, numpy.array, numpy.array)
        variant_mat (numpy.array)
            Variants presence absence matrix (n, k)
        lineage_effects (bool)
            Whether to fit lineage effects
        lineage_clusters (numpy.array)
            Population structure matrix or lineage association
            binary matrix (n, k)
        covariates (numpy.array)
            Covariates matrix (n, m)
        continuous (bool)
            Whether the phenotype is continuous
        filter_pvalue (float)
            Pre-filtering p-value threshold
        lrt_pvalue (float)
            Post-fitting p-value threshold

    Returns:
        all_variants (iterable)
            All variant objects fitted or filtered
    """
    all_variants = []
    filtered_variants = []
    for var_idx, variant in enumerate(variants):
        notes = set()
        var, p, k = variant
        if var.pattern is None or k is None:
            notes.add('af-filter')
            all_variants.append(var._replace(notes=notes,
                                             prefilter=True,
                                             filter=False))
            variant_mat[:, var_idx] = np.zeros(variant_mat.shape[0])
            continue
        if not continuous:
            prep, bad_chisq = pre_filtering(p, k, continuous)
        else:
            prep = 0
            bad_chisq = False
        if bad_chisq:
            notes.add('bad-chisq')
        if prep >= filter_pvalue or not np.isfinite(prep):
            notes.add('pre-filtering-failed')
            all_variants.append(var._replace(notes=notes,
                                             prep=prep,
                                             prefilter=True,
                                             filter=False))
            variant_mat[:, var_idx] = np.zeros(variant_mat.shape[0])
            continue
        var = var._replace(prep=prep,
                           notes=notes,
                           prefilter=False)
        filtered_variants.append(var)

    # remove empty rows from filtering
    variant_mat = variant_mat[:, ~np.all(variant_mat == 0, axis=0)]

    if variant_mat.shape[1] == 0:
        return all_variants

    # fit LMM to block
    res = fit_lmm_block(p, covar, K, variant_mat, continuous)
    assert len(res['p_values']) == len(filtered_variants), "length of LMM result does not match number of variants"

    passed_vars = []
    for lmm_result_idx, tested_variant in zip(range(len(res['p_values'])),
                                              filtered_variants):
        notes = tested_variant.notes
        if res['p_values'][lmm_result_idx] >= lrt_pvalue or not np.isfinite(
                                            res['p_values'][lmm_result_idx]):
            notes.add('lrt-filtering-failed')
            all_variants.append(tested_variant._replace(
                                      notes=notes,
                                      pvalue=res['p_values'][lmm_result_idx],
                                      filter=True))
        else:
            if lineage_effects:
                max_lineage = fit_lineage_effect(lineage_clusters,
                                                 covariates, k)
            else:
                max_lineage = None

            tested_variant = tested_variant._replace(
                    pvalue=res['p_values'][lmm_result_idx],
                    kbeta=res['beta'][lmm_result_idx],
                    bse=res['bse'][lmm_result_idx],
                    notes=notes,
                    filter=False,
                    max_lineage=max_lineage)

            all_variants.append(tested_variant)

    return all_variants

def fit_lmm_block(p, covar, K, variant_block, continuous=False):
    """Actually fits the LMM to numpy variant array

    Args:
        p (numpy.array)
            Phenotype vector from initialise_lmm
        covar (numpy.array)
            Covariates from initialise_lmm
        K (numpy.array)
            Kinship matrix from initialise_lmm
        variant_block (numpy.array)
            Variants presence absence matrix (n, k)
        continuous (bool)
            Whether the phenotype is continuous or binary

            [default = False]

    Returns:
        lmm_results (dict)
            LMM results for this variants block
    """
    if continuous:
        res = qtl_test_lmm(variant_block, p, K=K, covs=covar, test='lrt', verbose=False)
    else:
        res = qtl_test_glmm(variant_block, p, 'bernoulli', K, covs=covar, test='lrt', verbose=False)

    lmm_results = {'beta': res.getBetaSNP()[0, :],
                   'bse': res.getBetaSNPste()[0, :],
                   'p_values': res.getPv()[0, :]
    }

    return(lmm_results)
