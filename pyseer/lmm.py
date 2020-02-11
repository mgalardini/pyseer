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
from .fastlmm.lmm_cov import LMM as lmm_cov

from .model import pre_filtering
from .model import fit_lineage_effect


def initialise_lmm(p, cov, K_in, lmm_cache_in=None, lmm_cache_out=None, lineage_samples=None):
    """Initialises LMM using the similarity matrix
    see _internal_single in fastlmm.association.single_snp

    Args:
        p (pandas.Series)
            Phenotypes vector (n, 1)
        cov (pandas.DataFrame)
            Covariance matrix (n, m)
        K_in (str)
            Similarity matrix filename
        lmm_chache_in (str or None)
            Filename for an input LMM cache, None if it has to be computed
        lmm_chache_out (str or None)
            Filename to save the LMM cache, None otherwise.
        lineage_samples (list or None)
            Sample names used for lineage (must match K_in)

    Returns:
        p (pandas.Series)
            Phenotype vector with the samples present in the similarity matrix
        lmm (pyseer.fastlmm.lmm_cov.LMM)
            Initialised LMM model
        h2 (float)
            Trait's variance explained by covariates
    """
    if lmm_cache_in is not None and os.path.exists(lmm_cache_in):
        if len(p.index.intersection(cov.index)) == p.shape[0]:
            covar = np.c_[cov.loc[p.index].values, np.ones((p.shape[0], 1))]
        # empty covariate
        elif (cov.shape[0] == 0 and cov.shape[1] == 0) or len(cov.shape) == 0:
            covar = np.ones((p.shape[0], 1))
        # inform the user about missing covariate/phenotypes and exit
        # special case for LMM
        else:
            sys.stderr.write("Phenotype and covariate file should have "
                             "matching samples for LMM\n")
            sys.exit(1)
        y = np.reshape(p.values, (-1, 1))

        lmm = lmm_cov(X=covar, Y=y, G=None, K=None)
        with np.load(lmm_cache_in) as data:
            lmm.U = data['arr_0']
            lmm.S = data['arr_1']
            h2 = data['arr_2'][0]

            if (lmm.U.shape[0] != len(p)):
                sys.stderr.write("Phenotype different length from cache file\n")
                sys.exit(1)
    else:
        # read and normalise K
        K = pd.read_csv(K_in,
                        index_col=0,
                        sep='\t')
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
        if len(p.index.intersection(cov.index)) == p.shape[0]:
            covar = np.c_[cov.loc[p.index].values, np.ones((p.shape[0], 1))]
        # empty covariate
        elif (cov.shape[0] == 0 and cov.shape[1] == 0) or len(cov.shape) == 0:
            covar = np.ones((p.shape[0], 1))
        # inform the user about missing covariate/phenotypes and exit
        # special case for LMM
        else:
            sys.stderr.write("Phenotype and covariate file should have "
                             "matching samples for LMM\n")
            sys.exit(1)

        factor = float(len(p)) / np.diag(K.values).sum()
        if factor == math.inf:
            sys.stderr.write("Invalid similarity matrix. Did you use --calc-C?\n")
            sys.exit(1)
        elif abs(factor-1.0) > 1e-15:
            K *= factor

        lmm = lmm_cov(X=covar, Y=y, K=K.values, G=None, inplace=True)
        result = lmm.findH2()
        h2 = result['h2']

        if lmm_cache_out is not None and not os.path.exists(lmm_cache_out):
            lmm.getSU()
            np.savez(lmm_cache_out, lmm.U, lmm.S, np.array([h2]))

    return(p, lmm, h2)


def fit_lmm(lmm, h2, variants, variant_mat, lineage_effects,
            lineage_clusters, covariates, continuous,
            filter_pvalue, lrt_pvalue):
    """Fits LMM and returns LMM tuples for printing

    Args:
        lmm (pyseer.fastlmm.lmm_cov.LMM)
            Initialised LMM model
        h2 (float)
            Trait's variance explained by covariates
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

        prep, bad_chisq = pre_filtering(p, k, continuous)
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
    res = fit_lmm_block(lmm, h2, variant_mat)
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
                    frac_h2=res['frac_h2'][lmm_result_idx],
                    notes=notes,
                    filter=False,
                    max_lineage=max_lineage)

            all_variants.append(tested_variant)

    return all_variants

def fit_lmm_block(lmm, h2, variant_block):
    """Actually fits the LMM to numpy variant array
    see map/reduce section of _internal_single in fastlmm.association.single_snp

    Args:
        lmm (pyseer.fastlmm.lmm_cov.LMM)
            Initialised LMM model
        h2 (float)
            Trait's variance explained by covariates
        variant_block (numpy.array)
            Variants presence absence matrix (n, k)

    Returns:
        lmm_results (dict)
            LMM results for this variants block
    """
    res = lmm.nLLeval(h2=h2, dof=None, scale=1.0,
                      penalty=0.0, snps=variant_block)

    beta = res['beta']
    chi2stats = beta*beta/res['variance_beta']

    lmm_results = {}
    lmm_results['p_values'] = stats.f.sf(chi2stats,
                                         1,
                                         lmm.U.shape[0]-(lmm.linreg.D+1))[:, 0]
    lmm_results['beta'] = beta[:, 0]
    lmm_results['bse'] = np.sqrt(res['variance_beta'][:, 0])
    lmm_results['frac_h2'] = np.sqrt(
                        res['fraction_variance_explained_beta'][:, 0]
                                    )

    return(lmm_results)
