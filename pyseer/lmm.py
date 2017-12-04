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
from fastlmm.lmm_cov import LMM as lmm_cov

# Initialises LMM using K matrix
# see _internal_single in fastlmm.association.single_snp
def initialise_lmm(p, cov, K_in, lmm_cache_in = None, lmm_cache_out = None):
    if cov.shape[0] == p.shape[0]:
        covar = np.c_[cov.values,np.ones((p.shape[0], 1))]
    else:
        covar = np.ones((p.shape[0], 1))
    y = np.reshape(p.values, (-1,1))

    if lmm_cache_in is not None and os.path.exists(lmm_cache_in):
        lmm = lmm_cov(X=covar, Y=y, G=None, K=None)
        with np.load(lmm_cache) as data:
            lmm.U = data['arr_0']
            lmm.S = data['arr_1']
            h2 = data['arr_2']
    else:
        # read and normalise K
        K = pd.read_table(K_in,
                  index_col=0)
        K = K.loc[p.index, p.index]
        factor = float(len(p)) / np.diag(K.values).sum()
        if abs(factor-1.0)>1e-15:
            K.values *= factor

        lmm = lmm_cov(X=covar, Y=y, K=K.values, G=None, inplace=True)
        result = lmm.findH2()
        h2 = result['h2']

        if lmm_cache_out is not None and not os.path.exists(lmm_cache_out):
            lmm.getSU()
            np.savez(lmm_cache_out, lmm.U,lmm.S,np.array([h2]))

    return(lmm, h2)

# Fits LMM and returns LMM tuples for printing
def fit_lmm(lmm, h2, variants, variant_mat, lineage_effects, lineage_clusters, covariates, lrt_pvalue):

    # fit LMM to block
    res = fit_lmm_block(lmm, h2, variant_mat)
    assert len(res['p_values']) == len(variants), "length of LMM result does not match number of variants"

    passed_vars = []
    for lmm_result_idx, tested_variant in zip(range(len(res)), variants):
        if res['p_values'][lmm_result_idx] < lrt_pvalue:
            if lineage_effects:
                max_lineage = lineage_effect_fit(lin, covariates, k)
            else:
                max_lineage = None

            # Copy variant. Might be inefficient?
            return_var = var_obj.LMM(tested_variant.kmer,
                             tested_variant.af,
                             tested_variant.prep,
                             res['p_values'][lmm_result_idx],
                             res['beta'][lmm_result_idx],
                             res['bse'][lmm_result_idx],
                             res['frac_h2'][lmm_result_idx],
                             max_lineage,
                             tested_variant.kstrains,
                             tested_variant.nkstrains)

            passed_vars.append(return_var)

    return(passed_vars)

# Actually fits the LMM to numpy variant array
# see map/reduce section of _internal_single in fastlmm.association.single_snp
def fit_lmm_block(lmm, h2, variant_block):
    res = lmm.nLLeval(h2=h2, dof=None, scale=1.0, penalty=0.0, snps=variant_block)

    beta = res['beta']
    chi2stats = beta*beta/res['variance_beta']

    lmm_results = {}
    lmm_results['p_values'] = stats.f.sf(chi2stats,1,lmm.U.shape[0]-(lmm.linreg.D+1))[:,0]
    lmm_results['beta'] = beta[:,0]
    lmm_results['bse'] = np.sqrt(res['variance_beta'][:,0])
    lmm_results['frac_h2'] = np.sqrt(res['fraction_variance_explained_beta'][:,0])

    return(lmm_results)
