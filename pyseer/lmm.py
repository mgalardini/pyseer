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
from collections import namedtuple
import statsmodels.formula.api as smf

from fastlmm.lmm_cov import LMM as lmm_cov

# Initialises LMM using K matrix
# see _internal_single in fastlmm.association.single_snp
def initialise_lmm(p, cov, K_in, lmm_cache_in = None, lmm_cache_out = None):
    covar = np.c_[cov.values,np.ones((p.shape[0], 1))]
    y = p.values

    if lmm_cache is not None and os.path.exists(lmm_cache):
        lmm = lmm_cov(X=covar, Y=y, G=None, K=None)
        with np.load(lmm_cache) as data:
            lmm.U = data['arr_0']
            lmm.S = data['arr_1']
            h2 = data['arr_2']
    else:
        K = pd.read_table(K_in,
                  index_col=0)
        K = K.loc[p.index, p.index]
        lmm = lmm_cov(X=covar, Y=y, K=K.values, G=None, inplace=True)

        result = lmm.findH2()
        h2 = result['h2']

        if lmm_cache_out is not None and not os.path.exists(lmm_cache_out):
            lmm.getSU()
            np.savez(lmm_cache_out, lmm.U,lmm.S,np.array([h2]))

    return(lmm, h2)

# Fits LMM and returns LMM tuples for printing
def fit_lmm(lmm, variants, variant_mat, lineage_effects, lineage_clusters, covariates, lrt_pvalue):

    # fit LMM to block
    res = fit_lmm_block(lmm, variant_mat)

    passed_vars = []
    for lmm_result, tested_variant in zip(res, variants):
        if lmm_result['p_values'] < lrt_pvalue:
            tested_variant.pvalue = lmm_result['p_values']
            tested_variant.kbeta = lmm_result['beta']
            tested_variant.bse = lmm_result['bse']
            tested_variant.frac_h2 = lmm_result['frac_h2']

            if lineage_effects:
                tested_variant.max_lineage = lineage_effect_fit(lin, covariates, k)

            return_vars.append(tested_variant)

    return(passed_vars)

# Actually fits the LMM to numpy variant array
# see map/reduce section of _internal_single in fastlmm.association.single_snp
def fit_lmm_block(lmm, variant_block):
    res = lmm.nLLeval(h2=h2, dof=None, scale=1.0, penalty=0.0, snps=variant_block)

    beta = res['beta']
    chi2stats = beta*beta/res['variance_beta']

    lmm_results = {}
    lmm_results['p_values'] = stats.f.sf(chi2stats,1,lmm.U.shape[0]-(lmm.linreg.D+1))[:,0]
    lmm_results['beta'] = beta[:,0]
    lmm_results['bse'] = np.sqrt(res['variance_beta'][:,0])
    lmm_results['frac_h2'] = np.sqrt(res['fraction_variance_explained_beta'][:,0])

    return(lmm_results)
