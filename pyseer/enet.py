# Copyright 2018 Marco Galardini and John Lees

'''Elastic net model implementations'''

import os
import sys
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, vstack
import math
import pandas as pd
from decimal import Decimal

import glmnet_python
from cvglmnet import cvglmnet
from cvglmnetCoef import cvglmnetCoef

import pyseer.classes as var_obj
from .input import read_variant
from .model import pre_filtering
from .model import fit_lineage_effect

# Loads all variants into memory for use with elastic net
def load_all_vars(var_type, p, burden, burden_regions, infile,
                   all_strains, sample_order, min_af, max_af,
                   uncompressed):
    """Make in iterable to load blocks of variants for LMM

    Args:
        var_type (str)
            Variants type (one of: kmers, vcf or Rtab)
        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        burden (bool)
            Whether to slice a vcf file by burden regions
        burden_regions (collections.deque)
            Burden regions to slice the vcf with
        infile (opened file)
            Handle to opened variant file
        all_strains (set-like)
            All sample labels that should be present
        sample_order
            Sampes order to interpret each Rtab line
        min_af (float)
            Minimum allele frequency (inclusive)
        max_af (bool)
            maximum allele frequency (inclusive)
        uncompressed (bool)
            Whether the kmers file is uncompressed
        block_size (int)
            How many variants to be loaded at once

    Returns:
        variants (iterable)
            A collection of pyseer.classes.LMM objects describing the
            loaded variants (n,)
        variant_mat (numpy.array)
            Variant bloack presence/absence matrix (n, block_size)
        eof (bool)
            Whether we are at the end of the file
    """
    # For building sparse matrix
    data = []
    indices = []
    indptr = [0]
    selected_vars = []
    var_idx = 0

    # For correlation calculation
    correlations = []
    while True:
        eof, k, var_name, kstrains, nkstrains, af = read_variant(
                                        infile, p, var_type,
                                        burden, burden_regions,
                                        uncompressed, all_strains,
                                        sample_order)

        # check for EOF
        if eof:
            break

        if k is not None and af > min_af and af < max_af:
            # Calculate correlation
            cor_a = k - np.mean(k)
            correlations.append(cor_a)

            # Minor allele encoding - most efficient use of sparse structure
            if af > 0.5:
                pres = 0
            else:
                pres = 1

            for idx, obs in enumerate(k):
                if obs == pres:
                    indices.append(idx)
                    data.append(1)
            indptr.append(len(indices))

            selected_vars.append(var_idx)

        var_idx += 1

    # construct sparse matrix, then filter out correlations
    variants = csr_matrix((data, indices, indptr), dtype=float)

    return(variants, selected_vars, var_idx, correlations)

def fit_enet(p, variants, covariates, continuous, alpha, n_folds = 10, n_cpus = 1):
    if continuous:
        regression_type = 'gaussian'
    else:
        regression_type = 'binomial'

    if covariates.shape[0] > 0:
        variants = vstack(csc_matrix(covariates), variants)

    enet_fit = cvglmnet(x = variants, y = p.values.astype('float64'), family = regression_type,
                        nfolds = n_folds, alpha = alpha, parallel = n_cpus)
    betas = cvglmnetCoef(enet_fit, s = 'lambda_min')
    best_lambda_idx = np.argmin(enet_fit['cvm'])

    # R^2 = 1 - sum((yi_obs - yi_predicted)^2) /sum((yi_obs - yi_mean)^2)
    RSS = np.sum(np.square(p.values - np.mean(p.values)))
    R2 = 1 - (enet_fit['cvm'][best_lambda_idx]/RSS)
    R2_err = enet_fit['cvsd'][best_lambda_idx]/RSS
    sys.stderr.write("Best penalty from cross-validation: " + '%.2E' % Decimal(enet_fit['lambda_min'][0]) + "\n")
    sys.stderr.write("Best R^2 from cross-validation: " + '%.3f' % Decimal(R2) + " ± " + '%.2E' % Decimal(R2_err) + "\n")

    return(betas.reshape(-1,))

def correlation_filter(p, cor_a, quantile_filter = 0.25):

    b = p.values - np.mean(p.values)
    sum_b_squared = np.sum(np.power(b, 2))

    correlations = []
    for a in cor_a:
        cor = np.abs(np.dot(a, b) / np.sqrt(np.sum(np.power(a, 2)) * sum_b_squared))
        correlations.append(cor)

    cor_filter = np.nonzero(correlations > np.percentile(correlations, quantile_filter*100))[0]
    return(cor_filter)


def find_enet_selected(enet_betas, var_indices, p, c, var_type, burden,
                       burden_regions, infile, all_strains, sample_order,
                       continuous, find_lineage, lin, uncompressed):

    # skip covariates
    if c.shape[1] > 0:
        enet_betas = enet_betas[c.shape[1]:,:]
        covar_betas = enet_betas[0:c.shape[1],:]
        for beta, covariate in zip(covar_betas, c):
            sys.stderr.write("Kept covariate '" + covariate + "', slope: " + '%.2E' % Decimal(beta) + "\n")

    current_var = 0
    for beta, var_idx in zip(enet_betas, var_indices):
        # Only need to process selected variants
        if beta == 0:
            continue
        while current_var < var_idx:
            eof, k, var_name, kstrains, nkstrains, af = read_variant(
                                        infile, p, var_type,
                                        burden, burden_regions,
                                        uncompressed, all_strains,
                                        sample_order, noparse=True)
            current_var += 1

        eof, k, var_name, kstrains, nkstrains, af = read_variant(
                                        infile, p, var_type,
                                        burden, burden_regions,
                                        uncompressed, all_strains,
                                        sample_order)
        current_var += 1
        if af > 0.5:
            beta *= -1

        notes = []

        # find unadjusted p-value and lineage
        (pval, bad) = pre_filtering(p, k, continuous)
        if bad:
            notes.append("bad-chisq")
        if find_lineage:
            max_lineage = fit_lineage_effect(lin, c, k)
        else:
            max_lineage = None

        yield var_obj.Enet(var_name, af, pval, beta, max_lineage, kstrains, nkstrains, notes)
