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
from scipy.sparse import csr_matrix, csc_matrix, hstack
import math
import pandas as pd
from decimal import Decimal
from tqdm import tqdm

import glmnet_python
from cvglmnet import cvglmnet
from cvglmnetCoef import cvglmnetCoef
from cvglmnetPredict import cvglmnetPredict

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

    pbar = tqdm(unit="variants")
    while True:
        eof, k, var_name, kstrains, nkstrains, af = read_variant(
                                        infile, p, var_type,
                                        burden, burden_regions,
                                        uncompressed, all_strains,
                                        sample_order)

        # check for EOF
        if eof:
            pbar.close()
            break
        else:
            pbar.update(1)

        if k is not None and af > min_af and af < max_af:
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

    # construct sparse matrix
    variants = csr_matrix((data, indices, indptr), dtype=float)

    return(variants, selected_vars, var_idx)

def fit_enet(p, variants, covariates, continuous, alpha, n_folds = 10, n_cpus = 1):
    if continuous:
        regression_type = 'gaussian'
    else:
        regression_type = 'binomial'

    if covariates.shape[0] > 0:
        variants = hstack([csc_matrix(covariates.values), variants])

    # Run model fit
    enet_fit = cvglmnet(x = variants, y = p.values.astype('float64'), family = regression_type,
                        nfolds = n_folds, alpha = alpha, parallel = n_cpus)
    betas = cvglmnetCoef(enet_fit, s = 'lambda_min')

    # Extract best lambda and predict class labels/values
    best_lambda_idx = np.argmin(enet_fit['cvm'])
    if continuous:
        preds = cvglmnetPredict(enet_fit, newx=variants, s='lambda_min', ptype='link')
    else:
        preds = cvglmnetPredict(enet_fit, newx=variants, s='lambda_min', ptype='class')

    # Write some summary stats
    # R^2 = 1 - sum((yi_obs - yi_predicted)^2) /sum((yi_obs - yi_mean)^2)
    SStot = np.sum(np.square(p.values - np.mean(p.values)))
    SSerr = np.sum(np.square(p.values.reshape(-1, 1) - preds))
    R2 = 1 - (SSerr/SStot)
    sys.stderr.write("Best penalty (lambda) from cross-validation: " + '%.2E' % Decimal(enet_fit['lambda_min'][0]) + "\n")
    if not continuous:
        sys.stderr.write("Best model deviance from cross-validation: " + '%.3f' % Decimal(enet_fit['cvm'][best_lambda_idx]) +
                         " ± " + '%.2E' % Decimal(enet_fit['cvsd'][best_lambda_idx]) + "\n")
    sys.stderr.write("Best R^2 from cross-validation: " + '%.3f' % Decimal(R2) + "\n")

    return(betas.reshape(-1,))

def correlation_filter(p, all_vars, quantile_filter = 0.25):

    # a = snp - mean(snp)
    # b = y - mean(y)
    # cor = abs(a%*%b / sqrt(sum(a^2)*sum(b^2)) )
    b = p.values - np.mean(p.values)
    sum_b_squared = np.sum(np.power(b, 2))

    #TODO multithread
    correlations = []
    for row_idx in tqdm(range(all_vars.shape[0]), unit="variants"):
        k = all_vars.getrow(row_idx)
        k_mean = csr_matrix.mean(k)

        ab = k.dot(b) - np.sum(k_mean * b)
        sum_a_squared = k.dot(k.transpose()).data[0] - 2*k_mean*csr_matrix.sum(k) + pow(k_mean, 2) * all_vars.shape[1]
        cor = np.abs(ab / np.sqrt(sum_a_squared * sum_b_squared))
        correlations.append(cor)

    cor_filter = np.nonzero(correlations > np.percentile(correlations, quantile_filter*100))[0]
    return(cor_filter)


def find_enet_selected(enet_betas, var_indices, p, c, var_type, burden,
                       burden_regions, infile, all_strains, sample_order,
                       continuous, find_lineage, lin, uncompressed):

    # skip intercept and covariates
    enet_betas = enet_betas[c.shape[1]+1:]

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