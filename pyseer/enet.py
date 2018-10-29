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
from .model import fixed_effects_regression

# Loads all variants into memory for use with elastic net
def load_all_vars(var_type, p, burden, burden_regions, infile,
                   all_strains, sample_order, min_af, max_af,
                   uncompressed):
    """Load all variants in the input file into a sparse
    matrix representation

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
        variants (scipy.sparse.csr_matrix)
            A sparse matrix representation of all variants in
            the input
        selected_vars (list)
            0-Indices of variants in the input file in variants
            (which passed AF filtering)
        var_idx (int)
            The number of variants passing AF filtering (number
            of rows of variants)
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
    """Fit an elastic net model to a set of variants. Prints
    information about model fit and prediction quality to STDERR

    Args:

        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        variants (scipy.sparse.csc_matrix)
            Wide sparse matrix representation of all variants to fit to
            (rows = samples, columns = variants)
        covariates (pandas.DataFrame)
            Covariate matrix (n, j)
        continuous (bool)
            If True fit a Gaussian error model, otherwise Bionomial error
        alpha (float)
            Between 0-1, sets the mix between ridge regression and lasso
            regression
        n_folds (int)
            Number of folds in cross-validation

            [default = 10]
        n_cpus (int)
            Number of processes to use in cross-validation
            Set to -1 to use all available

            [default = 1]

    Returns:
        betas (numpy.array)
            The fitted betas (slopes) for each variant
    """
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
    """Calculates correlations between phenotype and variants,
    giving those that are above the specified quantile

    Args:
        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        all_vars (scipy.sparse.csr_matrix)
            Narrow sparse matrix representation of all variants to fit to
            (rows = variants, columns = samples)
        quantile_filter (float)
            The quantile to discard at e.g. 0.25, retain top 75%

            [default = 0.25]

    Returns:
        cor_filter (numpy.array)
            The indices of variants passing the filter
    """
    # a = snp - mean(snp)
    # b = y - mean(y)
    # cor = abs(a%*%b / sqrt(sum(a^2)*sum(b^2)) )
    b = p.values - np.mean(p.values)
    sum_b_squared = np.sum(np.power(b, 2))

    # NOTE: I couldn't get this to multithread efficiently using sparse matrices...
    # might work if the matrix was divided into chunks of rows first, but maybe not
    # worth it as it's pretty quick anyway
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


def find_enet_selected(enet_betas, var_indices, p, c, var_type, fit_seer, burden,
                       burden_regions, infile, all_strains, sample_order,
                       continuous, find_lineage, lin, uncompressed):
    """Read through the variant input file again, yielding just those variants
    which had a non-zero slope for printing

    Args:
        enet_betas (numpy.array)
            Fitted slopes of intercept, covariants and variants
            from elastic net
        var_indices (list)
            The 0-indexed locations (in the original file) of
            variants represented in enet_betas
        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        c (numpy.array)
            Covariate matrix (n, j)
        var_type (str)
            Variants type (one of: kmers, vcf or Rtab)
        fit_seer (tuple: m, null_model, null_firth)
            Distance projection and null models required to
            fit fixed effect regression
        burden (bool)
            Whether to slice a vcf file by burden regions
        burden_regions (collections.deque)
            Burden regions to slice the vcf with
        infile (opened file)
            Handle to opened variant file
        all_strains (set-like)
            All sample labels that should be present
        sample_order
            Sample order to interpret each Rtab line
        continuous (bool)
            Is phenotype/fit continuous?
        lineage_effects (bool)
            Whether to fit lineages or not
        lin (numpy.array)
            Lineages matrix (n, k)
        uncompressed (bool)
            Whether the kmers file is uncompressed

    Returns:
        variant (var_obj.Enet)
            Iterable of processed variants for printing
    """
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

        # find pvalues and lineages
        if fit_seer != None:
            m, null_res, null_firth = fit_seer
            seer_fit = fixed_effects_regression(var_name, p.values, k, m, c, af, None,
                             find_lineage, lin,
                             1, 1, null_res, null_firth,
                             kstrains, nkstrains, continuous)
            pval = seer_fit.prep
            adj_pval = seer_fit.pvalue
            max_lineage = seer_fit.max_lineage
            notes = seer_fit.notes
        # find just unadjusted p-value and lineage
        else:
            (pval, bad) = pre_filtering(p.values, k, continuous)
            adj_pval = math.nan
            if bad:
                notes.append("bad-chisq")
            if find_lineage:
                max_lineage = fit_lineage_effect(lin, c, k)
            else:
                max_lineage = None

        yield var_obj.Enet(var_name, af, pval, adj_pval, beta, max_lineage, kstrains, nkstrains, notes)
