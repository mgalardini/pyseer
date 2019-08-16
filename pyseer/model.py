# Copyright 2017 Marco Galardini and John Lees

'''Original SEER model (fixed effects) implementations'''

import os
import sys
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import math
import statsmodels
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
# handle different versions of statsmodels
try:
    smf.OLS
except AttributeError:
    smf.OLS = statsmodels.regression.linear_model.OLS
try:
    smf.Logit
except AttributeError:
    smf.Logit = statsmodels.discrete.discrete_model.Logit

import pyseer.classes as var_obj


def pre_filtering(p, k, continuous):
    """Calculate a naive p-value from a chisq test (binary phenotype)
    or a t-test (continuous phenotype) which is not adjusted for population
    structure

    Args:
        p (numpy.array)
            Phenotypes vector (n, 1)
        k (numpy.array)
            Variant presence-absence vector (n, 1)
        continous (bool)
            Whether phenotypes are continuous or binary

    Returns:
        prep (float)
            Naive p-value
        bad_chisq (boolean)
            Whether the chisq test had small values in the
            contingency table
    """
    bad_chisq = False
    if continuous:
        prep = stats.ttest_ind(p[k == 1],
                               p[k == 0],
                               equal_var=False)[1]
    else:
        t = np.concatenate((p.reshape(-1, 1), k.reshape(-1, 1)), axis=1).T
        table = [[t[0][(t[0] == 1) & (t[1] == 1)].shape[0],
                  t[0][(t[0] == 1) & (t[1] == 0)].shape[0]],
                 [t[0][(t[0] == 0) & (t[1] == 1)].shape[0],
                  t[0][(t[0] == 0) & (t[1] == 0)].shape[0]]]

        # check for small values
        table = np.array(table)
        if table[table <= 1].shape[0] > 0 or table[table <= 5].shape[0] > 1:
            bad_chisq = True

        prep = stats.chi2_contingency(table, correction=False)[1]

    return(prep, bad_chisq)


def fit_null(p, m, cov, continuous, firth=False):
    """Fit the null model i.e. regression without k-mer

    `y ~ Wa`

    Returns log-likelihood

    Args:
        p (numpy.array)
            Phenotypes vector (n, 1)
        m (numpy.array)
            Population structure matrix (n, k)
        cov (pandas.DataFrame)
            Covariants dataframe (n, j)
        continous (bool)
            Whether phenotypes are continuous or binary
        firth (bool)
            For binary phenotypes whether to use firth regression

    Returns:
        null_res (statsmodels.regression.linear_model.RegressionResultsWrapper or float or None)
            Fitted model or log-likelihood (if firth) or
            None if could not fit
    """
    v = np.ones(p.shape[0]).reshape(-1, 1)
    if m.shape[1] > 0:
        v = np.concatenate((v, m), axis=1)
    if cov.shape[1] > 0:
        v = np.concatenate((v, cov.values), axis=1)

    if continuous:
        null_mod = smf.OLS(p, v)
    else:
        start_vec = np.zeros(v.shape[1])
        start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))
        null_mod = smf.Logit(p, v)

    try:
        if continuous:
            null_res = null_mod.fit(disp=False)
        else:
            if firth:
                firth_res = fit_firth(null_mod, start_vec, v, p)
                if firth_res is None:
                    sys.stderr.write('Firth regression did not converge for null model\n')
                    return None
                (intercept, kbeta, beta, bse, fitll) = firth_res
                null_res = fitll
            else:
                try:
                    null_res = null_mod.fit(start_params=start_vec,
                                            method='newton',
                                            disp=False)
                # Null fit with default optimiser may fail, Powell
                # optimizer might work
                except np.linalg.linalg.LinAlgError:
                    null_res = null_mod.fit(start_params=start_vec,
                                            method='powell',
                                            disp=False)
    except np.linalg.linalg.LinAlgError:
        sys.stderr.write('Matrix inversion error for null model\n')
        return None
    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
        sys.stderr.write('Perfectly separable data error for null model\n')
        return None

    return null_res


def fit_lineage_effect(lin, c, k):
    """Fits the model `k ~ Wa` using binomial error with logit link.
    W are the lineages (either a projection of samples, or cluster indicators)
    and covariates.
    Returns the index of the most significant lineage

    Args:
        lin (numpy.array)
            Population structure matrix or lineage association
            binary matrix (n, k)
        c (numpy.array)
            Covariants matrix (n, j)
        k (numpy.array)
            Variant presence-absence vector (n, 1)

    Returns:
        max_lineage (int or None)
            Index of the most significant lineage
            or None is could not fit
    """
    if c.shape[0] == lin.shape[0]:
        X = np.concatenate((np.ones(lin.shape[0]).reshape(-1, 1),
                            lin,
                            c),
                           axis=1)
    else:
        X = np.concatenate((np.ones(lin.shape[0]).reshape(-1, 1),
                            lin),
                           axis=1)

    lineage_mod = smf.Logit(k, X)
    try:
        lineage_res = lineage_mod.fit(method='newton', disp=False)

        wald_test = np.divide(np.absolute(lineage_res.params), lineage_res.bse)
        # excluding intercept and covariates
        max_lineage = np.argmax(wald_test[1:lin.shape[1]+1])
    # In case regression fails
    except (statsmodels.tools.sm_exceptions.PerfectSeparationError,
            np.linalg.LinAlgError):
        max_lineage = None

    return max_lineage


def fixed_effects_regression(variant, p, k, m, c, af, pattern,
                             lineage_effects, lin,
                             pret, lrtt, null_res, null_firth,
                             kstrains, nkstrains, continuous):
    """Fits the model `y ~ Xb + Wa` using either binomial error with
    logit link (binary traits) or Gaussian error (continuous traits)

    * `y` is the phenotype
    * `X` is the variant presence/absence (fixed effects)
    * `W` are covariate fixed effects, including population structure
    * `a` and `b` are slopes to be fitted

    Args:
        variant (str)
            Variant identifier
        p (numpy.array)
            Phenotype vector (binary or continuous) (n, 1)
        k (numpy.array)
            Variant presence/absence vector (n, 1)
        m (numpy.array)
            Population structure matrix (n, m)
        c (numpy.array)
            Covariants matrix (n, j)
        af (float)
            Allele frequency
        pattern (str)
            Variant hashed pattern
        lineage_effects (bool)
            Whether to fit lineages or not
        lin (numpy.array)
            Lineages matrix (n, k)
        pret (float)
            Pre-filtering p-value threshold
        lrtt (float)
            Post-fitting p-value threshold
        null_res (float or statsmodels.regression.linear_model.RegressionResultsWrapper)
            Null-fit likelihood (binary) or model (continuous)
        null_firth (float)
            Firth regression likelihood
        kstrains (iterable)
            Sample labels with the variant
        nkstrains (iterable)
            Sample labels without the variant
        continuous (bool)
            Whether the phenotype is continuous or not

    Returns:
        result (pyseer.classes.Seer)
            Results container
    """
    notes = set()

    # was this af-filtered?
    if p is None:
        notes.add('af-filter')
        return var_obj.Seer(variant, pattern, af, np.nan, np.nan,
                            np.nan, np.nan, np.nan, np.array([]),
                            None, kstrains, nkstrains,
                            notes, True, False)

    # pre-filtering
    prep, bad_chisq = pre_filtering(p, k, continuous)
    if bad_chisq:
        notes.add('bad-chisq')
    if prep > pret or not np.isfinite(prep):
        notes.add('pre-filtering-failed')
        return var_obj.Seer(variant, pattern, af, prep, np.nan,
                            np.nan, np.nan, np.nan, np.array([]),
                            None, kstrains, nkstrains,
                            notes, True, False)

    # actual regression
    if m.shape[0] != k.shape[0]:
        # no distances
        if c.shape[0] == k.shape[0]:
            v = np.concatenate((np.ones(p.shape[0]).reshape(-1, 1),
                                k.reshape(-1, 1),
                                c),
                               axis=1)
        else:
            v = np.concatenate((np.ones(p.shape[0]).reshape(-1, 1),
                                k.reshape(-1, 1)),
                               axis=1)
    elif c.shape[0] == m.shape[0]:
        # covariates and distances
        v = np.concatenate((np.ones(m.shape[0]).reshape(-1, 1),
                            k.reshape(-1, 1),
                            m,
                            c),
                           axis=1)
    else:
        # no covariates
        v = np.concatenate((np.ones(m.shape[0]).reshape(-1, 1),
                            k.reshape(-1, 1),
                            m),
                           axis=1)
    try:
        if continuous:
            mod = smf.OLS(p, v)

            res = mod.fit()
            intercept = res.params[0]
            kbeta = res.params[1]
            beta = res.params[2:]
            bse = res.bse[1]
            lrt_pvalue = res.pvalues[1]
            #lrt_pvalue = res.compare_lr_test(null_res)[1]

        else:
            mod = smf.Logit(p, v)

            start_vec = np.zeros(v.shape[1])
            start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))

            if not bad_chisq:
                try:
                    res = mod.fit(start_params=start_vec,
                                  method='newton',
                                  disp=False)

                    if res.bse[1] > 3:
                        bad_chisq = True
                        notes.add('high-bse')
                    else:
                        lrstat = -2*(null_res - res.llf)
                        lrt_pvalue = 1
                        if lrstat > 0:  # non-convergence
                            lrt_pvalue = stats.chi2.sf(lrstat, 1)

                        intercept = res.params[0]
                        kbeta = res.params[1]
                        beta = res.params[2:]
                        bse = res.bse[1]
                except statsmodels.tools.sm_exceptions.PerfectSeparationError:
                    bad_chisq = True
                    notes.add('perfectly-separable-data')

            # Fit Firth regression with large SE, or nearly separable values
            if bad_chisq:
                firth_fit = fit_firth(mod, start_vec, v, p)
                if firth_fit is None:  # Firth failure
                    notes.add('firth-fail')
                    return var_obj.Seer(variant, pattern, af, prep, np.nan,
                                        np.nan, np.nan, np.nan, np.array([]),
                                        None, kstrains, nkstrains,
                                        notes, False, True)
                else:
                    intercept, kbeta, beta, bse, fitll = firth_fit
                    beta = np.array(beta)
                    lrstat = -2*(null_firth - fitll)
                    lrt_pvalue = 1
                    if lrstat > 0:  # check for non-convergence
                        lrt_pvalue = stats.chi2.sf(lrstat, 1)

    except np.linalg.linalg.LinAlgError:
        # singular matrix error
        notes.add('matrix-inversion-error')
        return var_obj.Seer(variant, pattern, af, prep, np.nan,
                            np.nan, np.nan, np.nan, np.array([]),
                            None, kstrains, nkstrains,
                            notes, False, True)

    if lineage_effects:
        max_lineage = fit_lineage_effect(lin, c, k)
    else:
        max_lineage = None

    if lrt_pvalue > lrtt or not np.isfinite(lrt_pvalue) or not np.isfinite(kbeta):
        notes.add('lrt-filtering-failed')
        return var_obj.Seer(variant, pattern, af, prep, lrt_pvalue,
                            kbeta, bse, intercept, beta,
                            max_lineage, kstrains, nkstrains,
                            notes, False, True)

    return var_obj.Seer(variant, pattern, af, prep, lrt_pvalue,
                        kbeta, bse, intercept, beta,
                        max_lineage, kstrains, nkstrains,
                        notes, False, False)


def firth_likelihood(beta, logit):
    """Convenience function to calculate likelihood of Firth regression

    Args:
        beta (numpy.array)
            (n, 1)
        logit (statsmodels.discrete.discrete_model.Logit)
            Logistic model

    Returns:
        likelihood (float)
            Firth likelihood
    """
    return -(logit.loglike(beta) +
             0.5*np.log(np.linalg.det(-logit.hessian(beta))))


def fit_firth(logit_model, start_vec, X, y,
              step_limit=1000, convergence_limit=0.0001):
    """Do firth regression

    Args:
        logit (statsmodels.discrete.discrete_model.Logit)
            Logistic model
        start_vec (numpy.array)
            Pre-initialized vector to speed-up convergence (n, 1)
        X (numpy.array)
            (n, m)
        y (numpy.array)
            (n, )
        step_limit (int)
            Maximum number of iterations
        convergence_limit (float)
            Convergence tolerance

    Returns:
        intercept (float)
           Intercept
        kbeta (float)
            Variant beta
        beta (iterable)
            Covariates betas (n-2)
        bse (float)
            Beta std-err
        fitll (float or None)
            Likelihood of fit or None if could not fit
    """

    beta_iterations = []
    beta_iterations.append(start_vec)
    for i in range(0, step_limit):
        pi = logit_model.predict(beta_iterations[i])
        W = np.diagflat(np.multiply(pi, 1-pi))
        var_covar_mat = np.linalg.pinv(
                        -logit_model.hessian(beta_iterations[i])
                        )

        # build hat matrix
        rootW = np.sqrt(W)
        H = np.dot(np.transpose(X), np.transpose(rootW))
        H = np.matmul(var_covar_mat, H)
        H = np.matmul(np.dot(rootW, X), H)

        # penalised score
        U = np.matmul(np.transpose(X),
                      y - pi + np.multiply(np.diagonal(H), 0.5 - pi))
        new_beta = beta_iterations[i] + np.matmul(var_covar_mat, U)

        # step halving
        j = 0
        while firth_likelihood(new_beta, logit_model) > firth_likelihood(
                                                        beta_iterations[i],
                                                        logit_model
                                                                        ):
            new_beta = beta_iterations[i] + 0.5*(new_beta - beta_iterations[i])
            j = j + 1
            if (j > step_limit):
                return None

        beta_iterations.append(new_beta)
        if i > 0 and (np.linalg.norm(beta_iterations[i] -
                      beta_iterations[i-1]) < convergence_limit):
            break

    return_fit = None
    if np.linalg.norm(beta_iterations[i] -
                      beta_iterations[i-1]) >= convergence_limit:
        pass
    else:
        # Calculate stats
        fitll = -firth_likelihood(beta_iterations[-1], logit_model)
        intercept = beta_iterations[-1][0]
        if len(beta_iterations[-1]) > 1:
            kbeta = beta_iterations[-1][1]
            bse = math.sqrt(-logit_model.hessian(beta_iterations[-1])[1, 1])
        else:
            # Encountered when fitting null without any distances/covariates
            kbeta = None
            bse = None

        if len(beta_iterations[-1]) > 2:
            beta = beta_iterations[-1][2:].tolist()
        else:
            beta = None

        return_fit = intercept, kbeta, beta, bse, fitll

    return return_fit
