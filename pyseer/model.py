# Copyright 2017 Marco Galardini and John Lees

'''Models implementations'''

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
from collections import namedtuple
import statsmodels.formula.api as smf


Seer = namedtuple('Seer', ['kmer',
                           'af', 'prep', 'lrt_pvalue',
                           'kbeta', 'bse',
                           'intercept', 'betas',
                           'kstrains', 'nkstrains',
                           'notes',
                           'prefilter', 'filter'])


def binary(kmer, p, k, m, c, af,
           pret, lrtt, null_res, null_firth,
           kstrains, nkstrains):
    notes = set()

    # was this af-filtered?
    if p is None:
        notes.add('af-filter')
        return Seer(kmer, af, np.nan, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, True, False)

    # pre-filtering
    t = np.concatenate((p.reshape(-1, 1), k.reshape(-1, 1)), axis=1).T
    table = [[t[0][(t[0] == 1) & (t[1] == 1)].shape[0],
              t[0][(t[0] == 1) & (t[1] == 0)].shape[0]],
             [t[0][(t[0] == 0) & (t[1] == 1)].shape[0],
              t[0][(t[0] == 0) & (t[1] == 0)].shape[0]]]

    # check for small values
    bad_chisq = 0
    bad_entries = 0
    for row in table:
        for entry in row:
            if entry <= 1:
                bad_chisq = True
                notes.add('bad-chisq')
            elif entry <= 5:
                bad_entries += 1
    if bad_entries > 1:
        bad_chisq = True
        notes.add('bad-chisq')

    prep = stats.chi2_contingency(table, correction=False)[1]
    if prep > pret or not np.isfinite(prep):
        notes.add('pre-filtering-failed')
        return Seer(kmer, af, prep, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, True, False)

    # actual logistic regression
    if c.shape[0] == m.shape[0]:
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
    mod = smf.Logit(p, v)

    start_vec = np.zeros(v.shape[1])
    start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))

    # suppress annoying stdout messages
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    try:
        if not bad_chisq:
            try:
                res = mod.fit(start_params=start_vec, method='newton')

                if res.bse[1] > 3:
                    bad_chisq = True
                    notes.add('high-bse')
                else:
                    lrstat = -2*(null_res - res.llf)
                    lrt_pvalue = 1
                    if lrstat > 0: # non-convergence
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
            firth_fit = fit_firth(mod, start_vec, kmer, v, p)
            if firth_fit is None: # Firth failure
                notes.add('firth-fail')
                return Seer(kmer, af, prep, np.nan,
                            np.nan, np.nan, np.nan, [],
                            kstrains, nkstrains,
                            notes, False, True)
            else:
                intercept, kbeta, beta, bse, fitll = firth_fit
                lrstat = -2*(null_firth - fitll)
                lrt_pvalue = 1
                if lrstat > 0: # non-convergence
                    lrt_pvalue = stats.chi2.sf(lrstat, 1)
    except np.linalg.linalg.LinAlgError:
        # singular matrix error
        notes.add('matrix-inversion-error')
        return Seer(kmer, af, prep, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, False, True)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    if lrt_pvalue > lrtt or not np.isfinite(lrt_pvalue) or not np.isfinite(kbeta):
        notes.add('lrt-filtering-failed')
        return Seer(kmer, af, prep, lrt_pvalue,
                    kbeta, bse, intercept, beta,
                    kstrains, nkstrains,
                    notes, False, True)

    return Seer(kmer, af, prep, lrt_pvalue,
                kbeta, bse, intercept, beta,
                kstrains, nkstrains,
                notes, False, False)


def firth_likelihood(beta, logit):
    return -(logit.loglike(beta) + 0.5*np.log(np.linalg.det(-logit.hessian(beta))))

# Do firth regression
# Note information = -hessian, for some reason available but not implemented in statsmodels
def fit_firth(logit_model, start_vec, kmer_name,
              X, y, step_limit=1000, convergence_limit=0.0001):

    beta_iterations = []
    beta_iterations.append(start_vec)
    for i in range(0, step_limit):
        pi = logit_model.predict(beta_iterations[i])
        W = np.diagflat(np.multiply(pi, 1-pi))
        var_covar_mat = np.linalg.pinv(-logit_model.hessian(beta_iterations[i]))

        # build hat matrix
        rootW = np.sqrt(W)
        H = np.dot(np.transpose(X), np.transpose(rootW))
        H = np.matmul(var_covar_mat, H)
        H = np.matmul(np.dot(rootW, X), H)

        # penalised score
        U = np.matmul(np.transpose(X), y - pi + np.multiply(np.diagonal(H), 0.5 - pi))
        new_beta = beta_iterations[i] + np.matmul(var_covar_mat, U)

        # step halving
        j = 0
        while firth_likelihood(new_beta, logit_model) > firth_likelihood(beta_iterations[i], logit_model):
            new_beta = beta_iterations[i] + 0.5*(new_beta - beta_iterations[i])
            j = j + 1
            if (j > step_limit):
                return None

        beta_iterations.append(new_beta)
        if i > 0 and (np.linalg.norm(beta_iterations[i] - beta_iterations[i-1]) < convergence_limit):
            break

    return_fit = None
    if np.linalg.norm(beta_iterations[i] - beta_iterations[i-1]) >= convergence_limit:
        pass
    else:
        # Calculate stats
        fitll = -firth_likelihood(beta_iterations[-1], logit_model)
        intercept = beta_iterations[-1][0]
        kbeta = beta_iterations[-1][1]
        beta = beta_iterations[-1][2:].tolist()
        bse = math.sqrt(-logit_model.hessian(beta_iterations[-1])[1,1])

        return_fit = intercept, kbeta, beta, bse, fitll

    return return_fit


def continuous(kmer, p, k, m, c, af,
               pret, lrtt, null_res, null_firth,
               kstrains, nkstrains):
    notes = set()

    # was this af-filtered?
    if p is None:
        notes.add('af-filter')
        return Seer(kmer, af, np.nan, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, True, False)

    # pre-filtering
    prep = stats.ttest_ind(p[k == 1],
                           p[k == 0],
                           equal_var=False)[1]
    if prep > pret or not np.isfinite(prep):
        notes.add('pre-filtering-failed')
        return Seer(kmer, af, prep, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, True, False)

    # actual linear regression
    if c.shape[0] == m.shape[0]:
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
    mod = smf.OLS(p, v)

    # suppress annoying stdout messages
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    try:
        res = mod.fit()
        intercept = res.params[0]
        kbeta = res.params[1]
        beta = res.params[2:]
        bse = res.bse[1]
    except np.linalg.linalg.LinAlgError:
        # singular matrix error
        # singular matrix error
        notes.add('matrix-inversion-error')
        return Seer(kmer, af, prep, np.nan,
                    np.nan, np.nan, np.nan, [],
                    kstrains, nkstrains,
                    notes, False, True)
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    lrt_pvalue = res.compare_lr_test(null_res)[1]
    if lrt_pvalue > lrtt or not np.isfinite(lrt_pvalue) or not np.isfinite(kbeta):
        notes.add('lrt-filtering-failed')
        return Seer(kmer, af, prep, lrt_pvalue,
                    kbeta, bse, intercept, beta,
                    kstrains, nkstrains,
                    notes, False, True)

    return Seer(kmer, af, prep, lrt_pvalue,
                kbeta, bse, intercept, beta,
                kstrains, nkstrains,
                notes, False, False)


# Fit the null model, regression without k-mer
def fit_null(p, m, cov, continuous, firth=False):
    if cov.shape[1] > 0:
        v = np.concatenate((np.ones(m.shape[0]).reshape(-1, 1),
                            m,
                            cov.values),
                           axis=1)
    else:
        # no covariates
        v = np.concatenate((np.ones(m.shape[0]).reshape(-1, 1),
                            m),
                           axis=1)

    if continuous:
        null_mod = mod = smf.OLS(p.values, v)
    else:
        start_vec = np.zeros(v.shape[1])
        start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))
        null_mod = smf.Logit(p, v)

    # suppress annoying stdout messages
    old_stdout = sys.stdout
    sys.stdout = open(os.devnull, "w")
    try:
        if continuous:
            null_res = null_mod.fit()
        else:
            if firth:
                (intercept, kbeta, beta, bse, fitll) = fit_firth(null_mod, start_vec, "null", v, p)
                null_res = fitll
            else:
                null_res = null_mod.fit(start_params=start_vec, method='newton')
                null_res = null_res.llf
    except np.linalg.linalg.LinAlgError:
        # singular matrix error
        sys.stderr.write('Matrix inversion error for null model\n')
        return None
    except statsmodels.tools.sm_exceptions.PerfectSeparationError:
        # singular matrix error
        sys.stderr.write('Perfetly separable data error for null model\n')
        return None
    finally:
        sys.stdout.close()
        sys.stdout = old_stdout

    return null_res

