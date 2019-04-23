# -*- coding: utf-8 -*-
# Copyright 2019 Marco Galardini and John Lees

'''Random forest model implementations'''

import os
import sys
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
from scipy.sparse import csr_matrix, csc_matrix, hstack
from decimal import Decimal

from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

def fit_rf(p, variants, covariates, weights, continuous, n_cpus = 1):
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
        weights (numpy.array)
            Sample weights (n, 1)
        continuous (bool)
            If True fit a regressor, otherwise classifier
        n_cpus (int)
            Number of processes to use in cross-validation
            Set to -1 to use all available

            [default = 1]

    Returns:
        betas (numpy.array)
            The fitted betas (slopes) for each variant
    """
    if continuous:
        clf = RandomForestRegressor(n_estimators=20, max_depth=None, min_samples_split=0.005,
                                    min_weight_fraction_leaf=0.002, max_leaf_nodes=None, n_jobs=n_cpus)
    else:
        clf = RandomForestClassifier(n_estimators=20, max_depth=None, min_samples_split=0.005,
                                     min_weight_fraction_leaf=0.002, max_leaf_nodes=None, n_jobs=n_cpus)

    if covariates.shape[0] > 0:
        variants = hstack([csc_matrix(covariates.values), variants])

    # Run model fit
    clf.fit(X = variants, y = p.values, sample_weight = weights)
    preds = clf.predict(variants)

    # Importances from sklearn
    # TODO should be using rfpimp? https://github.com/parrt/random-forest-importances
    betas = clf.feature_importances_

    # Write some summary stats
    # R^2 = 1 - sum((yi_obs - yi_predicted)^2) /sum((yi_obs - yi_mean)^2)
    SStot = np.sum(np.square(p.values - np.mean(p.values)))
    SSerr = np.sum(np.square(p.values.reshape(-1, 1) - preds))
    R2 = 1 - (SSerr/SStot)
    sys.stderr.write("Out of bag score: " + '%.3f' % Decimal(clf.oob_score_) + "\n")
    sys.stderr.write("Best R^2 from cross-validation: " + '%.3f' % Decimal(R2) + "\n")

    return(clf, betas.reshape(-1,))


