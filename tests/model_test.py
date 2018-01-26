import warnings
import unittest
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from pyseer.model import pre_filtering
from pyseer.model import fit_null
from pyseer.model import fit_firth


np.random.seed(42)


class TestPreFiltering(unittest.TestCase):
    def test_pre_filtering_binary(self):
        p = np.random.randint(2, size=100)
        k = np.random.randint(2, size=100)
        prep, bad_chisq = pre_filtering(p, k, False)
        self.assertEqual(prep, 0.8838393707778841)
        self.assertFalse(bad_chisq)
        # continous phenotype
        p = np.random.random(100)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            prep, bad_chisq = pre_filtering(p, k, False)
        self.assertTrue(np.isnan(prep))
        self.assertTrue(bad_chisq)
        # bad_chisq example
        p = np.concatenate((np.ones(50), np.zeros(50)))
        k = np.concatenate((np.ones(45), np.zeros(55)))
        prep, bad_chisq = pre_filtering(p, k, False)
        self.assertEqual(prep, 1.4919966396986922e-19)
        self.assertTrue(bad_chisq)

    def test_pre_filtering_continuous(self):
        p = np.random.random(100)
        k = np.random.randint(2, size=100)
        prep, bad_chisq = pre_filtering(p, k, True)
        self.assertEqual(prep, 0.54292660587236519)
        self.assertFalse(bad_chisq)
        # using a binary p Matrix
        p = np.concatenate((np.ones(50), np.zeros(50)))
        k = np.concatenate((np.ones(45), np.zeros(55)))
        prep, bad_chisq = pre_filtering(p, k, True)
        self.assertEqual(prep, 8.6308642007939013e-30)
        self.assertFalse(bad_chisq)


class TestFitNull(unittest.TestCase):
    def test_fit_null_binary(self):
        p = np.random.randint(2, size=100)
        m = np.random.random(size=(100, 10))
        cov = pd.DataFrame([])
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m),
                           axis=1)
        start_vec = np.zeros(v.shape[1])
        start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))
        null_mod = smf.Logit(p, v)
        null_test = null_mod.fit(start_params=start_vec,
                                 method='newton',
                                 disp=False)
        # no covariates
        null_res = fit_null(p, m, cov, False, firth=False)
        self.assertTrue(abs((null_test.params - null_res.params).max()) < 1E-15)
        self.assertAlmostEqual(null_test.llr_pvalue, null_res.llr_pvalue)
        # no covariates, firth regression
        (intercept, kbeta, beta, bse, fitll) = fit_firth(null_mod,
                                                         start_vec,
                                                         'null',
                                                         v, p)
        null_test = fitll
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_test, null_res)
        # covariates
        cov = pd.DataFrame(np.random.random(size=(100, 3)))
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m,
                            cov),
                           axis=1)
        start_vec = np.zeros(v.shape[1])
        start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))
        null_mod = smf.Logit(p, v)
        null_test = null_mod.fit(start_params=start_vec,
                                 method='newton',
                                 disp=False)
        null_res = fit_null(p, m, cov, False, firth=False)
        self.assertTrue(abs((null_test.params - null_res.params).max()) < 1E-15)
        self.assertAlmostEqual(null_test.llr_pvalue, null_res.llr_pvalue)
        # covariates, firth regression
        (intercept, kbeta, beta, bse, fitll) = fit_firth(null_mod,
                                                         start_vec,
                                                         'null',
                                                         v, p)
        null_test = fitll
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_test, null_res)


if __name__ == '__main__':
    unittest.main()
