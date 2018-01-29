import warnings
import unittest
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
from pyseer.model import pre_filtering
from pyseer.model import fit_null
from pyseer.model import fit_firth
from pyseer.model import fit_lineage_effect


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
        self.assertTrue(abs((null_test.params - null_res.params).max())
                        < 1E-15)
        # scipy >= 1.x.x removed a function used here
        # should not affect anything else in pyseer
        # so skipping this test for now
        # TODO: put it back when problem is fixed
        # self.assertAlmostEqual(null_test.llr_pvalue, null_res.llr_pvalue)
        #
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
        self.assertTrue(abs((null_test.params - null_res.params).max())
                        < 1E-15)
        # scipy >= 1.x.x removed a function used here
        # should not affect anything else in pyseer
        # so skipping this test for now
        # TODO: put it back when problem is fixed
        # self.assertAlmostEqual(null_test.llr_pvalue, null_res.llr_pvalue)
        #
        # covariates, firth regression
        (intercept, kbeta, beta, bse, fitll) = fit_firth(null_mod,
                                                         start_vec,
                                                         'null',
                                                         v, p)
        null_test = fitll
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_test, null_res)
        # perfectly separable data
        p = np.array([1]*10 + [0]*90)
        m = np.array([1]*10 + [0]*90).reshape(-1, 1)
        cov = pd.DataFrame([])
        self.assertEqual(fit_null(p, m, cov, False, False), None)

    def test_fit_null_continuous(self):
        p = np.random.random(size=100)
        m = np.random.random(size=(100, 10))
        cov = pd.DataFrame([])
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m),
                           axis=1)
        null_mod = smf.OLS(p, v)
        null_test = null_mod.fit(disp=False)
        # no covariates
        null_res = fit_null(p, m, cov, True, firth=False)
        self.assertTrue(abs((null_test.params - null_res.params).max())
                        < 1E-15)
        # covariates
        cov = pd.DataFrame(np.random.random(size=(100, 3)))
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m,
                            cov),
                           axis=1)
        null_mod = smf.OLS(p, v)
        null_test = null_mod.fit(disp=False)
        null_res = fit_null(p, m, cov, True, firth=False)
        self.assertTrue(abs((null_test.params - null_res.params).max())
                        < 1E-15)


class TestFitLineageEffect(unittest.TestCase):
    def test_fit_lineage_effect(self):
        k = np.random.randint(2, size=100)
        m = np.random.random(size=(100, 10))
        cov = pd.DataFrame([])
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m),
                           axis=1)
        lineage_mod = smf.Logit(k, v)
        lineage_res = lineage_mod.fit(method='newton',
                                      disp=False)
        wald_test = np.divide(np.absolute(lineage_res.params), lineage_res.bse)
        lineage_test = np.argmax(wald_test[1:m.shape[1]+1])
        # no covariates
        max_lineage = fit_lineage_effect(m, cov, k)
        self.assertEqual(lineage_test, max_lineage)
        # no covariates, binary lineage Matrix
        # as if the data was loaded from a file
        m = np.random.randint(2, size=(100, 10))
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m),
                           axis=1)
        lineage_mod = smf.Logit(k, v)
        lineage_res = lineage_mod.fit(method='newton',
                                      disp=False)
        wald_test = np.divide(np.absolute(lineage_res.params), lineage_res.bse)
        lineage_test = np.argmax(wald_test[1:m.shape[1]+1])
        max_lineage = fit_lineage_effect(m, cov, k)
        self.assertEqual(lineage_test, max_lineage)
        # covariates
        cov = np.random.random(size=(100, 3))
        v = np.concatenate((np.ones(100).reshape(-1, 1),
                            m,
                            cov),
                           axis=1)
        lineage_mod = smf.Logit(k, v)
        lineage_res = lineage_mod.fit(method='newton',
                                      disp=False)
        wald_test = np.divide(np.absolute(lineage_res.params), lineage_res.bse)
        lineage_test = np.argmax(wald_test[1:m.shape[1]+1])
        max_lineage = fit_lineage_effect(m, cov, k)
        self.assertEqual(lineage_test, max_lineage)
        # perfectly separable data
        k = np.array([1]*10 + [0]*90)
        m = np.array([1]*10 + [0]*90).reshape(-1, 1)
        cov = pd.DataFrame([])
        self.assertEqual(fit_lineage_effect(m, cov, k), None)


if __name__ == '__main__':
    unittest.main()
