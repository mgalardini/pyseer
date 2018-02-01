import os
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
from pyseer.model import firth_likelihood
from pyseer.model import fixed_effects_regression


DATA_DIR = 'tests/unit_tests_data'
P_BINARY = os.path.join(DATA_DIR, 'p_binary.txt')
P_CONT = os.path.join(DATA_DIR, 'p_continuous.txt')
K = os.path.join(DATA_DIR, 'k.txt')
M = os.path.join(DATA_DIR, 'm.txt')
COV = os.path.join(DATA_DIR, 'cov.txt')
LIN = os.path.join(DATA_DIR, 'lin.txt')
FIRTH_VARS = os.path.join(DATA_DIR, 'firth_vars.txt')


class TestPreFiltering(unittest.TestCase):
    def test_pre_filtering_binary(self):
        p = np.loadtxt(P_BINARY)
        k = np.loadtxt(K)
        prep, bad_chisq = pre_filtering(p, k, False)
        self.assertEqual(prep, 0.5365065578449575)
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
        p_cont = np.loadtxt(P_CONT)
        k = np.loadtxt(K)
        prep, bad_chisq = pre_filtering(p_cont, k, True)
        self.assertEqual(prep, 0.29623810011571716)
        self.assertFalse(bad_chisq)
        # using a binary p Matrix
        p = np.concatenate((np.ones(50), np.zeros(50)))
        k = np.concatenate((np.ones(45), np.zeros(55)))
        prep, bad_chisq = pre_filtering(p, k, True)
        self.assertEqual(prep, 8.6308642007939013e-30)
        self.assertFalse(bad_chisq)


class TestFitNull(unittest.TestCase):
    def test_fit_null_binary(self):
        p = np.loadtxt(P_BINARY)
        m = np.loadtxt(M)
        cov = pd.DataFrame([])
        # no covariates
        params = np.array([-1.41572498, 0.35847998, -0.03014792, 2.46252819, 0.96908425,
                           -0.20952455, -0.27988125, 0.36798503, -0.03278285, -1.34132024,
                           0.844149])
        null_res = fit_null(p, m, cov, False, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)
        # no covariates, firth regression
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_res, -57.884527394557985)
        # covariates
        cov = np.loadtxt(COV)
        cov = pd.DataFrame(cov)
        params = np.array([-0.87072948, 0.26456701, 0.03485904, 2.80243184, 1.086393,
                           -0.3882244, -0.46883396, 0.61387846, 0.09962477, -1.45376984,
                           0.93929299, 0.07927743, -1.54631396, 0.1098796])
        null_res = fit_null(p, m, cov, False, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)
        # covariates, firth regression
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_res, -55.60790630835098)
        # perfectly separable data
        p = np.array([1]*10 + [0]*90)
        m = np.array([1]*10 + [0]*90).reshape(-1, 1)
        cov = pd.DataFrame([])
        self.assertEqual(fit_null(p, m, cov, False, False), None)

    def test_fit_null_continuous(self):
        p_cont = np.loadtxt(P_CONT)
        m = np.loadtxt(M)
        # no covariates
        params = np.array([0.65572473, -0.16129649, 0.03417796, -0.08011702, 0.10902641,
                           0.00599514, -0.09081684, -0.13653787, 0.17798003, -0.16793408,
                           0.12959982])
        null_res = fit_null(p_cont, m, pd.DataFrame([]), True, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)
        # covariates
        cov = np.loadtxt(COV)
        cov = pd.DataFrame(cov)
        params = np.array([3.13564240e-01, 5.50424455e-02, 1.71810270e-04, 6.03266921e-01,
                           2.29686618e-01, -7.61535306e-02, -1.08730249e-01, 1.27931038e-01,
                           2.74059958e-02, -3.04257997e-01, 1.85369147e-01, 1.94048312e-02,
                           -3.16794302e-01, 1.64952764e-02])
        null_res = fit_null(p, m, cov, True, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)


class TestFitLineageEffect(unittest.TestCase):
    def test_fit_lineage_effect(self):
        k = np.loadtxt(K)
        m = np.loadtxt(M)
        lin = np.loadtxt(LIN)
        cov = np.loadtxt(COV)
        # no covariates
        max_lineage = fit_lineage_effect(m, pd.DataFrame([]), k)
        self.assertEqual(max_lineage, 2)
        # no covariates, binary lineage Matrix
        # as if the data was loaded from a file
        max_lineage = fit_lineage_effect(lin, pd.DataFrame([]), k)
        self.assertEqual(max_lineage, 2)
        # covariates
        max_lineage = fit_lineage_effect(m, cov, k)
        self.assertEqual(max_lineage, 2)
        # perfectly separable data
        k = np.array([1]*10 + [0]*90)
        m = np.array([1]*10 + [0]*90).reshape(-1, 1)
        cov = pd.DataFrame([])
        self.assertEqual(fit_lineage_effect(m, cov, k), None)


class TestFirthFit(unittest.TestCase):
    def test_firth_likelihood(self):
        p = np.loadtxt(P_BINARY)
        m = np.loadtxt(M)
        firth_vars = np.loadtxt(FIRTH_VARS)
        mod = smf.Logit(p, m)
        fll = firth_likelihood(firth_vars, mod)
        self.assertAlmostEqual(fll, 97.13375906431875)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fll = firth_likelihood(firth_vars + 100,
                                   mod)
        self.assertAlmostEqual(fll, np.inf)

    def test_fit_firth(self):
        p = np.loadtxt(P_BINARY)
        m = np.loadtxt(M)
        mod = smf.Logit(p, m)
        start_vec = np.zeros(m.shape[1])
        start_vec[0] = np.log(np.mean(p)/(1-np.mean(p)))
        (intercept, kbeta, beta, bse, fitll) = fit_firth(mod,
                                                         start_vec,
                                                         m, p)
        self.assertAlmostEqual(intercept, 0.13954805021495864)
        self.assertAlmostEqual(kbeta, -0.31901219992017243)
        tbeta = [1.9588025,  0.7251749, -0.5605268, -0.5396909,  0.0594742,
                 -0.2001795, -1.4873298,  0.5050208]
        self.assertTrue(abs((np.array(beta)- np.array(tbeta)).max())
                        < 1E-7)
        self.assertAlmostEqual(bse, 2.848207537910185)
        self.assertAlmostEqual(fitll, -58.249948818380204)
        fitll = fit_firth(mod,
                          start_vec,
                          m, p,
                          step_limit=10,
                          convergence_limit=1E-10)
        self.assertEqual(fitll, None)


class TestFixedEffectsRegression(unittest.TestCase):
    def test_fixed_effects_regression_binary(self):
        p = np.loadtxt(P_BINARY)
        k = np.loadtxt(K)
        m = np.loadtxt(M)
        lin = np.loadtxt(LIN)
        cov = pd.DataFrame([])
        variant = 'variant'
        af = 0.2
        pattern = 'test'
        lineage_effects = False
        lin = np.random.randint(2, size=(100, 4))
        null_fit = -9.9
        null_firth = -9.9
        kstrains = ['K%d' % i for i in range(k[k == 1].shape[0])]
        nkstrains = ['NK%d' % i for i in range(k[k == 0].shape[0])]
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, lineage_effects, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        # TODO: test var_obj, probe different exit statuses

    def test_fixed_effects_regression_continuos(self):
        p_cont = np.loadtxt(P_CONT)
        k = np.loadtxt(K)
        m = np.loadtxt(M)
        lin = np.loadtxt(LIN)
        cov = pd.DataFrame([])
        variant = 'variant'
        af = 0.2
        pattern = 'test'
        lineage_effects = False
        lin = np.random.randint(2, size=(100, 4))
        null_mod = smf.Logit(p, m)
        null_fit = null_mod.fit(disp=False)
        null_firth = -9.9
        kstrains = ['K%d' % i for i in range(k[k == 1].shape[0])]
        nkstrains = ['NK%d' % i for i in range(k[k == 0].shape[0])]
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, lineage_effects, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        # TODO: test var_obj, probe different exit statuses


if __name__ == '__main__':
    unittest.main()
