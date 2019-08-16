import os
import warnings
import unittest
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf
# handle different versions of statsmodels
try:
    smf.Logit
except AttributeError:
    smf.Logit = statsmodels.discrete.discrete_model.Logit
from pyseer.model import pre_filtering
from pyseer.model import fit_null
from pyseer.model import fit_firth
from pyseer.model import fit_lineage_effect
from pyseer.model import firth_likelihood
from pyseer.model import fixed_effects_regression
from pyseer.classes import Seer


DATA_DIR = 'tests/unit_tests_data'
P_BINARY = os.path.join(DATA_DIR, 'p_binary.txt')
P_CONT = os.path.join(DATA_DIR, 'p_continuous.txt')
K = os.path.join(DATA_DIR, 'k.txt')
M = os.path.join(DATA_DIR, 'm.txt')
COV = os.path.join(DATA_DIR, 'cov.txt')
LIN = os.path.join(DATA_DIR, 'lin.txt')
FIRTH_VARS = os.path.join(DATA_DIR, 'firth_vars.txt')


def eq_seer(s1, s2):
    """Test whether two Seer objects are the same"""
    diff = set()
    for p in ['kmer', 'pattern',
              'kstrains', 'nkstrains', 'notes',
              'prefilter', 'filter']:
        x = getattr(s1, p)
        y = getattr(s2, p)
        if x != y:
            diff.add(p)

    for p in ['af', 'prep', 'pvalue',
              'kbeta', 'bse', 'intercept']:
        x = getattr(s1, p)
        y = getattr(s2, p)
        if not np.isfinite(x) and not np.isfinite(y):
            continue
        if np.isfinite(x) and not np.isfinite(y):
            diff.add(p)
        if np.isfinite(y) and not np.isfinite(x):
            diff.add(p)
        if abs(x - y) > 1E-7:
            diff.add(p)

    if s1.max_lineage is not None and s2.max_lineage is not None:
        p = 'max_lineage'
        x = getattr(s1, p)
        y = getattr(s2, p)
        if not np.isfinite(x) and not np.isfinite(y):
            pass
        else:
            if np.isfinite(x) and not np.isfinite(y):
                diff.add(p)
            if np.isfinite(y) and not np.isfinite(x):
                diff.add(p)
            if x != y:
                diff.add(p)
    elif s1.max_lineage is None and s2.max_lineage is None:
        pass
    else:
        diff.add('max_lineage')

    if s1.betas.shape[0] > 0 and s2.betas.shape[0] > 0:
        if s1.betas.shape[0] != s2.betas.shape[0]:
            diff.add('betas')
        if abs((s1.betas - s2.betas).max()) > 1E-7:
            diff.add('betas')

    return diff


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
        params = np.array([-1.41572498, 0.35847998, -0.03014792, 2.46252819,
                           0.96908425, -0.20952455, -0.27988125, 0.36798503,
                           -0.03278285, -1.34132024, 0.844149])
        null_res = fit_null(p, m, cov, False, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)
        # no covariates, firth regression
        null_res = fit_null(p, m, cov, False, firth=True)
        self.assertAlmostEqual(null_res, -57.884527394557985)
        # covariates
        cov = np.loadtxt(COV)
        cov = pd.DataFrame(cov)
        params = np.array([-0.87072948, 0.26456701, 0.03485904, 2.80243184,
                           1.086393, -0.3882244, -0.46883396, 0.61387846,
                           0.09962477, -1.45376984, 0.93929299, 0.07927743,
                           -1.54631396, 0.1098796])
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
        params = np.array([0.65572473, -0.16129649, 0.03417796, -0.08011702,
                           0.10902641, 0.00599514, -0.09081684, -0.13653787,
                           0.17798003, -0.16793408, 0.12959982])
        null_res = fit_null(p_cont, m, pd.DataFrame([]), True, firth=False)
        self.assertTrue(abs((params - null_res.params).max())
                        < 1E-7)
        # covariates
        cov = np.loadtxt(COV)
        cov = pd.DataFrame(cov)
        params = np.array([0.49070237, -0.17284083, 0.00710691, -0.11784811,
                           0.07352861, 0.01219004, -0.04772721, -0.17089199,
                           0.18198025, -0.17141095, 0.11330439, 0.08887165,
                           0.20304982, 0.13802362])
        null_res = fit_null(p_cont, m, cov, True, firth=False)
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
        self.assertTrue(abs((np.array(beta)-np.array(tbeta)).max())
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
        null_fit = -9.9
        null_firth = -9.9
        kstrains = ['K%d' % i for i in range(k[k == 1].shape[0])]
        nkstrains = ['NK%d' % i for i in range(k[k == 0].shape[0])]
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, 0.5365065578449575, 1,
                     -0.668215625696782,
                     0.47087488598995186,
                     -1.29962042280822,
                     np.array([0.42265596, 0.10078512, 2.77587593,
                               0.94439244, -0.13846857, -0.14140035,
                               0.38328562, -0.1986484, -1.51779346,
                               0.94618541]),
                     None, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail af-filtering
        var_obj = fixed_effects_regression(variant, None, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, np.nan,
                     np.nan, np.nan, np.nan, np.nan,
                     np.array([]),
                     None, kstrains, nkstrains,
                     set(['af-filter']),
                     True, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail pre-filtering
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, False, None,
                                           0.05, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, 0.5365065578449575,
                     np.nan, np.nan, np.nan, np.nan,
                     np.array([]),
                     None, kstrains, nkstrains,
                     set(['pre-filtering-failed']),
                     True, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail filtering
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, False, None,
                                           1, 0.05, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, 0.5365065578449575, 1,
                     -0.668215625696782,
                     0.47087488598995186,
                     -1.29962042280822,
                     np.array([0.42265596, 0.10078512, 2.77587593,
                               0.94439244, -0.13846857, -0.14140035,
                               0.38328562, -0.1986484, -1.51779346,
                               0.94618541]),
                     None, kstrains, nkstrains,
                     set(['lrt-filtering-failed']),
                     False, True)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # bad-chisq
        p = np.array([1]*10 + [0]*90)
        k = np.array([1]*10 + [0]*90)
        m = np.array([1]*10 + [0]*90).reshape(-1, 1)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                               pattern, False, None,
                                               1, 1, null_fit, null_firth,
                                               kstrains, nkstrains,
                                               False)
        t_obj = Seer(variant, pattern,
                     af, 1.5239706048320995e-23, 1,
                     -88.72472137305186,
                     0.0,
                     890.8154121360252,
                     np.array([-88.72472137305188]),
                     None, kstrains, nkstrains,
                     set(['bad-chisq']),
                     False, False)
        #self.assertEqual(eq_seer(var_obj, t_obj), set())
        self.assertEqual(var_obj.notes, t_obj.notes)
        p = np.loadtxt(P_BINARY)
        k = np.loadtxt(K)
        m = np.loadtxt(M)
        # covariates
        cov = pd.DataFrame(np.loadtxt(COV))
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, 0.5365065578449575, 1,
                     -0.7082070719359966,
                     0.4852518061533321,
                     -0.809194818156449,
                     np.array([0.325464, 0.16147301, 3.17003634,
                               1.05383182, -0.31762591, -0.32545411,
                               0.65876263, -0.07939636, -1.61743885,
                               1.04396837, 0.13034889, -1.59225167,
                               0.1938934]),
                     None, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # lineage
        cov = pd.DataFrame([])
        var_obj = fixed_effects_regression(variant, p, k, m, cov, af,
                                           pattern, True, lin,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           False)
        t_obj = Seer(variant, pattern,
                     af, 0.5365065578449575, 1,
                     -0.668215625696782,
                     0.47087488598995186,
                     -1.29962042280822,
                     np.array([0.42265596, 0.10078512, 2.77587593,
                               0.94439244, -0.13846857, -0.14140035,
                               0.38328562, -0.1986484, -1.51779346,
                               0.94618541]),
                     2, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # TODO: probe different exit statuses:
        # 1- prefectly-separated-data (is it actually possible?)
        # 2- firth-fail
        # 3- matrix-inversion-error
        # 4- high-bse

    def test_fixed_effects_regression_continuous(self):
        p_cont = np.loadtxt(P_CONT)
        k = np.loadtxt(K)
        m = np.loadtxt(M)
        lin = np.loadtxt(LIN)
        cov = pd.DataFrame([])
        variant = 'variant'
        af = 0.2
        pattern = 'test'
        lineage_effects = False
        null_mod = smf.Logit(p_cont, m)
        null_fit = null_mod.fit(disp=False)
        null_firth = -9.9
        kstrains = ['K%d' % i for i in range(k[k == 1].shape[0])]
        nkstrains = ['NK%d' % i for i in range(k[k == 0].shape[0])]
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, lineage_effects, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af,
                     0.29623810011571716,
                     0.4694146479961355,
                     -0.043638262259610316,
                     0.06006023185402142,
                     0.6655803214920781,
                     np.array([-0.1560651, 0.04372272, -0.06398297,
                               0.10658197,  0.01046428, -0.08089156,
                               -0.13733075, 0.16774866, -0.17746121,
                               0.13386466]),
                     None, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail af-filtering
        var_obj = fixed_effects_regression(variant, None, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af, np.nan,
                     np.nan, np.nan, np.nan, np.nan,
                     np.array([]),
                     None, kstrains, nkstrains,
                     set(['af-filter']),
                     True, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail pre-filtering
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, False, None,
                                           0.05, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af, 0.29623810011571716,
                     np.nan, np.nan, np.nan, np.nan,
                     np.array([]),
                     None, kstrains, nkstrains,
                     set(['pre-filtering-failed']),
                     True, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # fail filtering
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1E-50, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af,
                     0.29623810011571716,
                     0.4694146479961355,
                     -0.043638262259610316,
                     0.06006023185402142,
                     0.6655803214920781,
                     np.array([-0.1560651, 0.04372272, -0.06398297,
                               0.10658197,  0.01046428, -0.08089156,
                               -0.13733075, 0.16774866, -0.17746121,
                               0.13386466]),
                     None, kstrains, nkstrains,
                     set(['lrt-filtering-failed']),
                     False, True)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # covariates
        cov = pd.DataFrame(np.loadtxt(COV))
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, False, None,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af,
                     0.29623810011571716,
                     0.4039092383440829,
                     -0.04946894010582922,
                     0.05897268709495734,
                     0.49957867277580303,
                     np.array([-0.16730353, 0.01750906, -0.09994545,
                               0.07018266, 0.01718979,
                               -0.03593312, -0.17211066,
                               0.17065225, -0.18230721, 0.11787759,
                               0.09058623, 0.20484901, 0.14072312]),
                     None, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # lineage
        cov = pd.DataFrame([])
        var_obj = fixed_effects_regression(variant, p_cont, k, m, cov, af,
                                           pattern, True, lin,
                                           1, 1, null_fit, null_firth,
                                           kstrains, nkstrains,
                                           True)
        t_obj = Seer(variant, pattern,
                     af,
                     0.29623810011571716,
                     0.4694146479961355,
                     -0.043638262259610316,
                     0.06006023185402142,
                     0.6655803214920781,
                     np.array([-0.1560651, 0.04372272, -0.06398297,
                               0.10658197,  0.01046428, -0.08089156,
                               -0.13733075, 0.16774866, -0.17746121,
                               0.13386466]),
                     2, kstrains, nkstrains,
                     set(),
                     False, False)
        self.assertEqual(eq_seer(var_obj, t_obj), set())
        # TODO: probe matrix-inversion


if __name__ == '__main__':
    unittest.main()
