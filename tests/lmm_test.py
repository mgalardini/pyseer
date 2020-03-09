import os
import warnings
import unittest
import numpy as np
import pandas as pd
from pyseer.lmm import initialise_lmm
from pyseer.lmm import fit_lmm
from pyseer.lmm import fit_lmm_block
from pyseer.classes import LMM


DATA_DIR = 'tests'
P_BINARY = os.path.join(DATA_DIR, 'subset.pheno')
S = os.path.join(DATA_DIR, 'similarity_subset.tsv.gz')
COV = os.path.join(DATA_DIR, 'covariates.txt')
C = os.path.join(DATA_DIR, 'lmm_cache.npz')
K = os.path.join(DATA_DIR, 'unit_tests_data', 'k.txt')
M = os.path.join(DATA_DIR, 'unit_tests_data', 'm.txt')


def eq_lmm(s1, s2):
    """Test whether two LMM objects are the same"""
    diff = set()
    for p in ['kmer', 'pattern',
              'kstrains', 'nkstrains', 'notes',
              'prefilter', 'filter']:
        x = getattr(s1, p)
        y = getattr(s2, p)
        if x != y:
            diff.add(p)

    for p in ['af', 'prep', 'pvalue',
              'kbeta', 'bse', 'frac_h2']:
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

    return diff


class TestInitialiseLmm(unittest.TestCase):
    def test_initialise_lmm(self):
        p = pd.read_csv(P_BINARY,
                        index_col=0,
                        sep='\t')['binary']
        cov = pd.DataFrame([])
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        self.assertEqual(x.shape[0], 50)
        self.assertAlmostEqual(y.findH2()['nLL'][0],
                               35.7033778)
        self.assertAlmostEqual(z, 0.0)
        # covariates
        cov = pd.read_csv(COV, index_col=0,
                          sep='\t')
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        self.assertEqual(x.shape[0], 50)
        self.assertAlmostEqual(y.findH2()['nLL'][0],
                               34.554038607321814)
        self.assertAlmostEqual(z, 0.0)
        # sample names not matching
        b = pd.Series(np.random.random(100),
                      index=['test_%d' % x for x in range(100)])
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            x, y, z = initialise_lmm(b, cov, S,
                                     lmm_cache_in=None,
                                     lmm_cache_out=None)
            self.assertEqual(x.shape[0], 0)
            self.assertTrue(not np.isfinite(y.findH2()['nLL'][0]))
            self.assertAlmostEqual(z, 0.0)
        # save cache
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=C)
        # load cache
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=C,
                                 lmm_cache_out=None)
        self.assertEqual(x.shape[0], 50)
        self.assertAlmostEqual(y.findH2()['nLL'][0],
                               34.554038607321814)
        self.assertAlmostEqual(z, 0.0)
        # different sizes
        b = pd.Series(np.random.random(10),
                      index=['test_%d' % x for x in range(10)])
        with self.assertRaises(SystemExit) as cm:
            initialise_lmm(b, cov, S,
                           lmm_cache_in=C,
                           lmm_cache_out=None)
            self.assertEqual(cm.exception.code, 1)
        # matching lineage samples
        cov = pd.DataFrame([])
        s = pd.read_csv(S, index_col=0,
                        sep='\t')
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None,
                                 lineage_samples=s.index)
        # non-matching lineage samples
        with self.assertRaises(SystemExit) as cm:
            x, y, z = initialise_lmm(p, cov, S,
                                     lmm_cache_in=None,
                                     lmm_cache_out=None,
                                     lineage_samples=s.index[:-1])


class TestFitLmm(unittest.TestCase):
    def test_fit_lmm(self):
        p = pd.read_csv(P_BINARY,
                        index_col=0,
                        sep='\t')['binary']
        cov = pd.DataFrame([])
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          False, 1, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.28252075514059294,
                            0.2920532220978148,
                            0.1513687600644123,
                            0.1420853593711293,
                            0.1519818397711344,
                            None,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(), False, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # af filtering
        var = LMM('variant',
                  None,
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          False, 1, 1)
        test_results = [LMM('variant', None, 0.2,
                            np.nan, np.nan, np.nan,
                            np.nan, np.nan, np.nan,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(['af-filter']), True, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # bad-chisq
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        bad_k = np.array([1]*5 + [0]*(p.shape[0]-5))
        variants = [(var, p.values, bad_k),]
        variant_mat = bad_k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          False, 1, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.2544505826463333,
                            0.263519965703956,
                            0.2666666666666663,
                            0.2357022603955158,
                            0.16116459280507586,
                            None,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(['bad-chisq']), False, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # pre-filtering
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          False, 0.05, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.28252075514059294,
                            np.nan, np.nan,
                            np.nan, np.nan, np.nan,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(['pre-filtering-failed']), True, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # lrt-filtering
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          False, 1, 0.05)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.28252075514059294,
                            0.2920532220978148,
                            np.nan, np.nan, np.nan, np.nan,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(['lrt-filtering-failed']), False, True),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # lineage fit
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        m = np.loadtxt(M)[:p.shape[0]]
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          True, m, cov,
                          False, 1, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.28252075514059294,
                            0.2920532220978148,
                            0.1513687600644123,
                            0.1420853593711293,
                            0.1519818397711344,
                            0,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(), False, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # lineage fit + covariates
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        m = np.loadtxt(M)[:p.shape[0]]
        cov = pd.read_csv(COV, index_col=0, header=None,
                          sep='\t').values
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          True, m, cov,
                          False, 1, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.28252075514059294,
                            0.2920532220978148,
                            0.1513687600644123,
                            0.1420853593711293,
                            0.1519818397711344,
                            0,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(), False, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())
        # continuous phenotype
        var = LMM('variant',
                  'pattern',
                  0.2,
                  np.nan, np.nan, np.nan,
                  np.nan, np.nan, np.nan,
                  ['k%d' % x
                   for x in range(p[p == 1].shape[0])],
                  ['nk%d' % x
                   for x in range(p[p == 0].shape[0])],
                  set(), True, True)
        k = np.loadtxt(K)[:p.shape[0]]
        variants = [(var, p.values, k),]
        variant_mat = k.reshape(-1, 1)
        results = fit_lmm(y, z,
                          variants, variant_mat,
                          False, [], cov,
                          True, 1, 1)
        test_results = [LMM('variant', 'pattern', 0.2,
                            0.2937152511367835,
                            0.2920532220978148,
                            0.1513687600644123,
                            0.1420853593711293,
                            0.1519818397711344,
                            None,
                            ['k%d' % x
                             for x in range(p[p == 1].shape[0])],
                            ['nk%d' % x
                             for x in range(p[p == 0].shape[0])],
                            set(), False, False),]
        for var, test_var in zip(results, test_results):
            self.assertEqual(eq_lmm(var, test_var), set())


class TestFitLmmBlock(unittest.TestCase):
    def test_fit_lmm_block(self):
        p = pd.read_csv(P_BINARY,
                        index_col=0,
                        sep='\t')['binary']
        cov = pd.DataFrame([])
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        k = np.loadtxt(K)[:p.shape[0]]
        variant_mat = k.reshape(-1, 1)
        result = fit_lmm_block(y, z, variant_mat)
        self.assertAlmostEqual(result['beta'][0],
                               0.15136876)
        self.assertAlmostEqual(result['bse'][0],
                               0.14208536)
        self.assertAlmostEqual(result['frac_h2'][0],
                               0.15198184)
        self.assertAlmostEqual(result['p_values'][0],
                               0.29205322)
        # impossibly high h2
        with self.assertRaises(KeyError):
            fit_lmm_block(y, 1, variant_mat)
        # shape mismatch
        with self.assertRaises(AssertionError):
            fit_lmm_block(y, z, variant_mat[:10])


if __name__ == '__main__':
    unittest.main()
