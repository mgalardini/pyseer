import os
import warnings
import unittest
import numpy as np
import pandas as pd
from pyseer.lmm import initialise_lmm


DATA_DIR = 'tests'
P_BINARY = os.path.join(DATA_DIR, 'subset.pheno')
S = os.path.join(DATA_DIR, 'similarity_subset.tsv.gz')
COV = os.path.join(DATA_DIR, 'covariates.txt')
C = os.path.join(DATA_DIR, 'lmm_cache.npz')


class TestInitialiseLmm(unittest.TestCase):
    def test_initialise_lmm(self):
        p = pd.read_table(P_BINARY,
                          index_col=0)['binary']
        cov = pd.DataFrame([])
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        self.assertEqual(x.shape[0], 50)
        self.assertAlmostEqual(y.findH2()['nLL'][0],
                               35.7033778)
        self.assertAlmostEqual(z, 0.0)
        # covariates
        cov = pd.read_table(COV, index_col=0,
                            header=None)
        x, y, z = initialise_lmm(p, cov, S,
                                 lmm_cache_in=None,
                                 lmm_cache_out=None)
        self.assertEqual(x.shape[0], 50)
        self.assertAlmostEqual(y.findH2()['nLL'][0],
                               34.55403861)
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
                               34.55403861)
        self.assertAlmostEqual(z, 0.0)
        # different sizes
        b = pd.Series(np.random.random(10),
                      index=['test_%d' % x for x in range(10)])
        with self.assertRaises(SystemExit) as cm:
            initialise_lmm(b, cov, S,
                           lmm_cache_in=C,
                           lmm_cache_out=None)
        self.assertEqual(cm.exception.code, 1)


if __name__ == '__main__':
    unittest.main()
