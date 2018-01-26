import warnings
import unittest
import numpy as np
from scipy import stats
from pyseer.model import pre_filtering


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

if __name__ == '__main__':
    unittest.main()
