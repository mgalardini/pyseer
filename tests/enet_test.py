import os
import sys
import gzip
import unittest
import numpy as np
import pandas as pd
from pysam import VariantFile
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from pyseer.enet import fit_enet
from pyseer.enet import load_all_vars
from pyseer.enet import correlation_filter
from pyseer.enet import find_enet_selected
from pyseer.model import fit_null


DATA_DIR = 'tests'
P = os.path.join(DATA_DIR, 'subset.pheno')
M = os.path.join(DATA_DIR, 'distances.tsv.gz')
KMER = os.path.join(DATA_DIR, 'kmers.gz')
PRES = os.path.join(DATA_DIR, 'presence_absence.Rtab.gz')
PRESSMALL = os.path.join(DATA_DIR, 'presence_absence_smaller.Rtab')
VCF = os.path.join(DATA_DIR, 'variants_smaller.vcf.gz')
VENET = os.path.join(DATA_DIR, 'unit_tests_data', 'enet_variants.txt')
PFIRTH = os.path.join(DATA_DIR, 'unit_tests_data', 'p_binary.txt')
PFIRTHC = os.path.join(DATA_DIR, 'unit_tests_data', 'p_continuous.txt')
MFIRTH = os.path.join(DATA_DIR, 'unit_tests_data', 'm.txt')


def open_rtab(fname, compressed=True):
    if compressed:
        infile = gzip.open(fname)
        header = infile.readline().decode().rstrip()
    else:
        infile = open(fname)
        header = infile.readline().rstrip()
    sample_order = header.split()[1:]
    return infile, sample_order


class TstFindEnetSelected(unittest.TestCase):
    def test_find_enet_selected_binary(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        b = np.array([5.60000000e-01, 2.16450216e-37, -1.81966726e-37,
                      1.99034682e-38, -1.34400134e-37, -2.94000294e-38,
                      -1.90213827e-37, 8.08080808e-38,
                      1.89393939e-37])
        idx = [10, 29, 39, 95, 110, 153, 156, 164]
        g = find_enet_selected(b, idx, p, np.array([[]]), 'vcf',
                               None, False, None,
                               VariantFile(VCF), set(p.index),
                               None, False, False, None, False)
        v = next(g)
        self.assertEqual(v.kmer, 'FM211187_83_G_A')
        self.assertEqual(v.af, 0.28)
        self.assertTrue(abs(v.prep - 0.17050715825327736) < 1E-7)
        self.assertTrue(np.isnan(v.pvalue))
        self.assertTrue(abs(v.kbeta - 2.164502164502162e-37) < 1E-7)
        self.assertEqual(v.max_lineage, None)
        self.assertEqual(v.kstrains,
                        [
                         'sample_10', 'sample_11', 'sample_13', 'sample_18',
                         'sample_19',
                         'sample_20', 'sample_23', 'sample_25', 'sample_26',
                         'sample_31', 'sample_34', 'sample_36', 'sample_40',
                         'sample_45'
                        ]
                        )
        self.assertEqual(v.nkstrains,
                        [
                        'sample_1', 'sample_12', 'sample_14', 'sample_15',
                        'sample_16', 'sample_17', 'sample_2', 'sample_21',
                        'sample_22', 'sample_24', 'sample_27', 'sample_28',
                        'sample_29', 'sample_3', 'sample_30', 'sample_32',
                        'sample_33', 'sample_35', 'sample_37', 'sample_38',
                        'sample_39', 'sample_4', 'sample_41', 'sample_42',
                        'sample_43', 'sample_44', 'sample_46', 'sample_47',
                        'sample_48', 'sample_49', 'sample_5', 'sample_50',
                        'sample_6', 'sample_7', 'sample_8', 'sample_9'
                        ]
                        )
        self.assertEqual(len(v.notes), 0)
        # read to exhaustion
        for v in g:
            pass
        self.assertEqual(v.kmer, 'FM211187_3592_G_A')
        # with fixed effects
        pf = np.loadtxt(PFIRTH)
        mf = np.loadtxt(MFIRTH)
        cov = pd.DataFrame([])
        null_res = fit_null(pf, mf, cov, False, firth=False).llr
        null_firth = fit_null(pf, mf, cov, False, firth=True)
        g = find_enet_selected(b, idx, p, np.array([[]]), 'vcf',
                               (mf, null_res, null_firth), False,
                               None, VariantFile(VCF), set(p.index),
                               None, False, False, None, False)
        v = next(g)
        self.assertEqual(v.kmer, 'FM211187_83_G_A')
        self.assertEqual(v.af, 0.28)
        self.assertTrue(abs(v.prep - 0.17050715825327736) < 1E-7)
        self.assertEqual(v.pvalue, 1)
        self.assertTrue(abs(v.kbeta - 2.164502164502162e-37) < 1E-7)
        self.assertEqual(v.max_lineage, None)
        self.assertEqual(v.kstrains,
                        [
                         'sample_10', 'sample_11', 'sample_13', 'sample_18',
                         'sample_19',
                         'sample_20', 'sample_23', 'sample_25', 'sample_26',
                         'sample_31', 'sample_34', 'sample_36', 'sample_40',
                         'sample_45'
                        ]
                        )
        self.assertEqual(v.nkstrains,
                        [
                        'sample_1', 'sample_12', 'sample_14', 'sample_15',
                        'sample_16', 'sample_17', 'sample_2', 'sample_21',
                        'sample_22', 'sample_24', 'sample_27', 'sample_28',
                        'sample_29', 'sample_3', 'sample_30', 'sample_32',
                        'sample_33', 'sample_35', 'sample_37', 'sample_38',
                        'sample_39', 'sample_4', 'sample_41', 'sample_42',
                        'sample_43', 'sample_44', 'sample_46', 'sample_47',
                        'sample_48', 'sample_49', 'sample_5', 'sample_50',
                        'sample_6', 'sample_7', 'sample_8', 'sample_9'
                        ]
                        )
        self.assertEqual(len(v.notes), 0)
        # read to exhaustion
        for v in g:
            pass
        self.assertEqual(v.kmer, 'FM211187_3592_G_A')


    def test_find_enet_selected_continuous(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['continuous']
        b = np.array([5.60000000e-01, 2.16450216e-37, -1.81966726e-37,
                      1.99034682e-38, -1.34400134e-37, -2.94000294e-38,
                      -1.90213827e-37, 8.08080808e-38,
                      1.89393939e-37])
        idx = [10, 29, 39, 95, 110, 153, 156, 164]
        g = find_enet_selected(b, idx, p, np.array([[]]), 'vcf',
                               None, False, None,
                               VariantFile(VCF), set(p.index),
                               None, True, False, None, False)
        v = next(g)
        self.assertEqual(v.kmer, 'FM211187_83_G_A')
        self.assertEqual(v.af, 0.28)
        self.assertTrue(abs(v.prep - 0.8807556966503836) < 1E-7)
        self.assertTrue(np.isnan(v.pvalue))
        self.assertTrue(abs(v.kbeta - 2.164502164502162e-37) < 1E-7)
        self.assertEqual(v.max_lineage, None)
        self.assertEqual(v.kstrains,
                        [
                         'sample_10', 'sample_11', 'sample_13', 'sample_18',
                         'sample_19',
                         'sample_20', 'sample_23', 'sample_25', 'sample_26',
                         'sample_31', 'sample_34', 'sample_36', 'sample_40',
                         'sample_45'
                        ]
                        )
        self.assertEqual(v.nkstrains,
                        [
                        'sample_1', 'sample_12', 'sample_14', 'sample_15',
                        'sample_16', 'sample_17', 'sample_2', 'sample_21',
                        'sample_22', 'sample_24', 'sample_27', 'sample_28',
                        'sample_29', 'sample_3', 'sample_30', 'sample_32',
                        'sample_33', 'sample_35', 'sample_37', 'sample_38',
                        'sample_39', 'sample_4', 'sample_41', 'sample_42',
                        'sample_43', 'sample_44', 'sample_46', 'sample_47',
                        'sample_48', 'sample_49', 'sample_5', 'sample_50',
                        'sample_6', 'sample_7', 'sample_8', 'sample_9'
                        ]
                        )
        self.assertEqual(len(v.notes), 0)
        # read to exhaustion
        for v in g:
            pass
        self.assertEqual(v.kmer, 'FM211187_3592_G_A')
        # binary
        with self.assertRaises(ValueError):
            next(find_enet_selected(b, idx, p, np.array([[]]), 'vcf',
                                    None, False, None,
                                    VariantFile(VCF), set(p.index),
                                    None, False, False, None, False))
        # with fixed effects
        pf = np.loadtxt(PFIRTHC)
        mf = np.loadtxt(MFIRTH)
        cov = pd.DataFrame([])
        null_res = fit_null(pf, mf, cov, True, firth=False)
        g = find_enet_selected(b, idx, p, np.array([[]]), 'vcf',
                               (mf, null_res, None), False,
                               None, VariantFile(VCF), set(p.index),
                               None, True, False, None, False)
        v = next(g)
        self.assertEqual(v.kmer, 'FM211187_83_G_A')
        self.assertEqual(v.af, 0.28)
        self.assertTrue(abs(v.prep - 0.8807556966503836) < 1E-7)
        self.assertTrue(abs(v.pvalue - 0.8984215870932599) < 1E-7)
        self.assertTrue(abs(v.kbeta - 2.164502164502162e-37) < 1E-7)
        self.assertEqual(v.max_lineage, None)
        self.assertEqual(v.kstrains,
                        [
                         'sample_10', 'sample_11', 'sample_13', 'sample_18',
                         'sample_19',
                         'sample_20', 'sample_23', 'sample_25', 'sample_26',
                         'sample_31', 'sample_34', 'sample_36', 'sample_40',
                         'sample_45'
                        ]
                        )
        self.assertEqual(v.nkstrains,
                        [
                        'sample_1', 'sample_12', 'sample_14', 'sample_15',
                        'sample_16', 'sample_17', 'sample_2', 'sample_21',
                        'sample_22', 'sample_24', 'sample_27', 'sample_28',
                        'sample_29', 'sample_3', 'sample_30', 'sample_32',
                        'sample_33', 'sample_35', 'sample_37', 'sample_38',
                        'sample_39', 'sample_4', 'sample_41', 'sample_42',
                        'sample_43', 'sample_44', 'sample_46', 'sample_47',
                        'sample_48', 'sample_49', 'sample_5', 'sample_50',
                        'sample_6', 'sample_7', 'sample_8', 'sample_9'
                        ]
                        )
        self.assertEqual(len(v.notes), 0)
        # read to exhaustion
        for v in g:
            pass
        self.assertEqual(v.kmer, 'FM211187_3592_G_A')

class TestCorrelationFilter(unittest.TestCase):
    def test_filter_binary(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        a = np.loadtxt(VENET)
        a = csr_matrix(a.T)
        f = correlation_filter(p, a, 0.75)
        self.assertTrue(abs(f - np.array([0, 5])).max() < 1E-7)
        # variant absent in all
        a = csr_matrix(np.zeros(a.T.shape))
        f = correlation_filter(p, a, 0.75)
        self.assertEqual(f.shape[0], 0)

    def test_filter_continuous(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['continuous']
        a = np.loadtxt(VENET)
        a = csr_matrix(a.T)
        f = correlation_filter(p, a, 0.75)
        self.assertTrue(abs(f - np.array([1, 2])).max() < 1E-7)
        # variant absent in all
        a = csr_matrix(np.zeros(a.T.shape))
        f = correlation_filter(p, a, 0.75)
        self.assertEqual(f.shape[0], 0)


class TestFitEnet(unittest.TestCase):
    def test_fit_binary(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        a = np.loadtxt(VENET)
        a = csc_matrix(a)
        weights = np.ones((p.shape[0], 1))
        # alpha = 1
        b = fit_enet(p, a, pd.DataFrame([]), weights, False, 1)
        self.assertTrue(b.sum() - 0.24116205681688876 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     0.24116206, 0., 0., 0., 0.,
                                     0., 0., 0., 0.
                                     ])).max() < 1E-7)
        # alpha = 0
        b = fit_enet(p, a, pd.DataFrame([]), weights, False, 0)
        self.assertTrue(b.sum() - 0.24116205681688876 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     2.41162057e-01, 2.16450216e-37,
                                     -1.81966726e-37, 1.99034682e-38,
                                     -1.34400134e-37, -2.94000294e-38,
                                     -1.90213827e-37, 8.08080808e-38,
                                     1.89393939e-37
                                     ])).max() < 1E-7)
        # alpha = 0.5
        b = fit_enet(p, a, pd.DataFrame([]), weights, False, 0.5)
        self.assertTrue(b.sum() - 0.24116205681688876 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     0.24116206, 0., 0., 0., 0.,
                                     0., 0., 0., 0.
                                     ])).max() < 1E-7)
        # alpha = 0.5, continuous
        b = fit_enet(p, a, pd.DataFrame([]), weights, True, 0.5)
        self.assertTrue(b.sum() - 0.5600000000000002 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     0.56, 0., 0., 0., 0., 0., 0., 0., 0.
                                     ])).max() < 1E-7)

    def test_fit_continuous(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['continuous']
        a = np.loadtxt(VENET)
        a = csc_matrix(a)
        weights = np.ones((p.shape[0], 1))
        # alpha = 1
        b = fit_enet(p, a, pd.DataFrame([]), weights, True, 1)
        self.assertTrue(b.sum() - 25.5 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     25.5, 0., 0., 0., 0.,
                                     0., 0., 0., 0.
                                     ])).max() < 1E-7)
        # alpha = 0
        b = fit_enet(p, a, pd.DataFrame([]), weights, True, 0)
        self.assertTrue(b.sum() - 25.5 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     2.55000000e+01, -6.01250601e-37,
                                     -5.19904932e-36, 5.43198819e-36,
                                     2.15250215e-36, 1.52250152e-36,
                                     -5.73921028e-37, 3.51515152e-36,
                                     1.76396316e-36
                                     ])).max() < 1E-7)
        # alpha = 0.5
        b = fit_enet(p, a, pd.DataFrame([]), weights, True, 0.5)
        self.assertTrue(b.sum() - 25.5 < 1E-7)
        self.assertTrue(abs(b -
                            np.array([
                                     25.5, 0., 0., 0., 0.,
                                     0., 0., 0., 0.
                                     ])).max() < 1E-7)
        # binary
        with self.assertRaises(ValueError):
            b = fit_enet(p, a, pd.DataFrame([]), weights, False, 0.5)


class TestLoadAllVars(unittest.TestCase):
    def test_unsupported(self):
        with self.assertRaises(ValueError):
            load_all_vars('test', None, None, None, None,
                          None, None, None, None, None, None)

    def test_no_file(self):
        with self.assertRaises(AttributeError):
            load_all_vars('kmers', None, None, None, None,
                          None, None, None, None, None, None)
        with self.assertRaises(TypeError):
            load_all_vars('vcf', None, None, None, None,
                          None, None, None, None, None, None)
        with self.assertRaises(AttributeError):
            load_all_vars('Rtab', None, None, None, None,
                          None, None, None, None, None, None)

    def test_load_all_vars_kmer(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        infile = gzip.open(KMER)
        variants, sidx, vidx  = load_all_vars('kmers', p, False, None,
                                              infile, set(p.index), None,
                                              0.45, 0.55, 1.0, False)
        self.assertEqual(variants.shape, (20, 50))
        self.assertEqual(variants.sum(), 474.0)
        self.assertTrue(abs(variants.toarray()[0] -
                            np.array([1., 1., 0., 1., 0., 0., 0., 0., 1.,
                                      0., 0., 1., 1., 0., 1., 1., 0.,
                                      1., 1., 1., 0., 0., 0., 0., 0., 1.,
                                      0., 0., 1., 0., 0., 1., 1., 0.,
                                      1., 1., 1., 0., 1., 1., 0., 0., 0.,
                                      1., 0., 1., 1.,
                                      1., 0., 1.])).max() < 1E-7)
        self.assertEqual(len(sidx), 20)
        self.assertEqual(sidx, [2, 6, 20, 32, 39, 54,
                                58, 60, 69, 89, 93,
                                123, 127, 134, 153,
                                156, 179, 180, 184, 194])
        self.assertEqual(vidx, 200)
        # not providing samples
        infile = gzip.open(KMER)
        with self.assertRaises(ZeroDivisionError):
            _  = load_all_vars('kmers', p, False, None,
                               infile, set([]), None,
                               0.45, 0.55, 1.0, False)
        # uncompressed option - only with python3+
        if sys.version_info[0] >= 3:
            infile = gzip.open(KMER)
            with self.assertRaises(TypeError):
                _  = load_all_vars('kmers', p, False, None,
                                   infile, set(p.index), None,
                                   0.45, 0.55, 1.0, True)
        # different type
        infile = gzip.open(KMER)
        with self.assertRaises(ValueError):
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set([]), None,
                               0.45, 0.55, 1.0, False)
        infile = gzip.open(KMER)
        with self.assertRaises(AttributeError):
            _  = load_all_vars('vcf', p, False, None,
                               infile, set([]), None,
                               0.45, 0.55, 1.0, False)
        # different file
        infile = gzip.open(PRES)
        with self.assertRaises(IndexError):
            _  = load_all_vars('kmers', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, False)
        infile = gzip.open(VCF)
        with self.assertRaises(IndexError):
            _  = load_all_vars('kmers', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, False)

    def test_load_all_vars_vcf(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        infile  = VariantFile(VCF)
        variants, sidx, vidx  = load_all_vars('vcf', p, False, None,
                                              infile, set(p.index),
                                              None,
                                              0.25, 0.75, 1.0, False)
        self.assertEqual(variants.shape, (8, 50))
        self.assertEqual(variants.sum(), 140)
        self.assertTrue(abs(variants.toarray()[0] -
                            np.array([
                                     0., 0., 0., 0., 0., 0., 0., 0., 0.,
                                     1., 1., 0., 1., 0., 0., 0., 0.,
                                     1., 1., 1., 0., 0., 1., 0., 1., 1.,
                                     0., 0., 0., 0., 1., 0., 0., 1.,
                                     0., 1., 0., 0., 0., 1., 0., 0., 0.,
                                     0., 1., 0., 0., 0., 0., 0.
                                     ])).max() < 1E-7)
        self.assertEqual(len(sidx), 8)
        self.assertEqual(sidx, [10, 29, 39, 95, 110, 153, 156, 164])
        self.assertEqual(vidx, 254)
        # not providing samples
        with self.assertRaises(ZeroDivisionError):
            infile  = VariantFile(VCF)
            _  = load_all_vars('vcf', p, False, None,
                               infile, set([]), None,
                               0.45, 0.55, 1.0, False)
        # different type
        with self.assertRaises(AttributeError):
            infile  = VariantFile(VCF)
            _  = load_all_vars('kmers', p, False, None,
                               infile, set(p.index), None,
                               0.01, 0.99, 1.0, True)
        with self.assertRaises(AttributeError):
            infile  = VariantFile(VCF)
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set(p.index), None,
                               0.01, 0.99, 1.0, False)
        # different file
        infile = gzip.open(KMER)
        with self.assertRaises(AttributeError):
            _  = load_all_vars('vcf', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, False)
        infile = gzip.open(PRES)
        with self.assertRaises(AttributeError):
            _  = load_all_vars('vcf', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, True)

    def test_load_all_vars_rtab(self):
        p = pd.read_csv(P,
                        index_col=0,
                        sep='\t')['binary']
        infile, sample_order = open_rtab(PRES)
        variants, sidx, vidx  = load_all_vars('Rtab', p, False, None,
                                              infile, set(p.index),
                                              sample_order,
                                              0.25, 0.75, 1.0, False)
        self.assertEqual(variants.shape, (7, 50))
        self.assertEqual(variants.sum(), 103.0)
        self.assertTrue(abs(variants.toarray()[0] -
                            np.array([0., 0., 1., 0., 0., 0., 1., 0., 1., 0.,
                                      0., 0., 1., 0., 0., 0., 0.,
                                      0., 1., 0., 0., 0., 0., 0., 0., 1., 0.,
                                      0., 1., 1., 0., 0., 0., 0.,
                                      0., 0., 0., 0., 1., 1., 0., 1., 0., 1.,
                                      0., 1., 0., 0., 1., 0.
                                      ])).max() < 1E-7)
        self.assertEqual(len(sidx), 7)
        self.assertEqual(sidx, [1426, 1436, 1463, 1484, 1492, 1496, 1498])
        self.assertEqual(vidx, 1499)
        # too few OGs
        with self.assertRaises(ValueError):
            infile, sample_order = open_rtab(PRESSMALL, compressed=False)
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set(p.index), sample_order,
                               0.01, 0.99, 1.0, False)
        # not providing samples
        with self.assertRaises(ValueError):
            infile, sample_order = open_rtab(PRESSMALL, compressed=False)
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set([]), [],
                               0.45, 0.55, 1.0, False)
        # different type
        with self.assertRaises(IndexError):
            infile, sample_orders = open_rtab(PRES)
            _  = load_all_vars('kmers', p, False, None,
                               infile, set(p.index), sample_order,
                               0.01, 0.99, 1.0, False)
        with self.assertRaises(AttributeError):
            infile, sample_orders = open_rtab(PRES)
            _  = load_all_vars('vcf', p, False, None,
                               infile, set(p.index), sample_order,
                               0.01, 0.99, 1.0, False)
        # different file
        infile = gzip.open(KMER)
        with self.assertRaises(ValueError):
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, False)
        infile = gzip.open(VCF)
        with self.assertRaises(ValueError):
            _  = load_all_vars('Rtab', p, False, None,
                               infile, set(p.index), None,
                               0.45, 0.55, 1.0, False)


if __name__ == '__main__':
    unittest.main()
