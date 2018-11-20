import os
import unittest
import numpy as np
from pyseer.classes import Seer, LMM
from pyseer.utils import set_env
from pyseer.utils import format_output


class TestEnvSet(unittest.TestCase):
    def test_set_env(self):
        old_environ = dict(os.environ)
        os.environ.update({'PYSEER_TEST': '1'})
        with set_env(PYSEER_TEST='2'):
            new_environ = dict(os.environ)
            self.assertEqual(new_environ['PYSEER_TEST'], '2')
        new_environ = dict(os.environ)
        self.assertEqual(new_environ['PYSEER_TEST'], '1')
        os.environ.update(old_environ)


class TestOutputFormatter(unittest.TestCase):
    variant = 'AAAAAAAAAAAGCATTTTACTATTTTA'
    pattern = 'fake_hash_for_testing'
    af = ('1.25E-01', 0.125)
    filter_p = ('9.14E-01', 0.914)
    lrt_p = ('3.24E-01', 0.324)
    beta = ('-5.93E-01', -0.593)
    beta_err = ('6.09E-01', 0.609)
    intercept = ('2.61E-01', 0.261)
    pc1 = ('-1.65E+00', -1.65)
    pc2 = ('-5.73E-01', -0.573)
    pc3 = ('2.15E+00', 2.15)
    pc4 = ('3.05E+00', 3.05)
    pc5 = ('-1.71E+00', -1.71)
    max_lineage = 0
    lineage_dict = ['MDS1', ]
    notes = ''
    kstrains = ('1,2,3,4', ['1', '2', '3', '4'])
    nkstrains = ('5,6,7,8', ['5', '6', '7', '8'])
    bprefilter = True
    bfilter = True
    variant_h2 = ('5.44E-02', 0.0544)
    # fixed effects output
    out1 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      intercept[0],
                      pc1[0], pc2[0], pc3[0], pc4[0], pc5[0],
                      notes))
    out1nan = '\t'.join((variant,
                         '',
                         filter_p[0],
                         lrt_p[0],
                         beta[0],
                         beta_err[0],
                         intercept[0],
                         '', pc2[0], pc3[0], pc4[0], pc5[0],
                         notes))
    # fixed effects with samples
    out2 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      intercept[0],
                      pc1[0], pc2[0], pc3[0], pc4[0], pc5[0],
                      kstrains[0],
                      nkstrains[0],
                      notes))
    # fixed effects with lineage
    out3 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      intercept[0],
                      pc1[0], pc2[0], pc3[0], pc4[0], pc5[0],
                      lineage_dict[max_lineage],
                      notes))
    # fixed effects with lineage and samples
    out4 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      intercept[0],
                      pc1[0], pc2[0], pc3[0], pc4[0], pc5[0],
                      lineage_dict[max_lineage],
                      kstrains[0],
                      nkstrains[0],
                      notes))
    # random effects
    out5 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      variant_h2[0],
                      notes))
    # random effects with nan
    out5nan = '\t'.join((variant,
                         '',
                         filter_p[0],
                         lrt_p[0],
                         beta[0],
                         beta_err[0],
                         '',
                         notes))
    # random effects with samples
    out6 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      variant_h2[0],
                      kstrains[0],
                      nkstrains[0],
                      notes))
    # random effects with lineage
    out7 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      variant_h2[0],
                      lineage_dict[max_lineage],
                      notes))
    # random effects with lineage and samples
    out8 = '\t'.join((variant,
                      af[0],
                      filter_p[0],
                      lrt_p[0],
                      beta[0],
                      beta_err[0],
                      variant_h2[0],
                      lineage_dict[max_lineage],
                      kstrains[0],
                      nkstrains[0],
                      notes))
    fixed = Seer(variant,
                 pattern,
                 af[1],
                 filter_p[1],
                 lrt_p[1],
                 beta[1],
                 beta_err[1],
                 intercept[1],
                 [pc1[1], pc2[1], pc3[1],
                  pc4[1], pc5[1]],
                 None,
                 kstrains[1],
                 nkstrains[1],
                 notes,
                 bprefilter,
                 bfilter)
    fixed_nan = Seer(variant,
                     pattern,
                     np.nan,
                     filter_p[1],
                     lrt_p[1],
                     beta[1],
                     beta_err[1],
                     intercept[1],
                     [np.nan, pc2[1], pc3[1],
                      pc4[1], pc5[1]],
                     None,
                     kstrains[1],
                     nkstrains[1],
                     notes,
                     bprefilter,
                     bfilter)
    fixed_lineage = Seer(variant,
                         pattern,
                         af[1],
                         filter_p[1],
                         lrt_p[1],
                         beta[1],
                         beta_err[1],
                         intercept[1],
                         [pc1[1], pc2[1], pc3[1],
                          pc4[1], pc5[1]],
                         max_lineage,
                         kstrains[1],
                         nkstrains[1],
                         notes,
                         bprefilter,
                         bfilter)
    bad_fixed = Seer(variant,
                     pattern,
                     'gotcha',
                     filter_p[1],
                     lrt_p[1],
                     beta[1],
                     beta_err[1],
                     intercept[1],
                     [pc1[1], pc2[1], pc3[1],
                      pc4[1], pc5[1]],
                     None,
                     kstrains[1],
                     nkstrains[1],
                     notes,
                     bprefilter,
                     bfilter)
    random = LMM(variant,
                 pattern,
                 af[1],
                 filter_p[1],
                 lrt_p[1],
                 beta[1],
                 beta_err[1],
                 variant_h2[1],
                 None,
                 kstrains[1],
                 nkstrains[1],
                 notes,
                 bprefilter,
                 bfilter)
    random_nan = LMM(variant,
                     pattern,
                     np.nan,
                     filter_p[1],
                     lrt_p[1],
                     beta[1],
                     beta_err[1],
                     np.nan,
                     None,
                     kstrains[1],
                     nkstrains[1],
                     notes,
                     bprefilter,
                     bfilter)
    random_lineage = LMM(variant,
                         pattern,
                         af[1],
                         filter_p[1],
                         lrt_p[1],
                         beta[1],
                         beta_err[1],
                         variant_h2[1],
                         max_lineage,
                         kstrains[1],
                         nkstrains[1],
                         notes,
                         bprefilter,
                         bfilter)
    bad_random = LMM(variant,
                     pattern,
                     'gotcha',
                     filter_p[1],
                     lrt_p[1],
                     beta[1],
                     beta_err[1],
                     variant_h2[1],
                     None,
                     kstrains[1],
                     nkstrains[1],
                     notes,
                     bprefilter,
                     bfilter)

    def test_formatting_fixed(self):
        self.assertEqual(format_output(self.fixed, lineage_dict=None,
                                       print_samples=False),
                         self.out1)
        self.assertEqual(format_output(self.fixed_nan, lineage_dict=None,
                                       print_samples=False),
                         self.out1nan)
        self.assertEqual(format_output(self.fixed, lineage_dict=None,
                                       print_samples=True),
                         self.out2)
        self.assertEqual(format_output(self.fixed_lineage, self.lineage_dict,
                                       print_samples=False),
                         self.out3)
        self.assertEqual(format_output(self.fixed_lineage, self.lineage_dict,
                                       print_samples=True),
                         self.out4)

        with self.assertRaises(TypeError):
            format_output(self.bad_fixed, [])

    def test_formatting_random(self):
        self.assertEqual(format_output(self.random, lineage_dict=None,
                                       model='lmm', print_samples=False),
                         self.out5)
        self.assertEqual(format_output(self.random_nan, lineage_dict=None,
                                       model='lmm', print_samples=False),
                         self.out5nan)
        self.assertEqual(format_output(self.random, lineage_dict=None,
                                       model='lmm', print_samples=True),
                         self.out6)
        self.assertEqual(format_output(self.random_lineage, self.lineage_dict,
                                       model='lmm', print_samples=False),
                         self.out7)
        self.assertEqual(format_output(self.random_lineage, self.lineage_dict,
                                       model='lmm', print_samples=True),
                         self.out8)

        with self.assertRaises(TypeError):
            format_output(self.bad_random, [])


if __name__ == '__main__':
    unittest.main()
