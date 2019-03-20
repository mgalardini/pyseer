# Copyright 2018 Marco Galardini and John Lees

'''Predict phenotypes using a fitted elastic net model'''

import os
import sys
from .utils import set_env
from collections import deque
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import pickle
import pandas as pd
from scipy.special import expit
from tqdm import tqdm

from .input import open_variant_file
from .input import load_covariates
from .input import read_variant

def get_options():
    import argparse

    description = 'Predict phenotypes using a fitted elastic net model'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('model',
                        help='Name of fitted model pickle file (.pkl)')
    parser.add_argument('samples',
                        help='File with samples to predict')
    parser.add_argument('--threshold',
                        help='Threshold to pick binary predictions',
                        type=float,
                        default=0.5)

    variants = parser.add_argument_group('Variants')
    variant_group = variants.add_mutually_exclusive_group(required=True)
    variant_group.add_argument('--kmers',
                               default=None,
                               help='Kmers file')
    variant_group.add_argument('--vcf',
                               default=None,
                               help='VCF file. Will filter any non '
                                    '\'PASS\' sites')
    variant_group.add_argument('--pres',
                               default=None,
                               help='Presence/absence .Rtab matrix as '
                                    'produced by roary and piggy')
    variants.add_argument('--burden',
                          help='VCF regions to group variants by for burden'
                          ' testing (requires --vcf). '
                          'Requires vcf to be indexed')
    variants.add_argument('--uncompressed',
                       action='store_true',
                       default=False,
                       help='Uncompressed kmers file [Default: gzipped]')

    covariates = parser.add_argument_group('Covariates')
    covariates.add_argument('--covariates',
                            default=None,
                            help='User-defined covariates file '
                                 '(tab-delimited, no header, '
                                 'first column contains sample names)')
    covariates.add_argument('--use-covariates',
                            default=None,
                            nargs='*',
                            help='Covariates to use. Format is "2 3q 4" '
                                 '(q for quantitative) '
                                 ' [Default: load covariates but don\'t use '
                                 'them]')

    return parser.parse_args()

def main():

    options = get_options()

    # Read in model pickle
    with open(options.model, 'rb') as pickle_obj:
        model_dict, continuous = pickle.load(pickle_obj)
    try:
        intercept = model_dict.pop('intercept')[1]
    except KeyError as e:
        sys.stderr.write("Intercept not found in model\n")
        intercept = 0

    # Read in samples, start building predictions
    samples = []
    with open(options.samples, 'r') as sample_file:
        for sample in sample_file:
            samples.append(sample.rstrip())

    p = pd.DataFrame(data=np.full(len(samples), intercept),
                     index=samples,
                     columns=['prediction'])
    predictions = np.array(p.values, dtype=np.float)

    # Read in covariates
    if options.covariates is not None:
        cov = load_covariates(options.covariates,
                              options.use_covariates,
                              p)
        if cov is None:
            sys.exit(1)
        else:
            for covariate in cov:
                pred_beta = model_dict.pop(covariate, (0, 0))
                if pred_beta[1] != 0:
                    predictions += (cov[covariate] * pred_beta[1]).reshape(-1, 1)

    # Open variant file
    sample_order = []
    all_strains = set(p.index)
    burden_regions = deque([])
    burden = False

    if options.kmers:
        var_type = "kmers"
        var_file = options.kmers
    elif options.vcf:
        var_type = "vcf"
        var_file = options.vcf
        if options.burden:
            burden = True
    else:
        var_type = "Rtab"
        var_file = options.pres

    infile, sample_order = open_variant_file(var_type, var_file, options.burden, burden_regions, options.uncompressed)

    # Read in all variants - only keep if matching name
    sys.stderr.write("Reading variants from input\n")
    pbar = tqdm(unit="variants")
    while True:
        eof, k, var_name, kstrains, nkstrains, af, missing = read_variant(
                                        infile, p, var_type,
                                        burden, burden_regions,
                                        options.uncompressed, all_strains,
                                        sample_order, keep_list = model_dict.keys())

        # check for EOF
        if eof:
            pbar.close()
            break
        else:
            pbar.update(1)

        # return 0 if not found. remove from dict if found
        pred_beta = model_dict.pop(var_name, (0, 0))
        if pred_beta[1] != 0:
            predictions += (k * pred_beta[1]).reshape(-1, 1)

    # Note those variants which did not appear - impute
    for missing in model_dict.keys():
        sys.stderr.write("Could not find covariate/variant " + missing + " in input file\n")
        predictions += model_dict[missing][0] * model_dict[missing][1]

    # apply link function
    if not continuous:
        link = expit(predictions)
    else:
        link = predictions

    # output
    if continuous:
        print("\t".join(['Sample','Link','Prediction']))
    else:
        print("\t".join(['Sample','Prediction','Link','Probability']))

    p = pd.DataFrame(data=np.hstack((predictions, link)),
                     index=samples,
                     columns=['link', 'prediction'])
    for row in p.itertuples():
        if not continuous:
            binary_prediction = 0
            if row[2] >= options.threshold:
                binary_prediction = 1
            print("\t".join([row[0], str(binary_prediction), str(row[1]), str(row[2])]))
        else:
            print("\t".join([row[0], str(row[1]), str(row[2])]))

