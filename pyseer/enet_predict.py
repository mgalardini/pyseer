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
from .input import load_lineage
from .input import load_phenotypes

from .enet import write_lineage_predictions

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
    parser.add_argument('--lineage-clusters',
                        help='Custom clusters to use as lineages '
                             'to report stratified accuracy')
    parser.add_argument('--true-values',
                        help='Pheno file with known phenotypes '
                             'to calculate accuracy',
                        default=None)
    parser.add_argument('--ignore-missing',
                        help='Treat missing values as REF/0 rather than '
                             'using the mean AF',
                        action='store_true',
                        default=False)

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

    # Read in lineages
    if options.lineage_clusters:
        lineage_clusters, lineage_dict = load_lineage(options.lineage_clusters, p)
        fold_ids = np.where(lineage_clusters == 1)[1]
    else:
        lineage_clusters, lineage_dict, fold_ids = (None, None, None)

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
        if eof or len(model_dict.keys()) == 0:
            pbar.close()
            break
        else:
            pbar.update(1)

        # return 0 if not found. remove from dict if found
        (pred_af, pred_beta) = model_dict.pop(var_name, (0, 0))
        if pred_beta != 0:
            # model is fitted to minor allele encoded variants so flip obs to make compatible
            if pred_af > 0.5:
                k = np.array(~np.array(k, dtype=bool), dtype=np.int64)
            predictions += (k * pred_beta).reshape(-1, 1)

    # Note those variants which did not appear - impute
    for missing in model_dict.keys():
        sys.stderr.write("Could not find covariate/variant " + missing + " in input file\n")
        if not options.ignore_missing:
            predictions += model_dict[missing][0] * model_dict[missing][1]

    # apply link function
    if not continuous:
        link = predictions
        predictions = expit(link)
        binary_predictions = np.zeros(predictions.shape[0])
        binary_predictions[np.where(predictions > options.threshold)[0]] = 1
    else:
        link = predictions

    # output
    if continuous:
        print("\t".join(['Sample','Link','Prediction']))
    else:
        print("\t".join(['Sample','Prediction','Link','Probability']))

    p = pd.DataFrame(data=np.hstack((link, predictions)),
                     index=samples,
                     columns=['link', 'prediction'])
    for row_idx, row in enumerate(p.itertuples()):
        if not continuous:
            print("\t".join([row[0], str(binary_predictions[row_idx]), str(row[1]), str(row[2])]))
        else:
            print("\t".join([row[0], str(row[1]), str(row[2])]))

    # report summary
    if options.true_values:
        y_true = load_phenotypes(options.true_values, None)
        intersecting_samples = p.index.intersection(y_true.index)
        y_true = y_true.loc[intersecting_samples]

        sys.stderr.write("Overall prediction accuracy\n")
        if not continuous:
            R2, confusion = write_lineage_predictions(y_true.values, binary_predictions, None,
                                        None, continuous, stderr_print=False)
            tn, fp, fn, tp = confusion[0]
            sys.stderr.write("R2: " + str(R2[0]) + "\n")
            sys.stderr.write("tn: " + str(tn) + "\n")
            sys.stderr.write("fp: " + str(fp) + "\n")
            sys.stderr.write("fn: " + str(fn) + "\n")
            sys.stderr.write("tp: " + str(tp) + "\n")
        else:
            R2, confusion = write_lineage_predictions(y_true.values, predictions, None,
                                None, continuous, stderr_print=False)
            sys.stderr.write("R2: " + str(R2[0]) + "\n")

        if fold_ids is not None:
            sys.stderr.write("Predictions within each lineage\n")
            if continuous:
                write_lineage_predictions(y_true.values, predictions, fold_ids,
                                        lineage_dict, continuous, stderr_print=True)
            else:
                write_lineage_predictions(y_true.values, binary_predictions, fold_ids,
                                        lineage_dict, continuous, stderr_print=True)
