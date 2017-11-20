# Copyright 2017 Marco Galardini and John Lees

'''Python reimplementation of SEER for bacterial GWAS'''

import os
import sys
import gzip
import warnings
import itertools
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import pandas as pd
from multiprocessing import Pool

from .__init__ import __version__

from .input import iter_kmers
from .input import load_structure
from .input import load_phenotypes
from .input import load_covariates

from .model import binary  
from .model import fit_null
from .model import continuous

from .utils import format_output

def get_options():
    import argparse

    description = 'SEER (doi: 10.1038/ncomms12797), reimplemented in python'
    parser = argparse.ArgumentParser(description=description,
                                     prog='pyseer')

    parser.add_argument('kmers',
                        help='Kmers file')
    parser.add_argument('phenotypes',
                        help='Phenotypes file')
    parser.add_argument('distances',
                        help='Strains distance square matrix')

    parser.add_argument('--continuous',
                        action='store_true',
                        default=False,
                        help='Force continuous phenotype [Default: binary auto-detect]')
    parser.add_argument('--print-samples',
                        action='store_true',
                        default=False,
                        help='Print sample lists [Default: hide samples]')
    parser.add_argument('--min-af',
                        type=float,
                        default=0.01,
                        help='Minimum AF [Default: 0.01]')
    parser.add_argument('--max-af',
                        type=float,
                        default=0.99,
                        help='Maximum AF [Default: 0.99]')
    parser.add_argument('--filter-pvalue',
                        type=float,
                        default=1,
                        help='Prefiltering t-test pvalue threshold [Default: 1]')
    parser.add_argument('--lrt-pvalue',
                        type=float,
                        default=1,
                        help='Likelihood ratio test pvalue threshold [Default: 1]')
    parser.add_argument('--max-dimensions',
                        type=int,
                        default=10,
                        help='Maximum number of dimensions to consider after MDS [Default: 10]')
    parser.add_argument('--covariates',
                        default=None,
                        help='User-defined covariates file (tab-delimited, no header, ' +
                             'first column contains sample names)')
    parser.add_argument('--use-covariates',
                        default=None,
                        nargs='*',
                        help='Covariates to use. Format is "2 3q 4" (q for quantitative)'
                             ' [Default: load covariates but don\'t use them]')
    parser.add_argument('--uncompressed',
                        action='store_true',
                        default=False,
                        help='Uncompressed kmers file [Default: gzipped]')
    parser.add_argument('--cpu',
                        type=int,
                        default=1,
                        help='Processes [Default: 1]')
    parser.add_argument('--save-m',
                        help='Prefix for saving matrix decomposition')
    scree_group = parser.add_mutually_exclusive_group()
    scree_group.add_argument('--load-m',
                        help='Load an existing matrix decomposition')
    scree_group.add_argument('--scree-plot',
                        action='store_true',
                        default=False,
                        help='Output MDS eigenvalues for a scree plot (exits without performing associations)')

    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()
    
    # check some arguments here
    if options.max_dimensions < 1:
        sys.stderr.write('Minimum number of dimensions after MDS scaling is 1\n')
        sys.exit(1)
    if options.cpu > 1 and sys.version_info[0] < 3:
        sys.stderr.write('pyseer requires python version 3 or above ' +
                         'unless the number of threads is 1\n')
        sys.exit(1)

    # silence warnings
    warnings.filterwarnings('ignore')
    #
    
    # reading phenotypes
    p = load_phenotypes(options.phenotypes)

    # Check whether any non 0/1 phenotypes
    if not options.continuous:
        if p.values[(p.values != 0) & (p.values != 1)].size > 0:
            options.continuous = True
            sys.stderr.write("Detected continuous phenotype\n")
        else:
            sys.stderr.write("Detected binary phenotype\n")

    # reading genome distances
    if options.load_m and os.path.isfile(options.load_m):
        m = pd.read_pickle(options.load_m)
        m = m.loc[p.index]
    else:
        m, evals = load_structure(options.distances, p)
        if options.save_m:
            m.to_pickle(options.save_m + ".pkl")

        if options.scree_plot:
            print('PC\teigenvalue')
            for i, e in enumerate(evals):
                print('PC%d\t%.5f' % (i+1, e))
            sys.stderr.write('Exiting without performing associations\n')
            sys.exit(0)
    if options.max_dimensions > m.shape[1]:
        sys.stderr.write('Population MDS scaling restricted to '+
                         '%d dimensions instead of requested %d\n' % (m.shape[1],
                                                                      options.max_dimensions))
        options.max_dimensions = m.shape[1]
    m = m.values[:, :options.max_dimensions]

    all_strains = set(p.index)

    # read covariates
    if options.covariates is not None:
        cov = load_covariates(options.covariates,
                              options.use_covariates,
                              p)
        if cov is None:
            sys.exit(1)
    else:
        cov = pd.DataFrame([])

    header = ['kmer', 'af', 'filter-pvalue',
              'lrt-pvalue', 'beta', 'beta-std-err',
              'intercept'] + ['PC%d' % i for i in range(1, options.max_dimensions+1)]
    if options.covariates is not None:
        header = header + [x for x in cov.columns]
    if options.print_samples:
        header = header + ['k-samples', 'nk-samples']
    header += ['notes']
    print('\t'.join(header))

    if options.uncompressed:
        infile = open(options.kmers)
    else:
        infile = gzip.open(options.kmers, 'r')

    # multiprocessing setup
    if options.cpu > 1:
        pool = Pool(options.cpu)

    # calculate null regressions once
    null_fit = fit_null(p, m, cov, options.continuous)
    if not options.continuous:
        firth_null = fit_null(p, m, cov, options.continuous, True)
    else:
        firth_null = True

    if null_fit is None or firth_null is None:
        sys.stderr.write('Could not fit null model, exiting\n')
        sys.exit(1)

    # iterator over each kmer
    # implements maf filtering
    k_iter = iter_kmers(p, m, cov,
                        infile, all_strains,
                        options.min_af, options.max_af,
                        options.filter_pvalue,
                        options.lrt_pvalue, null_fit, firth_null,
                        options.uncompressed)

    # keep track of the number of the total number of kmers and tests
    prefilter = 0
    tested = 0
    printed = 0

    # actual association test
    if options.cpu > 1:
        # multiprocessing proceeds 1000 kmers per core at a time
        if options.continuous:
            while True:
                ret = pool.starmap(continuous,
                                   itertools.islice(k_iter,
                                                    options.cpu*1000))
                if not ret:
                    break
                for x in ret:
                    if x.prefilter:
                        prefilter += 1
                        continue
                    tested += 1
                    if x.filter:
                        continue
                    printed += 1
                    print(format_output(x,
                                        options.print_samples))
        else:
            while True:
                ret = pool.starmap(binary,
                                   itertools.islice(k_iter,
                                                    options.cpu*1000))
                if not ret:
                    break
                for x in ret:
                    if x.prefilter:
                        prefilter += 1
                        continue
                    tested += 1
                    if x.filter:
                        continue
                    printed += 1
                    print(format_output(x,
                                        options.print_samples))
    else:
        for data in k_iter:
            if options.continuous:
                ret = continuous(*data)
            else:
                ret = binary(*data)
            if ret.prefilter:
                prefilter += 1
                continue
            tested += 1
            if ret.filter:
                continue
            printed += 1
            print(format_output(ret,
                                options.print_samples))

    sys.stderr.write('%d loaded kmers\n' % (prefilter + tested))
    sys.stderr.write('%d filtered kmers\n' % prefilter)
    sys.stderr.write('%d tested kmers\n' % tested)
    sys.stderr.write('%d printed kmers\n' % printed)


if __name__ == "__main__":
    main()
