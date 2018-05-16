# Copyright 2017 Marco Galardini and John Lees

'''Python reimplementation of SEER for bacterial GWAS'''

import os
import sys
import gzip
import warnings
import itertools
import operator
import re
from collections import deque
from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
from scipy.stats import norm
import pandas as pd
from sklearn import manifold
from multiprocessing import Pool
from pysam import VariantFile

from .__init__ import __version__

from .input import load_phenotypes
from .input import load_structure
from .input import load_lineage
from .input import load_covariates
from .input import load_burden
from .input import iter_variants
from .input import load_var_block
from .input import iter_variants_lmm

from .model import fixed_effects_regression
from .model import fit_null

from .lmm import initialise_lmm
from .lmm import fit_lmm

from .utils import format_output


def get_options():
    import argparse

    description = 'SEER (doi: 10.1038/ncomms12797), reimplemented in python'
    parser = argparse.ArgumentParser(description=description,
                                     prog='pyseer')

    phenotypes = parser.add_argument_group('Phenotype')
    phenotypes.add_argument('--phenotypes',
                            required=True,
                            help='Phenotypes file')
    phenotypes.add_argument('--phenotype-column',
                            default=None,
                            help='Phenotype file column to use '
                                 '[Default: last column]')

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

    distances = parser.add_argument_group('Distances')
    distance_group = distances.add_mutually_exclusive_group()
    distance_group.add_argument('--distances',
                                help='Strains distance square matrix '
                                     '(fixed or lineage effects)')
    distance_group.add_argument('--load-m',
                                help='Load an existing matrix decomposition')
    similarity_group = distances.add_mutually_exclusive_group()
    similarity_group.add_argument('--similarity',
                                  help='Strains similarity square matrix '
                                       '(for --lmm)')
    similarity_group.add_argument('--load-lmm',
                                  help='Load an existing lmm cache')
    distances.add_argument('--save-m',
                           help='Prefix for saving matrix decomposition')
    distances.add_argument('--save-lmm',
                           help='Prefix for saving LMM cache')
    distances.add_argument('--mds',
                           default="classic",
                           choices=['classic', 'metric', 'non-metric'],
                           help='Type of multidimensional scaling '
                                '[Default: classic]')
    distances.add_argument('--max-dimensions',
                           type=int,
                           default=10,
                           help='Maximum number of dimensions to consider '
                                'after MDS [Default: 10]')
    distances.add_argument('--no-distances',
                           action='store_true',
                           default=False,
                           help='Allow run without a distance '
                                'matrix')

    association = parser.add_argument_group('Association options')
    association.add_argument('--continuous',
                             action='store_true',
                             default=False,
                             help='Force continuous phenotype '
                                  '[Default: binary auto-detect]')
    association.add_argument('--lmm',
                             action='store_true',
                             default=False,
                             help='Use random instead of fixed effects '
                                  'to correct for population structure. '
                                  'Requires a similarity matrix')
    association.add_argument('--lineage',
                             action='store_true',
                             help='Report lineage effects')
    association.add_argument('--lineage-clusters',
                             help='Custom clusters to use as lineages '
                                  '[Default: MDS components]')
    association.add_argument('--lineage-file',
                             default="lineage_effects.txt",
                             help='File to write lineage association to '
                                  '[Default: lineage_effects.txt]')

    filtering = parser.add_argument_group('Filtering options')
    filtering.add_argument('--min-af',
                           type=float,
                           default=0.01,
                           help='Minimum AF [Default: 0.01]')
    filtering.add_argument('--max-af',
                           type=float,
                           default=0.99,
                           help='Maximum AF [Default: 0.99]')
    filtering.add_argument('--filter-pvalue',
                           type=float,
                           default=1,
                           help='Prefiltering t-test pvalue threshold '
                                '[Default: 1]')
    filtering.add_argument('--lrt-pvalue',
                           type=float,
                           default=1,
                           help='Likelihood ratio test pvalue threshold '
                                '[Default: 1]')

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

    other = parser.add_argument_group('Other')
    other.add_argument('--print-samples',
                       action='store_true',
                       default=False,
                       help='Print sample lists [Default: hide samples]')
    other.add_argument('--print-filtered',
                       action='store_true',
                       default=False,
                       help='Print filtered variants (i.e. fitting errors) [Default: hide them]')
    other.add_argument('--output-patterns',
                       default=False,
                       help='File to print patterns to, useful for finding '
                            'pvalue threshold')
    other.add_argument('--uncompressed',
                       action='store_true',
                       default=False,
                       help='Uncompressed kmers file [Default: gzipped]')
    other.add_argument('--cpu',
                       type=int,
                       default=1,
                       help='Processes [Default: 1]')
    other.add_argument('--block_size',
                       type=int,
                       default=3000,
                       help='Number of variants per core [Default: 3000]')

    other.add_argument('--version', action='version',
                       version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    # check some arguments here
    if options.max_dimensions < 1:
        sys.stderr.write('Minimum number of dimensions after MDS is 1\n')
        sys.exit(1)
    if options.cpu > 1 and sys.version_info[0] < 3:
        sys.stderr.write('pyseer requires python version 3 or above ' +
                         'unless the number of threads is 1\n')
        sys.exit(1)
    if options.burden and not options.vcf:
        sys.stderr.write('Burden test can only be performed with VCF input\n')
        sys.exit(1)
    if not options.no_distances:
        if (options.lmm and (options.distances or options.load_m) and not options.lineage) or (not options.lmm and (options.similarity or options.load_lmm)):
            sys.stderr.write('Must use distance matrix with fixed effects, or similarity matrix with random effects\n')
            sys.stderr.write('Unless performing a lineage analysis with random effects\n')
            sys.exit(1)
        if (options.lmm and not (options.distances or options.load_m) and options.lineage):
            sys.stderr.write('Must also provide a distance matrix to report lineage effects\n')
            sys.exit(1)
    else:
        if options.distances or options.load_m:
            sys.stderr.write('Cannot use --no-distances with --distances or --load-m\n')
            sys.exit(1)
        if options.lmm:
            sys.stderr.write('Cannot use --no-distances with --lmm\n')
            sys.exit(1)
    if (options.block_size < 1):
        sys.stderr.write('Block size must be at least 1\n')
        sys.exit(1)

    # silence warnings
    warnings.filterwarnings('ignore')
    #

    # reading phenotypes
    p = load_phenotypes(options.phenotypes, options.phenotype_column)
    sys.stderr.write("Read " + str(len(p)) + " phenotypes\n")

    # Check whether any non 0/1 phenotypes
    if not options.continuous:
        if p.values[(p.values != 0) & (p.values != 1)].size > 0:
            options.continuous = True
            sys.stderr.write("Detected continuous phenotype\n")
        else:
            sys.stderr.write("Detected binary phenotype\n")

    # read covariates
    if options.covariates is not None:
        cov = load_covariates(options.covariates,
                              options.use_covariates,
                              p)
        if cov is None:
            sys.exit(1)
    else:
        cov = pd.DataFrame([])

    # fixed effects or lineage effects require regressing p ~ m
    if (options.lineage and not options.lineage_clusters) or not options.lmm:
        # reading genome distances
        if not options.no_distances:
            if options.load_m and os.path.isfile(options.load_m):
                m = pd.read_pickle(options.load_m)
                sys.stderr.write("Loaded projection with dimension " + str(m.shape) + "\n")
            else:
                # see if we have setup a seed for non-classical mds
                # a bit of a ugly hack for testing
                seed = os.environ.get('PYSEERSEED', None)
                if seed is not None:
                    seed = int(seed)

                m = load_structure(options.distances, p, options.max_dimensions,
                                   options.mds, options.cpu, seed)
                if options.save_m:
                    m.to_pickle(options.save_m + ".pkl")

            if options.max_dimensions > m.shape[1]:
                sys.stderr.write('Population MDS scaling restricted to ' +
                                 '%d dimensions instead of requested %d\n' %
                                 (m.shape[1],
                                  options.max_dimensions))
                options.max_dimensions = m.shape[1]

            intersecting_samples = p.index.intersection(m.index)
            sys.stderr.write("Analysing " + str(len(intersecting_samples)) + " samples"
                             " found in both phenotype and structure matrix\n")
            p = p.loc[intersecting_samples]

            m = m.loc[p.index]
            m = m.values[:, :options.max_dimensions]
        else:
            m = np.empty(shape=(0, 0))

        if cov.shape[1] > 0:
            cov = cov.loc[p.index]

        # calculate null regressions once
        null_fit = fit_null(p.values, m, cov, options.continuous)
        if not options.continuous and not options.lmm:
            firth_null = fit_null(p.values, m, cov, options.continuous, True)
        else:
            firth_null = True

        if null_fit is None or firth_null is None:
            sys.stderr.write('Could not fit null model, exiting\n')
            sys.exit(1)

    # lineage effects using null model - read BAPS clusters and fit pheno ~ lineage
    lineage_clusters = None
    lineage_samples = None
    if options.lineage:
        lineage_samples = p.index # this is ensured in load_lineage

        lineage_dict = []
        lineage_wald = {}
        if options.lineage_clusters:
            lineage_clusters, lineage_dict = load_lineage(options.lineage_clusters, p)

            # The lineage design matrix is not full rank, so one lineage
            # predictor needs to be removed. Do multiple single variable linear
            # regressions (as lineages are orthogonal, same as multiple linear
            # regression) and remove the one least associated with the phenotype
            for lineage, lineage_design in zip(lineage_dict, lineage_clusters.T):
                lineage_fit = fit_null(p.values, lineage_design.reshape(-1, 1), cov,
                                   options.continuous)
                if lineage_fit is None:
                    sys.stderr.write('Could not fit lineage null model, exiting\n')
                    sys.exit(1)

                lineage_wald[lineage] = np.absolute(lineage_fit.params[1])/lineage_fit.bse[1]

            min_lineage = min(lineage_wald.items(), key=operator.itemgetter(1))[0]
            min_index = lineage_dict.index(min_lineage)
            lineage_clusters = np.delete(lineage_clusters, min_index, 1)
        else:
            lineage_dict = ["MDS" + str(i+1)
                            for i in range(options.max_dimensions)]
            lineage_clusters = m
            lineage_fit = null_fit

            # Calculate lineage effects
            for lineage, slope, se in zip(lineage_dict, lineage_fit.params[1:],
                                      lineage_fit.bse[1:]):
                lineage_wald[lineage] = np.absolute(slope)/se

        # sort and print lineage effects
        sys.stderr.write('Writing lineage effects to %s\n' %
                         options.lineage_file)
        with open(options.lineage_file, 'w') as lineage_out:
            lineage_out.write("\t".join(["lineage", "wald_test", "p-value"]) + "\n")
            for lineage, wald in sorted(lineage_wald.items(),
                                        key=operator.itemgetter(1),
                                        reverse=True):
                pval = 2 * (1 - norm.cdf(wald))
                lineage_out.write("\t".join([lineage, str(wald), str(pval)]) + "\n")
    else:
        lineage_dict = None

    # binary regression takes LLF as null, not full model fit
    if not options.continuous and not options.lmm:
        null_fit = null_fit.llf

    # LMM setup - see _internal_single in fastlmm.association.single_snp
    if options.lmm:
        sys.stderr.write("Setting up LMM\n")
        p, lmm, h2 = initialise_lmm(p, cov, options.similarity, options.load_lmm,
                                 options.save_lmm, lineage_samples)
        sys.stderr.write("h^2 = " + '{0:.2f}'.format(h2) + "\n")

    # Open variant file
    sample_order = []
    all_strains = set(p.index)
    burden_regions = deque([])
    burden = False

    if options.kmers:
        var_type = "kmers"
        if options.uncompressed:
            infile = open(options.kmers)
        else:
            infile = gzip.open(options.kmers, 'r')
    elif options.vcf:
        var_type = "vcf"
        infile = VariantFile(options.vcf)
        if options.burden:
            burden = True
            load_burden(options.burden, burden_regions)
    else:
        # Rtab files have a header, rather than sample names accessible by row
        var_type = "Rtab"
        infile = open(options.pres)
        header = infile.readline().rstrip()
        sample_order = header.split()[1:]

    # keep track of the number of the total number of kmers and tests
    prefilter = 0
    tested = 0
    printed = 0

    # open pattern file if specified
    if options.output_patterns:
        patterns = open(options.output_patterns, 'wb')

    # header fields
    header = ['variant', 'af', 'filter-pvalue',
              'lrt-pvalue', 'beta', 'beta-std-err']

    if not options.lmm:
        header = header + ['intercept']

        if not options.no_distances:
            header = header + ['PC%d' % i
                                for i in range(1, options.max_dimensions+1)]
        if options.covariates is not None:
            header = header + [x for x in cov.columns]
    else:
        header = header + ['variant_h2']

    if options.lineage:
        header = header + ['lineage']
    if options.print_samples:
        header = header + ['k-samples', 'nk-samples']
    header += ['notes']
    print('\t'.join(header))

    # multiprocessing setup
    if options.cpu > 1:
        pool = Pool(options.cpu)

    # actual association test
    if not options.lmm:
        # iterator over each variant
        # implements maf filtering
        v_iter = iter_variants(p, m, cov, var_type, burden, burden_regions,
                               infile, all_strains, sample_order,
                               options.lineage, lineage_clusters,
                               options.min_af, options.max_af,
                               options.filter_pvalue,
                               options.lrt_pvalue, null_fit, firth_null,
                               options.uncompressed, options.continuous)

        if options.cpu > 1:
            # multiprocessing proceeds X variants per core at a time
            while True:
                ret = pool.starmap(fixed_effects_regression,
                                   itertools.islice(
                                                v_iter,
                                                options.cpu*options.block_size))
                if not ret:
                    break
                for x in ret:
                    if x.prefilter:
                        prefilter += 1
                        continue
                    tested += 1
                    if options.output_patterns:
                        patterns.write(x.pattern)

                    if x.filter and not options.print_filtered:
                        continue
                    printed += 1
                    print(format_output(x,
                                        lineage_dict,
                                        options.lmm,
                                        options.print_samples))
        else:
            for data in v_iter:
                ret = fixed_effects_regression(*data)

                if ret.prefilter:
                    prefilter += 1
                    continue
                tested += 1
                if options.output_patterns:
                    patterns.write(ret.pattern)

                if ret.filter and not options.print_filtered:
                    continue
                printed += 1
                print(format_output(ret,
                                    lineage_dict,
                                    options.lmm,
                                    options.print_samples))
    else:
        v_iter = load_var_block(var_type, p, burden, burden_regions,
                                infile, all_strains, sample_order,
                                options.min_af, options.max_af,
                                options.uncompressed, options.block_size)
        lmm_iter = iter_variants_lmm(v_iter, lmm, h2,
                                     options.lineage, lineage_clusters,
                                     cov.values,
                                     options.continuous,
                                     options.filter_pvalue,
                                     options.lrt_pvalue)
        if options.cpu > 1:
            while True:
                ret = pool.starmap(fit_lmm,
                                   itertools.islice(
                                                lmm_iter,
                                                options.cpu))
                if not ret:
                    break
                for values in ret:
                    for x in values:
                        if x.prefilter:
                            prefilter += 1
                            continue
                        tested += 1
                        if options.output_patterns:
                            patterns.write(x.pattern)

                        if x.filter and not options.print_filtered:
                            continue
                        printed += 1
                        print(format_output(x,
                                            lineage_dict,
                                            options.lmm,
                                            options.print_samples))
        else:
            for data in lmm_iter:
                ret = fit_lmm(*data)

                for x in ret:
                    if x.prefilter:
                        prefilter += 1
                        continue
                    tested += 1
                    if options.output_patterns:
                        patterns.write(x.pattern)

                    if x.filter and not options.print_filtered:
                        continue
                    printed += 1
                    print(format_output(x,
                                        lineage_dict,
                                        options.lmm,
                                        options.print_samples))

    sys.stderr.write('%d loaded variants\n' % (prefilter + tested))
    sys.stderr.write('%d filtered variants\n' % prefilter)
    sys.stderr.write('%d tested variants\n' % tested)
    sys.stderr.write('%d printed variants\n' % printed)


if __name__ == "__main__":
    main()
