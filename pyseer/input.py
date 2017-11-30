# Copyright 2017 Marco Galardini and John Lees

'''Functions to read data into pyseer and iterate over instances'''

import sys
from .utils import set_env
import re
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import pandas as pd
from sklearn import manifold
from .cmdscale import cmdscale

LMM = namedtuple('LMM', ['kmer',
                           'af', 'prep', 'pvalue',
                           'kbeta', 'bse', 'frac_h2',
                           'max_lineage',
                           'kstrains', 'nkstrains'])

def load_phenotypes(infile, column):
    p = pd.Series([float(x.rstrip().split()[column-1])
                   for x in open(infile)],
                  index=[x.split()[0]
                         for x in open(infile)])
    # Remove missing values
    p = p.dropna()
    return p


def load_structure(infile, p, max_dimensions, mds_type="metric", n_cpus=1):
    m = pd.read_table(infile,
                      index_col=0)
    m = m.loc[p.index, p.index]

    # MDS
    if mds_type is "classic":
        projection, evals = cmdscale(m)
    else:
        metric_mds = True
        if mds_type is "non-metric":
            metric_mds = False
        elif mds_type is not "metric":
            sys.stderr.write("Unsupported mds type chosen. Assuming metric\n")

        mds = manifold.MDS(max_dimensions, metric_mds, n_jobs=n_cpus,
                           dissimilarity='precomputed')
        projection = mds.fit_transform(m.values)

    m = pd.DataFrame(projection,
                     index=m.index)
    for i in range(m.shape[1]):
        m[i] = m[i] / max(abs(m[i]))
    return m


# Loads custom cluster/lineage definitions
def load_lineage(infile, p):
    lin = pd.Series([x.rstrip().split()[1]
                   for x in open(infile)],
                  index=[x.split()[0]
                         for x in open(infile)])
    lin = lin.loc[p.index]
    lineages = set(lin.values)
    lineages.pop()

    lineage_design_mat = []
    lineage_assign = []
    for categ in lineages:
        lineage_design_mat.append(pd.Series([1 if x == categ
                                          else 0
                                          for x in lin.values],
                                         index=lin.index))
        lineage_assign.append(categ)
    lineage_design_mat = pd.concat(lineage_design_mat, axis=1)

    return(lineage_design_mat, lineage_assign)


def load_covariates(infile, covariates, p):
    c = pd.read_table(infile,
                      header=None,
                      index_col=0)
    c.columns = ['covariate%d' % (x+2) for x in range(c.shape[1])]
    c = c.loc[p.index]
    # which covariates to use?
    if covariates is None:
        cov = pd.DataFrame([])
    else:
        cov = []
        for col in covariates:
            cnum = int(col.rstrip('q'))
            if cnum == 1 or cnum > c.shape[1] + 1:
                sys.stderr.write('Covariates columns values should be '
                                 '> 1 and lower than total number of ' +
                                 'columns (%d)\n' % (c.shape[1] + 1))
                return None
            if col[-1] == 'q':
                # quantitative
                cov.append(c['covariate%d' % cnum])
            else:
                # categorical, dummy-encode it
                categories = set(c['covariate%d' % cnum])
                categories.pop()
                for i, categ in enumerate(categories):
                    cov.append(pd.Series([1 if x == categ
                                          else 0
                                          for x in
                                          c['covariate%d' % cnum].values],
                                         index=c.index,
                                         name='covariate%d_%d' % (cnum, i)))
        cov = pd.concat(cov, axis=1)
    return cov


def load_burden(infile, burden_regions):
    with open(infile, "r") as region_file:
        for region in region_file:
            (name, region) = region.rstrip().split()
            burden_regions.append((name, region))

# Parse depending on input file type. Need to end with a variant name
# and pres/abs dictionary
def read_variant(infile, line_in, var_type, uncompressed, all_strains, sample_order)

    d = {}
    if var_type == "kmers":
        if not uncompressed:
            line_in = line_in.decode()
        var_name, strains = (line_in.split()[0],
                             line_in.rstrip().split('|')[1].lstrip().split())

        d = {x.split(':')[0]: 1
             for x in strains}

    elif var_type == "vcf":
        if not burden:
            var_name = read_vcf_var(line_in, d)
            if var_name is None:
                return(None, None, None, None)
            else:
            # burden test. Regions are named contig:start-end.
            # Start is non-inclusive, so start one before to include
            (var_name, region) = line_in
            region = re.match('^(.+):(\d+)-(\d+)$', region)
            if region:
                # Adds presence to d for every variant
                # observation in region
                for variant in infile.fetch(region.group(1),
                                            int(region.group(2)) - 1,
                                            int(region.group(3))):
                    region_var_name = read_vcf_var(variant, d)
            else:  # stop trying to make 'fetch' happen
                sys.stderr.write("Could not parse region %s\n" %
                                 str(region))
                return (None, None, None, None)

    elif var_type == "Rtab":
        split_line = line_in.rstrip().split()
        var_name, strains = split_line[0], split_line[1:]
        for present, sample in zip(strains, sample_order):
            if present is not '0':
                d[sample] = 1

    # Use common dictionary to format design matrix etc
    kstrains = sorted(set(d.keys()).intersection(all_strains))
    nkstrains = sorted(all_strains.difference(set(kstrains)))

    # default for missing samples is absent kmer
    # currently up to user to be careful about matching pheno and var files
    for x in nkstrains:
        d[x] = 0

    af = float(len(kstrains)) / len(all_strains)
    # filter by AF
    if af < min_af or af > max_af:
        # pass it to the actual tests to keep track
        return(None, None, None, None)

    k = np.array([d[x] for x in p.index
                      if x in d])

    return(k, var_name, kstrains, nkstrains, af)


# Iterable to pass to fixed effects regression
def iter_variants(p, m, cov, var_type, burden, burden_regions, infile,
                  all_strains, sample_order, lineage_effects, lineage_clusters,
                  min_af, max_af, filter_pvalue, lrt_pvalue, null_fit,
                  firth_null, uncompressed, continuous):
    while True:
        if var_type is "vcf":
            # burden tests read through regions and slice vcf
            if burden:
                if len(burden_regions) > 0:
                    l = burden_regions.popleft()
                else:  # Last; to raise exception on next loop
                    l = None
            # read single vcf line
            else:
                l = next(infile)
        else:
            # kmers and Rtab plain text files
            l = infile.readline()

        # check for EOF
        if not l:
            raise StopIteration

        k, var_name, kstrains, nkstrains, af = read_variant(infile, l, var_type, uncompressed, all_strains, sample_order)
        if k is None:
            yield(None, None, None, None, None, None,
                  None, None, None, None, None, None,
                  None, None, None)

        v = p.values
        c = cov.values

        yield (var_name, v, k, m, c, af,
               lineage_effects, lineage_clusters,
               filter_pvalue, lrt_pvalue, null_fit, firth_null,
               kstrains, nkstrains, continuous)

# Loads a block of variants into memory for use with LMM
def load_var_block(var_type, burden, burden_regions, infile,
                  all_strains, sample_order, min_af, max_af, filter_pvalue,
                  uncompressed, continuous, block_size):
    prefilter = 0
    variants = []
    for var_idx in range(block_size):
        if var_type is "vcf":
            # burden tests read through regions and slice vcf
            if burden:
                if len(burden_regions) > 0:
                    l = burden_regions.popleft()
                else:  # Last; to raise exception on next loop
                    l = None
            # read single vcf line
            else:
                l = next(infile)
        else:
            # kmers and Rtab plain text files
            l = infile.readline()

        # check for EOF
        if not l:
            break

        k, var_name, kstrains, nkstrains, af = read_variant(infile, l, var_type, uncompressed, all_strains, sample_order)
        if k:
            if af > min_af and af < max_af:
                prep, bad_chisq = pre_filtering(p, k, continuous)
                if prep < filter_pvalue:
                    variants.append(LMM(var_name, af, prep, 0, 0, 0, 0, None, kstrains, nkstrains))
                    if variant_mat:
                        np.concatenate(variant_mat, k, axis=0)
                    else:
                        variant_mat = k
            else:
                prefilter += 1

    counts['prefilter'] = prefilter
    counts['tested'] = len(variants)

    return(variants, variant_mat, counts)

# Parses vcf variants from pysam. Returns None if filtered variant.
# Mutates passed dictionary d
def read_vcf_var(variant, d):

    var_name = "_".join([variant.contig, str(variant.pos)] +
                        [str(allele) for allele in variant.alleles])

    # Do not support multiple alleles. Use 'bcftools norm' to split these
    if len(variant.alts) > 1:
        sys.stderr.write("Multiple alleles at %s_%s. Skipping\n" %
                         (variant.contig, str(variant.pos)))
        var_name = None
    elif "PASS" not in variant.filter.keys() and "." not in variant.filter.keys():
        var_name = None
    else:
        for sample, call in variant.samples.items():
            # This is dominant encoding. Any instance of '1' will count as present
            # Could change to additive, summing instances, or reccessive only counting
            # when all instances are 1.
            # Shouldn't matter for bacteria, but some people call hets
            for haplotype in call['GT']:
                if str(haplotype) is not "." and haplotype != 0:
                    d[sample] = 1
                    break

    return(var_name)
