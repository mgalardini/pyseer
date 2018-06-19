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
import hashlib
import binascii

import pyseer.classes as var_obj
from .cmdscale import cmdscale


def load_phenotypes(infile, column):
    """Load phenotypes vector

    Args:
        infile (str)
            Matrix input file
        column (str or None)
            Phenotype column name or None to pick the last column

    Returns:
        p (pandas.Series)
            Phenotype vector (n, 1)
    """
    p = pd.read_table(infile, index_col=0)
    p.index = p.index.astype(str)
    if column is None:
        p = p[p.columns[-1]]
    else:
        p = p[column]
    # Remove missing values
    p = p.dropna()
    return p


def load_structure(infile, p, max_dimensions, mds_type="classic", n_cpus=1,
                   seed=None):
    """Load population structure and apply multidimensional scaling

    Args:
        infile (str)
            Population structure (distance matrix) input file
        p (pandas.Series)
            Phenotype vector (n, 1)
        max_dimensions (int)
            Maximum dimensions to consider when applying
            `metric` or `non-metric` MDS
        mds_type (str)
            MDS algorithm to apply. One of `classic`,
            `metric` or `non-metric`. Any other input will trigger
            the `metric` MDS
        n_cpus (int)
            Number of CPUs to be used for the `metric` or `non-metric`MDS
        seed (int or None)
            Random seed for `metric` or `non-metric` MDS, None if not required

    Returns:
        m (pandas.DataFrame)
            Population structure after MDS (n, m)
    """
    m = pd.read_table(infile,
                      index_col=0)
    m.index = m.index.astype(str)
    sys.stderr.write("Structure matrix has dimension " + str(m.shape) + "\n")

    # Also take intersection here, so that MDS isn't done using samples not present
    # in sample
    intersecting_samples = p.index.intersection(m.index)
    m = m.loc[intersecting_samples, intersecting_samples]

    if len(intersecting_samples) == 0:
        sys.stderr.write('None of the phenotyped samples were found in ' +
                         'population structure matrix\n')
        sys.exit(1)

    # MDS
    if mds_type == "classic":
        projection, evals = cmdscale(m)
    else:
        metric_mds = True
        if mds_type == "non-metric":
            metric_mds = False
        elif mds_type != "metric":
            sys.stderr.write("Unsupported mds type chosen. Assuming metric\n")

        mds = manifold.MDS(max_dimensions, metric_mds, n_jobs=n_cpus,
                           dissimilarity='precomputed',
                           random_state=seed)
        projection = mds.fit_transform(m.values)

    m = pd.DataFrame(projection,
                     index=m.index)
    for i in range(m.shape[1]):
        m[i] = m[i] / max(abs(m[i]))
    return m


def load_lineage(infile, p):
    """Load custom lineage clusters definitions

    Args:
        infile (str)
            Input file for lineage clusters
        p (pandas.Series)
            Phenotypes vector (n, 1)

    Returns:
        result (tuple of (numpy.array, list))
            Lineage binary matrix and cluster labels
    """
    lin = pd.Series([x.rstrip().split()[1]
                     for x in open(infile)],
                    index=[x.split()[0]
                           for x in open(infile)])

    if (len(p.index.difference(lin.index)) > 0):
        sys.stderr.write("All samples with a phenotype must be present in lineage file\n")
        sys.exit(1)
    else:
        lin = lin.loc[lin.index.intersection(p.index)]

    lineages = sorted(set(lin.values))

    lineage_design_mat = []
    lineage_assign = []
    for categ in lineages:
        lineage_design_mat.append(pd.Series([1 if x == categ
                                             else 0
                                             for x in lin.values],
                                            index=lin.index))
        lineage_assign.append(categ)
    lineage_design_mat = pd.concat(lineage_design_mat, axis=1)

    return(lineage_design_mat.values, lineage_assign)


def load_covariates(infile, covariates, p):
    """Load and encode a covariates matrix

    Args:
        infile (str)
            Input file for the covariates matrix
        covariates (iterable or None)
            List of string indicating which columns to use and their
            interpretation. Example: `2q` indicates that the second column
            from the file is a quantitative variable, `2` indicates that
            that same column is categorical. If None, the matrix is loaded
            but nothing is done with it.
        p (pandas.Series)
            Phenotypes vector (n, 1)

    Returns:
        cov (pandas.DataFrame)
            Covariance matrix (n, m)
    """
    c = pd.read_table(infile,
                      header=None,
                      index_col=0)
    c.index = c.index.astype(str)
    c.columns = ['covariate%d' % (x+2) for x in range(c.shape[1])]

    if (len(p.index.difference(c.index)) > 0):
        sys.stderr.write("All samples with a phenotype must be present in covariate file\n")
        sys.exit(1)
    else:
        c = c.loc[c.index.intersection(p.index)]

    # which covariates to use?
    if covariates is None:
        cov = pd.DataFrame([])
    else:
        cov = []
        for col in covariates:
            cnum = int(col.rstrip('q'))
            if cnum == 1 or cnum > c.shape[1] + 1:
                sys.stderr.write('Covariates columns values should be '
                                 '> 1 and less than or equal to total number of ' +
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
        if len(cov) > 0:
            cov = pd.concat(cov, axis=1)
        else:
            cov = pd.DataFrame([])
    return cov


def load_burden(infile, burden_regions):
    """Load burden regions for VCF analysis

    Args:
        infile (str)
            Input file for burden regions
        burden_regions (list)
            List to be filled in-place
    """
    with open(infile, "r") as region_file:
        for region in region_file:
            (name, region) = region.rstrip().split()
            burden_regions.append((name, region))


def read_variant(infile, p, var_type, burden, burden_regions,
                 uncompressed, all_strains, sample_order):
    """Read input line and parse depending on input file type

    Return a variant name and pres/abs vector

    Args:
        infile (opened file)
            Handle to opened variant file
        p (pandas.Series)
            Phenotypes vector (n, 1)
        var_type (str)
            Variants type (one of: kmers, vcf or Rtab)
        burden (bool)
            Whether to slice a vcf file by burden regions
        burden_regions (collections.deque)
            Burden regions to slice the vcf with
        uncompressed (bool)
            Whether the kmers file is uncompressed
        all_strains (set-like)
            All sample labels that should be present
        sample_order
            Sampes order to interpret each Rtab line

    Returns:
        eof (bool)
            Whether we are at the end of the file
        k (numpy.array)
            Variant presence/absence vector
        var_name (str)
            Variant name
        kstrains (list)
            Samples in which the variant is present
        nkstrains (list)
            Samples in which the variant is absent
        af (float)
            Allele frequency
    """
    if var_type not in {'kmers', 'vcf', 'Rtab'}:
        raise ValueError('Variants type not supported')
    if var_type is "vcf":
        # burden tests read through regions and slice vcf
        if burden:
            if len(burden_regions) > 0:
                line_in = burden_regions.popleft()
            else:  # Last; to raise exception on next loop
                line_in = None
        # read single vcf line
        else:
            try:
                line_in = next(infile)
            except StopIteration:
                line_in = None
    else:
        # kmers and Rtab plain text files
        line_in = infile.readline()

    if not line_in:
        eof = True
        return(eof, None, None, None, None, None)
    else:
        eof = False
        d = {}
        if var_type == "kmers":
            if not uncompressed:
                line_in = line_in.decode()
            var_name, strains = (line_in.split()[0],
                                 line_in.rstrip().split(
                                 '|')[1].lstrip().split())

            d = {x.split(':')[0]: 1
                 for x in strains}

        elif var_type == "vcf":
            if not burden:
                var_name = read_vcf_var(line_in, d)
                if var_name is None:
                    return (eof, None, None, None, None, None)
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
                        var_sub_name = read_vcf_var(variant, d)
                else:  # stop trying to make 'fetch' happen
                    sys.stderr.write("Could not parse region %s\n" %
                                     str(region))
                    return (eof, None, None, None, None, None)

        elif var_type == "Rtab":
            split_line = line_in.rstrip().split()
            var_name, strains = split_line[0], split_line[1:]
            # sanity check
            if len(strains) != len(sample_order):
                raise ValueError('Unexpected mismatch between header and data row')
            for present, sample in zip(strains, sample_order):
                # sanity check
                if present not in {'0', '1'}:
                    raise ValueError('Rtab file not binary')
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

        k = np.array([d[x] for x in p.index
                      if x in d])

    return (eof, k, var_name, kstrains, nkstrains, af)


def read_vcf_var(variant, d):
    """Parses vcf variants from pysam

    Returns None if filtered variant. Mutates passed dictionary d

    Args:
        variant (pysam.libcbcf.VariantRecord)
            Variant to be parsed
        d (dict)
            Dictionary to be populated in-place
    """
    var_name = "_".join([variant.contig, str(variant.pos)] +
                        [str(allele) for allele in variant.alleles])

    # Do not support multiple alleles. Use 'bcftools norm' to split these
    if variant.alts != None and len(variant.alts) > 1:
        sys.stderr.write("Multiple alleles at %s_%s. Skipping\n" %
                         (variant.contig, str(variant.pos)))
        var_name = None
    elif len(variant.filter.keys()) > 0 and "PASS" not in variant.filter.keys():
        var_name = None
    else:
        for sample, call in variant.samples.items():
            # This is dominant encoding. Any instance of '1' will count as present
            # Could change to additive, summing instances, or reccessive only counting
            # when all instances are 1.
            # Shouldn't matter for bacteria, but some people call hets
            for haplotype in call['GT']:
                if str(haplotype) is not "." and haplotype is not None and haplotype != 0:
                    d[sample] = 1
                    break

    return(var_name)


def iter_variants(p, m, cov, var_type, burden, burden_regions, infile,
                  all_strains, sample_order, lineage_effects, lineage_clusters,
                  min_af, max_af, filter_pvalue, lrt_pvalue, null_fit,
                  firth_null, uncompressed, continuous):
    """Make an iterable to pass single variants to fixed effects regression

    Args:
        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        m (numpy.array)
            Population structure matrix (n, k)
        cov (pandas.DataFrame)
            Covariates matrix (n, m)
        var_type (str)
            Variants type (one of: kmers, vcf or Rtab)
        burden (bool)
            Whether to slice a vcf file by burden regions
        burden_regions (collections.deque)
            Burden regions to slice the vcf with
        infile (opened file)
            Handle to opened variant file
        all_strains (set-like)
            All sample labels that should be present
        sample_order
            Sampes order to interpret each Rtab line
        lineage_effects (bool)
            Whether to fit lineage effects
        lineage clusters (list)
            Lineage clusters indexes
        min_af (float)
            Minimum allele frequency (inclusive)
        max_af (bool)
            maximum allele frequency (inclusive)
        filter_pvalue (float)
            Pre-filtering p-value threshold
        lrt_pvalue (float)
            Filtering p-value threshold
        null_fit (float or statsmodels.regression.linear_model.RegressionResultsWrapper)
            Null-fit likelihood (binary) or model (continuous)
        firth_null (float)
            Firth regression likelihood
        uncompressed (bool)
            Whether the kmers file is uncompressed
        continuous (bool)
            Whether the phenotype is continuous or not

    Returns:
        var_name (str)
            Variant name
        v (numpy.array)
            Phenotypes vector (n, 1)
        k (numpy.array)
            Variant presence/absence vector (n, 1)
        m (numpy.array)
            Population structure matrix (n, k)
        c (numpy.array)
            Covariates matrix (n, m)
        af (float)
            Allele frequency
        pattern (bytes)
            Variant hash
        lineage_effects (bool)
            Whether to fit lineage effects
        lineage clusters (list)
            Lineage clusters indexes
        filter_pvalue (float)
            Pre-filtering p-value threshold
        lrt_pvalue (float)
            Filtering p-value threshold
        null_fit (float or statsmodels.regression.linear_model.RegressionResultsWrapper)
            Null-fit likelihood (binary) or model (continuous)
        firth_null (float)
            Firth regression likelihood
        kstrains (iterable)
            Sample labels with the variant
        nkstrains (iterable)
            Sample labels without the variant
        continuous (bool)
            Whether the phenotype is continuous or not
    """
    while True:
        eof, k, var_name, kstrains, nkstrains, af = read_variant(infile,
                                                                 p,
                                                                 var_type,
                                                                 burden,
                                                                 burden_regions,
                                                                 uncompressed,
                                                                 all_strains,
                                                                 sample_order)

        # check for EOF
        if eof:
            raise StopIteration

        if (k is None) or not (min_af <= af <= max_af):
            yield (None, None, None, None, None, None,
                   None, None, None, None, None, None,
                   None, None, None, None)
        else:
            v = p.values
            c = cov.values
            pattern = hash_pattern(k)

            yield (var_name, v, k, m, c, af, pattern,
                   lineage_effects, lineage_clusters,
                   filter_pvalue, lrt_pvalue, null_fit, firth_null,
                   kstrains, nkstrains, continuous)


def iter_variants_lmm(variant_iter, lmm, h2,
                      lineage, lineage_clusters,
                      covariates, continuous, filter_pvalue, lrt_pvalue):
    """Make an iterable to pass single variants to fixed effects regression"""
    for variants, variant_mat, eof in variant_iter:
        if len(variants) == 0:
            break
        yield (lmm, h2, variants, variant_mat, lineage,
               lineage_clusters, covariates,
               continuous, filter_pvalue, lrt_pvalue)
        if eof:
            break


# Loads a block of variants into memory for use with LMM
def load_var_block(var_type, p, burden, burden_regions, infile,
                   all_strains, sample_order, min_af, max_af,
                   uncompressed, block_size):
    """Make in iterable to load blocks of variants for LMM

    Args:
        var_type (str)
            Variants type (one of: kmers, vcf or Rtab)
        p (pandas.DataFrame)
            Phenotype vector (n, 1)
        burden (bool)
            Whether to slice a vcf file by burden regions
        burden_regions (collections.deque)
            Burden regions to slice the vcf with
        infile (opened file)
            Handle to opened variant file
        all_strains (set-like)
            All sample labels that should be present
        sample_order
            Sampes order to interpret each Rtab line
        min_af (float)
            Minimum allele frequency (inclusive)
        max_af (bool)
            maximum allele frequency (inclusive)
        uncompressed (bool)
            Whether the kmers file is uncompressed
        block_size (int)
            How many variants to be loaded at once

    Returns:
        variants (iterable)
            A collection of pyseer.classes.LMM objects describing the
            loaded variants (n,)
        variant_mat (numpy.array)
            Variant bloack presence/absence matrix (n, block_size)
        eof (bool)
            Whether we are at the end of the file
    """
    while True:
        variants = []
        # pre-allocation of memory
        variant_mat = np.zeros((len(p), block_size))
        for var_idx in range(block_size):
            eof, k, var_name, kstrains, nkstrains, af = read_variant(
                                            infile, p, var_type,
                                            burden, burden_regions,
                                            uncompressed, all_strains,
                                            sample_order)

            # check for EOF
            if eof:
                break

            if k is None or af < min_af or af > max_af:
                pattern = None
            else:
                pattern = hash_pattern(k)
                variant_mat[:, var_idx] = k
            variants.append((var_obj.LMM(var_name, pattern, af, np.nan,
                                         np.nan, np.nan, np.nan, np.nan,
                                         np.nan, kstrains, nkstrains,
                                         set(), True, True),
                             p.values, k))
        yield (variants, variant_mat, eof)
        if eof:
            break

    yield None, None, True


def hash_pattern(k):
    """Calculates the hash of a presence/absence vector

    Args:
        k (numpy.array)
            Variant presence/absence binary vector (n, 1)

    Returns:
        hash (byte)
            Hashed pattern
    """
    pattern = k.view(np.uint8)
    hashed = hashlib.md5(pattern)
    return (binascii.b2a_base64(hashed.digest()))
