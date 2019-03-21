# Copyright 2017 Marco Galardini and John Lees

'''Calculate similariity matrix C from variants'''

import os
import sys
import gzip
import numpy as np
import pandas as pd
from pysam import VariantFile

from .__init__ import __version__

from .input import load_var_block

block_size = 1000


def get_options():
    import argparse

    description = 'Calculate a similarity matrix using variant presence/absence information'
    parser = argparse.ArgumentParser(description=description,
                                     prog='similarity')

    parser.add_argument("samples",
                        help="List of sample names to use")

    variant_group = parser.add_mutually_exclusive_group(required=True)
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

    parser.add_argument('--min-af',
                        type=float,
                        default=0.01,
                        help='Minimum AF [Default: 0.01]')
    parser.add_argument('--max-af',
                        type=float,
                        default=0.99,
                        help='Maximum AF [Default: 0.99]')
    parser.add_argument('--max-missing',
                           type=float,
                           default=0.05,
                           help='Maximum missing (vcf/Rtab) [Default: 0.05]')
    parser.add_argument('--uncompressed',
                        action='store_true',
                        default=False,
                        help='Uncompressed kmers file [Default: gzipped]')

    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    # Create dummy pheno object from sample list
    sample_list = []
    with open(options.samples, 'r') as sample_file:
        for sample in sample_file:
            sample_list.append(sample.rstrip())
    p = pd.Series(np.zeros(len(sample_list)),
                  index=sample_list)

    # Open variant file. Mostly copied from __main__
    sample_order = []
    all_strains = set(p.index)

    if options.kmers:
        var_type = "kmers"
        if options.uncompressed:
            infile = open(options.kmers)
        else:
            infile = gzip.open(options.kmers, 'r')
    elif options.vcf:
        var_type = "vcf"
        infile = VariantFile(options.vcf)
    else:
        # Rtab files have a header, rather than sample names accessible by row
        var_type = "Rtab"
        infile = open(options.pres)
        header = infile.readline().rstrip()
        sample_order = header.split()[1:]

    eof = 0
    # no copy of first variant_mat made. Reserve memory
    G = np.empty((len(p), block_size))
    sys.stderr.write("Reading in variants\n")
    v_iter = load_var_block(var_type, p, None, None, infile,
                            all_strains, sample_order,
                            options.min_af, options.max_af,
                            options.max_missing, options.uncompressed,
                            block_size)
    while not eof:
        variants, variant_mat, eof = next(v_iter)
        if G.shape[1] > block_size:
            G = np.concatenate(G, variant_mat)
        else:
            G = variant_mat

    sys.stderr.write("Calculating sample similarity\n")
    K = np.matmul(G, np.transpose(G))
    K_out = pd.DataFrame(K, index=p.index, columns=p.index)
    K_out.to_csv(sys.stdout,
                 sep='\t')


if __name__ == "__main__":
    main()
