#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Script to map kmer hits to reference'''

import sys
import os
import tempfile
from math import log10

from .bwa import bwa_index
from .bwa import bwa_iter

def get_options():
    import argparse

    description = 'Create Phandango Manhattan plot from kmer results'
    parser = argparse.ArgumentParser(description=description, prog="phandango_plot")

    parser.add_argument("kmers",
                        help="Kmers file, filtered output from SEER")
    parser.add_argument("reference",
                        help="Reference fasta file")
    parser.add_argument("output",
                        help="Output file")

    parser.add_argument("--bwa",
                        help="Location of bwa executable"
                        "[default=bwa]",
                        default="bwa")
    parser.add_argument("--tmp-prefix",
                        help="Directory to store temporary files"
                        "default=cwd",
                        default=os.getcwd())
    return parser.parse_args()


def main():
    options = get_options()

    # Open seer results
    # seer_remaining = seer_results
    seer_results = open(options.kmers, 'r')
    header_vals = seer_remaining.readline().rstrip.split("\t")
    p_val_col = 0
    for column in header_vals:
        if column == "lrt-pvalue":
            break
        p_val_col += 1

    tmp_fa = tempfile.NamedTemporaryFile(prefix=options.tmp_prefix + "/")
    with open(tmp_fa.name, 'w') as kmer_fa:
        for kmer in seer_results:
            kmer_fa.write(">" + str(kmers_remaining) + "\n")
            kmer_fa.write(kmer.split("\t")[0] + "\n")

    seer_results.seek(0)
    seer_results.readline()

    # index reference sequence
    bwa_index(options.reference)

    # run bwa mem -k 8
    with open(options.output, 'w') as outfile:
        outfile.write("\t".join(["SNP", "BP", "minLOG10(P)", "log10(p)", "r^2"]))

        mapped_kmers = bwa_iter(options.reference, remaining_fa_tmp, "mem")
        for mapping, kmer_line in zip(mapped_kmers, seer_remaining):
            p_val = float(kmer_line.split("\t")[p_val_col])
            if mapping.mapped and p_val > 0:
                log10p = -log10(p_val)
                for (contig, start, end, strand) in mapping.positions:
                    outfile.write("\t".join(["26", ".", str(start) + ".." + str(end), str(log10p), "0"]))

    # Clean up
    tmp_fa.close()

if __name__ == "__main__":
    main()
