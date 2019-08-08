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
    parser = argparse.ArgumentParser(description=description, prog="phandango_mapper")

    parser.add_argument("kmers",
                        help="Kmers file, filtered output from SEER")
    parser.add_argument("reference",
                        help="Reference fasta file")
    parser.add_argument("output",
                        help="Output file")

    parser.add_argument("--bwa",
                        help="Location of bwa executable "
                        "[Default: bwa]",
                        default="bwa")
    parser.add_argument("--tmp-prefix",
                        help="Directory to store temporary files "
                        "[Default: cwd]",
                        default=os.getcwd())
    parser.add_argument("--use-filter-p",
                        help="Plot the unadjusted p-value "
                        "[Default: lrt-p-value]",
                        action='store_true',
                        default=False)

    return parser.parse_args()


def main():
    options = get_options()

    # Open seer results
    # seer_remaining = seer_results
    seer_results = open(options.kmers, 'r')
    header_vals = seer_results.readline().rstrip().split("\t")
    lrt_idx = None
    lin_idx = None
    for idx, column in enumerate(header_vals):
        if column == "lrt-pvalue" and not options.use_filter_p:
            lrt_idx = idx
        elif column == "filter-pvalue" and options.use_filter_p:
            lrt_idx = idx
        elif column == "lineage":
            lin_idx = idx

    if not lrt_idx:
        sys.stderr.write("Could not find 'lrt-pvalue' field in header\n")
        sys.exit(1)

    tmp_fa = tempfile.NamedTemporaryFile(prefix=options.tmp_prefix + "/")
    kmer_idx = 0
    with open(tmp_fa.name, 'w') as kmer_fa:
        for kmer in seer_results:
            kmer_idx += 1
            kmer_fa.write(">" + str(kmer_idx) + "\n")
            kmer_fa.write(kmer.split("\t")[0] + "\n")

    seer_results.seek(0)
    seer_results.readline()

    # index reference sequence
    bwa_index(options.reference)

    # run bwa mem -k 8
    mapped = 0
    total = 0
    with open(options.output, 'w') as outfile, open(options.output + ".lineage", 'w') as linfile:
        outfile.write("\t".join(["SNP", "BP", "minLOG10(P)", "log10(p)", "r^2"]) + "\n")

        contigs = set()
        mapped_kmers = bwa_iter(options.reference, tmp_fa.name, "mem")
        for mapping, kmer_line in zip(mapped_kmers, seer_results):
            total += 1
            p_val = float(kmer_line.split("\t")[lrt_idx])
            if mapping.mapped and p_val > 0:
                mapped += 1
                log10p = -log10(p_val)
                for (contig, start, end, strand) in mapping.positions:
                    if (contig not in contigs):
                        contigs.add(contig)
                        if (len(contigs) > 2):
                            sys.stderr.write("WARNING: New chromosome/contig observed: " + contig + "\n")
                            sys.stderr.write("WARNING: Points will be overplotted - add previous contig length to BP\n")

                    outfile.write("\t".join([contig, ".", str(start) + ".." + str(end), str(log10p), "0"]) + "\n")
                    if lin_idx:
                        linfile.write(kmer_line.split("\t")[lin_idx] + "\n")

    # Clean up
    sys.stderr.write("Read " + str(total) + " k-mers\n")
    sys.stderr.write("Mapped " + str(mapped) + " k-mers\n")
    tmp_fa.close()
    if not lin_idx:
        os.remove(options.output + ".lineage")

if __name__ == "__main__":
    main()
