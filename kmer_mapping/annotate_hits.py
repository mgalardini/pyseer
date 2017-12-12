#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Script to annotate kmer hits'''

import sys
import os
import re
import tempfile
import subprocess
import pybedtools

from .bwa import bwa_index
from .bwa import bwa_iter

def get_options():
    import argparse

    description = 'Iteratively annotate significant kmers from SEER'
    parser = argparse.ArgumentParser(description=description, prog="annotate_hits")

    parser.add_argument("kmers",
                        help="Kmers file, filtered output from SEER")
    parser.add_argument("references",
                        help="File of reference annotations. "
                        "First column fasta sequence, second column gff annotation, "
                        "third column 'ref' or 'draft'")
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


# Return overlapping and closest annotations in GFF file
def create_annotation(gff_bed, contig, start, end, strand):
    start = str(start)
    end = str(end)

    query_interval = pybedtools.BedTool(' '.join([contig, start, end, 'kmer', '0', strand]), from_string=True)
    in_gene = extract_gene(query_interval.intersect(b=gff_bed, s=False, stream=True, wb=True))
    up_gene = extract_gene(query_interval.closest(b=gff_bed, s=False, D="ref", iu=True, stream=True))
    down_gene = extract_gene(query_interval.closest(b=gff_bed, s=False, D="ref", id=True, stream=True))
    pybedtools.cleanup() # delete the bed file

    # return contig:pos;gene_down;gene_in;gene_up
    annotation_string = ";".join([contig + ":" + start + "-" + end, down_gene, in_gene, up_gene])
    return annotation_string


# returns first overlapping feature with gene= annotation. Otherwise first feature ID
def extract_gene(bedtools_intervals):
    ID = None
    gene = None

    for match in bedtools_intervals.features():
        for tag in match.fields[15].split(";"):
            parse_tag = re.search('^(.+)="(.+)"$', tag)
            if parse_tag:
                if parse_tag.group(1) == "gene":
                    gene = parse_tag.group(2)
                    break
                elif parse_tag.group(1) == "ID" and ID is None:
                    ID = parse_tag.group(2)
        if gene is not None:
            break

    if gene is None:
        if ID is not None:
            gene = ID
        else:
            gene = ""

    return gene

def main():
    options = get_options()

    # tmp file locations
    remaining_tmp = options.tmp_prefix + "/remaining_kmers.txt"
    remaining_next_tmp = options.tmp_prefix + "/remaining_kmers_next.txt"
    remaining_fa_tmp = options.tmp_prefix + "/remaining_kmers.fa"
    remaining_fa_next_tmp = options.tmp_prefix + "/remaining_kmers_next.fa"
    pybedtools.helpers.set_tempdir(options.tmp_prefix)

    # read references and drafts into list
    references = []
    with open(options.references, 'r') as reference_files:
        for reference in reference_files:
            (fa, gff, ref) = reference.rstrip().split()
            references.append((fa, gff, ref))

    output_file = open(options.output, 'w')

    # Open seer results
    # seer_remaining = seer_results
    seer_remaining = open(options.kmers, 'r')
    header = seer_remaining.readline()

    # Write out kmer fasta file, keep track of count
    kmers_remaining = 0
    with open(remaining_fa_tmp, 'w') as kmer_fa:
        for kmer in seer_remaining:
            kmers_remaining += 1
            kmer_fa.write(">" + str(kmers_remaining) + "\n")
            kmer_fa.write(kmer.split("\t")[0] + "\n")

    seer_remaining.seek(0)
    seer_remaining.readline()

    # for each reference, then draft
    ref_id = 0
    for reference in references:
        (ref_fa, ref_gff, ref_type) = reference
        ref_id += 1

        # print number of kmers remaining. if zero, break
        if ref_type == "ref":
            sys.stderr.write("Reference " + str(ref_id) + "\n")
        else:
            sys.stderr.write("Draft reference " + str(ref_id) + "\n")
        sys.stderr.write(str(kmers_remaining) + " kmers remain\n")
        if kmers_remaining == 0:
            break

        # index reference sequence
        bwa_index(ref_fa)
        if ref_type == "ref":
            bwa_algorithm = "mem"
        elif ref_type == "draft":
            bwa_algorithm = "fastmap"
        else:
            bwa_algorithm = "fastmap"
            sys.stderr.write("Unknown reference type " + ref_type + " for " + ref_fa + ". Assuming draft\n")

        # Fix ref annotation
        tmp_bed = tempfile.NamedTemporaryFile(prefix=options.tmp_prefix + "/")
        subprocess.run("gff2bed < " + ref_gff + " > " + tmp_bed.name, shell=True, check=True)
        ref_annotation = pybedtools.BedTool(tmp_bed.name)
        filtered_ref = ref_annotation.filter(lambda x: True if x[7] == "CDS" else False).saveas('tmp_bed')

        next_seer_remaining = open(remaining_next_tmp, 'w')
        next_fasta_remaining = open(remaining_fa_next_tmp, 'w')

        # run bwa mem -k 8 for ref, bwa fastmap for draft of remaining.fa
        new_idx = 0
        mapped_kmers = bwa_iter(ref_fa, remaining_fa_tmp, bwa_algorithm)
        for mapping, kmer_line in zip(mapped_kmers, seer_remaining):
            if mapping.mapped:
                kmers_remaining -= 1
                annotations = []
                for (contig, start, end, strand) in mapping.positions:
                    annotations.append(create_annotation(filtered_ref, contig, start, end, strand))
                output_file.write("\t".join([kmer_line.rstrip(), ",".join(annotations)]) + "\n")
            else:
                # if unmapped write to seer_remaining and remaining.fa
                next_seer_remaining.write(kmer_line)

                new_idx += 1
                next_fasta_remaining.write(">" + str(new_idx) + "\n")
                next_fasta_remaining.write(kmer_line.split("\t")[0] + "\n")

        # Clean up
        seer_remaining.close()
        next_seer_remaining.close()
        next_fasta_remaining.close()
        tmp_bed.close()
        os.rename(remaining_next_tmp, remaining_tmp)
        os.rename(remaining_fa_next_tmp, remaining_fa_tmp)

        # Open next kmer file
        seer_remaining = open(remaining_tmp, 'r')

if __name__ == "__main__":
    main()
