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
                        help="Location of bwa executable "
                        "[default=bwa]",
                        default="bwa")
    parser.add_argument("--tmp-prefix",
                        help="Directory to store temporary files "
                        "[default=./]",
                        default=os.getcwd())
    return parser.parse_args()


# returns first overlapping feature with gene= annotation. Otherwise first feature ID
def extract_genes(bedtools_intervals):
    annotations = {}
    for match in bedtools_intervals.features():
        kmer_id, hit_id = match.fields[3].split("_")
        if annotations.get(int(kmer_id)) == None:
            annotations[int(kmer_id)] = {}

        ID = None
        gene = None
        for tag in match.fields[15].split(";"):
            parse_tag = re.search('^(.+)=(.+)$', tag)
            if parse_tag:
                if parse_tag.group(1) == "gene":
                    gene = parse_tag.group(2)
                    break
                elif parse_tag.group(1) == "ID" and ID is None:
                    ID = parse_tag.group(2)
        if gene is None:
            if ID is not None:
                gene = ID
            else:
                gene = ""

        if annotations[int(kmer_id)].get(int(hit_id)) == None or annotations[int(kmer_id)][int(hit_id)] == "":
            annotations[int(kmer_id)][int(hit_id)] = gene
        else:
            annotations[int(kmer_id)][int(hit_id)] += "|" + gene

    return annotations

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
        if kmers_remaining == 0:
            break
        sys.stderr.write(str(kmers_remaining) + " kmers remain\n")
        if ref_type == "ref":
            sys.stderr.write("Reference " + str(ref_id) + "\n")
        else:
            sys.stderr.write("Draft reference " + str(ref_id) + "\n")

        # index reference sequence
        bwa_index(ref_fa)
        if ref_type == "ref":
            bwa_algorithms = ["mem", "fastmap"]
        elif ref_type == "draft":
            bwa_algorithms = ["fastmap"]
        else:
            bwa_algorithms = ["fastmap"]
            sys.stderr.write("Unknown reference type " + ref_type + " for " + ref_fa + ". Assuming draft\n")

        # Fix ref annotation
        tmp_bed = tempfile.NamedTemporaryFile(prefix=options.tmp_prefix + "/")
        try:
            subprocess.run("gff2bed < " + ref_gff + " > " + tmp_bed.name, shell=True, check=True)
        except AttributeError:
            # python prior to 3.5
            subprocess.check_call("gff2bed < " + ref_gff + " > " + tmp_bed.name, shell=True)
        ref_annotation = pybedtools.BedTool(tmp_bed.name)
        filtered_ref = ref_annotation.filter(lambda x: True if x[7] == "CDS" else False).saveas('tmp_bed')
        ref_annotation = pybedtools.BedTool('tmp_bed')

        for bwa_algorithm in bwa_algorithms:
            next_seer_remaining = open(remaining_next_tmp, 'w')
            next_fasta_remaining = open(remaining_fa_next_tmp, 'w')

            # run bwa mem -k 8 for ref, bwa fastmap for draft of remaining.fa
            new_idx = 0
            kmer_lines = []
            map_pos = {}

            mapped_kmers = bwa_iter(ref_fa, remaining_fa_tmp, bwa_algorithm)
            with tempfile.NamedTemporaryFile('w', prefix=options.tmp_prefix + "/") as query_bed:
                kmer_idx = 0
                for mapping, kmer_line in zip(mapped_kmers, seer_remaining):
                    if mapping.mapped:
                        kmers_remaining -= 1
                        kmer_lines.append(kmer_line.rstrip())
                        map_pos[kmer_idx] = []
                        for hit_idx, (contig, start, end, strand) in enumerate(mapping.positions):
                            map_pos[kmer_idx].append(contig + ":" + str(start) + "-" + str(end))
                            query_bed.write('\t'.join([contig, str(start), str(end), str(kmer_idx) + "_" + str(hit_idx), '0', strand]) + "\n")
                        kmer_idx += 1
                    else:
                        # if unmapped write to seer_remaining and remaining.fa
                        next_seer_remaining.write(kmer_line)

                        new_idx += 1
                        next_fasta_remaining.write(">" + str(new_idx) + "\n")
                        next_fasta_remaining.write(kmer_line.split("\t")[0] + "\n")

                if kmer_idx > 0:
                    query_bed.flush()
                    query_interval = pybedtools.BedTool(query_bed.name)
                    sorted_query = query_interval.sort()

                    in_genes = extract_genes(query_interval.intersect(b=ref_annotation, s=False, stream=True, wb=True))
                    up_genes = extract_genes(sorted_query.closest(b=ref_annotation, s=False, D="ref", iu=True, stream=True))
                    down_genes = extract_genes(sorted_query.closest(b=ref_annotation, s=False, D="ref", id=True, stream=True))

                    for kmer_idx, kmer_line  in enumerate(kmer_lines):
                        annotations = []
                        for hit_idx, hit in enumerate(map_pos[kmer_idx]):
                            annotation = hit + ";"
                            if kmer_idx in down_genes and hit_idx in down_genes[kmer_idx]:
                                annotation += down_genes[kmer_idx][hit_idx]
                            annotation += ";"
                            if kmer_idx in in_genes and hit_idx in in_genes[kmer_idx]:
                                annotation += in_genes[kmer_idx][hit_idx]
                            annotation += ";"
                            if kmer_idx in up_genes and hit_idx in up_genes[kmer_idx]:
                                annotation += up_genes[kmer_idx][hit_idx]
                            annotations.append(annotation)

                        output_file.write("\t".join([kmer_line, ",".join(annotations)]) + "\n")
                else:
                    # something went wrong, write down remaining kmers
                    for kmer_line in seer_remaining:
                        # if unmapped write to seer_remaining and remaining.fa
                        next_seer_remaining.write(kmer_line)

                        new_idx += 1
                        next_fasta_remaining.write(">" + str(new_idx) + "\n")
                        next_fasta_remaining.write(kmer_line.split("\t")[0] + "\n")

                pybedtools.cleanup() # delete the bed file

            # Clean up
            seer_remaining.close()
            next_seer_remaining.close()
            next_fasta_remaining.close()
            os.rename(remaining_next_tmp, remaining_tmp)
            os.rename(remaining_fa_next_tmp, remaining_fa_tmp)

            # Open next kmer file
            seer_remaining = open(remaining_tmp, 'r')

        # Clean up
        tmp_bed.close()
        os.remove('tmp_bed')

    sys.stderr.write(str(kmers_remaining) + " kmers remain unannotated\n")


if __name__ == "__main__":
    main()
