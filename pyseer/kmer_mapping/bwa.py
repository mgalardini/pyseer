#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Functions to run and interface with bwa'''

import subprocess
import os
import sys

from collections import namedtuple

BWA = namedtuple('BWA', ['mapped', 'positions'])

# Creates a bwa index, if it does not already exist
def bwa_index(fasta_file):
    suffix_list = [".amb", ".ann", ".bwt", ".pac", ".sa"]

    create_index = 0
    for index_file in [fasta_file + suf for suf in suffix_list]:
        if not os.path.isfile(index_file):
            create_index = 1
            break

    if create_index:
        command = "bwa index " + fasta_file
        try:
            subprocess.run(command, shell=True, check=True, stderr=subprocess.DEVNULL)
        except AttributeError:
            # python prior to 3.5
            devnull = open(os.devnull, 'w')
            subprocess.check_call(command, shell=True, stderr=devnull)

# Runs bwa, iterates over results and parses them
def bwa_iter(reference, fasta, algorithm):

    if algorithm == "mem":
        command = "bwa mem -k 8 " + reference + " " + fasta
    elif algorithm == "fastmap":
        command = "bwa fastmap -l 9 " + reference + " " + fasta
    else:
        sys.stderr.write("Unknown algorithm type for bwa\n")
        raise ValueError(algorithm)

    try:
        devnull = stderr=subprocess.DEVNULL
    except AttributeError:
        # python 2.x
        devnull = open(os.devnull, 'w')
    bwa_p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=devnull, shell=True, universal_newlines=True)
    # read sam file from bwa mem
    if algorithm == "mem":
        # discard header
        bwa_p.stdout.readline()
        bwa_p.stdout.readline()

        for sam_line in bwa_p.stdout:
            sam_fields = sam_line.rstrip().split("\t")
            positions = []
            if int(sam_fields[1]) & 4 == 4:
                mapped = False
            else:
                mapped = True

                # primary mapping
                if int(sam_fields[1]) & 16 == 16:
                    strand = "-"
                else:
                    strand = "+"
                positions.append((sam_fields[2], sam_fields[3], int(sam_fields[3]) + len(sam_fields[9]) - 1, strand))

                # secondary mappings (as good as primary - same CIGAR string)
                if len(sam_fields) > 15:
                    secondary = sam_fields[15].split(":")
                    if secondary[0] == "XA" and secondary[1] == "Z":
                        for secondary_mapping in secondary[2].split(";"):
                            if secondary_mapping != '':
                                (contig, pos, cigar, edit_distance) = secondary_mapping.split(",")
                                if cigar == sam_fields[5]:
                                    strand = pos[0]
                                    positions.append((contig, pos[1:], int(pos[1:]) + len(sam_fields[9]) - 1, strand))

            yield(BWA(mapped, positions))
    # read bwa fastmap output
    else:
        mapped = False
        positions = []

        first_line = bwa_p.stdout.readline().rstrip().split("\t")
        (sq, idx, length) = first_line
        while True:
            fastmap_line = bwa_p.stdout.readline()
            fastmap_line = fastmap_line.rstrip()
            if fastmap_line == "//":
                next_line = bwa_p.stdout.readline().rstrip().split("\t")
                fastmap_hit = BWA(mapped, positions)
                if len(next_line) < 3:  # EOF reached
                    yield(fastmap_hit)
                    raise StopIteration
                else:
                    (sq, idx, length) = next_line
                    mapped = False
                    positions = []
                    yield(fastmap_hit)
            else:
                hits = []
                fastmap_fields = fastmap_line.split("\t")
                # in case a line is missing a few fields
                if len(fastmap_fields) < 5:
                    continue
                #
                if fastmap_fields[1] == '0' and fastmap_fields[2] == length: #  full hits only
                    mapped = True
                    for hit in fastmap_fields[4:]:
                        (contig, pos) = hit.split(":")
                        strand = pos[0]
                        positions.append((contig, int(pos[1:]), int(pos[1:]) + int(length) - 1, strand))
