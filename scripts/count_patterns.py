#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Count unique patterns'''

mem_adjust = 10

def get_options():
    import argparse

    description = 'Extract a distance matrix from a phylogeny'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('patterns',
                        help='File of patterns from pyseer')

    parser.add_argument('--cores',
                        default=1,
                        help='Number of cores to use')
    parser.add_argument('--memory',
                        default=1000,
                        help='Maximum memory to use (in Mb)')
    parser.add_argument('--temp',
                        default='/tmp',
                        help='Directory to write tmp files to')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import sys
    import subprocess

    command = ("LC_ALL=C sort -u " +
              "--parallel=" + str(options.cores) +
              " -S " + str(options.memory - mem_adjust) + "M" +
              " -T " + options.temp +
              " " + options.patterns +
              " | wc -l")
    p = subprocess.check_output(command, shell=True, encoding='utf-8')
    print(p.rstrip())
