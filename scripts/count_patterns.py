#!/usr/bin/env python
# Copyright 2017 Marco Galardini and John Lees

'''Count unique patterns'''

mem_adjust = 10


def get_options():
    import argparse

    description = 'Calculate p-value threshold using Bonferroni correction'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('patterns',
                        help='File of patterns from pyseer')

    parser.add_argument('--alpha',
                        default=0.05,
                        type=float,
                        help='Family-wise error rate')
    parser.add_argument('--cores',
                        default=1,
                        help='Number of cores to use')
    parser.add_argument('--memory',
                        default=1024,
                        help='Maximum memory to use (in Mb)')
    parser.add_argument('--temp',
                        default='/tmp',
                        help='Directory to write tmp files to')

    return parser.parse_args()


if __name__ == "__main__":
    options = get_options()

    import sys
    import subprocess
    from decimal import Decimal

    command = "LC_ALL=C sort -u "
    if options.cores > 1:
        command +=  "--parallel=" + str(options.cores)
    command += (" -S " + str(int(options.memory) - mem_adjust) + "M" +
               " -T " + options.temp +
               " " + options.patterns +
               " | wc -l")

    p = subprocess.check_output(command, shell=True, universal_newlines=True)
    print("Patterns:\t" + p.rstrip())
    print("Threshold:\t" + '%.2E' % Decimal(options.alpha/float(p.rstrip())))
