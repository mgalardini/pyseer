# Copyright 2017 Marco Galardini and John Lees

'''Square a Mash distance matrix'''

import os
import sys
import pandas as pd

from .__init__ import __version__


def get_options():
    import argparse

    description = 'Make a square matrix out of a mash stream'
    parser = argparse.ArgumentParser(description=description,
                                     prog='square_mash')

    parser.add_argument('--classic',
                        action='store_true',
                        default=False,
                        help='Output table in a format suitable for R_mds.pl')

    parser.add_argument('--version', action='version',
                        version='%(prog)s '+__version__)

    return parser.parse_args()


def main():
    options = get_options()

    d = {}
    for l in sys.stdin:
        g1, g2, dist = l.split()[:3]
        g1 = os.path.split(g1)[-1].split('.')[0]
        g2 = os.path.split(g2)[-1].split('.')[0]
        dist = float(dist)
        d[g1] = d.get(g1, {})
        d[g1][g2] = dist
        d[g2] = d.get(g2, {})
        d[g2][g1] = dist

    m = pd.DataFrame(d)
    if not options.classic:
        m.to_csv(sys.stdout,
                 sep='\t')
    else:
        m.to_csv(sys.stdout,
                 index=False,
                 header=False,
                 sep=',')


if __name__ == "__main__":
    main()
