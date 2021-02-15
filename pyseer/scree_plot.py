# Copyright 2017 Marco Galardini and John Lees

'''Draw a scree plot'''

from .cmdscale import cmdscale


def get_options():
    import argparse

    description = 'Draw a scree-plot from MDS eigenvalues'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('distances',
                        help='Strains distance square matrix')

    parser.add_argument('--max-dimensions',
                        type=int,
                        default=30,
                        help='Maximum dimensions to plot [Default: 30]')
    parser.add_argument('--output',
                        default='scree_plot.png',
                        help='Plot filename [Default: scree_plot.png]')

    return parser.parse_args()


def main():
    options = get_options()

    import sys
    import numpy as np
    import pandas as pd

    import matplotlib
    try:
        matplotlib.use("TkAgg")
        import matplotlib.pyplot as plt
    except ImportError as e:
        matplotlib.use("AGG")
        import matplotlib.pyplot as plt


    m = pd.read_csv(options.distances,
                    index_col=0,
                    sep='\t')
    # metric MDS scaling
    projection, evals = cmdscale(m)

    if evals.shape[0] > options.max_dimensions:
        sys.stderr.write('Plotting only the first %d eigenvalues out of %d\n' %
                         (options.max_dimensions,
                          evals.shape[0]))
        evals = evals[:options.max_dimensions]

    plt.figure(figsize=(0.25*evals.shape[0], 3))

    plt.plot(range(evals.shape[0]),
             evals,
             'ko-')
    plt.ylabel('eigenvalue')
    plt.xlabel('PCs')
    plt.xticks(range(evals.shape[0]),
               range(1, evals.shape[0]+1),
               rotation=90)

    plt.xlim(-0.25, evals.shape[0]-0.75)
    if options.max_dimensions >= 5:
        plt.tight_layout()

    plt.savefig(options.output, dpi=150)


if __name__ == "__main__":
    main()
