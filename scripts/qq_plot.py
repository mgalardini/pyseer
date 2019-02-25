# Copyright 2017 Marco Galardini and John Lees

'''Draw a qq-plot'''


def get_options():
    import argparse

    description = 'Draw a QQ-plot from pyseer lrt-pvalue results'
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('table',
                        help='Pyseer output')

    parser.add_argument('--output',
                        default='qq_plot.png',
                        help='Plot filename [Default: qq_plot.png]')

    return parser.parse_args()


def main():
    options = get_options()

    import sys
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    import matplotlib.pyplot as plt

    m = pd.read_csv(options.table,
                    usecols=['lrt-pvalue'],
                    sep='\t')['lrt-pvalue']

    plt.figure(figsize=(4, 3.75))
    ax = plt.subplot(111)

    y = -np.log10(m)
    x = -np.log10(np.random.uniform(0, 1, m.shape[0]))

    fig = sm.qqplot_2samples(y,
                             x,
                             xlabel='Expected $-log_{10}(pvalue)$',
                             ylabel='Observed $-log_{10}(pvalue)$',
                             line='45',
                             ax=ax)

    ax = fig.axes[0]
    ax.lines[0].set_color('k')
    ax.lines[0].set_alpha(0.3)

    ax.set_xlim(-0.5, x.max()+0.5)
    ax.set_ylim(-0.5, y.max()+0.5)

    plt.tight_layout()

    plt.savefig(options.output, dpi=150)


if __name__ == "__main__":
    main()
