# Copyright 2017 Marco Galardini and John Lees

'''Utilities'''

import os
import contextlib
from decimal import Decimal

# thanks to Laurent LAPORTE on SO
@contextlib.contextmanager
def set_env(**environ):
    """
    Temporarily set the process environment variables.

    >>> with set_env(PLUGINS_DIR=u'test/plugins'):
    ...   "PLUGINS_DIR" in os.environ
    True

    >>> "PLUGINS_DIR" in os.environ
    False
    """
    old_environ = dict(os.environ)
    os.environ.update(environ)
    try:
        yield
    finally:
        os.environ.clear()
        os.environ.update(old_environ)


# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np


def format_output(item, print_samples=False):
    out = '%s' % item.kmer
    out += '\t' + '\t'.join(['%.2E' % Decimal(x)
                             if np.isfinite(x)
                             else ''
                             for x in (item.af,
                                       item.prep,
                                       item.lrt_pvalue,
                                       item.kbeta,
                                       item.bse,
                                       item.intercept)])
    out += '\t' + '\t'.join(['%.2E' % x
                             if np.isfinite(x)
                             else ''
                             for x in item.betas])
    if print_samples:
        out += '\t' + '\t'.join((','.join(item.kstrains),
                                 ','.join(item.nkstrains)))
    out += '\t%s' % ','.join(item.notes)
    return out

