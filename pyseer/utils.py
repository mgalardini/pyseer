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


def format_output(item, lineage_dict=None, model='seer', print_samples=False):
    """Format results for a variant for stdout printing

    Args:
        item (pyseer.classes.Seer or pyseer.classes.LMM)
            Variant results container
        lineage_dict (list)
            Lineage labels
        model (str)
            The model used
        print_samples (bool)
            Whether to add the samples list to the putput

    Returns:
        out (str)
            Tab-delimited string to be printed
    """
    out = '%s' % item.kmer

    if model == "enet" or model == "rf":
        out += '\t' + '\t'.join(['%.2E' % Decimal(x)
                                 if np.isfinite(x)
                                 else ''
                                 for x in (item.af,
                                           item.prep,
                                           item.pvalue,
                                           item.kbeta)])
    else:
        out += '\t' + '\t'.join(['%.2E' % Decimal(x)
                                 if np.isfinite(x)
                                 else ''
                                 for x in (item.af,
                                           item.prep,
                                           item.pvalue,
                                           item.kbeta,
                                           item.bse)])
        if model == 'lmm':
            if np.isfinite(item.frac_h2):
                frac_h2 = '%.2E' % Decimal(item.frac_h2)
            else:
                frac_h2 = ''
            out += '\t' + frac_h2
        else:
            if np.isfinite(item.intercept):
                intercept = '%.2E' % Decimal(item.intercept)
            else:
                intercept = ''
            out += '\t' + intercept + '\t'
            # No distances may not have further betas
            if not np.all(np.equal(item.betas, None)):
                out += '\t'.join(['%.2E' % Decimal(x)
                                 if np.isfinite(x)
                                 else ''
                                 for x in item.betas])

    if lineage_dict is not None:
        if item.max_lineage is not None and np.isfinite(item.max_lineage):
            out += '\t' + lineage_dict[item.max_lineage]
        else:
            out += '\tNA'
    if print_samples:
        out += '\t' + '\t'.join((','.join(item.kstrains),
                                 ','.join(item.nkstrains)))
    out += '\t%s' % ','.join(item.notes)

    return out
