from collections import namedtuple

LMM = namedtuple('LMM', ['kmer',
                           'af', 'prep', 'pvalue',
                           'kbeta', 'bse', 'frac_h2',
                           'max_lineage',
                           'kstrains', 'nkstrains'])

Seer = namedtuple('Seer', ['kmer',
                           'af', 'prep', 'pvalue',
                           'kbeta', 'bse',
                           'intercept', 'betas',
                           'max_lineage',
                           'kstrains', 'nkstrains',
                           'notes',
                           'prefilter', 'filter'])
