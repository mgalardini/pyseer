# Copyright 2017 Marco Galardini and John Lees

'''Functions to read data into pyseer and iterate over instances'''

from .utils import set_env
# avoid numpy taking up more than one thread
with set_env(MKL_NUM_THREADS='1',
             NUMEXPR_NUM_THREADS='1',
             OMP_NUM_THREADS='1'):
    import numpy as np
import pandas as pd


def load_phenotypes(infile):
    p = pd.Series([float(x.rstrip().split()[-1])
                   for x in open(infile)],
                  index=[x.split()[0]
                         for x in open(infile)])
    return p


# thanks to Francis Song for this function
# source: http://www.nervouscomputer.com/hfs/cmdscale-in-python/
def cmdscale(D):
    """
    Classical multidimensional scaling (MDS)

    Parameters
    ----------
    D : (n, n) array
        Symmetric distance matrix.

    Returns
    -------
    Y : (n, p) array
        Configuration matrix. Each column represents a dimension. Only the
        p dimensions corresponding to positive eigenvalues of B are returned.
        Note that each dimension is only determined up to an overall sign,
        corresponding to a reflection.

    e : (n,) array
        Eigenvalues of B.
    """
    # Number of points
    n = len(D)

    # Centering matrix
    H = np.eye(n) - np.ones((n, n))/n

    # YY^T
    B = -H.dot(D**2).dot(H)/2

    # Diagonalize
    evals, evecs = np.linalg.eigh(B)

    # Sort by eigenvalue in descending order
    idx   = np.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]

    # Compute the coordinates using positive-eigenvalued components only
    w, = np.where(evals > 0)
    L  = np.diag(np.sqrt(evals[w]))
    V  = evecs[:,w]
    Y  = V.dot(L)

    return Y, evals[evals > 0]


def load_structure(infile, p):
    m = pd.read_table(infile,
                      index_col=0)
    m = m.loc[p.index, p.index]
    # metric MDS scaling
    projection, evals = cmdscale(m)
    m = pd.DataFrame(projection,
                     index=m.index)
    for i in range(m.shape[1]):
        m[i] = m[i] / max(abs(m[i]))
    return m, evals


def load_covariates(infile, covariates, p):
    c = pd.read_table(infile,
                      header=None,
                      index_col=0)
    c.columns = ['covariate%d' % (x+2) for x in range(c.shape[1])]
    c = c.loc[p.index]
    # which covariates to use?
    if covariates is None:
        cov = pd.DataFrame([])
    else:
        cov = []
        for col in covariates:
            cnum = int(col.rstrip('q'))
            if cnum == 1 or cnum > c.shape[1] + 1:
                sys.stderr.write('Covariates columns values should be > 1 and lower ' +
                                 'than total number of columns (%d)\n' % (c.shape[1] + 1))
                return None
            if col[-1] == 'q':
                # quantitative
                cov.append(c['covariate%d' % cnum])
            else:
                # categorical, dummy-encode it
                categories = set(c['covariate%d' % cnum])
                categories.pop()
                for i, categ in enumerate(categories):
                    cov.append(pd.Series([1 if x == categ
                                          else 0
                                          for x in c['covariate%d' % cnum].values],
                                         index=c.index,
                                         name='covariate%d_%d' % (cnum, i)))
        cov = pd.concat(cov, axis=1)
    return cov


def iter_kmers(p, m, cov, infile, all_strains,
               min_af, max_af,
               filter_pvalue, lrt_pvalue, null_fit, firth_null,
               uncompressed):
    for l in infile:
        if not uncompressed:
            l = l.decode()
        kmer, strains = l.split()[0], l.rstrip().split('|')[1].lstrip().split()

        d = {x.split(':')[0]: 1
             for x in strains}
        kstrains = sorted(set(d.keys()).intersection(all_strains))
        nkstrains = sorted(all_strains.difference(
                           all_strains.intersection(
                           {x.split(':')[0] for x in strains})))
        for x in nkstrains:
            d[x] = 0

        af = float(len(kstrains)) / len(all_strains)
        # filter by AF
        if af < min_af or af > max_af:
            # pass it to the actual tests to keep track
            yield (kmer, None, None, None, None, af,
                   None, None, None, None,
                   kstrains, nkstrains)
            continue

        v = p.values
        k = np.array([d[x] for x in p.index
                      if x in d])
        c = cov.values

        yield (kmer, v, k, m, c, af,
               filter_pvalue, lrt_pvalue, null_fit, firth_null,
               kstrains, nkstrains)

