# Copyright 2014 Microsoft Corporation

# Notice of changes made to this file as part of their integration
# into pyseer
'''FaST-LMM 1D fit for h2. Modified to python3 syntax and removed logging'''
#

import scipy as SP
import scipy.optimize as opt


def minimize1D(f, evalgrid = None, nGrid=10, minval=0.0, maxval = 0.99999, verbose=False, brent=True,check_boundaries = True, resultgrid=None, return_grid=False):
    '''
    minimize a function f(x) in the grid between minval and maxval.
    The function will be evaluated on a grid and then all triplets,
    where the inner value is smaller than the two outer values are optimized by
    Brent's algorithm.
    --------------------------------------------------------------------------
    Input:
    f(x)    : callable target function
    evalgrid: 1-D array prespecified grid of x-values
    nGrid   : number of x-grid points to evaluate f(x)
    minval  : minimum x-value for optimization of f(x)
    maxval  : maximum x-value for optimization of f(x)
    brent   : boolean indicator whether to do Brent search or not.
              (default: True)
    --------------------------------------------------------------------------
    Output list:
    [xopt, f(xopt)]
    xopt    : x-value at the optimum
    f(xopt) : function value at the optimum
    --------------------------------------------------------------------------
    '''
    #evaluate the target function on a grid:
    if verbose: print("evaluating target function on a grid")
    if evalgrid is not None and brent:# if brent we need to sort the input values
        i_sort = evalgrid.argsort()
        evalgrid = evalgrid[i_sort]
    if resultgrid is None:
        [evalgrid,resultgrid] = evalgrid1D(f, evalgrid = evalgrid, nGrid=nGrid, minval=minval, maxval = maxval  )

    i_currentmin=resultgrid.argmin()
    minglobal = (evalgrid[i_currentmin],resultgrid[i_currentmin])
    if brent:#do Brent search in addition to rest?
        if check_boundaries:
            if verbose: print("checking grid point boundaries to see if further search is required")
            if resultgrid[0]<resultgrid[1]:#if the outer boundary point is a local optimum expand search bounded between the grid points
                if verbose: print("resultgrid[0]<resultgrid[1]--> outer boundary point is a local optimum expand search bounded between the grid points")
                minlocal = opt.fminbound(f,evalgrid[0],evalgrid[1],full_output=True)
                if minlocal[1]<minglobal[1]:
                    if verbose: print("found a new minimum during grid search")
                    minglobal=minlocal[0:2]
            if resultgrid[-1]<resultgrid[-2]:#if the outer boundary point is a local optimum expand search bounded between the grid points
                if verbose: print("resultgrid[-1]<resultgrid[-2]-->outer boundary point is a local optimum expand search bounded between the grid points")
                minlocal = opt.fminbound(f,evalgrid[-2],evalgrid[-1],full_output=True)
                if minlocal[1]<minglobal[1]:
                    if verbose: print("found a new minimum during grid search")
                    minglobal=minlocal[0:2]
        if verbose: print("exploring triplets with brent search")
        onebrent=False
        for i in range(resultgrid.shape[0]-2):#if any triplet is found, where the inner point is a local optimum expand search
            if (resultgrid[i+1]<resultgrid[i+2]) and (resultgrid[i+1]<resultgrid[i]):
                onebrent=True
                if verbose: print("found triplet to explore")
                minlocal = opt.brent(f,brack = (evalgrid[i],evalgrid[i+1],evalgrid[i+2]),full_output=True)
                if minlocal[1]<minglobal[1]:
                    minglobal=minlocal[0:2]
                    if verbose: print("found new minimum from brent search")
    if return_grid:
        return (minglobal[0], minglobal[1], evalgrid, resultgrid)
    else:
        return minglobal


def evalgrid1D(f, evalgrid = None, nGrid=10, minval=0.0, maxval = 0.99999, dimF=0):
    '''
    evaluate a function f(x) on all values of a grid.
    --------------------------------------------------------------------------
    Input:
    f(x)    : callable target function
    evalgrid: 1-D array prespecified grid of x-values
    nGrid   : number of x-grid points to evaluate f(x)
    minval  : minimum x-value for optimization of f(x)
    maxval  : maximum x-value for optimization of f(x)
    --------------------------------------------------------------------------
    Output:
    evalgrid    : x-values
    resultgrid  : f(x)-values
    --------------------------------------------------------------------------
    '''
    if evalgrid is None:
        step = (maxval-minval)/(nGrid)
        evalgrid = SP.arange(minval,maxval+step,step)
    if dimF:
        resultgrid = SP.ones((evalgrid.shape[0],dimF))*9999999999999.0
    else:
        resultgrid = SP.ones(evalgrid.shape[0])*9999999999999.0
    for i in range(evalgrid.shape[0]):
        fevalgrid = f(evalgrid[i])
        assert SP.isreal(fevalgrid).all(),"function returned imaginary value"
        resultgrid[i] = fevalgrid
    return (evalgrid,resultgrid)
