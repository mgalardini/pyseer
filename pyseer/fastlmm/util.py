import numpy as np
import warnings
import sys
import time

def thin_results_file(myfile,dup_postfix="v2"):
    '''
    Used in score vs lrt to remove any lines in the results
    ending with "v2", as these were replicate gene set entries.
    '''
    sets = np.loadtxt(myfile,dtype=str,comments=None)
    nodup_ind = []
    dup_ind = []

    #indexes of non-duplicates, as indicated by dup_postfix
    for i in range(0,sets.shape[0]):
        tmpset=sets[i,0]
        if tmpset[-2:]!=dup_postfix:
            nodup_ind.append(i)
        else:
            dup_ind.append(i)

    sets_nodup = sets[nodup_ind]
    print("%i reps, and %i non-reps" % (len(dup_ind),len(nodup_ind)))
    return sets_nodup

def compare_files(file1,file2,tol=1e-8,delimiter="\t"):
    '''
    Given two files, compare the contents, including numbers up to absolute tolerance, tol
    Returns: val,msg
    where val is True/False (true means files to compare to each other) and a msg for the failure.
    '''
    dat1=np.loadtxt(file1,dtype='str',delimiter=delimiter,comments=None)
    dat2=np.loadtxt(file2,dtype='str',delimiter=delimiter,comments=None)

    ncol1=dat1[0].size
    ncol2=dat2[0].size

    if ncol1!=ncol2:
        return False,"num columns do not match up"

    try:
        head1=dat1[0,:]
        head2=dat2[0,:]
    except:
        #file contains just a single column.
        return np.all(dat1==dat2), "single column result doesn't match exactly ('{0}')".format(file1)
    #logging.warn("DO headers match up? (file='{0}', '{1}' =?= '{2}')".format(file1, head1,head2))
    if not np.all(head1==head2):
        return False, "headers do not match up (file='{0}', '{1}' =?= '{2}')".format(file1, head1,head2)

    for c in range(ncol1):
        checked=False
        col1=dat1[1:,c]
        col2=dat2[1:,c]
        try:
            #if it is numeric
            col1=np.array(col1,dtype='float64')
            col2=np.array(col2,dtype='float64')
        except Exception:
            # if it is a string
            pass
            if not np.all(col1==col2):
                return False, "string column %s does not match" % head1[c]
            checked=True

        #if it is numeric
        if not checked:
            absdiff=np.absolute(col1-col2)
            if np.any(absdiff>tol):
                try:
                    return False, "numeric column %s does diff of %e not match within tolerance %e" % (head1[c],max(absdiff),  tol)
                except:
                    return False, "Error trying to print error message while comparing '{0}' and '{1}'".format(file1,file2)

    return True, "files are comparable within abs tolerance=%e" % tol

def compare_mixed_files(file1,file2,tol=1e-8,delimiter="\t"):
    '''
    Given two files, compare the contents, including numbers up to absolute tolerance, tol
    Returns: val,msg
    where val is True/False (true means files to compare to each other) and a msg for the failure.
    '''
    dat1=np.loadtxt(file1,dtype='str',delimiter=delimiter,comments=None)
    dat2=np.loadtxt(file2,dtype='str',delimiter=delimiter,comments=None)

    ncol1=dat1[0].size
    ncol2=dat2[0].size

    if ncol1!=ncol2:
        return False,"num columns do not match up"

    try:
        r_count = dat1.shape[0]
        c_count = dat1.shape[1]
    except:
        #file contains just a single column.
        return np.all(dat1==dat2), "single column result doesn't match exactly ('{0}')".format(file1)

    for r in range(r_count):
        for c in range(c_count):
            val1 = dat1[r,c]
            val2 = dat2[r,c]
            if val1!=val2:
                try:
                    f1 = float(val1)
                    f2 = float(val2)
                except:
                    return False, "Values do not match up (file='{0}', '{1}' =?= '{2}')".format(file1, val1, val2)
                if abs(f1-f2) > tol:
                    return False, "Values too different (file='{0}', '{1}' =?= '{2}')".format(file1, val1, val2)
    return True, "files are comparable within abs tolerance=%e" % tol

#could make this more efficient by reading in blocks of SNPs, as in
#FastLmmSet.py:KfromAltSnps()
def write_kernel(iid,K,fileout):
    '''
    writes out kernel
    assumes that iid contains a list of the ids, or else a list of [famid, personid] which
    then get merged with a space in between
    '''
    nInd = K.shape[0]
    header = 'var'
    iid_merged = []

    # first line contains iids
    for i in range(nInd):
        if iid.ndim==1 or iid.shape[1]==1:
            header += '\t%s'%(iid[i])
            iid_merged.append('%s'%(iid[i]))
        else:
            header += '\t%s %s'%(iid[i,0],iid[i,1])
            iid_merged.append('%s %s'%(iid[i,0],iid[i,1]))

    # each row of the matrix is one line
    f = open(fileout,'w')
    f.write(header+'\n')
    for i in range(nInd):
        row = ['\t%.4f'%x for x in K[i,:]]
        f.write('%s%s\n'%(iid_merged[i],''.join(row)))
    f.close()

def write_plink_covariates(iid,X,fileout):
    '''
    writes out plink-style covariates/phen file
    assuming that X is [N,M] for N individuals and M features
    assumes that iid contains a list of [famid, personid]
    '''
    [nInd,M] = X.shape

    # each row of the matrix is one line
    f = open(fileout,'w')
    for i in range(nInd):
        row = ['\t%.4f'%x for x in X[i,:]]
        f.write('%s\t%s%s\n'%(iid[i,0],iid[i,1],''.join(row)))
    f.close()


def combineseeds(seed1,seed2):
    import hashlib
    import sys
    seed=int(hashlib.md5(str(seed1) + "_" + str(seed2)).hexdigest()[-8:], 16)    #as of numpy 1.9, seeds must be 32-bit, so keep only the 8 right-most hex digits
    return seed


def standardize_col(dat,meanonly=False):
    '''
    Mean impute each columns of an array.
    '''
    colmean=np.nanmean(dat,axis=0)
    if ~meanonly:
        colstd=np.nanstd(dat,axis=0)
    else:
        colstd=None
    ncol=dat.shape[1]
    nmissing=np.zeros((ncol))
    datimp=np.empty_like(dat); datimp[:]=dat
    for c in np.arange(0,ncol):
        datimp[np.isnan(datimp[:,c]),c]=colmean[c]
        datimp[:,c]=datimp[:,c]-colmean[c]
        if not meanonly:
            if colstd[c]>1e-6:
                datimp[:,c]=datimp[:,c]/colstd[c]
            else:
                print("warning: colstd=" + colstd[c] + " during normalization")
        nmissing[c]=float(np.isnan(dat[:,c]).sum())
    fracmissing=nmissing/dat.shape[0]
    return datimp,fracmissing

def extractcols(filein,colnameset=None,dtypeset=None):
    if colnameset is None: raise Exception("must specify column names to read")
    import pandas as pd
    data=pd.read_csv(filein,delimiter = '\t',dtype=dtypeset,usecols=colnameset)
    r={}
    for j in np.arange(0,len(colnameset)):
        name=colnameset.pop()
        r[name]=(data[name].values)
    return r


def argintersect_left(a, b):
    """
    find indices in a, whose corresponding values are in b
    ----------------------------------------------------------------------
    Input:
    a        : array, for which indices are returned that are in the intersect with b
    b        : array to be intersected with a
    ----------------------------------------------------------------------
    Output:
    the indices of elements of a, which are in intersect of a and b
    ----------------------------------------------------------------------
    """
    return np.arange(a.shape[0])[np.in1d(a,b)]


def intersect_ids(idslist,sep="Q_Q"):
    '''
    Takes a list of 2d string arrays of family and individual ids.
    These are intersected.
    "sep" is used to concatenate the family and individual ids into one unique string
    Returns: indarr, an array of size N x L, where N is the number of
             individuals in the intersection, and L is the number of lists in idslist, and which
             contains the index to use (in order) such that all people will be identical and in order
             across all data sets.
    If one of the lists=None, it is ignored (but still has values reported in indarr, all equal to -1),
    but the first list must not be None.
    '''
    #!!warnings.warn("This intersect_ids is deprecated. Pysnptools includes newer versions of intersect_ids", DeprecationWarning)
    id2ind={}
    L=len(idslist)
    observed=np.zeros(L,dtype='bool')

    for l, id_list in enumerate(idslist):
            if id_list is not None:
                observed[l]=1
                if l==0:
                    if ~observed[l]:
                        raise Exception("first list must be non-empty")
                    else:
                        for i in range(id_list.shape[0]):
                            id=id_list[i,0] +sep+ id_list[i,1]
                            entry=np.zeros(L)*np.nan #id_list to contain the index for this id, for all lists provided
                            entry[l]=i                 #index for the first one
                            id2ind[id]=entry
                elif observed[l]:
                    for i in range(id_list.shape[0]):
                        id=id_list[i,0] +sep+ id_list[i,1]
                        if id in id2ind:
                            id2ind[id][l]=i

    indarr=np.array(id2ind.values(),dtype='float')  #need float because may contain NaNs
    indarr[:,~observed]=-1                          #replace all Nan's from empty lists to -1
    inan = np.isnan(indarr).any(1)                  #find any rows that contain at least one Nan
    indarr=indarr[~inan]                            #keep only rows that are not NaN
    indarr=np.array(indarr,dtype='int')             #convert to int so can slice
    return indarr

def indof_constfeatures(X,axis=0):
    '''
    Assumes features are columns (by default, but can do rows), and checks to see if all features are simply constants,
    such that it is equivalent to a bias and nothing else
    '''
    featvar=np.var(X,axis=axis)
    badind = np.nonzero(featvar==0)[0]
    return badind

def constfeatures(X,axis=0):
    '''
    Assumes features are columns (by default, but can do rows), and checks to see if all features are simply constants,
    such that it is equivalent to a bias and nothing else
    '''
    featmeans=np.mean(X,axis=axis)
    return (X-featmeans==0).all()


def appendtofilename(filename,midfix,sep="."):
        import os
        dir, fileext = os.path.split(filename)
        file, extension = os.path.splitext(fileext)
        infofilename = dir + os.path.sep + file + sep + midfix + extension
        return infofilename

def datestamp(appendrandom=False):
    import datetime
    now = datetime.datetime.now()
    s = str(now)[:19].replace(" ","_").replace(":","_")
    if appendrandom:
        import random
        s += "_" + str(random.random())[2:]
    return s



#not needed, just use the sp RandomState.permutation
#def permute(numbersamples):
#    perm = np.random.permutation(numbersamples)
#    return perm

#Not needed because enumerate is built in to the language
#def appendindex(iter):
#    index = -1;
#    for item in iter:
#        index += 1
#        yield item, index

def create_directory_if_necessary(name, isfile=True, robust=False):
    import os
    if isfile:
        directory_name = os.path.dirname(name)
    else:
        directory_name = name

    if directory_name != "":
        if not robust:
            try:
                os.makedirs(directory_name)
            except OSError:
                if not os.path.isdir(directory_name):
                    raise Exception("not valid path: '{0}'. (Working directory is '{1}'".format(directory_name,os.getcwd()))
        else:
            is_ok = False
            time_to_sleep = 10.0
            for i in range(25):
                try:
                    os.makedirs(directory_name)
                    is_ok = True
                    break
                except OSError:
                    if not os.path.isdir(directory_name):
                        time_to_sleep *= 1.1
                        warnings.warn("creating directory robust=True, try#{0},time={3} error: not valid path: '{1}'. (Working directory is '{2}'".format(i, directory_name,os.getcwd(),int(time_to_sleep)))
                        time.sleep(int(time_to_sleep))  #make random?
                    else:
                        is_ok = True
                        break
            if not is_ok:
                raise Exception("not valid path: '{0}'. (Working directory is '{1}'".format(directory_name,os.getcwd()))


def which(vec):
    '''
    find the True from the index 0 with bool vector vec
    ----------------------------------------------------------------------
    Input:
    vec        : vector of bool
    ----------------------------------------------------------------------
    Output:
    index of the first True from the bool vector vec
    ----------------------------------------------------------------------
    '''
    for i, item in enumerate(vec):
        if (item):
            return(i)
    return(-1)

def which_opposite(vec):
    '''
    find the True from the index 0 with bool vector vec
    ----------------------------------------------------------------------
    Input:
    vec        : vector of bool
    ----------------------------------------------------------------------
    Output:
    index of the last True from the bool vector vec
    ----------------------------------------------------------------------
    '''
    for i in reversed(range(len(vec))):
        item = vec[i]
        if (item):
            return(i)
    return(-1)


def generatePermutation(numbersamples,randomSeedOrState):
    from numpy.random import RandomState

    if isinstance(randomSeedOrState,RandomState):
        randomstate = randomSeedOrState
    else:
        randomstate = RandomState(int(randomSeedOrState % sys.maxint))

    perm = randomstate.permutation(numbersamples)
    return perm

def excludeinds(pos0, pos1, mindist = 10.0,idist = 2):
    '''
    get the indices of SNPs that have to be excluded from the set of null SNPs when testing alternative SNPs to correct for proximal contamination.
    --------------------------------------------------------------------------
    Input:
    pos0        : [S0*3] array of null-model SNP positions
    pos1        : [S0*3] array of alternative-model SNP positions
    idist       : index in pos array that the exclusion is based on.
                  (1=genetic distance, 2=basepair distance)
    --------------------------------------------------------------------------
    Output:
    i_exclude   : [S] 1-D boolean array indicating excluson of SNPs
                  (True: exclude, False: do not exclude)
    --------------------------------------------------------------------------
    '''
    chromosomes1 = np.unique(pos1[:,0])
    i_exclude = np.zeros(pos0[:,0].shape[0],dtype = 'bool')
    if (mindist>=0.0):
        for ichr in range(chromosomes1.shape[0]):
            i_SNPs1_chr=pos1[:,0] == chromosomes1[ichr]
            i_SNPs0_chr=pos0[:,0] == chromosomes1[ichr]
            pos1_ = pos1[i_SNPs1_chr,idist]
            pos0_ = pos0[i_SNPs0_chr,idist]
            distmatrix = pos1_[np.newaxis,:] - pos0_[:,np.newaxis]
            i_exclude[i_SNPs0_chr] = (np.absolute(distmatrix)<=mindist).any(1)
    return i_exclude


def dotDotRange(dotDotString):
    '''
    A method for generating integers.
    For example:

> for i in util.dotDotRange("1..4,100,-1..1"): print i
1
2
3
4
100
-1
0
1

    '''
    for intervalString in dotDotString.split(","):
        parts = intervalString.split("..")
        if len(parts) > 2 : raise Exception("Expect at most one '..' between commas. (see {0})".format(intervalString))
        start = int(parts[0])
        if len(parts) == 1:
            yield start
        else:
            lastInclusive = int(parts[1])
            for i in range(start,lastInclusive+1):
                yield i


def _run_length_encode(seq):
    count = 0
    previous = None
    for item in seq:
        if count == 0:
            count = 1
            previous = item
        elif item == previous:
            count += 1
        else:
            yield previous, count
            previous = item
            count =1
    if count > 0:
        yield previous, count

def _rel_to_midpoint(rle):
    previous_count = 0
    for item, count in rle:
        yield previous_count + count // 2
        previous_count += count

def _color_list(chr_list,rle):
    chr_to_index = dict((chr,index) for index,(chr,count) in enumerate(rle))
    index_to_color = {0:"b",1:"g"}
    result = [index_to_color[chr_to_index[chr]%len(index_to_color)] for chr in chr_list]
    return result

def manhattan_plot(chr_pos_pvalue_array,pvalue_line=None,plot_threshold=1.0,vline_significant=False,marker="o", chromosome_starts=None, xaxis_unit_bp=True, alpha=0.5):
    """
    Function to create a Manhattan plot.  See http://en.wikipedia.org/wiki/Manhattan_plot.

    Args:
        chr_pos_pvalue_array:   an n x 3 numpy array. The three columns are the chrom number
                                (as a number), the position, and pvalue.
                                :type chr_pos_pvalue_array: numpy array
        pvalue_line:            (Default: None). If given, draws a line at that PValue.
                                :type pvalue_line: a 'pheno dictionary' or a string
        plot_threshold:         plot only SNPs that achieve a P-value smaller than pvalue_threshold
                                to speed up plotting
        vline_significant:      boolean. Draw a vertical line at each significant Pvalue?
                                :rtype: none, but changes the global current figure.
        marker:                 marker for the scatter plot. default: "o"
        chromosome_starts:      [Nchrom x 3] ndarray: chromosome, cumulative start position, cumulative stop position
                                cumulative chromosome starts, for plotting. If None (default), this is estimated from data
        xaxis_unit_bp:          plot cumulative position in basepair units on x axis? If False, only
                                use rank of SNP positions. (default: True)
        alpha:                  alpha (opaquness) for P-value markers in scatterplot (default 0.5)

    Returns:
        chromosome_starts       [Nchrom x 3] ndarray: chromosome, cumulative start position, cumulative stop position
                                cumulative chromosome starts used in plotting.

    :Example:

    >>> from fastlmm.association import single_snp
    >>> from pysnptools.snpreader import Bed
    >>> import matplotlib.pyplot as plt
    >>> import fastlmm.util.util as flutil
    >>> pheno_fn = "../feature_selection/examples/toydata.phe"
    >>> results_dataframe = single_snp(test_snps="../feature_selection/examples/toydata.5chrom", pheno=pheno_fn, h2=.2)
    >>> #chromosome_starts = flutil.manhattan_plot(results_dataframe.as_matrix(["Chr", "ChrPos", "PValue"]),pvalue_line=1e-7)
    >>> #plt.show()

    """
    import matplotlib
    matplotlib.use('Agg') #This lets it work even on machines without graphics displays
    import matplotlib.pyplot as plt

    # create a copy of the data and sort it by chrom and then position
    array = np.array(chr_pos_pvalue_array)
    if plot_threshold:
        array = array[array[:,2]<=plot_threshold]
    else:
        plot_threshold = 1.0
    array=array[np.argsort(array[:,1]),:] #sort by ChrPos
    array=array[np.argsort(array[:,0],kind='mergesort'),:] #Finally, sort by Chr (but keep ChrPos in case of ties)
    rle = list(_run_length_encode(array[:,0]))

    if xaxis_unit_bp:   #compute and use cumulative basepair positions for x-axis
        if chromosome_starts is None:
            chromosome_starts = _compute_x_positions_chrom(array)
        chr_pos_list = _compute_x_positions_snps(array, chromosome_starts)
        plt.xlim([0,chromosome_starts[-1,2]+1])
        plt.xticks(chromosome_starts[:,1:3].mean(1),chromosome_starts[:,0])
    else:               #use rank indices for x-axis
        chr_pos_list = np.arange(array.shape[0])
        xTickMarks = [str(int(item)) for item,count in rle]
        plt.xlim([0,array.shape[0]])
        plt.xticks(list(_rel_to_midpoint(rle)), xTickMarks)
    y = -np.log10(array[:,2])
    max_y = y.max()

    if pvalue_line and vline_significant:   #mark significant associations (ones that pass the pvalue_line) by a red vertical line:
        idx_significant = array[:,2]<pvalue_line
        if np.any(idx_significant):
            y_significant = y[idx_significant]
            chr_pos_list_significant = chr_pos_list[idx_significant]
            for i in range(len(chr_pos_list_significant)):
                plt.axvline(x=chr_pos_list_significant[i],ymin = 0.0, ymax = y_significant[i], color = 'r',alpha=0.8)

    plt.scatter(chr_pos_list,y,marker=marker,c=_color_list(array[:,0],rle),edgecolor='none',s=y/max_y*20+0.5, alpha=alpha)
    plt.xlabel("chromosome")
    plt.ylabel("-log10(P value)")

    if pvalue_line:
        plt.axhline(-np.log10(pvalue_line),linestyle="--",color='gray')
    plt.ylim([-np.log10(plot_threshold),None])
    return chromosome_starts

def _compute_x_positions_chrom(positions, offset=1e5):
    chromosomes = np.unique(positions[:,0])
    chromosomes.sort()
    chromosome_starts = np.zeros((chromosomes.shape[0],3),dtype="object")
    chr_start_next = 0
    for i, chromosome in enumerate(chromosomes):
        pos_chr = positions[positions[:,0]==chromosome]
        chromosome_starts[i,0] = chromosome                     #the chromosome
        chromosome_starts[i,1] = chr_start_next                 #start of the chromosome
        chromosome_starts[i,2] = chr_start_next + pos_chr.max() #end of the chromosome
        chr_start_next = chromosome_starts[i,2] + offset
    return chromosome_starts

def _compute_x_positions_snps(positions, chromosome_starts):
    cumulative_pos = np.zeros(positions.shape[0])
    for i, chromosome_start in enumerate(chromosome_starts):
        idx_chr = positions[:,0]==chromosome_start[0]
        cumulative_pos[idx_chr] = positions[idx_chr][:,1] + chromosome_start[1]
    return cumulative_pos

if __name__ == "__main__":

    #logging.basicConfig(level=logging.INFO)

    import doctest
    doctest.testmod()
