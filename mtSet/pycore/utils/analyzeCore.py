import sys
sys.path.append('./../../..')
import os
import subprocess
import pdb
import sys
import csv
import numpy as NP
from optparse import OptionParser
import time
import mtSet.pycore.modules.multiTraitSetTest as MTST
from mtSet.pycore.utils.read_utils import readNullModelFile
from mtSet.pycore.utils.read_utils import readWindowsFile
from mtSet.pycore.utils.read_utils import readCovarianceMatrixFile
from mtSet.pycore.utils.read_utils import readCovariatesFile
from mtSet.pycore.utils.read_utils import readPhenoFile
from mtSet.pycore.external.limix import plink_reader
import scipy as SP
import warnings

def scan(bfile,Y,cov,null,wnds,minSnps,i0,i1,perm_i,resfile,F):

    if perm_i is not None:
        print 'Generating permutation (permutation %d)'%perm_i
        NP.random.seed(perm_i)
        perm = NP.random.permutation(Y.shape[0])

    mtSet = MTST.MultiTraitSetTest(Y,S_XX=cov['eval'],U_XX=cov['evec'],F=F)
    mtSet.setNull(null)
    bim = plink_reader.readBIM(bfile,usecols=(0,1,2,3))
    fam = plink_reader.readFAM(bfile,usecols=(0,1))
   
    print 'fitting model'
    wnd_file = csv.writer(open(resfile,'wb'),delimiter='\t')
    for wnd_i in range(i0,i1):
        print '.. window %d - (%d, %d-%d) - %d snps'%(wnd_i,int(wnds[wnd_i,1]),int(wnds[wnd_i,2]),int(wnds[wnd_i,3]),int(wnds[wnd_i,-1]))
        if int(wnds[wnd_i,-1])<minSnps:
            print 'SKIPPED: number of snps lower than minSnps'
            continue
        #RV = bed.read(PositionRange(int(wnds[wnd_i,-2]),int(wnds[wnd_i,-1])))
        RV = plink_reader.readBED(bfile, useMAFencoding=True, blocksize = 1, start = int(wnds[wnd_i,4]), nSNPs = int(wnds[wnd_i,5]), order  = 'F',standardizeSNPs=False,ipos = 2,bim=bim,fam=fam)
        
        Xr = RV['snps']
        if perm_i is not None:
            Xr = Xr[perm,:]
        rv = mtSet.optimize(Xr)
        line = NP.concatenate([wnds[wnd_i,:],rv['LLR']])
        wnd_file.writerow(line)
    pass

def analyze(options):

    # load data
    print 'import data'
    if options.cfile is None:
        cov = {'eval':None,'evec':None}
        warnings.warn('warning: cfile not specifed, a one variance compoenent model will be considered')
    else:
        cov = readCovarianceMatrixFile(options.cfile,readCov=False)
    Y = readPhenoFile(options.pfile,idx=options.trait_idx)
    null = readNullModelFile(options.nfile)
    wnds = readWindowsFile(options.wfile)

    F = None
    if options.ffile:
        F = readCovariatesFile(options.ffile)
        #null['params_mean'] = SP.loadtxt(options.nfile + '.f0')
        

    if F is not None: assert Y.shape[0]==F.shape[0], 'dimensions mismatch'

            
    if options.i0 is None: options.i0 = 1
    if options.i1 is None: options.i1 = wnds.shape[0]

    # name of output file
    if options.perm_i is not None:
        res_dir = os.path.join(options.resdir,'perm%d'%options.perm_i)
    else:
        res_dir = os.path.join(options.resdir,'test')
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    n_digits = len(str(wnds.shape[0]))
    fname = str(options.i0).zfill(n_digits)
    fname+= '_'+str(options.i1).zfill(n_digits)+'.res'
    resfile = os.path.join(res_dir,fname)

    # analysis
    t0 = time.time()
    scan(options.bfile,Y,cov,null,wnds,options.minSnps,options.i0,options.i1,options.perm_i,resfile,F)
    t1 = time.time()
    print '... finished in %s seconds'%(t1-t0)

