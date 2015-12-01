
import os
import subprocess
import pdb
import sys
import numpy as np
import numpy.linalg as la
import scipy as sp
import warnings
from optparse import OptionParser
import time
import scipy as sp
import warnings
import scipy.sparse.linalg as ssla
import pdb
import glob

import GNetLMM.pycore.io.bedReader as bedReader
import GNetLMM.pycore.io.phenoReader as phenoReader
import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.io.writer as writer

import gnetlmm
import assoc_results


def merge_files(fns_in, fn_out):
    """
    writing all files fns_ins into fn_out
    """
    with open(fn_out, 'w') as fout:
        for fn_in in fns_in:
            with open(fn_in) as fin:
                for line in fin:
                    fout.write(line)
                    
def concatenate(files):
    fns_in = glob.glob(files + "*.csv")
    fn_out = files + '.csv'
    merge_files(fns_in, fn_out)
    

def scan(bfile,pfile,cfile,ffile,vfile,assocfile,startTraitIdx,nTraits):
    """
    running association scan

    input:
    bfile      :   basefilename of plink file
    pfile      :   phenotype file
    cfile      :   covariance file
    ffile      :   fixed effects file
    vfile      :   file containing vstructures
    assocfile  :   file for saving results
    """
    K = None
    if cfile is not None:
        K = np.loadtxt(cfile)

    Covs = None
    if ffile is not None:
        Covs = np.loadtxt(ffile)
        
    preader = phenoReader.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader, Covs=Covs,K=K)
    model.load_vstructures(vfile+".csv")
    model.update_associations(startTraitIdx, nTraits)
    model.save_updates(assocfile)

    

def find_vstructures(bfile, pfile,gfile,cisfile, assoc0file,window,vfile,startTraitIdx,nTraits):
    """
    running association scan

    input:
    pfile      :   phenotype file
    cfile      :   covariance file
    ffile      :   fixed effects file
    cisfile    :   file containing cis anchors
    assoc0file :   file containing the results from the initial association scan
    vfile      :   file containing v-structures
    """
    preader = phenoReader.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    
    model = gnetlmm.GNetLMM(preader, greader, window=window)

    genecorr_reader = reader.FileReader(gfile + '.pv')
    model.set_genecorr_reader(genecorr_reader)
    
    assoc0Reader = reader.FileReader(assoc0file + '.pv')
    model.set_assoc0_reader(assoc0Reader)
    
    model.load_cis_anchors(cisfile)
    model.find_vstructures(startTraitIdx, nTraits)
    model.save_vstructures(vfile+'.csv')


def merge_results(assocfile, assoc0file):
    """
    merging association updates with lmm-scan
    """
    results = assoc_results.AssocResultsList()
    results.load_csv(assocfile + '.csv')
    results.save_matrix(assoc0file, assocfile)


def analyse(options):

    """ finding v-structures """
    if options.find_vstructures:
        t0 = time.time()
        print 'Finding v-structures'
        assert options.bfile!=None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify the phenotype file'
        assert options.gfile is not None, 'Please specify the gene-gene correlation file'
        assert options.cisfile is not None, 'Please specify the cis-anchor file'
        assert options.assoc0file is not None, 'Please specify the assoc0 file'
        assert options.vfile is not None, 'Please specify output file for vstructures'
        find_vstructures(options.bfile,options.pfile,options.gfile,options.cisfile,options.assoc0file, options.window, options.vfile,options.startTraitIdx,options.nTraits)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
    """ updating associations """
    if options.update_assoc:
        t0 = time.time()
        print 'Updating associations'
        assert options.bfile!=None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify the phenotype file'
        assert options.vfile is not None, 'Please specify vstructure file'
        assert options.assocfile is not None, 'Please specify an output file'
        scan(options.bfile,options.pfile,options.cfile,options.ffile,options.vfile, options.assocfile,options.startTraitIdx,options.nTraits)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
    """ merging results """
    if options.merge_assoc:
        t0 = time.time()
        print 'Merging associations'
        assert options.assocfile is not None, 'Please specify an output file'
        assert options.assoc0file is not None, 'Please specify assoc0 results'
        merge_results(options.assocfile, options.assoc0file)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
        
    """ concatenate files """
    if options.concatenate:
        t0 = time.time()
        print 'Concatenating files'
        assert options.files is not None, 'Please specify file start'
        concatenate(options.files)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
