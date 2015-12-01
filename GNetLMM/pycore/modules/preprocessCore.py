
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
import glob
import limix.io.plink as plink

import gnetlmm

import GNetLMM.pycore.io.bedReader as bedReader
import GNetLMM.pycore.io.phenoReader as phenoReader
import GNetLMM.pycore.io.reader as reader
import GNetLMM.pycore.io.writer as writer


def merge_files(fns_in, fn_out):
    """
    writing all files fns_ins into fn_out
    """
    with open(fn_out, 'w') as fout:
        for fn_in in fns_in:
            with open(fn_in) as fin:
                for line in fin:
                    fout.write(line)


def merge_assoc0_scan(assoc0file):
    """
    merging associations files
    """
    fn_beta = glob.glob(assoc0file + '*.beta.matrix')
    fn_pv = glob.glob(assoc0file + '*.pv.matrix')
    merge_files(fn_beta, assoc0file + '.beta.matrix')
    merge_files(fn_pv, assoc0file + '.pv.matrix')
    

    
def computeCovarianceMatrixPlink(plink_path,out_dir,bfile,cfile,sim_type='RRM'):
    """
    computing the covariance matrix via plink
    """
    
    print "Using plink to create covariance matrix"
    cmd = '%s --bfile %s '%(plink_path,bfile)

    if sim_type=='RRM':
        # using variance standardization
        cmd += '--make-rel square '
    else:
        raise Exception('sim_type %s is not known'%sim_type)

    cmd+= '--out %s'%(os.path.join(out_dir,'plink'))
    
    subprocess.call(cmd,shell=True)

    # move file to specified file
    if sim_type=='RRM':
        old_fn = os.path.join(out_dir, 'plink.rel')
        os.rename(old_fn,cfile+'.cov')
        
        old_fn = os.path.join(out_dir, 'plink.rel.id')
        os.rename(old_fn,cfile+'.cov.id')

    if sim_type=='IBS':
        old_fn = os.path.join(out_dir, 'plink.mibs')
        os.rename(old_fn,cfile+'.cov')

        old_fn = os.path.join(out_dir, 'plink.mibs.id')
        os.rename(old_fn,cfile+'.cov.id')

    os.remove(os.path.join(out_dir, 'plink.nosex'))
    os.remove(os.path.join(out_dir, 'plink.log'))                 
    

    
def computeCovarianceMatrixPython(out_dir,bfile,cfile,sim_type='RRM'):
    print "Using python to create covariance matrix. This might be slow. We recommend using plink instead."

    if sim_type is not 'RRM':
        raise Exception('sim_type %s is not known'%sim_type)

    """ loading data """
    data = plink.readBED(bfile)
    iid  = data['iid']
    X = data['snps']
    N = X.shape[1]
    print '%d variants loaded.'%N
    print '%d people loaded.'%X.shape[0]
    
    """ normalizing markers """
    print 'Normalizing SNPs...'
    p_ref = X.mean(axis=0)/2.
    X -= 2*p_ref

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        X /= sp.sqrt(2*p_ref*(1-p_ref))
        
    hasNan = sp.any(sp.isnan(X),axis=0)
    print '%d SNPs have a nan entry. Exluding them for computing the covariance matrix.'%hasNan.sum()

    """ computing covariance matrix """
    print 'Computing relationship matrix...'
    K = sp.dot(X[:,~hasNan],X[:,~hasNan].T)
    K/= 1.*N
    print 'Relationship matrix calculation complete'
    print 'Relationship matrix written to %s.cov.'%cfile
    print 'IDs written to %s.cov.id.'%cfile

    """ saving to output """
    np.savetxt(cfile + '.cov', K, delimiter='\t',fmt='%.6f')
    np.savetxt(cfile + '.cov.id', iid, delimiter=' ',fmt='%s')
    





def eighCovarianceMatrix(cfile):
    """
    compute similarity matrix using plink

    Input:
    cfile        :   the covariance matrix will be read from cfile.cov while the eigenvalues and the eigenverctors will
                        be written to cfile.cov.eval and cfile.cov.evec respectively
    """
    # precompute eigenvalue decomposition
    K = np.loadtxt(cfile+'.cov')
    K+= 1e-4*sp.eye(K.shape[0])
    S,U = la.eigh(K); S=S[::-1]; U=U[:,::-1]
    np.savetxt(cfile+'.cov.eval',S,fmt='%.6f')
    np.savetxt(cfile+'.cov.evec',U,fmt='%.6f')


def computeCovarianceMatrix(plink_path,bfile,cfile,sim_type='RRM'):
    """
    compute similarity matrix using plink

    Input:
    plink_path   :   plink path
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    cfile        :   the covariance matrix will be written to cfile.cov and the corresponding identifiers
                         to cfile.cov.id. If not specified, the covariance matrix will be written to cfile.cov and
                         the individuals to cfile.cov.id in the current folder.
    sim_type     :   {IBS/RRM} are supported
    """
    try:
        output    = subprocess.check_output('%s --version --noweb'%plink_path,shell=True)
        use_plink = float(output.split(' ')[1][1:-3])>=1.9
    except:
        use_plink = False

    assert bfile!=None, 'Path to bed-file is missing.'
    assert os.path.exists(bfile+'.bed'), '%s.bed is missing.'%bfile
    assert os.path.exists(bfile+'.bim'), '%s.bim is missing.'%bfile
    assert os.path.exists(bfile+'.fam'), '%s.fam is missing.'%bfile

    # create dir if it does not exist
    out_dir = os.path.split(cfile)[0]
    if out_dir!='' and (not os.path.exists(out_dir)):
        os.makedirs(out_dir)


    if use_plink:
        computeCovarianceMatrixPlink(plink_path,out_dir,bfile,cfile,sim_type=sim_type)
    else:
        computeCovarianceMatrixPython(out_dir,bfile,cfile,sim_type=sim_type)



def initial_scan(bfile, pfile, cfile, ffile, assoc0file, startSnpIdx=0, nSnps=np.inf,
                 memory_efficient = False):
    """
    running initial scan using a standard linear mixed model

    Input:
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    pfile        :   phenotype file
    cfile        :   covariance matrix file
    ffile        :   covariates file

    assoc0file   :   basename of output file 
    """
    K = None
    if cfile is not None:
        K = np.loadtxt(cfile)

    Covs = None
    if ffile is not None:
        Covs = np.loadtxt(ffile)

    preader = phenoReader.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader,K=K,Covs=Covs)
    beta0, pv0 = model.initial_scan(startSnpIdx, nSnps,memory_efficient)

    write = writer.Writer(assoc0file+'.pv')
    write.writeMatrix(pv0, fmt='%.4e')

    write = writer.Writer(assoc0file+'.beta')
    write.writeMatrix(beta0, fmt='%.4f')

   

def marginal_genecorr(bfile,pfile, gfile):
    """
    running marginal gene-gene correlations

    Input:
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    pfile        :   phenotype file
    cfile        :   covariance matrix file
    ffile        :   covariates file

    gfile        :   basename of output file 
    """
    preader = phenoReader.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader)
    corr, pv = model.marginal_gene_correlations()
    write = writer.Writer(gfile+'.pv')
    write.writeMatrix(pv, fmt='%.4e')
    write = writer.Writer(gfile+'.corr')
    write.writeMatrix(corr, fmt='%.4f')
    
  
def gene_has_cis_anchor(bfile, pfile, assoc0, cis_thresh, cis_window, cisfile):
    """
    tests if a gene has a cis anchor

    input:
    bfile        :   binary bed file (bfile.bed, bfile.bim and bfile.fam are required)
    pfile        :   phenotype file
    assoc0file   :   basefilename of initial association scan
    cis_thresh   :   thrshold for cis-association
    cis_window   :   maximal distance between cis-snp and gene
    cis_file     :   filename for saving cis assocaitions
    """
    preader = phenoReader.PhenoReaderFile(pfile)
    greader =  bedReader.BedReader(bfile)
    model = gnetlmm.GNetLMM(preader,greader, window=cis_window)

    assoc0Reader = reader.FileReader(assoc0 + '.pv')
    model.set_assoc0_reader(assoc0Reader)
    model.gene_has_cis_anchor(cis_thresh)
    model.save_cis_anchors(cisfile)

    

def preprocess(options):

    """ computing the covariance matrix """
    if options.compute_cov:
       assert options.bfile!=None, 'Please specify a bfile.'
       assert options.cfile is not None, 'Specify covariance matrix basename'
       print 'Computing covariance matrix'
       t0 = time.time()
       computeCovarianceMatrix(options.plink_path,options.bfile,options.cfile,options.sim_type)
       t1 = time.time()
       print '... finished in %s seconds'%(t1-t0)
       #print 'Computing eigenvalue decomposition'
       #t0 = time.time()
       #eighCovarianceMatrix(options.cfile) 
       #t1 = time.time()
       #print '... finished in %s seconds'%(t1-t0)

    """ running initial scan """
    if options.initial_scan:
        assert options.bfile is not None, 'Please specify a bfile.'
        assert options.pfile is not None, 'Please specify a phenotypic file.'
        assert options.assoc0file is not None, 'Please specify an output file for saving the associations results.'
        t0 = time.time()
        print 'Running initial scan'
        initial_scan(options.bfile,options.pfile,options.cfile,options.ffile, options.assoc0file,startSnpIdx=options.startSnpIdx, nSnps=options.nSnps,
                     memory_efficient=options.memory_efficient)
        t1 = time.time()
        print '... finished in %s seconds'%(t1-t0)
        
    """ computing marginal gene-gene correlations """
    if options.gene_corr:
        assert options.bfile is not None, 'Please specify a bfile.'
        assert options.pfile!=None, "Please specify a phenotypic file."
        assert options.gfile!=None, "Please specify an output file for saving the gene-gene correlations"
        t0 = time.time()
        print "Computing marginal gene-gene correlations"
        marginal_genecorr(options.bfile, options.pfile,options.gfile)
        t1=time.time()
        print '... finished in %s seconds'%(t1-t0)

    """ determining cis anchors """
    if options.cis_anchor:
        assert options.bfile is not None, 'Please specify a bfile.'
        assert options.pfile!=None, "Please specify a phenotypic file."
        assert options.assoc0file is not None, "Please specifz assoc0 basefilename."
        assert options.cisfile is not None, 'Please specify file for saving cis anchors.'
        t0 = time.time()
        print "Determining genes having a cis anchor"
        gene_has_cis_anchor(options.bfile,options.pfile, options.assoc0file, options.cis_thresh, options.cis_window,options.cisfile)
        t1=time.time()
        print '... finished in %s seconds'%(t1-t0)

    """ merging assoc0 files """
    if options.merge_assoc0_scan:
        assert options.assoc0file is not None, 'Please specify assoc0 basefilename.'
        t0 = time.time()
        print 'Merge assoc0 files'
        merge_assoc0_scan(options.assoc0file)
        t1 = time.time()
        print '.... finished in %s seconds'%(t1-t0)
