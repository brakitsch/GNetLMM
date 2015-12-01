import GNetLMM.pycore.modules.utils as utils
import reader
import numpy as np

class GenoReaderMatrix(reader.MatrixReader):
    def __init__(self, M, pos, chrom, ids):
        """
        constructor

        input:
        M       :   genotype matrix [FxN]
        pos     :   positions [F]
        chrom   :   chromosomal positions [F]
        ids     :   snp identifiers [TF]
        """
        self.M = M
        self.pos = pos
        self.chrom = chrom
        self.ids = ids
        

    def getSnpPos(self):
        return self.pos

    def getSnpChrom(self):
        return self.chrom

    def getSnpIds(self):
        return self.ids

    def get_nrows(self):
        return self.M.shape[0]

    def get_ncols(self):
        return self.M.shape[1]

    def compute_Kbg_excluded_chrom(self,chrom_i):
        idx_bg = self.chrom!=chrom_i
        M = self.M[idx_bg]
        Kbg = utils.computeLinearKernel(M.T)
        return Kbg

    def get_chrom_idx(self,chrom_i):
        idx = np.nonzero(self.chrom==chrom_i)[0]
        start = idx[0]
        nSnps = len(idx)
        return start, nSnps
        
        return self.chrom==chrom_i

    def loadSnpBlock(self,start = 0, nSNPs = np.inf):
      """
      read [basefilename].bed,[basefilename].bim,[basefilename].fam
      --------------------------------------------------------------------------
      Input:
    
      start           : index of the first SNP to be loaded from the .bed-file (default 0)
      nSNPs           : load nSNPs from the .bed file (default SP.inf, meaning all)
      """
      return self.M[start:start+nSNPs,:]
      

        
