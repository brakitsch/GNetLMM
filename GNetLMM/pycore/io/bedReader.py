


import numpy as np
import pdb




class BedReader():
    def __init__(self, basefilename):
        """
        constructor

        input:
        basefilename   : name of basefile
        """
        self.basefilename = basefilename
        self.bim = self.readBIM()
        self.fam = self.readFAM()

    def get_ncols(self):
        return self.fam.shape[0]

    def get_nrows(self):
        return self.bim.shape[0]


    def getRowInfo(self):
        return {'rs':self.bim[:,1], 'chrom':self.bim[:,0], 'pos':self.bim[:,3]}
 
        
    def getSnpPos(self):
        try:
            return np.array(self.bim[:,3], dtype=int)
        except:
            pass
        return self.bim[:,3]

    def getSnpChrom(self):
        try:
            return np.array(self.bim[:,0], dtype=int)
        except:
            pass
        return self.bim[:,0]

    def getSnpIds(self):
        return self.bim[:,1]

    def readBIM(self,usecols=None):
        """
        helper method for speeding up read BED
        """
        bim = self.basefilename+ '.bim'
        bim = np.loadtxt(bim,dtype = 'str',usecols=usecols)
        return bim
    

    def readFAM(self,usecols=None):
        """
        helper method for speeding up read FAM
        """
        fam = self.basefilename+'.fam'
        fam = np.loadtxt(fam,dtype = 'str',usecols=usecols)
        return fam

    def loadSnpBlock(self,start=0, nSNPs=np.inf):
        return self.readBED(start=start, nSNPs=nSNPs)['snps']


    def compute_Kbg_excluded_chrom(self,chrom_i):
        raise Exception('not implemented yet')
    
    def readBED(self,blocksize = 1, start = 0, nSNPs = np.inf, startpos = None, endpos = None, order  = 'F',ipos = 2):
      """
      read [basefilename].bed,[basefilename].bim,[basefilename].fam
      --------------------------------------------------------------------------
      Input:
      
      blocksize       : load blocksize SNPs at a time (default 1)
      start           : index of the first SNP to be loaded from the .bed-file (default 0)
      nSNPs           : load nSNPs from the .bed file (default SP.inf, meaning all)
      startpos        : starting position of the loaded genomic region[chr,bpdist]
      endpos          : end-position of the loaded genomic region     [chr,bpdist]
      order           : memory layout of the returned SNP array (default 'F')
                        'F'   : Fortran-style column-major array (SNP-major)
                        'C'   : C-style row-major array (individual-major)
      standardizeSNPs : bool indeicator if the resulting SNP array is supposed to 
      be zero-mean and unit-vatiance with mean imputed missing
      values (default False)
      ipos            : the index of the position index to use (default 2)
                        1 : genomic distance
                        2 : base-pair distance
      --------------------------------------------------------------------------
      Output dictionary:
      'rs'     : [S] array rs-numbers
      'pos'    : [S*3] array of positions [chromosome, genetic dist, basepair dist]
      'snps'   : [N*S] array of snp-data
      'iid'    : [N*2] array of family IDs and individual IDs
       --------------------------------------------------------------------------
      """
      rs = self.bim[:,1]
      pos = np.array(self.bim[:,(0,2,3)],dtype = 'float')

      if startpos is not None:
          #pdb.set_trace()
          i_c = pos[:,0]==startpos[0]
          i_largerbp = pos[:,ipos]>=startpos[ipos]
          start = which(i_c * i_largerbp)
          while (start-1 >= 0 and pos[start-1,ipos] == startpos[ipos]):
              start = start -1
              i_c = pos[:,0]==endpos[0]
              i_smallerbp = pos[:,ipos]>=endpos[ipos]
              end = which(i_c * i_smallerbp)
              while (end+1 < pos.shape[0] and pos[end+1,ipos] == endpos[ipos]):
                  end = end + 1
                  nSNPs = end - start
                  if (nSNPs<=0) or (end==0) or (start<=0):
                    ret = {
                        'pos':np.zeros((0,3)),
                        'rs':np.zeros((0)),
                        'iid':self.fam,
                        'snps':np.zeros((self.fam.shape[0],0))
                        }
                    return ret
        
      N = self.fam.shape[0]
      S = self.bim.shape[0]
      S_res = min(S,start + nSNPs)
      nSNPs = min(S-start,nSNPs)
  
      if nSNPs<=0:
          ret = {
              'rs'     :rs[start:start],
              'pos'    :pos[start:start,:],
              #'snps'   :SNPs[0:N,start:start],
              'snps'   :np.zeros((N,0)),
              'iid'    :self.fam
              }
          return ret
       
      SNPs = np.zeros(((np.ceil(0.25*N)*4),nSNPs),order=order)
      bed = self.basefilename + '.bed'
      with open(bed, "rb") as f:
          mode = f.read(2)
          if mode != 'l\x1b':
              raise Exception('No valid binary PED file')
          mode = f.read(1) #\x01 = SNP major \x00 = individual major
          if mode != '\x01':
              raise Exception('only SNP-major is implemented')
          startbit = np.ceil(0.25*N)*start+3
          f.seek(startbit)
     
          for blockStart in np.arange(0,nSNPs,blocksize):
              blockEnd = min(S,blockStart+blocksize)
              Sblock = min(nSNPs-blockStart,blocksize)
              nbyte = int(np.ceil(0.25*N)*Sblock)
              bytes = np.array(bytearray(f.read(nbyte))).reshape((np.ceil(0.25*N),Sblock),order='F')
            
              SNPs[3::4,blockStart:blockEnd][bytes>=64]=np.nan
              SNPs[3::4,blockStart:blockEnd][bytes>=128]=1
              SNPs[3::4,blockStart:blockEnd][bytes>=192]=2
              bytes=np.mod(bytes,64)
              SNPs[2::4,blockStart:blockEnd][bytes>=16]=np.nan
              SNPs[2::4,blockStart:blockEnd][bytes>=32]=1
              SNPs[2::4,blockStart:blockEnd][bytes>=48]=2
              bytes=np.mod(bytes,16)
              SNPs[1::4,blockStart:blockEnd][bytes>=4]=np.nan
              SNPs[1::4,blockStart:blockEnd][bytes>=8]=1
              SNPs[1::4,blockStart:blockEnd][bytes>=12]=2
              bytes=np.mod(bytes,4)
              SNPs[0::4,blockStart:blockEnd][bytes>=1]=np.nan
              SNPs[0::4,blockStart:blockEnd][bytes>=2]=1
              SNPs[0::4,blockStart:blockEnd][bytes>=3]=2
  
      snps = SNPs[0:N,:]

 
      ret = {
           'rs'     :rs[start:S_res],
           'pos'    :pos[start:S_res,:],
           'snps'   :snps.T,
           'iid'    :self.fam
           }
      return ret


   
