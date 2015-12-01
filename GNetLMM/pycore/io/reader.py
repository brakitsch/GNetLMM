import pdb
import numpy as np
import os
import linecache

class MatrixReader:
    """
    container for holding matrix or hdf object
    """
    def __init__(self,M):
        """
        constructor

        input:
        M   :   data matrix/hdf object
        """
        self.M = M

    def getMatrix(self):
        """
        returns complete matrix
        """
        return self.M[:]

    def getRows(self,idx):
        """
        returns the rows M[idx] of the matrix

        input:
        idx   :   row indices
        """
        return self.M[idx,:]


class FileReader:
    """
    container for holding phenotype files of the following format
    - phenotypes are saved in basefile.matrix
    - column information is saved in basefile.cols
    - row information is saved in basefile.rows
    """
    
    def __init__(self, basefile, load_rowinfo=True, load_colinfo=False):
        """
        constructor

        input:
        basefile      : name of basefile
        load_rowinfo  : if True, row information is loaded (default: True)
        load_colinfo  : if True, column information is loaded (default: False)
        """
        self.basefile = basefile
        if load_rowinfo: self.row_info = self.getInfo('rows')
        if load_colinfo: self.col_info = self.getInfo('cols')
        
    def getMatrix(self):
        """
        returns complete matrix
        """
        M = np.loadtxt(self.basefile + '.matrix')
        return M

    def getRows(self,idx):
        """
        returns the rows M[idx] of the matrix

        input:
        idx   :   row indices
        """
        RV = []
        for i in idx:
            line = linecache.getline(self.basefile + '.matrix', i+1)
            line = line.split(' ')
            RV.append(line)

        RV = np.array(RV,dtype=float)
        return RV
        
    def getInfo(self, which):
        """
        loads the row information into memory
        """
        fn = self.basefile+'.' + which
        if not(os.path.exists(fn)): return None
        M = np.loadtxt(fn,dtype=str)
        if M.ndim==1: M = M[:,np.newaxis]
        header = M[0]
        data = {}
        for ikey,key in enumerate(header):
            arr = M[:,ikey][1:]
            try:
                arr = np.array(arr, dtype=int)
            except:
                pass
            data[key] = arr
            
        return data

 