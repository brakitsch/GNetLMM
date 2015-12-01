
import reader
import pdb
import GNetLMM.pycore.modules.utils as utils


class PhenoReaderFile(reader.FileReader):

    def getGeneStart(self):
        return self.row_info['gene_start']

    def getGeneEnd(self):
        if 'gene_end' in self.row_info:
            gene_end =  self.row_info['gene_end']
        elif 'gene_stop' in self.row_info:
            gene_end =  self.row_info['gene_stop']
        else:
            gene_end =  self.row_info['gene_start']
        return gene_end
        
    def getGeneChrom(self):
        gene_chrom = self.row_info['gene_chrom']
        return gene_chrom

    def getGeneIds(self):
        if 'gene_ids' in self.row_info:
            return self.row_info['gene_ids']

    def get_nrows(self):
        keys = self.row_info.keys()
        return len(self.row_info[keys[0]])


class PhenoReaderMatrix(reader.MatrixReader):

    def __init__(self, M, gene_start, gene_stop, chrom,ids):
        """
        constructor

        input:
        M       :   phenotype matrix [TxN]
        pos     :   positions [T]
        chrom   :   chromosomal positions [T]
        ids     :   snp identifiers [T]
        """
        self.M = M
        self.gene_start = gene_start
        self.gene_stop  = gene_stop
        self.chrom = chrom
        self.ids = ids

    def getGeneStart(self):
        return self.gene_start

    def getGeneEnd(self):
        return self.gene_stop
        
    def getGeneChrom(self):
        return self.chrom

    def getGeneIds(self):
        return self.ids

    def get_nrows(self):
        return self.M.shape[0]

    def get_ncols(self):
        return self.M.shape[1]
 
