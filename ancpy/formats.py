from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import sys
import itertools as it
import pandas as pd
import numpy as np
import pysam


class Genotypes(object):
    """A class for storing genotype data which at its core 
    consists of a dense p x n genotype matrix, a dataframe 
    of metadata associated which each snp, and a list of 
    individual ids. We also provide computation options for 
    setting up dataset for commonly using analyses

    Arguments
    ---------
    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    Attributes
    ----------
    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int 
        number of snps (features)

    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, y=None, snp_df=None, inds=None):

        # p x n genotype matrix
        self.y = y

        # DataFrame of SNP level information
        self.snp_df = snp_df
        
        # list of individual ids
        self.inds = inds

        if ((type(y) != type(None)) and 
            (type(snp_df) != type(None)) and 
            (type(inds) != type(None))):
            # number of individuals
            self.n = len(self.inds)

            # number of snps
            self.p = self.snp_df.shape[0]
        else:
            # number of individuals
            self.n = None

            # number of snps
            self.p = None

        # below is adapted from 
        # https://github.com/bulik/ldsc/blob/master/munge_sumstats.py

        # complement of each base
        self.base_complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

        # bases 
        self.bases = self.base_complement.keys()

        # dict of allele combinations that are strand ambiguous
        self.strand_ambiguous_alleles = {''.join(x): x[0] == self.base_complement[x[1]]
                                         for x in it.product(self.bases, self.bases)
                                         if x[0] != x[1]} 

    def remove_chrom_prefix(self):
        """Remove the chromsome prefix from the snp_df
        this helps for making datasets consistent before 
        merging
        """
        self.snp_df["CHROM"] = self.snp_df["CHROM"].apply(lambda x: int(x[3:]))

    def keep_autosomes(self):
        """Keep snps only present on the autosomes
        """
        # filter out snps not on the autosome
        self.snp_df = self.snp_df[self.snp_df["CHROM"] <= 22]

    def remove_rare_common(self, eps):
        """Removes SNPs that are too rare or 
        too common given some threshold eps

        Arguments
        ---------
        eps : float
            frequency threshold eps
        """
        # compute allele frequencies at each SNP
        f = np.nansum(self.y, axis=1) / (2 * self.n)

        # get the indicies of SNPs that dont pass the thresh
        idx = np.where((f >= eps) & (f <= (1. - eps)))[0]
        
        # keep only these SNPs 
        self.snp_df = self.snp_df.iloc[idx]
    
    def remove_strand_ambiguous_alleles(self):
        """Removes any snps that have alleles that
        are strand ambiguous. 

        TODO: speed this up its currently a bit slow
        also add removal of incosistent alleles
        """
        # lambda function for applying to each snp
        strand_unamb_func = lambda row: not self.strand_ambiguous_alleles[row["A1"] + row["A2"]]

        # snps that have unambigous strand
        idx = self.snp_df.apply(strand_unamb_func, axis=1)

        # index the the snps
        self.snp_df = self.snp_df[idx]
    
    def reindex(self):
        """Reset the index column in the pandas dataframe
        """
        # get indicies of snps to keep
        idx = self.snp_df["idx"].tolist()

        # reindex data matrix
        self.y = self.y[idx, :] 

        # re-compute number of snps
        self.p = self.snp_df.shape[0]

        # reset index column
        self.snp_df["idx"] = np.arange(0, self.p, 1)
    
    def binarize(self):
        """Convert genotype data to binary
        by mapping 2->1, 0->0 and randomly 
        selecting a 0 or 1 for heterozygotes.
        This emulates the read sampled data found
        in ancient DNA studies
        """
        # randomly sample hets
        p_het = np.sum(self.y[self.y == 1.])
        if p_het > 0:
            samps = np.random.binomial(1, .5, int(p_het)).astype(np.float32) 
            np.place(self.y, self.y == 1., samps)

        # set 2s to 1s
        self.y[self.y == 2.] = 1.

    def normalize(self, scale, impute):
        """Mean center the data so each snp has mean 0. 
        If scale true then scale each SNP to have var 1. 
        If impute is true it sets missing sites to 0. i.e. 
        the approx mean. Note for all computations we 
        ignore nans.

        Arguments
        ---------
        scale : bool
            scale the data to have variance 1

        impute : bool
            predict missing sites with the expected mean 0
        """
        if scale:
            # compute mean genotype for every SNP
            mu = np.nanmean(self.y, axis=1)

            # compute std dev for every SNP
            std = np.nanstd(self.y, axis=1)

            # center and scale
            zt = (self.y.T - mu) / std
            self.z = zt.T
        else:
            # compute mean genotype for every SNP
            mu = np.nanmean(self.y, axis=1)

            # center
            zt = self.y.T - mu
            self.z = zt.T

        if impute:
            # set nans to zero
            self.z = np.nan_to_num(self.z)


class AncestryMap(Genotypes):
    """A class for the eigenstrat genotype format
    which consists of a geno, snp, and ind file

    Arguments
    ---------
    prefix : str
        prefix of  path to data files

    Attributes
    ----------
    prefix : str
        prefix of  path to data files

    y : np.array
        n x p genotype matrix

    snp_df : pd.DataFrame

        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int
        number of snps (features)

    geno_path : str
        path to eigenstrat geno file

    ind_path : str
        path to eigenstrat ind file

    snp_path : str
        path to eigenstart snp file
    
    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, prefix):
    
        # inherit attributes from Genotypes
        super().__init__()
    
        # prefix to paths
        self.prefix = prefix
        
        # path to eigenstart geno file
        self.geno_path = "{}.geno".format(self.prefix)
        
        # path to eigenstrat ind file
        self.ind_path = "{}.ind".format(self.prefix)
        
        # path to eigenstrat snp file
        self.snp_path = "{}.snp".format(self.prefix)

    def _get_snp_dataframe(self):
        """Gets snp level dataframe stored as 
        pandas DataFrame
        """
        # read into pandas df
        self.snp_df = pd.read_table(self.snp_path, header=None,
                                    delim_whitespace=True)
        
        # check if the dataframe has 6 columns
        if self.snp_df.shape[1] < 6:
            raise ValueError("{}.snp must have 6 cols".format(self.snp_path))
        
        # name columns
        self.snp_df.columns = ["SNP", "CHROM", "CM_POS", "POS", "A1", "A2"]
        
        # number of snps
        self.p = self.snp_df.shape[0]
        
        # add rownumbers
        self.snp_df["idx"] = np.arange(0, self.p, 1)

    def _get_inds(self):
        """Get list of individual ids
        """
        ind_df = pd.read_table(self.ind_path, header=None,
                               delim_whitespace=True, delimiter="\t")
        
        # check if the dataframe has 6 columns
        if ind_df.shape[1] < 3:
            raise ValueError("{}.ind must have 3 cols".format(self.ind_path))

        ind_df.columns = ["IND", "SEX", "CLST"]
        self.inds = ind_df["IND"].tolist()
        
        # number of inds
        self.n = len(self.inds)


class PackedAncestryMap(AncestryMap):
    """Class for packed ancestry map eigenstrat format. 
    The packed format's geno is binary so it requires a 
    couple steps processing before being loaded into a 
    numpy array. We modify some code from
    
    https://github.com/mathii/pyEigenstrat
    
    Arguments
    ---------
    prefix : str
        prefix of  path to data files

    Attributes
    ----------
    prefix : str
        prefix of  path to data files

    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int
        number of snps (features)

    geno_path : str
        path to eigenstrat geno file

    ind_path : str
        path to eigenstrat ind file

    snp_path : str
        path to eigenstart snp file
    
    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, prefix):

        # inhert attributes from AncestryMap
        super().__init__(prefix)

        # get the snp dataframe
        super()._get_snp_dataframe()

        # get list of individuals
        super()._get_inds()

        # get the genotype matrix
        self._get_genotype_matrix()

    def _get_genotype_matrix(self):
        """Gets the genotype matrix stored as a
        numpy array
        """
        rlen = max(48, int(np.ceil(self.n * 2 / 8)))

        # read in binary 
        self.y = np.fromfile(self.geno_path, dtype='uint8')[rlen:]
        self.y.shape = (self.p, rlen)

        # unpack 
        self.y = np.unpackbits(self.y, axis=1)[:, :(2 * self.n)]
        self.y = 2 * self.y[:, ::2] + self.y[:, 1::2]

        # convert to float to allow missing data stored as nan
        self.y = self.y.astype(np.float32)

        # set missing data to nan
        self.y[self.y == 3] = np.nan


class UnpackedAncestryMap(AncestryMap):
    """Class for unpacked ancestry map eigenstrat format. The unpacked format's
    geno is not binary so it doesn't require the same processing before 
    being loaded into a numpy array. 
    
    Arguments
    ---------
    prefix : str
        prefix of  path to data files

    Attributes
    ----------
    prefix : str
        prefix of  path to data files

    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int
        number of snps (features)

    geno_path : str
        path to eigenstrat geno file

    ind_path : str
        path to eigenstrat ind file

    snp_path : str
        path to eigenstart snp file
    
    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, prefix):

        # inhert attributes from AncestryMap
        super().__init__(prefix)

        # get the snp dataframe
        super()._get_snp_dataframe()

        # get list of individuals
        super()._get_inds()

        # get the genotype matrix
        self._get_genotype_matrix()

    def _get_genotype_matrix(self):
        """Gets the genotype matrix stored as a
        numpy array
        """
        # read the geno file
        with open(self.geno_path, "r") as f:
            matstr = f.read().replace("\n", "")

        # convert to long array
        ys = np.fromstring(matstr, dtype=np.int8) - 48

        # reshape to p x n matrix
        self.y = ys.reshape(self.p, self.n).astype(np.float32)

        # replace 9s with nans
        self.y[self.y == 9] = np.nan


class PlinkTraw(Genotypes):
    """Class for plink traw format

    TODO: double check COUNTED vs ALT
    
    Arguments
    ---------
    prefix : str
        prefix of  path to data files

    Attributes
    ----------
    prefix : str
        prefix of  path to data files

    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int
        number of snps (features)

    traw_path : str
        path to eigenstrat geno file
    
    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, prefix):

        # inherit attributes from Genotypes
        super().__init__()

        # prefix to paths
        self.prefix = prefix

        # path to traw file
        self.traw_path = "{}.traw".format(self.prefix)

        # traw dataframe
        traw_df = pd.read_table(self.traw_path, header=0, delimiter="\t")

        self._get_snp_dataframe(traw_df)
        self._get_inds(traw_df)
        self._get_genotype_matrix(traw_df)

    def _get_snp_dataframe(self, traw_df):
        """Gets snp level dataframe stored as 
        pandas DataFrame
        """
        # create snp dataframe
        self.snp_df = pd.DataFrame({"SNP": traw_df["SNP"].tolist(),
                                    "CHROM": traw_df["CHR"].tolist(),
                                    "CM_POS": traw_df["(C)M"].tolist(),
                                    "POS": traw_df["POS"].tolist(),
                                    "A1": traw_df["COUNTED"].tolist(),
                                    "A2": traw_df["ALT"].tolist()
                                    })

        # add number of snps
        self.p = self.snp_df.shape[0]

        # add rownumbers
        self.snp_df["idx"] = np.arange(0, self.p, 1)

    def _get_inds(self, traw_df):
        """Get list of individual ids
        """
        # get ind list
        self.inds = list(map(lambda x: x.split("_")[0], traw_df.columns.tolist()[6:]))

        # add number of individuals
        self.n = len(self.inds)

    def _get_genotype_matrix(self, traw_df):
        """Gets the genotype matrix stored as a
        numpy array
        """
        # genotype matrix
        drop_cols = ["SNP", "CHR", "(C)M", "POS", "ALT", "COUNTED"] 
        self.y = traw_df.drop(drop_cols, axis=1).as_matrix()
        self.y = self.y.astype(np.float32)


class ReadSampledVcf(Genotypes):
    """Read in variant call format for read sampled dataset
    with the allele depth (AD) field. For example this could be
    ouput by https://github.com/mathii/gdc3/blob/master/apulldown.py

    Arguments
    ---------
    prefix : str
        prefix of  path to data files
    
    p : int
        number of snps (features)

    Attributes
    ----------
    prefix : str
        prefix of  path to data files

    y : np.array
        p x n genotype matrix 

    snp_df : pd.DataFrame
        Dataframe storing snp level meta data

    inds : list
        list of individual ids

    n : int
        number of individuals (samples)

    p : int
        number of snps (features)

    vcf_path : str
        path to read sampled vcf file
    
    base_complement : dict
        dictionary where keys are bases and values
        are the bases complement

    bases : list
        list of bases

    strand_ambiguous_alleles : dict
        dictionary of strand ambiguous alleles 
        combinations as keys and bools as values
    """
    def __init__(self, prefix, p):

        # inherit attributes from Genotypes
        super().__init__()

        # prefix to paths
        self.prefix = prefix

        self.p = p

        # path to traw file
        self.vcf_path = "{}.vcf.gz".format(self.prefix)

        # get genotype matrix, SNP DataFrame, and inviduals all at once!
        self._get_data()

    def _get_data(self):
        """Get a DataFrame with columns as individuals
        """
        # read vcf
        vcf = pysam.VariantFile(self.vcf_path)

        # list of sample ids
        self.inds = str(vcf.header).split("\n")[-2].split("\t")[9:]
        self.n = len(self.inds) # number of samples
        self.y = np.empty((self.p, self.n)) # genotype matrix

        # SNP DataFrame lists
        snps = []
        chroms = []
        poss = []
        a1s = []
        a2s = []

        # fill up genotypes
        for i,rec in enumerate(vcf.fetch()):

            # variant level data
            chroms.append(rec.chrom)
            poss.append(rec.pos)
            snps.append(rec.id)
            a1s.append(rec.ref[0])
            a2s.append(rec.alts[0])            

            # sample level data
            for j,sample in enumerate(self.inds):
                #ad = rec.samples[sample]["AD"]
                gt = rec.samples[sample]["GT"]
                if gt[0] != None and gt[1] != None:
                    self.y[i, j] = 2 - (gt[0] + gt[1])

                # if coverage is 0
                #if ad[0] == 0 and ad[1] == 0:
                #    self.y[i, j] = np.nan
                #else:
                #    q = ad[0] / (ad[0] + ad[1]) 
                #    self.y[i, j] = np.random.choice([0, 2], p=[1.-q, q])
                # if we observe at least single alt allele
                #elif ad[0] > 0:
                #    self.y[i, j] = 1
                # if we dont observe any alt alleles
                #else:
                #    self.y[i, j] == 0

        # create SNP DataFrame
        cm_poss = np.empty(len(snps))
        cm_poss[:] = np.nan
        self.snp_df = pd.DataFrame({"SNP": snps,
                                    "CHROM": chroms,
                                    "CM_POS": cm_poss,
                                    "POS": poss,
                                    "A1": a1s,
                                    "A2": a2s
                                   })

        # add number of snps
        self.p = self.snp_df.shape[0]

        # add rownumbers
        self.snp_df["idx"] = np.arange(0, self.p, 1)



