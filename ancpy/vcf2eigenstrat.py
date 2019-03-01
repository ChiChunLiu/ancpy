#!/usr/bin/python3

from pandas_plink import read_plink
import argparse
import pandas as pd
import numpy as np
import utils as utl




if __name__== "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--input', type = str, default = "", help = "prefix for plink files")
    parser.add_argument('-o', '--output', type = str, default = "", help = "prefix for outputs")
    args = parser.parse_args()
    
    while True:
        line = a.readline()
        if not line: 
            break


    
    # write files
    out = utl.eigenstrat(geno = gt,
                         snp = snp,
                         ind = sample)
    out.write_eigenstrat(args.output)
