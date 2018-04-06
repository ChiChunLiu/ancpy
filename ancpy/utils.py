from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import pandas as pd
import numpy as np


class Merger(object)
    """
    """
    def __init__(self, genotypes_0, genotypes_1):
        
        self.genotypes_0 = genotypes_0
        self.genotypes_1 = genotypes_1
        
    """
    intersect snp dataframe to find overlapping snps
    index genotpye matrix 0 and 1 by overlapping snps
    find snps that need to be flipped

    """
        