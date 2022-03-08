""" loads a file containing pvalues from some
independent test
"""

from __future__ import print_function, division

def load_pvalues(path):
    """ load mutations into dict indexed by HGNC ID.
    
    Args:
        path: path to file containing de novo data. This should have five tab-
            separated columns e.g.
    GENE   PVALUE
    Returns:
        dictionary of gene and pvalue pairs
    """
    
    pvalues = {}
    with open(path, "r") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            line = line.rstrip().split("\t")
            gene = line[0]
            pvalue = line[1]

            # trim out variants that are missing data
            if gene == "" or gene == "." or gene == "NA" or pvalue == "." or pvalue == "" or pvalue=="NA":
                continue
            pvalues[gene] = float(pvalue)
    return pvalues
        
        
        
    
    
    
    
