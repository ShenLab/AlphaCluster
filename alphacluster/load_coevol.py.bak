""" functions to load coevol strength
"""

import asyncio

async def load_basic_coevol_strength(path_dir, gene):
    """ Get the coevol strength
    
    Args:
    path_dir: where coevol are located
    gene:  HGNC symbol for gene
        
    Returns:
    array of coevol strength
    """
    path = path_dir + "/" + gene + ".coevol.txt"
    #    coevol = {}
    coevol = []
    with open(path, "r") as f:
        for line in f:
            coevol += [float(x) for x in line.split()]
    print(coevol[1])
    print(len(coevol))
    return coevol
