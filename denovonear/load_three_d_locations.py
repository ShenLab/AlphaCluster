""" functions to load genes, and identify transcripts containing de novos.
"""

import asyncio

async def load_three_d_locations(path_dir, gene):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        list of Transcript objects for gene, including genomic ranges and sequences
    """
    path = path_dir + "/" + gene + ".xyz"
    #    three_d_locations = {}
    three_d_locations = []
    with open(path, "r") as handle:
        aa_num = 0
        for line in handle:
            line = line.rstrip().split("\t")
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])
            #            three_d_locations[aa_num]= {x,y,z}
            three_d_locations.append([x,y,z])
            aa_num += 1
    return three_d_locations
