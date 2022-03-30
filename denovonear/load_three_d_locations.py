""" functions to load genes, and identify transcripts containing de novos.
"""

import asyncio

def load_three_d_locations(path_dir, gene):
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

def load_three_d_multimer(path_dir, multimer):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        list of Transcript objects for gene, including genomic ranges and sequences
    """
    #for gene in sort(genes):
    #    multimer += str(gene) + ","
    #multimer = multimer[:-1]
    path = path_dir + "/" + multimer+ ".multimer"
    #    three_d_locations = {}
    three_d_locations = []
    with open(path, "r") as handle:
        aa_num = 0
        for line in handle:
            line = line.rstrip().split("\t")
            chain = str(line[0])
            residue_number = int(line[1])
            x = float(line[2])
            y = float(line[3])
            z = float(line[4])
            #            three_d_locations[aa_num]= {x,y,z}
            three_d_locations.append([chain, residue_number,x,y,z])
            aa_num += 1
    return three_d_locations

d = {'CYS': 'C',
     'ASP': 'D',
     'SER': 'S',
     'GLN': 'Q',
     'LYS': 'K',
     'ILE': 'I',
     'PRO': 'P',
     'THR': 'T',
     'PHE': 'F',
     'ASN': 'N', 
     'GLY': 'G',
     'HIS': 'H',
     'LEU': 'L',
     'ARG': 'R',
     'TRP': 'W', 
     'ALA': 'A',
     'VAL':'V',
     'GLU': 'E',
     'TYR': 'Y',
     'MET': 'M'}

def load_three_d_locations_from_pdb(path_dir, gene):
    """ sort out all the necessary sequences and positions for a gene
    
    Args:
        ensembl: EnsemblRequest object to request data from ensembl
        gene_id: HGNC symbol for gene
        de_novos: list of de novo positions, so we can check they all fit in
            the gene transcript
        
    Returns:
        list of Transcript objects for gene, including genomic ranges and sequences
    """
    path = path_dir + "/" + gene + ".pdb"
    #    three_d_locations = {}
    three_d_locations = []
    with open(path, "r") as handle:
        for line in handle:
            list = line.rstrip().split()
            id = list[0]
            if id == 'ATOM':
                type = list[2]
                if type == 'CA':
                    residue = d[list[3]]
                    chain = list[4]
                    atom_count = int(list[5])
                    x = float(list[6])
                    y = float(list[7])
                    z = float(list[8])
                    three_d_locations.append([x,y,z,residue,atom_count,chain])

    return three_d_locations
