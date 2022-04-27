""" functions to load genes, and identify transcripts containing de novos.
"""

import asyncio
import gzip


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
    path = path_dir + "/" + gene + ".pdb.gz"
    #    three_d_locations = {}
    three_d_locations = []
    with gzip.open(path, "r") as handle:
        for line in handle:
            list = line.rstrip().decode("utf-8").split()
            line = line.decode("utf-8")
            #id = list[0]
            id = line[0:4].strip()
            #print(id)
            if id == 'ATOM':
                #type = list[2]
                type = line[12:16].strip()
                if type == 'CA':
                    #print(type)
                    #residue = d[list[3]]
                    #chain = list[4]
                    #atom_count = int(list[5])
                    #x = float(list[6])
                    #y = float(list[7])
                    #z = float(list[8])
                    residue = d[line[17:20].strip()]
                    chain = line[20:22].strip()
                    atom_count = line[22:27].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())

                    three_d_locations.append([x,y,z,residue,atom_count,chain])

    return three_d_locations
