""" loads a file containing known de novo mutations
"""

from __future__ import print_function, division

missense = ["missense_variant", "missense"]
#missense = ["missense_variant", "stop_lost", "inframe_deletion",
#    "inframe_insertion", "coding_sequence_variant", "protein_altering_variant"]

lof = ["stop_gained", "splice_acceptor_variant",
    "splice_donor_variant", "frameshift_variant", "initiator_codon_variant",
    "start_lost", "conserved_exon_terminus_variant"]
    
synonymous = ["synonymous_variant"]

def load_de_novos(path, exclude_indels=True, inherited_variants=False):
    """ load mutations into dict indexed by HGNC ID.
    
    Args:
        path: path to file containing de novo data. This should have five tab-
            separated columns e.g.
            hgnc  chr  pos      consequence         var_type    aa_num
            CTC1  17   8139190  missense_variant    snv         3
            CTC1  17   8139191  frameshift_variant  indel       6
        exclude_indels: True/False for whether we want to exclude indels. If we
            are testing clustering of de novos, we can only use SNVs, but if we
            are determining mutation rates for a gene, then we include indels,
            in order to identify the transcripts that the de novos lie within.
    
    Returns:
        dictionary of missense and lof counts for each gene, indexed by HGNC symbols
    """
    
    genes = {}
    with open(path, "r") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            
            line = line.rstrip().split("\t")
            gene = line[0]
            position = int(line[2])
            #position = int(line[2]) - 1
            ref = line[3]
            alt = str(line[4])
            consequence = line[5]
            var_type = line[6]
            snp_of_indel = line[7]
            #score = line[8]
            #consequence = line[3]
            #var_type = line[4]

            if inherited_variants and not "inherited" in var_type.lower():
                continue
                
            # ignore indels (some splice_acceptor_variants (in the
            # functional_consequences) are indels
            if exclude_indels and "indel" in var_type.lower():
                continue
            
            # trim out variants that are missing data
            if gene == "" or gene == "." or position == "NA":
                continue
            
            if gene not in genes:
                genes[gene] = {"missense": [], "nonsense": []}
            
            if consequence in missense:
                genes[gene]["missense"].append((position,alt))
            elif consequence in lof:
                genes[gene]["nonsense"].append((position,alt))
        
    return genes
        

def load_de_novos_chrom_pos_alt(path, exclude_indels=True, inherited_variants=False):
    """ load mutations into dict indexed by HGNC ID.
    
    Args:
        path: path to file containing de novo data. This should have five tab-
            separated columns e.g.
            hgnc  chr  pos      consequence         var_type    aa_num
            CTC1  17   8139190  missense_variant    snv         3
            CTC1  17   8139191  frameshift_variant  indel       6
        exclude_indels: True/False for whether we want to exclude indels. If we
            are testing clustering of de novos, we can only use SNVs, but if we
            are determining mutation rates for a gene, then we include indels,
            in order to identify the transcripts that the de novos lie within.
    
    Returns:
        dictionary of missense and lof counts for each gene, indexed by HGNC symbols
    """
    
    genes = {}
    with open(path, "r") as handle:
        header = handle.readline().strip().split("\t")
        for line in handle:
            
            line = line.rstrip().split("\t")
            gene = line[0]
            chrom = line[1]
            position = int(line[2])
            #position = int(line[2]) - 1
            ref = line[3]
            alt = str(line[4])
            consequence = line[5]
            var_type = line[6]
            snp_of_indel = line[7]
            #score = line[8]
            #consequence = line[3]
            #var_type = line[4]

            if inherited_variants and not "inherited" in var_type.lower():
                continue
                
            # ignore indels (some splice_acceptor_variants (in the
            # functional_consequences) are indels
            if exclude_indels and "indel" in var_type.lower():
                continue
            
            # trim out variants that are missing data
            if gene == "" or gene == "." or position == "NA":
                continue
            
            if gene not in genes:
                genes[gene] = {"missense": [], "nonsense": []}
            
            if consequence in missense:
                genes[gene]["missense"].append((chrom, position,alt))
            elif consequence in lof:
                genes[gene]["nonsense"].append((chrom, position,alt))
        
    return genes
        

        
    
    
    
    
