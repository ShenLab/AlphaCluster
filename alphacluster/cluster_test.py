from math import log, isnan, floor
import os
from scipy.stats import chi2
from difflib import SequenceMatcher

from alphacluster.load_inherited import load_inherited
from alphacluster.load_gene import load_gene, get_de_novos_in_transcript, minimise_transcripts
from alphacluster.load_mutation_rates import load_mutation_rates
from alphacluster.load_de_novos import load_de_novos
from alphacluster.site_specific_rates import SiteRates
#from alphacluster.genotype_rates import GenotypeRates
from alphacluster.simulate import get_p_value, get_p_value_1d, get_p_value_coevol, get_p_value_multi

def lcs_func(s1, s2):
    matrix = [["" for x in range(len(s2))] for x in range(len(s1))]
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                if i == 0 or j == 0:
                    matrix[i][j] = s1[i]
                else:
                    matrix[i][j] = matrix[i-1][j-1] + s1[i]
            else:
                matrix[i][j] = max(matrix[i-1][j], matrix[i][j-1], key=len)

    cs = matrix[-1][-1]

    return len(cs), cs

def create_windows(variants, window_size = 5000):
    windows = []
    for v in variants:
        start = v - window_size/2
        end = v + window_size/2
        if len([d for d in variants if d >= start and d <= end]) > 1:
            windows += [[start,end]]
    return windows

def fishers_method(values):
    """ function to combine p values, using Fisher's method
    
    We occasionally have multiple P values for a mutation type, obtained from
    alternative transcripts for the gene. If we have only one P value for the
    gene for the mutation type, we just use that value, if we don't have any
    data, we use "NA", otherwise combine p-values using Fisher's method.
    
    Args:
        x: list of P-values for a gene
    
    Returns:
        combined P-value
    """
    values = [ x for x in values if not isnan(x) ]
    # use Fisher's combined method to estimate the P value from multiple
    # P-values. The chi square statistic is -2*sum(ln(P-values))
    if len(values)>1:
        if min(values)>1E-20:
            return chi2.sf(-2 * sum(map(log, values)), 2 * len(values))
        else :
            return min(values)
    else:
        return values

def cluster_de_novos(symbol,
                     variants,
                     three_d_locations,
                     gene,
                     scale,
                     threshold = None,
                     p = 0,
                     dist_file_output = "",
                     iterations=1000000,
                     mut_dict=None,
                     scores_file = None,
                     score_col = 5,
                     dbNSFP_score_obj = None,
                     annotator =None,
                     pvalues_in = None,
                     variant_type = "de_novo",
                     inherited_controls = {}):

    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        symbol: HGNC symbol for a gene
        variants: dictionary of de novo positions for the HGNC gene,
        indexed by functional type
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    if variant_type == "de_novos" and  mut_dict is None:
        mut_dict = load_mutation_rates()

    missense = [pos for pos,alt in variants["missense"]]
    nonsense = [pos for pos,alt in variants["nonsense"]]

    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = gene.transcripts
        #minimized = minimise_transcripts(transcripts, missense + nonsense)
        #transcripts = [x for x in transcripts if x.get_name() in minimized]
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    counts = {"miss_count":[], "nons_count": []}
    results = {"missense":[], "nonsense": []}
    iteration = 0
    transcripts = [gene.canonical]

    for transcript in transcripts:
        try:
            t_aa = "".join(list(transcript.translate(transcript.get_cds_sequence())))
            p_aa = "".join([row[3] for row in three_d_locations])

            #Find longest common substring
            lcs_len,lcs = lcs_func(t_aa,p_aa)
            start_offset = t_aa.find(lcs)
            end_offset = t_aa[::-1].find(lcs[::-1])-1

            if(start_offset != 0):
                transcript.set_cds_start(transcript.get_cds_start(),start_offset*3)
            if(end_offset != 0):
                transcript.set_cds_end(transcript.get_cds_end(),end_offset*-3)                
        except:
            print("AlphaFold and transcript disagree for transcript " + transcript.get_name())
            continue
            
        iteration += 1
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        
        scores = {}
        if dbNSFP_score_obj is not None:
            dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
            if annotator in dbNSFP_header:
                anno_idx = dbNSFP_header.index(annotator)
                alt_idx = dbNSFP_header.index('alt')
                pos_idx = dbNSFP_header.index('pos(1-based)')
                for line in dbNSFP_score_obj.fetch(transcript.get_chrom()[3:],
                                                   transcript.get_start()-1,
                                                   transcript.get_end()):
                    values = line.split('\t')
                    pos = values[pos_idx]
                    alt = values[alt_idx]
                    if values[anno_idx] == '.' :
                        continue
                    score = float(values[anno_idx])
                    if int(pos) not in scores:
                        scores[int(pos)] = {}
                    scores[int(pos)][ord(alt)] = float(score)
                    
            try:
                missense_scores = [scores[int(pos)][ord(alt)] for pos,alt in variants["missense"] if pos in missense_events and int(pos) in scores and ord(alt) in scores[int(pos)]]
            except KeyError:
                print("No score exists for this position")
                return(None)
            except:
                print("Error while loading scores for missense positions")
                return(None)
        elif scores_file is not None:
            for line in scores_file.fetch(transcript.get_chrom()[3:],
                                          transcript.get_start()-1,
                                          transcript.get_end()):
                #_, pos, _, alt, _, score = line.split('\t')
                pos = line.split('\t')[1]
                alt = line.split('\t')[3]
                #score = line.split('\t')[5]
                score = line.split('\t')[int(score_col)]                
                if int(pos) not in scores:
                    scores[int(pos)] = {}
                scores[int(pos)][ord(alt)] = float(score)
            missense_scores = []
            for pos,alt in variants["missense"]:
                if pos in missense_events and int(pos) in scores:
                    if ord(alt) in scores[int(pos)]:
                        missense_scores.append(scores[int(pos)][ord(alt)])
        else:
            missense_scores = [-1 for x in missense_events]

        # Handling thresholding of distribution
        if scale:
            rates = SiteRates(transcript, mut_dict, scores)
        else:
            rates = SiteRates(transcript, mut_dict, {})

        # Move to no longer need transcript in next layer
        miss_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in missense_events ]
        nons_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]        

        miss_cds_positions = [ x for x in miss_cds_positions if x>=0]
        nons_cds_positions = [ x for x in nons_cds_positions if x>=0]        

        missense_dist_file_output=""
        nonsense_dist_file_output=""
        if(dist_file_output != None):
            missense_dist_file_output = dist_file_output + symbol + ".missense.3d.dist"
            nonsense_dist_file_output = dist_file_output + symbol + ".nonsense.3d.dist" + dist_file_output

        (miss_dist, miss_prob, miss_count) = get_p_value(rates,
                                             three_d_locations,
                                             iterations,
                                             "missense",
                                             miss_cds_positions,
                                             missense_scores,
                                             scale,
                                             threshold,
                                             p,
                                             missense_dist_file_output.encode('utf-8'))


        dists["miss_dist"].append(miss_dist)
        counts["miss_count"].append(miss_count)
        probs["miss_prob"].append(miss_prob)
        results["missense"] += [[transcript.get_start(), transcript.get_end(), len(missense_events), miss_dist, miss_prob]]

        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])

    if pvalues_in != None:
        dists["miss_dist_pois"] = dists["miss_dist"]
        dists["miss_dist_comb"] = dists["miss_dist"]

        probs["miss_prob_pois"] = [pvalues_in[symbol]]
        
        probs["miss_prob_comb"] = probs["miss_prob"].copy()
        probs["miss_prob_comb"].append(pvalues_in[symbol])
        
    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    probs.update(counts)
    return probs, results


def cluster_de_novos_multi(multimer_name,
                           all_chains,
                           proteins,
                           proteins_to_chains,
                           variants,
                           three_d_locations,
                           genes,#ensembl,
                           scale,
                           threshold = None,
                           p = 0,
                           dist_file_output = "",
                           iterations=1000000,
                           mut_dict=None,
                           scores_file = None,
                           score_col = 5,
                           pvalues_in = None):
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        symbol: HGNC symbol for a gene
        variants: dictionary of de novo positions for the HGNC gene,
        indexed by functional type
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    results = {"missense":[], "nonsense": []}

    # for rates we need 1) transcripts 2) mut rates 3) scores 4)sites
    # so we have to load each of these
    rates = []
    all_missense = []
    
    # load mutation background rates
    if mut_dict is None:
        mut_dict = load_mutation_rates()

    # Get list of valid amino acids
    residues = {}
    for x,y,z,residue,residue_number,chain in three_d_locations:
        if chain not in residues:
            residues[chain] = set()
        residues[chain].add(int(residue_number))
        
    # Make three_d_locations into list of dictionary
    three_d_list_of_dict = [{}]*len(all_chains)
    for x,y,z,residue,residue_number,chain in three_d_locations:
        if chain in all_chains:
            chain_number = all_chains.index(chain)
            three_d_list_of_dict[int(chain_number)][int(residue_number)] = [float(x),float(y),float(z)]

    missense_cds_positions = {}
    all_missense_locations = []
    all_missense_scores = []
    de_novo_count_per_chain = [None]*len(all_chains)
    
    for protein, chains in proteins_to_chains.items():
        missense = [int(pos) for pos,alt in variants[protein]["missense"]]
    
        try:
            #Currently use the canonical transcript
            transcript = genes[protein].canonical
            print("BEFORE")
            print(transcript.get_cds_start())
            print(transcript.get_cds_end())                       
        except IndexError as e:
            print(e)
            return None

        # Collect residues in the PDB model
        these_residues = set(residues[chains[0]])
        for chain in chains:
            these_residues = set.intersection(residues[chain], these_residues)
        
        # Confirm that transcript and protein align well
        t_aa = "".join(list(transcript.translate(transcript.get_cds_sequence())))
        p_aa = "".join([row[3] for row in three_d_locations if row[5] == chain])
            
        #Find longest common substring
        lcs_len,lcs = lcs_func(t_aa,p_aa)
            
        #find offsets
        start_offset = t_aa.find(lcs)
        end_offset = t_aa[::-1].find(lcs[::-1])-1

        start_codon = start_offset
        end_codon = len(t_aa)-end_offset-2

        # Load missense variants for this protein
        missense_records = {}
        missense_records_map = {}
        missense_records[protein] = [(pos, int(transcript.get_coding_distance(pos)['pos']), alt) for pos,alt in variants[protein]["missense"]]

        # Trim down to missense variants in the PDB model
        missense_records_map[protein] = {pos:alt for pos, cds_pos, alt in missense_records[protein] if floor(float(cds_pos)/3) in these_residues}

        # Import scores
        scores = {}
        scores[protein] = {}
        if scores_file is not None:
            if scores[protein] == {}:
                for line in scores_file.fetch(transcript.get_chrom()[3:],
                                              transcript.get_start()-1,
                                              transcript.get_end()):
                    pos = line.split('\t')[1]
                    alt = line.split('\t')[3]
                    score = line.split('\t')[int(score_col)]                
                    if int(pos) not in scores[protein]:
                        scores[protein][int(pos)] = {}
                    scores[protein][int(pos)][ord(alt)] = float(score)

        # Remove missense variants below threshold
        if threshold is not None:
            missense_records[protein] = [(pos, cds_pos,alt) for pos, cds_pos,alt in missense_records[protein] if floor(float(cds_pos)/3) in these_residues and float(scores[protein][int(pos)][ord(alt)]) > float(threshold)]
            missense_records_map[protein] = {pos:alt for pos, cds_pos, alt in missense_records[protein] }
            
        # Mark how many de novos to expect on chains corresponding to this protein
        for chain in chains:
            de_novo_count_per_chain[all_chains.index(chain)] = len(missense_records[protein])

        # Collect the missense for this protein
        for chain in chains:
            #for cds_pos in missense_cds_positions[protein]:
            for pos,cds_pos,alt in missense_records[protein]:
                # Collect positions
                print(cds_pos)
                print(floor(float(cds_pos/3))+1)
                print(three_d_list_of_dict[all_chains.index(chain)][floor(float(cds_pos/3))+1])
                all_missense_locations.append(three_d_list_of_dict[all_chains.index(chain)][floor(float(cds_pos/3))+1])

                # Collect scores
                if scores[protein] != {}:
                    all_missense_scores.append(scores[protein][pos][ord(alt)])
                    #alt = missense_records_map[protein][transcript.get_position_on_chrom(cds_pos,0)]
                    #all_missense_scores.append(scores[protein][transcript.get_position_on_chrom(cds_pos,0)][ord(alt)])
                else:
                    all_missense_scores.append(-1)

        #Create and append the rate object for protein
        if len(missense_records[protein]) > 0:
            rates.append({"protein" : protein,
                              "rate" : SiteRates(transcript,
                                                 mut_dict,
                                                 scores[protein],
                                                 float(-2),
                                                 residues = list(these_residues))})


    missense_dist_file_output=""
    nonsense_dist_file_output=""
    if(dist_file_output != None):
        missense_dist_file_output = dist_file_output + multimer_name+ ".missense.3d.dist"
        nonsense_dist_file_output = dist_file_output + multimer_name + ".nonsense.3d.dist" 


    (miss_dist, miss_prob) = get_p_value_multi(all_chains,
                                               proteins,
                                               de_novo_count_per_chain,
                                               rates,
                                               three_d_list_of_dict,
                                               iterations,
                                               "missense",
                                               all_missense_locations,
                                               all_missense_scores,
                                               p,
                                               missense_dist_file_output.encode('utf-8'))

    dists["miss_dist"].append(miss_dist)
    probs["miss_prob"].append(miss_prob)

    results["missense"] += [[transcript.get_start(), transcript.get_end(), len(all_missense_locations), miss_dist, miss_prob]]

    # remove the de novos analysed in the current transcript, so that
    # analysis of subsequent transcripts uses independent events. NOTE THAT
    # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
    # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
    #missense = [x for x in missense if x not in missense_events]
    #nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])

    if pvalues_in != None:
        dists["miss_dist_comb"] = dists["miss_dist"]
        probs["miss_prob_comb"] = probs["miss_prob"].copy()
        dists["miss_dist_pois"] = dists["miss_dist"]
        probs["miss_prob_pois"] = []
        for protein in proteins:
            probs["miss_prob_comb"].append(pvalues_in[protein])
            probs["miss_prob_pois"].append(pvalues_in[protein])
        
    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    return probs, results



def cluster_de_novos_coevol(symbol, variants, coevol, gene, p = 0, dist_file_output = "", iterations=1000000, mut_dict=None,scores = {}, pvalues_in=None,
                           variant_type = "de_novo"):
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        symbol: HGNC symbol for a gene
        variants: dictionary of de novo positions for the HGNC gene,
        indexed by functional type
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    if variant_type == "de_novos" and  mut_dict is None:
        mut_dict = load_mutation_rates()

    missense = [pos for pos,alt in variants["missense"]]
    nonsense = [pos for pos,alt in variants["nonsense"]]

    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = gene.transcripts
        minimized = minimise_transcripts(transcripts, missense + nonsense)
        transcripts = [x for x in transcripts if x.get_name() in minimized]
        #transcripts = await load_gene(ensembl, symbol, missense + nonsense)
    except IndexError as e:
        print(e)
        return None

    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    #try:
    #    transcripts = await load_gene(ensembl, symbol, missense + nonsense)
    #except IndexError as e:
    #    print(e)
    #    return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    iteration = 0
    for transcript in transcripts:
        iteration+=1
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)

        if variant_type == "de_novo":
            rates = SiteRates(transcript, mut_dict, {})
        else:
            return

        missense_dist_file_output=""
        nonsense_dist_file_output=""
        if(dist_file_output != None):
            missense_dist_file_output = dist_file_output + symbol + ".missense.3d.dist"
            nonsense_dist_file_output = dist_file_output + symbol + ".nonsense.3d.dist"

        (miss_dist, miss_prob) = get_p_value_coevol(transcript, rates, coevol, iterations, "missense", missense_events, p, missense_dist_file_output.encode('utf-8'))
#        if iteration == 1:
#            os.rename('current_distribution.txt', symbol + ".missense.coev.dist.txt")
        (nons_dist, nons_prob) = get_p_value_coevol(transcript, rates, coevol, iterations, "lof", nonsense_events, p, nonsense_dist_file_output.encode('utf-8'))
#        if iteration == 1:
#            os.rename('current_distribution.txt', symbol + ".nonsense.coevol.dist.txt")
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        
        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])

    if pvalues_in != None:
        dists["miss_dist_comb"] = dists["miss_dist"]
        probs["miss_prob_comb"] = probs["miss_prob"].copy()
        probs["miss_prob_comb"].append(pvalues_in[symbol])

        probs["miss_prob_pois"] = [pvalues_in[symbol]]
        dists["miss_dist_pois"] = dists["miss_dist"]


    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    
    return probs




def cluster_de_novos_1d(symbol,
                        variants,
                        gene,
                        scale,
                        threshold,
                        p = 0,
                        dist_file_output = "",
                        iterations=1000000,
                        mut_dict=None,
                        scores_file = None,
                        score_col = 5,
                        pvalues_in = None):
    """ analysis proximity cluster of de novos in a single gene
    
    Args:
        symbol: HGNC symbol for a gene
        de_novos: dictionary of de novo positions for the HGNC gene,
        indexed by functional type
        iterations: number of simulations to run
        ensembl: EnsemblRequest object, for obtaing info from ensembl
        mut_dict: dictionary of mutation rates, indexed by trinuclotide sequence
    
    Returns:
        a dictionary containing P values, and distances for missense, nonsense,
        and synonymous de novos events. Missing data is represented by "NA".
    """
    
    if mut_dict is None:
        mut_dict = load_mutation_rates()
    
    missense = [pos for pos,alt in variants["missense"]]
    nonsense = [pos for pos,alt in variants["nonsense"]]
    
    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = gene.transcripts
        minimized = minimise_transcripts(transcripts, missense + nonsense)
        transcripts = [x for x in transcripts if x.get_name() in minimized]
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}

    iteration = 0
    for transcript in transcripts:
        iteration += 1
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)

        scores = {}
        if scores_file is not None:
            for line in scores_file.fetch(transcript.get_chrom()[3:],
                                          transcript.get_start()-1,
                                          transcript.get_end()):
                #_, pos, _, alt, _, score = line.split('\t')
                pos = line.split('\t')[1]
                alt = line.split('\t')[3]
                score = line.split('\t')[int(score_col)]                
                if int(pos) not in scores:
                    scores[int(pos)] = {}
                scores[int(pos)][ord(alt)] = float(score)                
            missense_scores = [scores[int(pos)][ord(alt)] for pos,alt in variants["missense"]]
            nonsense_scores = [1 for x in nonsense_events]
        else:
            missense_scores = [1 for x in missense_events]
            nonsense_scores = [1 for x in nonsense_events]

        # Handling thresholding of distribution
        if scale:
            rates = SiteRates(transcript, mut_dict, scores)
        else:
            rates = SiteRates(transcript, mut_dict, {})
        #if threshold == None:
        #    rates = SiteRates(transcript, mut_dict, scores)
        #else:
        #    rates = SiteRates(transcript, mut_dict, scores, float(threshold))

        #rates = SiteRates(transcript, mut_dict, scores)

        missense_dist_file_output=""
        nonsense_dist_file_output=""
        if(dist_file_output != None):
            missense_dist_file_output = dist_file_output + symbol + ".missense.3d.dist"
            nonsense_dist_file_output = dist_file_output + symbol + ".nonsense.3d.dist"
        
        (miss_dist, miss_prob) = get_p_value_1d(transcript,
                                                rates,
                                                iterations,
                                                "missense",
                                                missense_events,
                                                missense_scores,
                                                scale,
                                                threshold,
                                                p,
                                                missense_dist_file_output.encode('utf-8'))
#        if iteration == 1:
#            os.rename('current_distribution.txt', symbol + ".missense.1d.dist.txt")
        (nons_dist, nons_prob) = get_p_value_1d(transcript,
                                                rates,
                                                iterations,
                                                "lof",
                                                nonsense_events,
                                                nonsense_scores,
                                                p,
                                                nonsense_dist_file_output.encode('utf-8'))
        
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        
        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])

    if pvalues_in != None:
        dists["miss_dist_comb"] = dists["miss_dist"]
        probs["miss_prob_comb"] = probs["miss_prob"].copy()
        probs["miss_prob_comb"].append(pvalues_in[symbol])
        probs["miss_prob_pois"] = [pvalues_in[symbol]]
        dists["miss_dist_pois"] = dists["miss_dist"]
    
    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    
    return probs
