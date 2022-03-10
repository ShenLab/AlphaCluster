from math import log, isnan, floor
import os
from scipy.stats import chi2
from load_inherited import load_inherited
from load_gene import load_gene, get_de_novos_in_transcript, best_transcript, minimise_transcripts_2
from denovonear.load_mutation_rates import load_mutation_rates
from denovonear.load_de_novos import load_de_novos
from denovonear.site_specific_rates import SiteRates
#from denovonear.genotype_rates import GenotypeRates
from simulate import get_p_value, get_p_value_1d, get_p_value_coevol, get_p_value_multi, get_p_value_entropy

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
    if len(values)>0:
        if min(values)>1E-20:
            return chi2.sf(-2 * sum(map(log, values)), 2 * len(values))
        else :
            return min(values)
    else:
        return values

async def cluster_analysis(symbol,
                           de_novos,
                           inherited,
                           inherited_controls,
                           three_d_locations,
                           coevolution,
                           ensembl,
                           degree = 0,
                           iterations=1000000,
                           mut_dict=None,
                           scores_file = None,
                           variant_type = ["de_novo"]):
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
    de_novo = False
    missense = []
    nonsense = []
    inherited_missense = []
    inherited_nonsense = []
    inherited = False
    # De novos were loaded
    if len(de_novos) > 0:
        de_novo = True
        if mut_dict is None:
            mut_dict = load_mutation_rates()
        missense = [pos for pos,alt in de_novos["missense"]]
        nonsense = [pos for pos,alt in de_novos["nonsense"]]

    # Inherited were loaded
    if len(inherited) > 0:
        inherited = True
        inherited_missense = [pos for pos,alt in inherited["missense"]]
        inherited_nonsense = [pos for pos,alt in inherited["nonsense"]]        
    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = await load_gene(ensembl, symbol, missense + nonsense + inherited_missense + inherited_nonsense)
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    results = {"missense":[], "nonsense": []}

    iteration = 0
    for transcript in transcripts:
        iteration += 1
        # Get variants location in transcript
        if de_novo:
            missense_events = get_de_novos_in_transcript(transcript, missense)
            nonsense_events = get_de_novos_in_transcript(transcript, nonsense)

        if inherited:
            inherited_missense_events = get_de_novos_in_trnascript(transcript, inherited_missense)
            inherited_nonsense_events = get_de_novos_in_trnascript(transcript, inherited_nonsense)            
            inherited_controls["missense"] = get_de_novos_in_transcript(transcript,
                                                                        inherited_controls["missense"])
            inherited_controls["missense"] = [ transcript.get_coding_distance(x)['pos'] for x in inherited_controls["missense"]]
            rates.set_inherited_controls(inherited_controls["missense"],
                                        "missense")
            inherited_controls["nonsense"] = get_de_novos_in_transcript(transcript,
                                                                        inherited_controls["nonsense"])
            inherited_controls["nonsense"] = [ transcript.get_coding_distance(x)['pos'] for x in inherited_controls["nonsense"]]
            rates.set_inherited_controls(inherited_controls["nonsense"],
                                        "nonsense")

        #         mode = "k-clusters"
        #         mode = ""
        #         if mode == "k-clusters":
        #             windows = create_windows(missense_events)
        # #            windows = list(set(windows))
        #             for w in windows:
        #                 dists["miss_dist_" + str(w[0]) + "-" + str(w[1])] = []
        #                 probs["miss_prob_" + str(w[0]) + "-" + str(w[1])] = []                
        #             for w in windows:
        #                 print("NEW WINDOW " + str(w[0]) + " to " + str(w[1]))
        #                 window_variants = [ x for x in missense_events if (x >= w[0] and x <= w[1])]
        #                 print(window_variants)
        #                 rates = SiteRates(transcript, mut_dict, w[0], w[1])
        #                 (miss_dist, miss_prob) = get_p_value_1d(transcript,
        #                                                         rates,
        #                                                         #three_d_locations,
        #                                                         iterations,
        #                                                         "missense",
        #                                                         window_variants,
        #                                                         p)
        #                 dists["miss_dist_" + str(w[0]) + "-" + str(w[1])].append(miss_dist)
        #                 probs["miss_prob_" + str(w[0]) + "-" + str(w[1])].append(miss_prob)
        #                 results["missense"] += [[w[0], w[1], len(window_variants), miss_dist, miss_prob]]
        #             windows = create_windows(nonsense_events)
        # #            windows = list(set(windows))            
        #             for w in windows:
        #                 dists["nons_dist_" + str(w[0]) + "-" + str(w[1])] = []
        #                 probs["nons_prob_" + str(w[0]) + "-" + str(w[1])] = []                
        #             for w in windows:
        #                 print("NEW NONSENSE WINDOW " + str(w[0]) + " to " + str(w[1]))
        #                 window_variants = [ x for x in nonsense_events if (x >= w[0] and x <= w[1])]
        #                 print(window_variants)                
        #                 rates = SiteRates(transcript, mut_dict, w[0], w[1])
        #                 (nons_dist, nons_prob) = get_p_value_1d(transcript,
        #                                                         rates,
        #                                                         #three_d_locations,
        #                                                         iterations,
        #                                                         "lof",
        #                                                         window_variants,
        #                                                         p)
        #                 dists["nons_dist_" + str(w[0]) + "-" + str(w[1])].append(nons_dist)
        #                 probs["nons_prob_" + str(w[0]) + "-" + str(w[1])].append(nons_prob)
        #                 results["nonsense"] += [[w[0], w[1], len(window_variants), nons_dist, nons_prob]]


        rates = SiteRates(transcript, mut_dict, scores)

        if len(scores) > 1:
            missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
            nonsense_scores = [1 for x in nonsense_events]            
        else:
            missense_scores = [1 for x in missense_events]
            nonsense_scores = [1 for x in nonsense_events]            
        (miss_dist, miss_prob) = get_p_value(transcript,
                                             rates,
                                             three_d_locations,
                                             iterations,
                                             "missense",
                                             missense_events,
                                             missense_scores,
                                             p)
        (nons_dist, nons_prob) = get_p_value(transcript,
                                             SiteRates(transcript,mut_dict, {}),
                                             three_d_locations,
                                             iterations,
                                             "lof",
                                             nonsense_events,
                                             nonsense_scores,
                                             p)
        
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        results["missense"] += [[transcript.get_start(), transcript.get_end(), len(missense), miss_dist, miss_prob]]
        results["nonsense"] += [[transcript.get_start(), transcript.get_end(), len(nonsense), nons_dist, nons_prob]]        
        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])
    
    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    return probs, results


async def cluster_de_novos(symbol, variants, three_d_locations, gene, threshold = None, p = 0, dist_file_output = "", iterations=1000000, mut_dict=None, scores_file = None, dbNSFP_score_obj = None, annotator =None, pvalues_in = None,variant_type = "de_novo", inherited_controls = {}):
#async def cluster_de_novos(symbol, variants, three_d_locations, ensembl, p = 0, iterations=1000000, mut_dict=None, scores_file = None, pvalues_in = None,variant_type = "de_novo", inherited_controls = {}):    
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
        minimized = minimise_transcripts_2(transcripts, missense + nonsense)
        transcripts = [x for x in transcripts if x.get_name() in minimized]
        #transcripts = await load_gene(ensembl, symbol, missense + nonsense)
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    results = {"missense":[], "nonsense": []}
    iteration = 0
    for transcript in transcripts:
        iteration += 1
        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)

        if variant_type == "inherited":
            inherited_controls["missense"] = get_de_novos_in_transcript(transcript,
                                                                    inherited_controls["missense"])
            inherited_controls["missense"] = [ transcript.get_coding_distance(x)['pos'] for x in inherited_controls["missense"]]
            rates.set_inherited_controls(inherited_controls["missense"],
                                         "missense")
            # NOTE: other possible approaches
            # rates = SiteRates(transcript, mut_dict, inherited_controls["missense"])
            # rates["missense"].set_inherited_choices(inherited_controls["missense"])

            # print(inherited_controls["missense"])
        mode = "k-clusters"
        mode = ""
        if mode == "k-clusters":
            windows = create_windows(missense_events)
#            windows = list(set(windows))
            for w in windows:
                dists["miss_dist_" + str(w[0]) + "-" + str(w[1])] = []
                probs["miss_prob_" + str(w[0]) + "-" + str(w[1])] = []                
            for w in windows:
                print("NEW WINDOW " + str(w[0]) + " to " + str(w[1]))
                window_variants = [ x for x in missense_events if (x >= w[0] and x <= w[1])]
                print(window_variants)
                rates = SiteRates(transcript, mut_dict, w[0], w[1])
                (miss_dist, miss_prob) = get_p_value_1d(transcript,
                                                        rates,
                                                        #three_d_locations,
                                                        iterations,
                                                        "missense",
                                                        window_variants,
                                                        p)
                dists["miss_dist_" + str(w[0]) + "-" + str(w[1])].append(miss_dist)
                probs["miss_prob_" + str(w[0]) + "-" + str(w[1])].append(miss_prob)
                results["missense"] += [[w[0], w[1], len(window_variants), miss_dist, miss_prob]]
            windows = create_windows(nonsense_events)
#            windows = list(set(windows))            
            for w in windows:
                dists["nons_dist_" + str(w[0]) + "-" + str(w[1])] = []
                probs["nons_prob_" + str(w[0]) + "-" + str(w[1])] = []                
            for w in windows:
                print("NEW NONSENSE WINDOW " + str(w[0]) + " to " + str(w[1]))
                window_variants = [ x for x in nonsense_events if (x >= w[0] and x <= w[1])]
                print(window_variants)                
                rates = SiteRates(transcript, mut_dict, w[0], w[1])
                (nons_dist, nons_prob) = get_p_value_1d(transcript,
                                                        rates,
                                                        #three_d_locations,
                                                        iterations,
                                                        "lof",
                                                        window_variants,
                                                        p)


                dists["nons_dist_" + str(w[0]) + "-" + str(w[1])].append(nons_dist)
                probs["nons_prob_" + str(w[0]) + "-" + str(w[1])].append(nons_prob)
                results["nonsense"] += [[w[0], w[1], len(window_variants), nons_dist, nons_prob]]



        scores = {}
        if dbNSFP_score_obj is not None:
            dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
            if annotator in dbNSFP_header:
                anno_idx = dbNSFP_header.index(annotator)
                alt_idx = dbNSFP_header.index('alt')
                pos_idx = dbNSFP_header.index('pos(1-based)')
                for line in dbNSFP_score_obj.fetch(transcript.get_chrom(),
                                                   transcript.get_start()-1,
                                                   transcript.get_end()):
                    values = line.split('\t')
                    pos = values[pos_idx]
                    alt = values[alt_idx]
                    if values[anno_idx] == '.' :
                        #print(pos)
                        #print(alt)
                        continue
                    score = float(values[anno_idx])
                    if int(pos) not in scores:
                        scores[int(pos)] = {}
                    scores[int(pos)][alt.encode('utf-8')] = float(score)
                    
            missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
            nonsense_scores = [-1 for x in nonsense_events]                        
        elif scores_file is not None:
            for line in scores_file.fetch(transcript.get_chrom(),
                                          transcript.get_start()-1,
                                          transcript.get_end()):
                #_, pos, _, alt, _, score = line.split('\t')
                pos = line.split('\t')[1]
                alt = line.split('\t')[3]
                score = line.split('\t')[13]                
                if int(pos) not in scores:
                    scores[int(pos)] = {}
                scores[int(pos)][alt.encode('utf-8')] = float(score)
            missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
            nonsense_scores = [-1 for x in nonsense_events]            
        else:
            missense_scores = [-1 for x in missense_events]
            nonsense_scores = [-1 for x in nonsense_events]
        print(transcript.get_name())

        scores = {}

        # Handling thresholding of distribution
        rates = SiteRates(transcript, mut_dict, scores)
        #if threshold == None:
        #    rates = SiteRates(transcript, mut_dict, scores)
        #else:
        #    rates = SiteRates(transcript, mut_dict, scores, float(threshold))

        # Move to no longer need transcript in next layer
        miss_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in missense_events ]
        nons_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]        

        missense_dist_file_output=""
        nonsense_dist_file_output=""
        if(dist_file_output != None):
            missense_dist_file_output = dist_file_output + symbol + ".missense.3d.dist"
            nonsense_dist_file_output = dist_file_output + symbol + ".nonsense.3d.dist" + dist_file_output
            
        (miss_dist, miss_prob) = get_p_value(rates,
                                             three_d_locations,
                                             iterations,
                                             "missense",
                                             miss_cds_positions,
                                             missense_scores,
                                             threshold,
                                             p,
                                             missense_dist_file_output.encode('utf-8'))
        #if iteration == 1:
        #    os.rename('current_distribution.txt', symbol + ".missense.3d.dist.txt")

        (nons_dist, nons_prob) = get_p_value(SiteRates(transcript,mut_dict, {}),
                                             three_d_locations,
                                             iterations,
                                             "lof",
                                             nons_cds_positions,
                                             nonsense_scores,
                                             p,
                                             nonsense_dist_file_output.encode('utf-8'))
        #if iteration == 1:
        #    os.rename('current_distribution.txt', symbol + ".nonsense.3d.dist.txt")

        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        results["missense"] += [[transcript.get_start(), transcript.get_end(), len(missense_events), miss_dist, miss_prob]]
        results["nonsense"] += [[transcript.get_start(), transcript.get_end(), len(nonsense_events), nons_dist, nons_prob]]        
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
    return probs, results

async def de_novos_entropy(symbol, variants, gene, p = 0, dist_file_output = "", iterations=1000000, mut_dict=None, scores_file = None, dbNSFP_score_obj = None, annotator =None):
#async def cluster_de_novos(symbol, variants, three_d_locations, ensembl, p = 0, iterations=1000000, mut_dict=None, scores_file = None, pvalues_in = None,variant_type = "de_novo", inherited_controls = {}):    
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
    mut_dict = load_mutation_rates()
    missense = [pos for pos,alt in variants["missense"]]
    nonsense = [pos for pos,alt in variants["nonsense"]]

    # load the set of transcripts that are the  minimum set of transcripts
    # required to contain all the de novos, unless we can't find any coding
    # transcripts that contain the de novos.
    try:
        transcripts = gene.transcripts
        minimized = minimise_transcripts_2(transcripts, missense + nonsense)
        transcripts = [x for x in transcripts if x.get_name() in minimized]
        #transcripts = await load_gene(ensembl, symbol, missense + nonsense)
    except IndexError as e:
        print(e)
        return None
    
    probs = {"miss_prob": [], "nons_prob": []}
    dists = {"miss_dist": [], "nons_dist": []}
    results = {"missense":[], "nonsense": []}

    iteration = 0
    for transcript in transcripts:
        iteration += 1

        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        scores = {}
        # if dbNSFP_score_obj is not None:
        #     dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
        #     if annotator in dbNSFP_header:
        #         anno_idx = dbNSFP_header.index(annotator)
        #         alt_idx = dbNSFP_header.index('alt')
        #         pos_idx = dbNSFP_header.index('pos(1-based)')
        #         for line in dbNSFP_score_obj.fetch(transcript.get_chrom(),
        #                                            transcript.get_start()-1,
        #                                            transcript.get_end()):
        #             values = line.split('\t')
        #             pos = values[pos_idx]
        #             alt = values[alt_idx]
        #             if values[anno_idx] == '.' :
        #                 #print(pos)
        #                 #print(alt)
        #                 continue
        #             score = float(values[anno_idx])
        #             if int(pos) not in scores:
        #                 scores[int(pos)] = {}
        #             scores[int(pos)][alt.encode('utf-8')] = float(score)
                    
        #     missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
        #     nonsense_scores = [-1 for x in nonsense_events]                        
        # elif scores_file is not None:
        #     for line in scores_file.fetch(transcript.get_chrom(),
        #                                   transcript.get_start()-1,
        #                                   transcript.get_end()):
        #         #_, pos, _, alt, _, score = line.split('\t')
        #         pos = line.split('\t')[1]
        #         alt = line.split('\t')[3]
        #         score = line.split('\t')[13]                
        #         if int(pos) not in scores:
        #             scores[int(pos)] = {}
        #         scores[int(pos)][alt.encode('utf-8')] = float(score)
        #     missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
        #     nonsense_scores = [-1 for x in nonsense_events]            
        # else:
        missense_scores = [-1 for x in missense_events]
        nonsense_scores = [-1 for x in nonsense_events]
        # print(transcript.get_name())

        rates = SiteRates(transcript, mut_dict, scores)


        # Move to no longer need transcript in next layer
        miss_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in missense_events ]
        nons_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]        

        missense_dist_file_output=""
        nonsense_dist_file_output=""
        if(dist_file_output != None):
            missense_dist_file_output = dist_file_output + symbol + ".missense.3d.dist"
            nonsense_dist_file_output = dist_file_output + symbol + ".nonsense.3d.dist"
        
        (miss_dist, miss_prob) = get_p_value_entropy(rates,
                                                     iterations,
                                                     "missense",
                                                     miss_cds_positions,
                                                     missense_scores,
                                                     p,
                                                     missense_dist_file_output.encode('utf-8'))
        (nons_dist, nons_prob) = get_p_value_entropy(SiteRates(transcript,mut_dict, {}),
                                                     iterations,
                                                     "lof",
                                                     nons_cds_positions,
                                                     nonsense_scores,
                                                     p,
                                                     nonsense_dist_file_output.encode('utf-8'))
        
        dists["miss_dist"].append(miss_dist)
        dists["nons_dist"].append(nons_dist)
        probs["miss_prob"].append(miss_prob)
        probs["nons_prob"].append(nons_prob)
        results["missense"] += [[transcript.get_start(), transcript.get_end(), len(missense_events), miss_dist, miss_prob]]
        results["nonsense"] += [[transcript.get_start(), transcript.get_end(), len(nonsense_events), nons_dist, nons_prob]]        
        # remove the de novos analysed in the current transcript, so that
        # analysis of subsequent transcripts uses independent events. NOTE THAT
        # THIS MIGHT MISS SOME CLUSTERING ACROSS MUTUALLY EXCLUSIVE TRANSCRIPTS
        # IF THE DE NOVO EVENTS ARE NEAR THE TRANSCRIPT DIVERGENCE.
        missense = [x for x in missense if x not in missense_events]
        nonsense = [x for x in nonsense if x not in  nonsense_events]
        
    for key in dists:
        dists[key] = ",".join([ str(x) for x in dists[key] ])

    probs = {k: fishers_method(probs[k]) for k in probs}
    probs.update(dists)
    return probs, results



async def cluster_de_novos_multi(multimer_name,
                                 all_chains,
                                 proteins,
                                 proteins_to_chains,
                                 variants,
                                 three_d_locations,
                                 ensembl,
                                 p = 0,
                                 dist_file_output = "",
                                 iterations=1000000,
                                 mut_dict=None,
                                 scores_file = None,
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
    for chain,residue_number,x,y,z in three_d_locations:
        if chain not in residues:
            residues[chain] = set()
        residues[chain].add(int(residue_number))
        

    # Make three_d_locations into list of dictionary
    three_d_list_of_dict = [{}]*len(all_chains)
    for chain,residue_number,x,y,z in three_d_locations:
        if chain in all_chains:
            chain_number = all_chains.index(chain)
            three_d_list_of_dict[int(chain_number)][int(residue_number)] = [float(x),float(y),float(z)]

    missense_cds_positions = {}
    all_missense_locations = []
    all_missense_scores = []
    de_novo_count_per_chain = [None]*len(all_chains)
    
    for protein, chains in proteins_to_chains.items():
        missense = [int(pos) for pos,alt in variants[protein]["missense"]]
#        nonsense = [pos for pos,alt in variants[protein]["nonsense"]]
    
        try:
            transcript = await best_transcript(ensembl,protein,missense)
        except IndexError as e:
            print(e)
            return None
        print(protein)
        print(missense)
        missense_events = get_de_novos_in_transcript(transcript, missense)
        print(missense_events)
        #        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)
        
        missense_cds_positions[protein] = [ int(transcript.get_coding_distance(x)['pos']) for x in missense_events ]
        #        nonsense_cds_positions[protein] = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]
        print(missense_cds_positions[protein])

        these_residues = set(residues[chains[0]])
        for chain in chains:
            these_residues = set.intersection(residues[chain], these_residues)
        
        # Trim down de novos to those in the 3d multimer model valid amino acids

        print([floor(float(cds_pos)/3) for cds_pos in missense_cds_positions[protein]])
        for cds_pos in missense_cds_positions[protein]:
            print(cds_pos)
            print(floor(float(cds_pos)/3))
            print(floor(float(cds_pos)/3) in these_residues)
            print(transcript.get_position_on_chrom(cds_pos,0))
        missense_cds_positions[protein] = [cds_pos for cds_pos in missense_cds_positions[protein] if floor(float(cds_pos)/3) in these_residues]

        missense_records = {}
        missense_records_map = {}
        missense_records[protein] = [(pos, int(transcript.get_coding_distance(pos)['pos']), alt) for pos,alt in variants[protein]["missense"]]        
        missense_records_map[protein] = {pos:alt for pos, cds_pos, alt in missense_records[protein] if floor(float(cds_pos)/3) in these_residues}
        print(missense_records[protein])
        print(missense_records_map[protein])
        #        nonsense_cds_positions[protein] = [cds_pos in nonsense_cds_positions if these_residues.contains(floor(cds_pos/3))]
        # Mark how many de novos to expect on chains corresponding to this protein
        for chain in chains:
            de_novo_count_per_chain[all_chains.index(chain)] = len(missense_cds_positions[protein])

        # Import scores
        scores = {}
        scores[protein] = {}
        if scores_file is not None:
            if scores[protein] == {}:
                for line in scores_file.fetch(transcript.get_chrom(),
                                              transcript.get_start()-1,
                                              transcript.get_end()):
                    pos = line.split('\t')[1]
                    alt = line.split('\t')[3]
                    score = line.split('\t')[13]                
                    if int(pos) not in scores[protein]:
                        scores[protein][int(pos)] = {}
                    scores[protein][int(pos)][alt.encode('utf-8')] = float(score)


        # Collect the missense for this protein
        for chain in chains:
            for cds_pos in missense_cds_positions[protein]:
                all_missense_locations.append(three_d_list_of_dict[all_chains.index(chain)][floor(cds_pos/3)])
                if scores[protein] != {}:
                    print(transcript.get_position_on_chrom(cds_pos,0))
                    print(scores[protein][transcript.get_position_on_chrom(cds_pos,0)])
                    alt = missense_records_map[protein][transcript.get_position_on_chrom(cds_pos,0)]
                    all_missense_scores.append(scores[protein][transcript.get_position_on_chrom(cds_pos,0)][alt.encode('utf-8')])
                else:
                    all_missense_scores.append(-1)

        
        #create and append the rate object for protein
        rates.append({"protein" : protein, "rate" : SiteRates(transcript, mut_dict, scores[protein], these_residues)})

        for rate in rates:
            print(rate["rate"]["missense"].get_summed_rate())
            print(len(rate["rate"]["missense"]))

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

#    (nons_dist, nons_prob) = get_p_value_multi(SiteRates(transcript,mut_dict, {}),
#                                               three_d_locations,
#                                               iterations,
#                                               "lof",
#                                               nonsense_cds_positions,
#                                               nonsense_scores,
#                                               p)
        
    dists["miss_dist"].append(miss_dist)
#    dists["nons_dist"].append(nons_dist)
    probs["miss_prob"].append(miss_prob)
#    probs["nons_prob"].append(nons_prob)
    results["missense"] += [[transcript.get_start(), transcript.get_end(), len(missense_events), miss_dist, miss_prob]]
#    results["nonsense"] += [[transcript.get_start(), transcript.get_end(), len(nonsense_events), nons_dist, nons_prob]]        
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



async def cluster_de_novos_coevol(symbol, variants, coevol, gene, p = 0, dist_file_output = "", iterations=1000000, mut_dict=None,scores = {}, pvalues_in=None,
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
        minimized = minimise_transcripts_2(transcripts, missense + nonsense)
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




async def cluster_de_novos_1d(symbol,
                              variants,
                              ensembl,
                              threshold,
                              p = 0,
                              dist_file_output = "",
                              iterations=1000000,
                              mut_dict=None,
                              scores_file = None,
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
        transcripts = await load_gene(ensembl, symbol, missense + nonsense)
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
            for line in scores_file.fetch(transcript.get_chrom(),
                                          transcript.get_start()-1,
                                          transcript.get_end()):
                #_, pos, _, alt, _, score = line.split('\t')
                pos = line.split('\t')[1]
                alt = line.split('\t')[3]
                score = line.split('\t')[13]                
                if int(pos) not in scores:
                    scores[int(pos)] = {}
                scores[int(pos)][alt.encode('utf-8')] = float(score)                
            missense_scores = [scores[int(pos)][alt.encode('utf-8')] for pos,alt in variants["missense"]]
            nonsense_scores = [1 for x in nonsense_events]
        else:
            missense_scores = [1 for x in missense_events]
            nonsense_scores = [1 for x in nonsense_events]

        # Handling thresholding of distribution
        rates = SiteRates(transcript, mut_dict, scores)
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
