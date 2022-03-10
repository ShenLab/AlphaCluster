""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from denovonear.weights import geomean, get_distances, get_distances_1d, get_distances_coevol, scale_distances, analyse_de_novos, analyse_de_novos_1d, analyse_de_novos_entropy, analyse_de_novos_coevol, analyse_de_novos_multi, WeightedChoices

from west_weights import get_pred_count

from math import floor, log

def get_p_value(rates, three_d_locations, iterations, consequence, cds_positions, de_novos_scores, threshold = None, p = 0,dist_file_output=""):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        three_d_locations: dictionary of three_d location for each base-pair location
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(cds_positions) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    weights = rates[consequence]
    #cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    #    codon_numbers = [ transcript.get_codon_info(x)['codon_number'] for x in de_novos]
    #    distances = get_distances(cds_positions)
    print(cds_positions)
    for cds_pos in cds_positions:
        print(floor(float(cds_pos)/3)+1)

    #Handle thresholding of input
    if threshold != None:
        print("THRESHOLDING")
        print(threshold)
        temp_pos = []        
        temp_sco = []
        for i, val in enumerate(de_novos_scores):
            if val >= float(threshold):
                temp_pos.append(cds_positions[i])
                temp_sco.append(de_novos_scores[i])
        cds_positions = temp_pos
        de_novos_scores = temp_sco
        print(cds_positions)
        print(de_novos_scores)
        if len(cds_positions) < 2:
            return (float('nan'), float('nan'))
        
    distances = get_distances(cds_positions,
                              three_d_locations)

    if len(de_novos_scores) > 0 and de_novos_scores[0] != -1 and threshold == None:
        print("SCORE SCALING")
        distances = scale_distances(distances,
                                    de_novos_scores)

    print(distances)
    observed = geomean(distances, p, 3.5)
    print("observed mean = " + str(observed))
    # call a cython wrapped C++ library to handle the simulations

    print("calling analyse_de_novos")
    #sys.stderr.write(transcript)
    print(len(cds_positions))
    print(observed)
    sim_prob = analyse_de_novos(weights, three_d_locations, iterations, len(cds_positions), observed, p, dist_file_output)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)

def get_p_value_entropy(rates, iterations, consequence, cds_positions, de_novos_scores, p = 0,dist_file_output=""):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        three_d_locations: dictionary of three_d location for each base-pair location
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(cds_positions) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    weights = rates[consequence]
    #cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    #    codon_numbers = [ transcript.get_codon_info(x)['codon_number'] for x in de_novos]
    #    distances = get_distances(cds_positions)
    print(cds_positions)
    for cds_pos in cds_positions:
        print(floor(float(cds_pos)/3)+1)

    n_total = len(cds_positions)
    n_dict = {}
    obs_entropy = 0
    for cds_pos in cds_positions:
        if cds_pos not in n_dict:
            n_dict[cds_pos] = 0
        n_dict[cds_pos] += 1
    for key,n_i in n_dict.items():
        obs_entropy += -1*n_i/n_total * log(n_i/n_total)
    
    print("observed entropy = " + str(obs_entropy))
    print("observed entropy normalized = " + str(obs_entropy/log(n_total)))    
    # call a cython wrapped C++ library to handle the simulations

    print("calling analyse_de_novos")
    #sys.stderr.write(transcript)
    sim_prob = analyse_de_novos_entropy(weights, iterations, len(cds_positions), obs_entropy, p, dist_file_output)
    observed = "{0:0.1f}".format(obs_entropy)
    
    return (observed, sim_prob)


def get_p_value_west(N_males, N_females, rates, three_d_locations, iterations, consequence, cds_positions, de_novos_scores, p = 0):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        three_d_locations: dictionary of three_d location for each base-pair location
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(cds_positions) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    weights = rates[consequence]
    #cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    #    codon_numbers = [ transcript.get_codon_info(x)['codon_number'] for x in de_novos]
    #    distances = get_distances(cds_positions)
    print(three_d_locations[0])
    print(cds_positions)
    for cds_pos in cds_positions:
        print(floor(float(cds_pos)/3))

    distances = get_distances(cds_positions, three_d_locations)
    print(distances[0])
    observed = sum(de_novos_scores)
    #observed = geomean(distances, p)
    print(observed)

    #if de_novos_scores[0] != -1 :
    #    distances = scale_distances(distances, de_novos_scores)
    #observed = geomean(distances, p)
    #print(observed)
    # call a cython wrapped C++ library to handle the simulations

    lambda_ = get_pred_count(weights.get_summed_rate(),
                             N_male,
                             N_female,
                             transcript.get_chrom())
    
    print("calling analyse_de_novos_west")
    #sys.stderr.write(transcript)
    sim_prob = analyse_de_novos_west(weights,
                                     three_d_locations,
                                     iterations,
                                     lambda_,
                                     len(cds_positions),
                                     observed,
                                     p)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)


def get_p_value_multi(chains,
                      proteins,
                      de_novo_count_per_chain,
                      rates_vector,
                      three_d_list_of_dict,
                      iterations,
                      consequence,
                      de_novo_locations,
                      de_novo_scores,
                      p = 0,
                      dist_file_output=""):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        three_d_locations: dictionary of three_d location for each base-pair location
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(de_novo_locations) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]

    #    cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    #    codon_numbers = [ transcript.get_codon_info(x)['codon_number'] for x in de_novos]
    #    distances = get_distances(cds_positions)

    #distances = get_distances_multi(cds_positions, three_d_locations)
    distances = []
    idx = 0
    for i in range(len(de_novo_locations)):
        for j in range(i+1, len(de_novo_locations)):
            distance = pow(
                           pow(de_novo_locations[i][0] - de_novo_locations[j][0],2) +
                           pow(de_novo_locations[i][1] - de_novo_locations[j][1],2) +
                           pow(de_novo_locations[i][2] - de_novo_locations[j][2],2)
                           ,0.5)
            distances.append(distance)

    assert(len(distances) == len(de_novo_locations)*(len(de_novo_locations)-1)/2)
            
    observed = geomean(distances, p, 3.5)
    print(observed)

    if de_novo_scores[0] != -1 :
        distances = scale_distances(distances, de_novo_scores)
    observed = geomean(distances, p, 3.5)
    print(observed)
    # call a cython wrapped C++ library to handle the simulations

    print("calling analyse_de_novos")
    #sys.stderr.write(transcript)
    print(de_novo_count_per_chain)
    weights_vector = []
    for rate in rates_vector:
        print(type(rate["rate"][consequence]))
        weights_vector.append({"protein": rate["protein"].encode('utf-8'), "rate" : rate["rate"][consequence]})

    chains = [a.encode('utf-8') for a in chains]
    proteins = [a.encode('utf-8') for a in proteins]
    print(chains)
    sim_prob = analyse_de_novos_multi(chains,
                                      proteins,
                                      de_novo_count_per_chain,
                                      weights_vector,
                                      three_d_list_of_dict,
                                      iterations,
                                      len(de_novo_locations),
                                      observed,
                                      p,
                                      dist_file_output)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)



def get_p_value_coevol(transcript, rates, coevol, iterations, consequence, de_novos, p = 0,dist_file_output=""):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        coevol: coevolutionary strengths
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(de_novos) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    
    weights = rates[consequence]
    
    cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
#    codon_numbers = [ transcript.get_codon_info(x)['codon_number'] for x in de_novos]
#    distances = get_distances(cds_positions)
    distances = get_distances_coevol(cds_positions, coevol)
    print(distances)
    observed = geomean(distances, p, 1)
    print(observed)
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos_coevol(weights, coevol, iterations, len(de_novos), observed, p, dist_file_output)
    
    observed = "{0:0.0001f}".format(observed)
    
    return (observed, sim_prob)



def get_p_value_1d(transcript, rates, iterations, consequence, de_novos, de_novos_scores = [], threshold = None, p = 0,dist_file_output=""):
    """ find the probability of getting de novos with a mean conservation
    
    The probability is the number of simulations where the mean conservation
    between simulated de novos is less than the observed conservation.
    
    Args:
        transcript: Transcript object for the current gene.
        rates: SiteRates object, which contains WeightedChoice entries for
            different consequence categories.
        iterations: number of simulations to perform
        consequence: string to indicate the consequence type e.g. "missense, or
            "lof", "synonymous" etc. The full list is "missense", "nonsense",
            "synonymous", "lof", "loss_of_function", "splice_lof",
            "splice_region".
        de_novos: list of de novos within a gene
    
    Returns:
        tuple of mean proximity for the observed de novos and probability of
        obtaining a value less than or equal to the observed proximity from the
        null distribution.
    """
    
    if len(de_novos) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    
    weights = rates[consequence]
    
    cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]

    if threshold != None:
        print("THRESHOLDING")
        print(threshold)
        temp_pos = []        
        temp_sco = []
        for i, val in enumerate(de_novos_scores):
            if val >= float(threshold):
                temp_pos.append(cds_positions[i])
                temp_sco.append(de_novos_scores[i])
        cds_positions = temp_pos
        de_novos_scores = temp_sco
        print(cds_positions)
        print(de_novos_scores)

    distances = get_distances_1d(cds_positions)
    if de_novos_scores[0] != -1 and threshold == None:
        print("SCORE SCALING")
        distances = scale_distances(distances, de_novos_scores)
    observed = geomean(distances, p, 1)
    print(observed)
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos_1d(weights, iterations, len(cds_positions), observed, p, dist_file_output)
    
    observed = "{0:0.1f}".format(observed)
    
    return (observed, sim_prob)
