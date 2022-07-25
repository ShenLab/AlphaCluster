""" class to analyse clustering of known de novos in genes according to their
distances apart within the gene, and compare that to simulated de novo events
within the same gene.
"""

from alphacluster.weights import geomean, get_distances, get_distances_1d, get_distances_coevol, scale_distances, analyse_de_novos, analyse_de_novos_1d,  analyse_de_novos_coevol, analyse_de_novos_multi, WeightedChoices

from alphacluster.west_weights import get_pred_count

from math import floor, log

def get_p_value(rates, three_d_locations, iterations, consequence, cds_positions, de_novos_scores, scale = False, threshold = None, p = 0,dist_file_output=""):
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
        return (float('nan'), float('nan'), len(cds_positions))
    
    weights = rates[consequence]
    three_d_locations_xyz = [[row[0],row[1], row[2]] for row in three_d_locations]
    three_d_locations_xyz = [[0,0,0]] + three_d_locations_xyz

    #Handle thresholding of input
    if threshold != None:
        temp_pos = []        
        temp_sco = []
        for i, val in enumerate(de_novos_scores):
            print(i)
            print(val)
            if val >= float(threshold):
                temp_pos.append(cds_positions[i])
                temp_sco.append(de_novos_scores[i])
        cds_positions = temp_pos
        de_novos_scores = temp_sco
        if len(cds_positions) < 2:
            return (float('nan'), float('nan'), len(cds_positions))

    # Calculate observed mean
    distances = get_distances(cds_positions, three_d_locations_xyz)
    if len(de_novos_scores) > 0 and de_novos_scores[0] != -1 and scale:
        distances = scale_distances(distances, de_novos_scores)
    observed = geomean(distances, p, 3.5)

    # Simulate to obtain distribution of means and p-value of observed
    sim_prob = analyse_de_novos(weights, three_d_locations_xyz, iterations, len(cds_positions), observed, p, dist_file_output)
    observed = "{0:0.1f}".format(observed)
    return (observed, sim_prob, len(cds_positions))


def get_p_value_west(N_males, N_females, rates, three_d_locations, iterations, consequence, cds_positions, de_novos_scores, p = 0):
    """ find the probability of getting de novos with a mean weight
    
    The probability is the number of simulations where the mean weight
    between simulated de novos is less than the observed weight.
    """
    
    if len(cds_positions) < 2:
        return (float('nan'), float('nan'))
    
    rename = {"lof": "loss_of_function"}
    if consequence in rename:
        consequence = rename[consequence]
    weights = rates[consequence]

    distances = get_distances(cds_positions, three_d_locations)
    observed = sum(de_novos_scores)
    
    lambda_ = get_pred_count(weights.get_summed_rate(),
                             N_male,
                             N_female,
                             transcript.get_chrom())
    
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

    # Calculate observed mean
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
    if de_novo_scores[0] != -1 :
        distances = scale_distances(distances, de_novo_scores)
    observed = geomean(distances, p, 3.5)

    # Create weights vector for each protein in chain
    weights_vector = []
    for rate in rates_vector:
        print(type(rate["rate"][consequence]))
        weights_vector.append({"protein": rate["protein"].encode('utf-8'), "rate" : rate["rate"][consequence]})

    chains = [a.encode('utf-8') for a in chains]
    proteins = [a.encode('utf-8') for a in proteins]

    # Run simulation for expected mean distribution and p-value
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
    
    weights = rates[consequence]
    cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in de_novos ]
    distances = get_distances_coevol(cds_positions, coevol)
    observed = geomean(distances, p, 1)
    
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos_coevol(weights, coevol, iterations, len(de_novos), observed, p, dist_file_output)

    observed = "{0:0.0001f}".format(observed)        
    return (observed, sim_prob)



def get_p_value_1d(transcript, rates, iterations, consequence, de_novos, de_novos_scores = [], scale = False, threshold = None, p = 0,dist_file_output=""):
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
        temp_pos = []        
        temp_sco = []
        for i, val in enumerate(de_novos_scores):
            if val >= float(threshold):
                temp_pos.append(cds_positions[i])
                temp_sco.append(de_novos_scores[i])
        cds_positions = temp_pos
        de_novos_scores = temp_sco
        if len(cds_positions) < 2:
            return (float('nan'), float('nan'))

    distances = get_distances_1d(cds_positions)
    if de_novos_scores[0] != -1 and scale:
        distances = scale_distances(distances, de_novos_scores)
    observed = geomean(distances, p, 1)
    # call a cython wrapped C++ library to handle the simulations
    sim_prob = analyse_de_novos_1d(weights, iterations, len(cds_positions), observed, p, dist_file_output)
    
    observed = "{0:0.1f}".format(observed)
    return (observed, sim_prob)
