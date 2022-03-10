""" Script to investigate the probability of multiple mutations clustering
within a single gene.
"""

import os
import sys
import asyncio
import argparse
import logging
import pysam
import random
from scipy.stats import poisson

from rate_limiter import RateLimiter
from load_mutation_rates import load_mutation_rates
from load_de_novos import load_de_novos, load_de_novos_chrom_pos_alt
from load_inherited import load_inherited, load_inherited_controls
from west_weights import generate_expected_buckets, generate_observed_buckets, generate_scores_from_buckets, get_pred_count
from cluster_test import cluster_de_novos, cluster_de_novos_1d, de_novos_entropy, cluster_de_novos_coevol, cluster_de_novos_multi, fishers_method
from load_three_d_locations import load_three_d_locations, load_three_d_multimer
from load_coevol import load_basic_coevol_strength
from load_pvalues import load_pvalues
from load_gene import (construct_gene_object,
                                  count_de_novos_per_transcript, minimise_transcripts, minimise_transcripts_2, get_transcript_ids,  get_de_novos_in_transcript)
from simulate import get_p_value_west
from denovonear.site_specific_rates import SiteRates
from denovonear.frameshift_rate import include_frameshift_rates
from denovonear.log_transform_rates import log_transform
from denovonear.gencode import Gencode
import operator as op
from functools import reduce

def ncr(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer // denom  # or / in Python 2

async def _load_gencode(symbols):
    ''' load gene coords and sequence via ensembl
    '''
    gencode = Gencode()
    async with RateLimiter(per_second=15) as ensembl:
        tasks = [load_gene(ensembl, symbol) for symbol in symbols]
        genes = await asyncio.gather(*tasks)
        for gene in genes:
            gencode.add_gene(gene)
        return gencode

def load_gencode(symbols, gencode=None, fasta=None):
    ''' load genes from gencode annotations file, with ensembl as backup
    '''
    if gencode and fasta:
        return Gencode(gencode, fasta)
    
    # use ensembl as backup if gencode file not available. This restricts the
    # asynchronous calls to within one section, and ensures they are called together
    return asyncio.get_event_loop().run_until_complete(_load_gencode(symbols))


async def west_obs_buckets_creation(ensembl, mut_dict, output, args):
    de_novos = load_de_novos_chrom_pos_alt(args.input)
    bucket_def = {}

    for idx in range(len(args.bucket_ann)):
        bucket_def[args.bucket_ann[idx]] = [float(x) for x in (args.bucket_cutoffs[idx]).strip("[]").split(',')]

    annotation_files = {"gMVP_rankscore" : args.gMVP_file ,
                        "dbNSFP": args.dbNSFP_file,
                        "dbNSFP_gene": args.dbNSFP_gene_file,                        
                        "MCR_region" : args.MCR_file}
    await generate_observed_buckets(args.output_dir,
                                    de_novos,
                                    bucket_def,
                                    annotation_files)



async def west_exp_buckets_creation(ensembl, mut_dict, output, args):
    N_males = args.N_males
    N_females = args.N_females
    bucket_def = {}

    for idx in range(len(args.bucket_ann)):
        bucket_def[args.bucket_ann[idx]] = [float(x) for x in (args.bucket_cutoffs[idx]).strip("[]").split(',')]

    annotation_files = {"gMVP_rankscore" : args.gMVP_file ,
                        "dbNSFP": args.dbNSFP_file,
                        "dbNSFP_gene": args.dbNSFP_gene_file,                        
                        "MCR_region" : args.MCR_file}

    await generate_expected_buckets(output,
                                    N_males,
                                    N_females,
                                    args.gene,
                                    load_gencode(args.gene,
                                                 args.gencode,
                                                 args.fasta),
                                    mut_dict,
                                    bucket_def,
                                    annotation_files)


async def west_scores_creation(ensembl, mut_dict, output, args):
    bucket_def = {}

    for idx in range(len(args.bucket_ann)):
        bucket_def[args.bucket_ann[idx]] = [float(x) for x in (args.bucket_cutoffs[idx]).strip("[]").split(',')]

    print("BEGIN")

    genes = load_gencode(None, args.gencode, args.fasta)            
    await generate_scores_from_buckets(args.output_dir, 
                                       genes,
                                       args.obs_buckets, 
                                       args.buckets_dir,
                                       bucket_def)
 

'''
async def west_weight_creation(ensembl, mut_dict, output, args):
    de_novos = load_de_novos_chrom_pos_alt(args.input)
    #print(de_novos)
    N_males = args.N_males
    N_females = args.N_females
    bucket_def = {}

    for idx in range(len(args.bucket_ann)):
        bucket_def[args.bucket_ann[idx]] = [float(x) for x in (args.bucket_cutoffs[idx]).strip("[]").split(',')]
        #print(args.bucket_ann[idx])
        #print(args.bucket_cutoffs[idx])

    #print(bucket_def)
    #bucket_def = {"gMVP_ranksscore" : [.25, .5, .75, .9,1.0], 
    #              "CADD_rankscore"  : [.25, .5, .75, .9,1.0], 
    #              "MCR_region" : [True, False],               
    #              "gnomad_pLI" : [.1,.9,1.0]                  
    #              }                                          

    annotation_files = {"gMVP_rankscore" : args.gMVP_file ,
                        "dbNSFP": args.dbNSFP_file,
                        "dbNSFP_gene": args.dbNSFP_gene_file,                        
                        "MCR_region" : args.MCR_file}

    #with open(args.all_genes) as f:
    #    all_genes = f.read().splitlines()

    if(args.gene == "observed"):
        await generate_observed_buckets(output,
                                        de_novos,
                                        bucket_def,
                                        annotation_files)
    else:
        await generate_buckets(output,
                               N_males,
                               N_females,
                               args.gene,
                               load_gencode(args.gene,
                                            args.gencode,
                                            args.fasta),
                               mut_dict,
                               bucket_def,
                               annotation_files)

'''
async def mvp3d(ensembl, mut_dict, output, args):
    # Load variants, both de_novo and inherited/ultrarare
    de_novos = {}
    inherited = {}
    inherited_controls = {}
    if "de_novos" in args and args.de_novos == True:
        de_novos = load_de_novos(args.input)
        de_novos = {protein : de_novos[protein] for protein in args.proteins.split(",")}
    if "inherited" in args and args.inherited == True:
        inherited = []
        for protein in args.proteins:
            transcripts = await get_transcript_ids(ensembl, protein)
            min_start = float("inf")
            max_end = 0
            for transcript in transcripts:
                gene = await construct_gene_object(ensembl, transcript)            
                chrom = gene.get_chrom()
                if min_start > gene.get_start():
                    min_start = gene.get_start()
                    if max_end < gene.get_end():
                        max_end = gene.get_end()
                        
                        # TODO: once UKBB is done annotating, make inherited controls
                        # capable of picking only missense and/or LoF
                        inherited_controls = {protein : {"missense": load_inherited_controls(chrom,
                                                                                             min_start,
                                                                                             max_end)
                                              }}
        
                        inherited = {protein : {'missense' : load_inherited(chrom,
                                                                            min_start,
                                                                            max_end,
                                                                            ["missense_variant"]),
                                                'nonsense' : load_inherited(chrom,
                                                                            min_start,
                                                                            max_end,
                                                                            ["stop_gained",
                                                                             "stop_lost",
                                                                             "start_lost",
                                                                             "frameshift_variant",
                                                                             "splice_acceptor_variant",
                                                                             "splice_donor_variant",
                                                                             "transcript_ablation"])}
                                     }


    #Handle score loading
    scores = {}
    if args.scores is not None:
        score_obj = pysam.TabixFile(args.scores)
        for line in scores_obj.fetch(transcript.get_chrom(),
                                     transcript.get_start()-1,
                                     transcript.get_end()):
            _, pos, _, alt, _, score = line.split('\t')
            if int(pos) not in scores:
                scores[int(pos)] = {}
            scores[int(pos)][alt.encode('utf-8')] = float(score)

    # Start clustering analysis
    iterations = 1000000
    output_file_has_header_line = False

    if args.mode == "multimer":
        three_d_locations = await load_three_d_multimer(args.protein_dir, args.proteins)

        
        
    else:
        for symbol in args.proteins:
            if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
                continue

            three_d_locations = []
            coev_relations = []

            if args.mode == "3d":
                three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
            elif args.mode == "coevol":
                coev_relations = await load_basic_coevol_strength(args.coev_dir, symbol)    
        
            # (probs,results) = await cluster_analysis(symbol,
            #                                          de_novos[symbol],
            #                                          inherited[symbol],
            #                                          inherited_controls[symbol],
            #                                          three_d_locations_list,
            #                                          coev_relations,
            #                                          ensembl,
            #                                          float(args.degree),
            #                                          iterations,
            #                                          mut_dict,
            #                                          scores)
        
            if probs is None:
                continue
            with open(output, "x") as output_file:
                # output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
                output_file.write("gene_id\tmutation_category\tstart \tend \t events_n \tdist \tprobability \n")                
                output_file_has_header_line = True
                for result in results["missense"]:
                    output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                                            "missense",
                                                                            result[0],
                                                                            result[1],
                                                                            result[2],
                                                                            result[3],
                                                                            result[4]))
                for result in results["nonsense"]:
                    output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                                            "nonsense",
                                                                            result[0],
                                                                            result[1],
                                                                            result[2],
                                                                            result[3],
                                                                            result[4]))

            # output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            # len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            # output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            # len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

        for symbol in sorted(inherited):
            if len(inherited[symbol]["missense"] + inherited[symbol]["nonsense"]) < 2:
                continue
            three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)

            genes = load_gencode(None, args.gencode, args.fasta)            
            probs = await cluster_de_novos(symbol,
                                           inherited[symbol],
                                           three_d_locations_list,
                                           genes[symbol],
                                           args.degree,
                                           args.dist_file_output,
                                           args.runs,#iterations,
                                           mut_dict,
                                           "inherited",
                                           inherited_controls)
        
            if probs is None:
                continue
            with open(output, "x") as output_file:
                if output_file_has_header_line == False:
                    output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")    
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "inherited_missense",
                                                                len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "inherited_nonsense",
                                                                len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

async def type_II_testing(ensembl, mut_dict, output, args):
    # Get gene set
#    genes_to_test = get_genes_for_test(args.min_length,
#                               args.max_length,
#                               args.min_af_err,
#                               args.max_af_err,
#                               args.constraint,
#                               args.min_avg_gmvp,
#                               args.max_avg_gmpv,
#                               args.min_avg_revel,
#                               args.max_avg_revel)

    symbol = args.protein
    genes = load_gencode(None, args.gencode, args.fasta)
    mut_dict = load_mutation_rates()
    try:
        gene = genes[symbol]
        transcripts = gene.transcripts
    except IndexError as e:
        print(e)
        return None

    transcript=transcripts[0]
    
    #     scores = {}
    #     if dbNSFP_score_obj is not None:
    #         dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
    #         annotator = args.dbNSFP_annotator
    #         if annotator in dbNSFP_header:
    #             anno_idx = dbNSFP_header.index(annotator)
    #             alt_idx = dbNSFP_header.index('alt')
    #             pos_idx = dbNSFP_header.index('pos(1-based)')
    #             for line in dbNSFP_score_obj.fetch(transcript.get_chrom(),
    #                                                    transcript.get_start()-1,
    #                                                    transcript.get_end()):
    #                 values = line.split('\t')
    #                 pos = values[pos_idx]
    #                 alt = values[alt_idx]
    #                 if values[anno_idx] == '.' :
    #                     print(pos)
    #                     print(alt)
    #                     continue
    #                 score = float(values[anno_idx])
    #                 if int(pos) not in scores:
    #                     scores[int(pos)] = {}
    #                 scores[int(pos)][alt.encode('utf-8')] = float(score)
    #     elif scores_file is not None:
    #         for line in scores_file.fetch(transcript.get_chrom(),
    #                                       transcript.get_start()-1,
    #                                       transcript.get_end()):
    #             pos = line.split('\t')[1]
    #             alt = line.split('\t')[3]
    #             score = line.split('\t')[13]                
    #             if int(pos) not in scores:
    #                 scores[int(pos)] = {}
    #             scores[int(pos)][alt.encode('utf-8')] = float(score)


    rates = SiteRates(transcript, mut_dict, {})
    total_mut = rates["missense"].get_summed_rate()
    de_novos = load_de_novos(args.input)
    if "protein" in args and args.protein in de_novos:
        de_novos = {args.protein : de_novos[args.protein]}
        symbol = args.protein
    else:
        return

    if len(de_novos[symbol]["missense"])  < 2:
        return

    three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
        
    if args.scores is not None:
        scores = pysam.TabixFile(args.scores)
    else:
        scores = None
           
    if args.dbNSFP is not None:
        dbNSFP_score_obj = pysam.TabixFile(args.dbNSFP)
    else:
        dbNSFP_score_obj = None
    
    if args.pvalues_in is not None:
        pvalues_in = load_pvalues(args.pvalues_in)
        if not symbol in pvalues_in:
            pvalues_in = None
    else:
        pvalues_in = None

    cohort = load_cohort()

        
    for iterations in range(1):#args.iterations):
        # num_of_dnvs = poisson.rvs(2*args.cohort_size*total_mut)#args.cohort_size)
        # # sample dnvs
        # de_novos[symbol]['missense'] = random.sample(de_novos[symbol]['missense'], num_of_dnvs)

        # de_novos[symbol]['nonsense'] = []
        # print(de_novos)
        # if len(de_novos[symbol]["missense"])  < 2:
        #     return


        #
        cohort = random.sample(cohort,args.cohort_size)
        de_novos[symbol]['nonsense'] = []
        de_novos[symbol]['missense'] = [(pos, alt) for pos, alt, iid  in de_novos_with_idd if iid in cohort]
        # do burden analysis
           
           
        #print(num_of_dnvs)
        #three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
        #for i in range(num_of_dnvs):
        #x = rates["missense"].choice_with_alleles()
        #de_novos[symbol] += [(transcript.get_position_on_chrom(x["pos"],x["offset"]), x["alt"])]

        #print(de_novos)
        #missense_scores = [scores[int(pos)][alt] for pos,alt in de_novos[symbol]]            
        #calculate pvalue_in from burden test, denovoWEST, or Poisson test
        
        #1d
        probs_1d = await cluster_de_novos_1d(symbol,
                                             de_novos[symbol],
                                             ensembl,
                                             args.degree,
                                             args.runs,#iterations,
                                             mut_dict)
        
        #1d gMVP
        probs_1d_gMVP = await cluster_de_novos_1d(symbol,
                                                  de_novos[symbol],
                                                  ensembl,
                                                  args.degree,
                                                  args.runs,#iterations,
                                                  mut_dict,
                                                  scores)

        #3d
        (probs_3d,results_3d) = await cluster_de_novos(symbol,
                                                       de_novos[symbol],
                                                       three_d_locations_list,
                                                       genes[symbol],
                                                       float(args.degree),
                                                       args.runs,#iterations,
                                                       mut_dict)
        #3d gMVP
        (probs_3d_gMVP,results_3d_gMVP) = await cluster_de_novos(symbol,
                                                                 de_novos[symbol],
                                                                 three_d_locations_list,
                                                                 genes[symbol],
                                                                 float(args.degree),
                                                                 args.runs,#iterations,
                                                                 mut_dict,
                                                                 scores)
        
        #3d REVEL
        (probs_3d_REVEL,results_3d_REVEL) = await cluster_de_novos(symbol,
                                                                   de_novos[symbol],
                                                                   three_d_locations_list,
                                                                   genes[symbol],
                                                                   float(args.degree),
                                                                   args.runs,#iterations,
                                                                   mut_dict,
                                                                   None,
                                                                   dbNSFP_score_obj,
                                                                   "REVEL_rankscore")

        #3d CADD
        (probs_3d_CADD,results_3d_CADD) = await cluster_de_novos(symbol,
                                                                 de_novos[symbol],
                                                                 three_d_locations_list,
                                                                 genes[symbol],
                                                                 float(args.degree),
                                                                 args.runs,#iterations,
                                                                 mut_dict,
                                                                 None,
                                                                 dbNSFP_score_obj,
                                                                 "CADD_raw_rankscore",
                                                                 None)

    
        for test_type, p in [("1d", probs_1d),("1d.gMVP", probs_1d_gMVP),("3d", probs_3d), ("3d.gMVP", probs_3d_gMVP), ("3d.REVEL", probs_3d_REVEL), ("3d.CADD", probs_3d_CADD)]:
            with open(symbol+"." + test_type+".out", "x") as output_file:
                output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), p["miss_dist"], p["miss_prob"]))
                if("miss_prob_comb" in p):
                    output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",len(de_novos[symbol]["missense"]), p["miss_dist_comb"], p["miss_prob_comb"]))
                    output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",len(de_novos[symbol]["missense"]), p["miss_dist_pois"], p["miss_prob_pois"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense", len(de_novos[symbol]["nonsense"]), p["nons_dist"], p["nons_prob"]))

            #p_cluster = probs["miss_prob"]
            #p_combine = probs["miss_prob_comb"]
            #alpha = 2.5E-6
            #if p_cluster < alpha:
            #    cluster_type_II += 1
            #if p_combine < alpha :
            #    combine_type_II += 1

    #print("type I of cluster test is " + str(cluster_type_II/len(genes)/args.iterations))
    #print("type I of combine test is " + str(combine_type_II/len(genes)/args.iterations))    
    
async def west_run(ensembl, mut_dict, output, args):

    # TODO: make load_de_denovos take args.protein as input
    de_novos = load_de_novos(args.input)
    iterations = 1000000
    # if protein specific run, limit the input variants to this protein/gene
    if "protein" in args and args.protein in de_novos:
        de_novos = {args.protein : de_novos[args.protein]}
    else:
        return
    print(de_novos)

    genes = load_gencode(args.protein,
                         args.gencode,
                         args.fasta)
    output_file_has_header_line = False
    for symbol in sorted(de_novos):
        gene = genes[symbol]
        transcript = gene.canonical
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            continue
        three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)

        if args.scores is not None:
            scores_file = pysam.TabixFile(args.scores)
        else:
            scores_file = None

        if args.pvalues_in is not None:
            pvalues_in = load_pvalues(args.pvalues_in)
            if not symbol in pvalues_in:
                pvalues_in = None
        else:
            pvalues_in = None

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
                nonsense_scores = [-1 for x in nonsense_events]            
        else:
            missense_scores = [-1 for x in missense_events]
            nonsense_scores = [-1 for x in nonsense_events]
        rates = SiteRates(transcript, mut_dict, scores)


        # Move to no longer need transcript in next layer
        miss_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in missense_events ]
        nons_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]        

        (miss_dist, miss_prob) = get_p_value_west(args.N_males,
                                                  args.N_females,
                                                  rates,
                                                  three_d_locations,
                                                  iterations,
                                                  "missense",
                                                  miss_cds_positions,
                                                  missense_scores,
                                                  0)
       
        if probs is None:
            continue

        with open(output, "x") as output_file:
            output_file_has_header_line = True
            # output_file.write("gene_id\tmutation_category\tstart \tend \t events_n \tdist \tprobability \n")                
            # for result in results["missense"]:
            #     output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
            #                                                             "missense",
            #                                                             result[0],
            #                                                             result[1],
            #                                                             result[2],
            #                                                             result[3],
            #                                                             result[4]))
            # for result in results["nonsense"]:
            #     output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
            #                                                             "nonsense",
            #                                                             result[0],
            #                                                             result[1],
            #                                                             result[2],
            #                                                             result[3],
            #                                                             result[4]))
            output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            if("miss_prob_comb" in probs):
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",
                len(de_novos[symbol]["missense"]), probs["miss_dist_comb"], probs["miss_prob_comb"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",
                len(de_novos[symbol]["missense"]), probs["miss_dist_pois"], probs["miss_prob_pois"]))
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

async def recurrent_test_run(ensembl, mut_dict, output, args):
    N_males = args.N_males
    N_females = args.N_females
    genes = load_gencode(None, args.gencode, args.fasta)
    gene_list = list(genes.get_genes().keys())
    random_genes = random.choices(gene_list,k= args.n)
    #random_genes = random.choices(genes,k=int(args.n))

    de_novos = {}
    for symbol in random_genes:
        try:
            print(symbol)
            gene = genes[symbol]
            transcript = gene.canonical
            rates = SiteRates(transcript, mut_dict, {})
            weights = rates["missense"]

            
            #generate random choice
            choice = weights.choice_with_alleles()
            de_novos[symbol] = {"missense":[], "nonsense":[]}
            de_novos[symbol]["missense"] = [(transcript.get_position_on_chrom(choice["pos"]), choice["alt"])] * args.count
            
            # Recurrent result
            lambda_1 = get_pred_count(choice["prob"],
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            
            p_recurrent = poisson.pmf(args.count, lambda_1)
            lambda_2 = get_pred_count(weights.get_summed_rate(),
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            
            p_count = poisson.pmf(args.count, lambda_2)
            print("p_recurrent = " + str(p_recurrent/p_count))
            ratio =  pow(lambda_1/lambda_2, args.count)
            print("test = "+ str(ratio))

            # Clustering result
            pvalues_in = None
            three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
            (probs,results) = await cluster_de_novos(symbol,
                                                 de_novos[symbol],
                                                 three_d_locations_list,
                                                 genes[symbol],
                                                 float(0),
                                                 1E7-1,#iterations,
                                                 mut_dict,
                                                 None,
                                                 pvalues_in)
            print("p_miss_add_one = " + str(probs["miss_prob"]))


            
            # Clustering result
            pvalues_in = None
            three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
            (probs,results) = await cluster_de_novos(symbol,
                                                 de_novos[symbol],
                                                 three_d_locations_list,
                                                 genes[symbol],
                                                 float(.001),
                                                 1E7-1,#iterations,
                                                 mut_dict,
                                                 None,
                                                 pvalues_in)
            print("p_miss_gener = " + str(probs["miss_prob"]))
            codon = transcript.get_codon_number_for_cds_position(choice["pos"])                                             
            prob = weights.get_prob(codon*3) + weights.get_prob(codon*3+1) + weights.get_prob(codon*3+2)
            
            lambda_1 = get_pred_count(prob,
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            
            lambda_2 = get_pred_count(weights.get_summed_rate(),
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            ratio = pow(prob/weights.get_summed_rate(), args.count)
            
            print("p_mis_comb = " + str(ratio))


            # Add more choices
            
            for i in range(5):
                other_choice = weights.choice_with_alleles()
                de_novos[symbol]["missense"].append((transcript.get_position_on_chrom(other_choice["pos"]), other_choice["alt"]))

            # Clustering result, all sites
            (probs,results) = await cluster_de_novos(symbol,
                                                 de_novos[symbol],
                                                 three_d_locations_list,
                                                 genes[symbol],
                                                 float(.001),
                                                 1E7-1,#iterations,
                                                 mut_dict,
                                                 None,
                                                 pvalues_in)
            print("p_miss = " + str(probs["miss_prob"]))

            
            # Clustering results, unique sites
            # 1 remove duplicates
            # 2 calculate duplicates probability
            miss = de_novos[symbol]["missense"]
            de_novos_counts = {i:miss.count(i) for i in miss}
            de_novos[symbol]["missense"] = list(set(de_novos[symbol]["missense"]))
            
            (probs,results) = await cluster_de_novos(symbol,
                                                 de_novos[symbol],
                                                 three_d_locations_list,
                                                 genes[symbol],
                                                 0,
                                                 1E8-1,#iterations,
                                                 mut_dict,
                                                 None,
                                                 pvalues_in)

            codon = transcript.get_codon_number_for_cds_position(choice["pos"])                                             
            prob = weights.get_prob(codon*3) + weights.get_prob(codon*3+1) + weights.get_prob(codon*3+2)
            
            lambda_1 = get_pred_count(prob,
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            
            lambda_2 = get_pred_count(weights.get_summed_rate(),
                                     N_males,
                                     N_females,
                                     transcript.get_chrom())
            ratio = pow(prob/weights.get_summed_rate(), args.count)
            
            print("p_miss_unique_sites = " + str(probs["miss_prob"]))
            print("ratio = " + str(ratio))
            print("p_miss_adj = " + str(probs["miss_prob"]*ratio))            
        except Exception as e:
            print(e)
            
async def poisson_test(ensembl, mut_dict, output, args):    
    de_novos = load_de_novos(args.input)
    if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
    else:
        return
    print(de_novos)    

    for symbol in sorted(de_novos):
        genes = load_gencode(None, args.gencode, args.fasta)
        gene = genes[symbol]

        transcript = gene.canonical
        print(transcript)
        # transcripts = gene.transcripts
        # minimized =  minimise_transcripts_2(transcripts,de_novos = dnvs)
        # transcripts = [x for x in transcripts if x.get_name() in minimized]

        dnvs = [i[0] for i in de_novos[symbol]["missense"]] + [i[0] for i in de_novos[symbol]["nonsense"]]
        missense = [x[0] for x in de_novos[symbol]["missense"]]
        nonsense = [x[0] for x in de_novos[symbol]["nonsense"]]

        if args.dbNSFP is not None:
            dbNSFP_score_obj = pysam.TabixFile(args.dbNSFP)
        else:
            dbNSFP_score_obj = None

        if args.scores is not None:
            scores_file = pysam.TabixFile(args.scores)
        else:
            scores_file = None
            
        scores = {}
        has_score = 0
        has_no_score = 0
        above_threshold = 0
        below_threshold = 0            
        if dbNSFP_score_obj is not None:
            annotator = args.dbNSFP_annotator
            dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
            if annotator in dbNSFP_header:
                anno_idx = dbNSFP_header.index(annotator)
                alt_idx = dbNSFP_header.index('alt')
                ref_idx = dbNSFP_header.index('ref')                
                pos_idx = dbNSFP_header.index('pos(1-based)')
                print(transcript.get_chrom())
                print(transcript.get_start())
                print(transcript.get_end())
                with open('position_scores_level.txt', 'a') as the_file:
                    for line in dbNSFP_score_obj.fetch(transcript.get_chrom(),
                                                       transcript.get_start()-1,
                                                       transcript.get_end()):
                        values = line.split('\t')
                        pos = values[pos_idx]
                        ref = values[ref_idx]                        
                        alt = values[alt_idx]
                        if values[anno_idx] == '.' :
                            has_no_score = has_no_score + 1
                            #print(pos)
                            #print(alt)
                            continue
                        else:
                            has_score = has_score + 1
                            score = float(values[anno_idx])
                        if int(pos) not in scores:
                            scores[int(pos)] = {}
                        scores[int(pos)][alt.encode('utf-8')] = float(score)

                        if(args.threshold != None):
                            if(score > float(args.threshold)):
                                above_threshold = above_threshold + 1
                            else:
                                below_threshold = below_threshold + 1
                        the_file.write(str(pos) + "\t" + str(ref) + "\t" + str(score) + "\n")                        
                        
                print("has score = " + str(has_score))
                print("has no score = " + str(has_no_score))
                print("above threshold = " + str(above_threshold))
                print("below threshold = " + str(below_threshold))                                                
        elif scores_file is not None:
            with open('position_scores_level.txt', 'w') as the_file:
                for line in scores_file.fetch(transcript.get_chrom(),
                                              transcript.get_start()-1,
                                              transcript.get_end()):
                    pos = line.split('\t')[1]
                    ref = line.split('\t')[2]
                    alt = line.split('\t')[3]
                    score = line.split('\t')[5]                

                    if score == '.' :
                        has_no_score = has_no_score + 1
                        continue
                    else:
                        has_score = has_score + 1
                        score = float(score)
                        if int(pos) not in scores:
                            scores[int(pos)] = {}
                        scores[int(pos)][alt.encode('utf-8')] = float(score)

                        if(args.threshold != None):
                            if(score > float(args.threshold)):
                                above_threshold = above_threshold + 1
                            else:
                                below_threshold = below_threshold + 1
                        the_file.write(str(pos) + "\t" + str(ref) + "\t" + str(score) + "\n")
                print("has score = " + str(has_score))
                print("has no score = " + str(has_no_score))
                print("above threshold = " + str(above_threshold))
                print("below threshold = " + str(below_threshold))                                                


        # Just take top transcript
        for transcript in [transcript]: 
            if scores != {} and args.threshold != None:
                rates = SiteRates(transcript, mut_dict, scores, float(args.threshold))
                missense = [pos for pos,alt in de_novos[symbol]["missense"] if scores[int(pos)][alt.encode('utf-8')] >= float(args.threshold) ]
                print(len(rates["missense"]))
            else:
                rates = SiteRates(transcript, mut_dict, {})


            print("Summed rates = " + str(rates["missense"].get_summed_rate()))
            missense_events = get_de_novos_in_transcript(transcript, missense)
            weights = rates["missense"]

            miss_codons = [transcript.get_codon_number(i) for i in missense_events]
            count = len(missense_events)

            #prob_aa = 0
            #for codon in miss_codons:
            #    prob_aa += weights.get_prob(codon*3) + weights.get_prob(codon*3+1) + weights.get_prob(codon*3+2)

            prob_aa = weights.get_summed_rate()
            
            pois_rate = poisson.pmf(count,
                                     get_pred_count(prob_aa,
                                                    args.N/2.0,
                                                    args.N/2.0,
                                                    transcript.get_chrom())) 
            
            
            print("Poisson test = "+str(pois_rate))
                
        with open(output, "w") as output_file:
            output_file.write("gene_id\tpoisson_p_value\n")
            output_file.write("{}\t{}\n".format(symbol, pois_rate)) 

async def clustering(ensembl, mut_dict, output, args):

    # TODO: make load_de_denovos take args.protein as input
    de_novos = load_de_novos(args.input)
    iterations = 1000000
    # if protein specific run, limit the input variants to this protein/gene
    if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
    else:
        return
    print(de_novos)
    output_file_has_header_line = False
    for symbol in sorted(de_novos):
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
        #if len(de_novos[symbol]["missense"]) < 2:
            with open(output, "x") as output_file:
                output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense", len(de_novos[symbol]["missense"]), "N/A", "N/A"))
            
            continue
        three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)

        if args.scores is not None:
            scores = pysam.TabixFile(args.scores)
        else:
            scores = None

        if args.dbNSFP is not None:
            dbNSFP_score_obj = pysam.TabixFile(args.dbNSFP)
        else:
            dbNSFP_score_obj = None

            
        if args.pvalues_in is not None:
            pvalues_in = load_pvalues(args.pvalues_in)
            if not symbol in pvalues_in:
                pvalues_in = None
        else:
            pvalues_in = None

        genes = load_gencode(None, args.gencode, args.fasta)
        (probs,results) = await cluster_de_novos(symbol,
                                                 de_novos[symbol],
                                                 three_d_locations_list,
                                                 genes[symbol],
                                                 args.threshold,
                                                 float(args.degree),
                                                 args.dist_file_output,
                                                 args.runs,#iterations,
                                                 mut_dict,
                                                 scores,
                                                 dbNSFP_score_obj,
                                                 args.dbNSFP_annotator,
                                                 pvalues_in)
        
        if probs is None:
            continue

        with open(output, "x") as output_file:
            output_file_has_header_line = True
            # output_file.write("gene_id\tmutation_category\tstart \tend \t events_n \tdist \tprobability \n")                
            # for result in results["missense"]:
            #     output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
            #                                                             "missense",
            #                                                             result[0],
            #                                                             result[1],
            #                                                             result[2],
            #                                                             result[3],
            #                                                             result[4]))
            # for result in results["nonsense"]:
            #     output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
            #                                                             "nonsense",
            #                                                             result[0],
            #                                                             result[1],
            #                                                             result[2],
            #                                                             result[3],
            #                                                             result[4]))
            output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            if("miss_prob_comb" in probs):
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",
                len(de_novos[symbol]["missense"]), probs["miss_dist_comb"], probs["miss_prob_comb"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",
                len(de_novos[symbol]["missense"]), probs["miss_dist_pois"], probs["miss_prob_pois"]))
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

    if "inherited" in args and args.inherited == True:
        inherited = []

        # TODO: figure out if I want to have bulk protein runs
        # and if not, just delete the following commented code
        # if "protein" in args:
        #     if args.protein in inherited and args.protein in de_novos:
        #         inherited = {args.protein : inherited[args.protein]}
        #         de_novos = {args.protein : de_novos[args.protein]}
        #     elif args.protein in inherited:
        #         inherited = {args.protein : inherited[args.protein]}
        #     elif args.protein in de_novos:
        #         de_novos = {args.protein : de_novos[args.protein]}
        #     else:
        #         return
        # else:
        #     return

        #Load inherited variants
        transcripts = await get_transcript_ids(ensembl, args.protein)
        min_start = float("inf")
        max_end = 0
        for transcript in transcripts:
            gene = await construct_gene_object(ensembl, transcript)            
            chrom = gene.get_chrom()
            if min_start > gene.get_start():
                min_start = gene.get_start()
            
            if max_end < gene.get_end():
                max_end = gene.get_end()

        # TODO: once UKBB is done annotating, make inherited controls
        # capable of picking only missense and/or LoF
        inherited_controls = {"missense": load_inherited_controls(chrom,
                                                                   min_start,
                                                                   max_end)
                              }
        
        inherited = {args.protein : {'missense' : load_inherited(chrom,
                                                                 min_start,
                                                                 max_end,
                                                                 ["missense_variant"]),
                                     'nonsense' : load_inherited(chrom,
                                                            min_start,
                                                            max_end,
                                                            ["stop_gained",
                                                             "stop_lost",
                                                             "start_lost",
                                                             "frameshift_variant",
                                                             "splice_acceptor_variant",
                                                             "splice_donor_variant",
                                                             "transcript_ablation"])}
                         }

    
        for symbol in sorted(inherited):
            if len(inherited[symbol]["missense"] + inherited[symbol]["nonsense"]) < 2:
                continue
            three_d_locations_list = await load_three_d_locations(args.protein_dir, symbol)
            genes = load_gencode(None, args.gencode, args.fasta)            
            probs = await cluster_de_novos(symbol,
                                           inherited[symbol],
                                           three_d_locations_list,
                                           genes[symbol],
                                           args.degree,
                                           args.dist_file_output,
                                           args.runs,#iterations,
                                           mut_dict,
                                           "inherited",
                                           inherited_controls)
        
            if probs is None:
                continue
            with open(output, "x") as output_file:
                if output_file_has_header_line == False:
                    output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")    
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "inherited_missense",
                                                                len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "inherited_nonsense",
                                                                len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

async def entropy(ensembl, mut_dict, output, args):
    # TODO: make load_de_denovos take args.protein as input
    de_novos = load_de_novos(args.input)
    iterations = 1000000
    # if protein specific run, limit the input variants to this protein/gene
    if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
    else:
        return
    print(de_novos)
    output_file_has_header_line = False
    for symbol in sorted(de_novos):
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            with open(output, "x") as output_file:
                output_file.write("gene_id\tmutation_category\tevents_n\tentropy\tprobability\n")
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense", len(de_novos[symbol]["missense"]), "N/A", "N/A"))
            continue

        genes = load_gencode(None, args.gencode, args.fasta)
        (probs,results) = await de_novos_entropy(symbol,
                                                 de_novos[symbol],
                                                 genes[symbol],
                                                 float(args.degree),
                                                 args.dist_file_output,
                                                 args.runs,#iterations,
                                                 mut_dict)
        if probs is None:
            continue

        with open(output, "x") as output_file:
            output_file_has_header_line = True
            output_file.write("gene_id\tmutation_category\tevents_n\tentropy\tprobability\n")
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))


async def clustering_multi(ensembl, mut_dict, output, args):

    assert len(args.chains) == len(args.proteins)
    
    # TODO: make load_de_denovos take args.protein as input
    all_de_novos = load_de_novos(args.input)
    de_novos = {}
    three_d_locations = []
    proteins_to_chains = {}
    iterations = args.runs

    three_d_locations += await load_three_d_multimer(args.multimer_dir, args.multimer)
    
    if "proteins" in args:
        for i in range(len(args.proteins)):
            if args.proteins[i] not in proteins_to_chains:
                proteins_to_chains[args.proteins[i]] = []
            proteins_to_chains[args.proteins[i]] += args.chains[i]
            if args.proteins[i] in all_de_novos:
                de_novos[args.proteins[i]] = all_de_novos[args.proteins[i]]
    else:
        return

    if args.scores is not None:
        scores = pysam.TabixFile(args.scores)
    else:
        scores = None
    if args.pvalues_in is not None:
        pvalues_in = load_pvalues(args.pvalues_in)
        for protein in proteins:
            if not protein in pvalues_in:
                pvalues_in = None
                continue
    else:
        pvalues_in = None

    print(proteins_to_chains)
    (probs,results) = await cluster_de_novos_multi(args.multimer,
                                                   args.chains,
                                                   args.proteins,
                                                   proteins_to_chains,
                                                   de_novos,
                                                   three_d_locations,
                                                   ensembl,
                                                   float(args.degree),
                                                   args.dist_file_output,
                                                   args.runs,#iterations,
                                                   mut_dict,
                                                   scores,
                                                   pvalues_in)
    os.rename('current_distribution.txt',  args.multimer + ".missense.multi.dist.txt")
    if probs is None:
        return

    with open(output, "x") as output_file:
        output_file_has_header_line = True
        # output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")


        symbol = ""
        for p in args.proteins:
            symbol += p+"_"
        output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
        output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
                                                        len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
        if("miss_prob_comb" in probs):
                    output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",
                                                        len(de_novos[symbol]["missense"]), probs["miss_dist_comb"], probs["miss_prob_comb"]))
                    output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",
                len(de_novos[symbol]["missense"]), probs["miss_dist_pois"], probs["miss_prob_pois"]))
        output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
                                                        len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))
        # output_file.write("gene_id\tmutation_category\tstart \tend \t events_n \tdist \tprobability \n")
        # for result in results["missense"]:
        #     output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
        #                                                             "missense",
        #                                                             result[0],
        #                                                             result[1],
        #                                                             result[2],
        #                                                             result[3],
        #                                                             result[4]))
        #     for result in results["nonsense"]:
        #         output_file.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(symbol,
        #                                                                 "nonsense",
        #                                                                 result[0],
        #                                                                 result[1],
        #                                                                 result[2],
        #                                                                 result[3],
        #                                                                 result[4]))


                
async def clustering_coev(ensembl, mut_dict, output, args):
    
    de_novos = load_de_novos(args.input)
    
    
    iterations = 1000000

    if "protein" in args:
        if args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
        else:
            return
    else:
        return
    
    for symbol in sorted(de_novos):
        if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
            continue
        
        coev_relations = await load_basic_coevol_strength(args.coev_dir, symbol)
        genes = load_gencode(None, args.gencode, args.fasta)
        probs = await cluster_de_novos_coevol(symbol,
                                              de_novos[symbol],
                                              coev_relations,
                                              genes[symbol],#ensembl,
                                              args.degree,
                                              args.dist_file_output,
                                              args.runs,#iterations,
                                              mut_dict)
        
        if probs is None:
            continue

        with open(output, "x") as output_file:
            output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")    
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            if("miss_dist_comb" in probs):
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",
                len(de_novos[symbol]["missense"]), probs["miss_dist_comb"], probs["miss_prob_comb"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",
                len(de_novos[symbol]["missense"]), probs["miss_dist_pois"], probs["miss_prob_pois"]))
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))

            
async def clustering_1d(ensembl, mut_dict, output, args):
    
    de_novos = load_de_novos(args.input)
    
    iterations = 1000000

    if "protein" in args:
        if args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
        else:
            return
    else:
        return
    
    for symbol in sorted(de_novos):
        #if len(de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]) < 2:
        if len(de_novos[symbol]["missense"]) < 2:
            with open(output, "x") as output_file:
                output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense", len(de_novos[symbol]["missense"]), "N/A", "N/A"))
            continue

        if args.scores is not None:
            scores = pysam.TabixFile(args.scores)
        else:
            scores = None
        if args.pvalues_in is not None:
            pvalues_in = load_pvalues(args.pvalues_in)
            if not symbol in pvalues_in:
                pvalues_in = None
        else:
            pvalues_in = None
        
        probs = await cluster_de_novos_1d(symbol,
                                          de_novos[symbol],
                                          ensembl,
                                          args.threshold,
                                          args.degree,
                                          args.dist_file_output,
                                          args.runs,#iterations,
                                          mut_dict,
                                          scores,
                                          pvalues_in)
        
        if probs is None:
            continue

        
        with open(output, "x") as output_file:
            output_file.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")    
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense",
            len(de_novos[symbol]["missense"]), probs["miss_dist"], probs["miss_prob"]))
            if("miss_dist_comb" in probs):
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_comb",
                len(de_novos[symbol]["missense"]), probs["miss_dist_comb"], probs["miss_prob_comb"]))
                output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense_pois",
                len(de_novos[symbol]["missense"]), probs["miss_dist_pois"], probs["miss_prob_pois"]))
            output_file.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "nonsense",
            len(de_novos[symbol]["nonsense"]), probs["nons_dist"], probs["nons_prob"]))


async def find_transcripts(ensembl, mut_dict, output, args):
    
    de_novos = load_de_novos(args.de_novos)
    
    output.write("hgnc_symbol\ttranscript_id\tlength\tde_novos\n")
    
    for symbol in sorted(de_novos):
        print(symbol)
        func_events = de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]
        
        # find the counts per transcript, depending on whether we want to count
        # for all transcripts containing one or more de novos, or to find the
        # minimum set of transcripts to contain the de novos
        try:
            if args.all_transcripts:
                counts = await count_de_novos_per_transcript(ensembl, symbol, func_events)
            elif args.minimal_transcripts:
                counts = await minimise_transcripts(ensembl, symbol, func_events)
        except (ValueError, IndexError):
            print("error occured with {0}".format(symbol))
            continue
        
        # write the transcript details to a file
        for key in counts:
            line = "{}\t{}\t{}\t{}\n".format(symbol, key, counts[key]["len"],
                counts[key]["n"])
            output.write(line)

def load_genes(path):
    """ load a file listing gene and transcript IDs
    
    Args:
        path: path to file containing gene IDs and transcript IDs e.g.
            gene_1    transcript_1.1    length_1    denovo_count
            gene_2    transcript_2.1    length_3    denovo_count
    
    Returns:
        dict of transcripts eg {'CTC1': ["ENST00000315684", "ENST00000485511"]}
    """
    
    with open(path, 'rt') as f:
        lines = [ x.split('\t')[:2] for x in f if not x.startswith('hgnc') ]
    
    transcripts = {}
    for symbol, tx in lines:
        if symbol not in transcripts:
            transcripts[symbol] = []
        transcripts[symbol].append(tx)
    
    return transcripts

async def get_mutation_rates(transcripts, mut_dict, ensembl):
    """ determines mutation rates per functional category for transcripts
    
    Args:
        transcripts: list of transcript IDs for a gene
        mut_dict: dictionary of local sequence context mutation rates
        ensembl: EnsemblRequest object, to retrieve information from Ensembl.
    
    Returns:
        tuple of (rates, merged transcript, and transcript CDS length)
    """
    
    rates = {'missense': 0, 'nonsense': 0, 'splice_lof': 0,
        'splice_region': 0, 'synonymous': 0}
    combined = None
    
    for tx_id in transcripts:
        try:
            tx = await construct_gene_object(ensembl, tx_id)
        except ValueError:
            continue
        
        if len(tx.get_cds_sequence()) % 3 != 0:
            raise ValueError("anomalous_coding_sequence")
        
        # ignore mitochondrial genes
        if tx.get_chrom() == "MT":
            continue
        
        sites = SiteRates(tx, mut_dict, masked_sites=combined)
        combined = tx + combined
        
        for cq in ['missense', 'nonsense', 'splice_lof', 'splice_region', 'synonymous']:
            rates[cq] += sites[cq].get_summed_rate()
    
    if combined is None:
        raise ValueError('no tx found')
    
    length = combined.get_coding_distance(combined.get_cds_end())['pos']
    
    return rates, combined, length

async def gene_rates(ensembl, mut_dict, output, args):
    
    transcripts = load_genes(args.genes)
    
    header = ['transcript_id', 'chrom', 'length', 'missense_rate', 'nonsense_rate',
        'splice_lof_rate', 'splice_region_rate', 'synonymous_rate']
    output.write('\t'.join(header) + "\n")
    
    for symbol in sorted(transcripts):
        print(symbol)
        try:
            rates, tx, length = await get_mutation_rates(transcripts[symbol],
                mut_dict, ensembl)
            # log transform rates, for consistency with Samocha et al.
            line = "{}\t{}\t{}\t{}".format(symbol, tx.get_chrom(), length, log_transform(rates))
        except (ValueError, KeyError) as error:
            print("{}\t{}\n".format(symbol, error))
            line = "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA".format(symbol, tx.get_chrom())
        
        output.write(line + '\n')
    
    output.close()
    include_frameshift_rates(args.out)

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description='denovonear cli interface')
    
    ############################################################################
    # CLI options in common
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("--out", help="output filename")
    parent.add_argument("--rates",
        help="Path to file containing sequence context-based mutation rates.")
    parent.add_argument("--gencode",
        help="optional path to gencode annotations file. If not provided, gene " \
                        "coordinates will be obtained via the Ensembl REST API.")
    parent.add_argument("--fasta",
        help="optional path to genome fasta file. If not provided, gene " \
            "coordinates will be obtained via the Ensembl REST API.")    
    parent.add_argument("--genome-build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch37", help="Genome build "
        "that the de novo coordinates are based on (GRCh37 or GRCh38")
    parent.add_argument("--log", default='ensembl_requests.log', help="where to write log files")
    
    subparsers = parser.add_subparsers()

    ############################################################################
    # CLI options for weight creation
    recurrent_test = subparsers.add_parser('recurrent_test', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    recurrent_test.add_argument("--males", dest="N_males", type = int, required=True)
    recurrent_test.add_argument("--females", dest="N_females", type = int, required=True)
    recurrent_test.add_argument("--n", dest="n", type = int, required=True)
    recurrent_test.add_argument("--count", dest="count", type = int, required=True)
    recurrent_test.add_argument("--protein_dir", dest="protein_dir",  required=True)            

    recurrent_test.set_defaults(func=recurrent_test_run)

    ############################################################################
    # CLI options for weight creation
    west_exp_buckets = subparsers.add_parser('west_exp_buckets', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    west_exp_buckets.add_argument("--out_dir", dest="output_dir")
    west_exp_buckets.add_argument("--males", dest="N_males", type = int, required=True)
    west_exp_buckets.add_argument("--females", dest="N_females", type = int, required=True)    
    west_exp_buckets.add_argument("--gMVP_file", dest="gMVP_file")
    west_exp_buckets.add_argument("--dbNSFP_file", dest="dbNSFP_file")
    west_exp_buckets.add_argument("--MCR_file", dest="MCR_file")
    west_exp_buckets.add_argument("--dbNSFP_gene_file", dest="dbNSFP_gene_file")
    west_exp_buckets.add_argument("--gene", dest="gene")    
    west_exp_buckets.add_argument("--bucket_ann",  nargs='+', dest="bucket_ann")
    west_exp_buckets.add_argument("--bucket_cutoffs",  nargs='+', dest="bucket_cutoffs")
    west_exp_buckets.set_defaults(func=west_exp_buckets_creation)
    
    ############################################################################
    west_obs_buckets = subparsers.add_parser('west_obs_buckets', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    west_obs_buckets.add_argument("--in", dest="input", required=True, help="Path to file listing known mutations in genes. See example file in data folder for format.")
    west_obs_buckets.add_argument("--out_dir", dest="output_dir")
    west_obs_buckets.add_argument("--gMVP_file", dest="gMVP_file")
    west_obs_buckets.add_argument("--dbNSFP_file", dest="dbNSFP_file")
    west_obs_buckets.add_argument("--MCR_file", dest="MCR_file")
    west_obs_buckets.add_argument("--dbNSFP_gene_file", dest="dbNSFP_gene_file")

    west_obs_buckets.add_argument("--bucket_ann",  nargs='+', dest="bucket_ann")
    west_obs_buckets.add_argument("--bucket_cutoffs",  nargs='+', dest="bucket_cutoffs")
    west_obs_buckets.set_defaults(func=west_obs_buckets_creation)
    ############################################################################
    west_scores = subparsers.add_parser('west_scores', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    west_scores.add_argument("--obs_buckets", dest="obs_buckets", required=True, help="Path to file listing known mutations in genes. See example file in data folder for format.")

    west_scores.add_argument("--buckets_dir", dest="buckets_dir", required=True)
    west_scores.add_argument("--out_dir", dest="output_dir")

    west_scores.add_argument("--bucket_ann",  nargs='+', dest="bucket_ann", required=True)
    west_scores.add_argument("--bucket_cutoffs",  nargs='+', dest="bucket_cutoffs", required=True)
    west_scores.set_defaults(func=west_scores_creation)
    ############################################################################
    # CLI options for west
    west = subparsers.add_parser('west', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    west.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    west.add_argument("--protein_dir", dest="protein_dir", required=True)
    west.add_argument("--scores", dest="scores")
    west.add_argument("--runs", dest="runs", type=float, default = 1E6)    
    west.add_argument("--protein", dest="protein", required=False)
    west.add_argument("--inherited", dest="inherited", required=False)
    west.add_argument("--males", dest="N_males", type = int, required=True)
    west.add_argument("--females", dest="N_females", type = int, required=True)    
    west.add_argument("--pvalues_in", dest="pvalues_in", required=False)
    west.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    west.set_defaults(func=west_run)

    ############################################################################
    # CLI options for type_II_testing
    type_II_test = subparsers.add_parser('type_II_test', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    type_II_test.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    type_II_test.add_argument("--protein_dir", dest="protein_dir", required=True)
    type_II_test.add_argument("--scores", dest="scores")
    type_II_test.add_argument("--dbNSFP", dest="dbNSFP")
    type_II_test.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")    
    type_II_test.add_argument("--runs", dest="runs", type=float, default = 1E6)
    type_II_test.add_argument("--cohort_size", dest="cohort_size", type=float)        
    type_II_test.add_argument("--protein", dest="protein", required=False)
    type_II_test.add_argument("--pvalues_in", dest="pvalues_in", required=False)
    type_II_test.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    type_II_test.set_defaults(func=type_II_testing)

    
    ############################################################################
    # CLI options for clustering
    cluster = subparsers.add_parser('cluster', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    cluster.add_argument("--protein_dir", dest="protein_dir", required=True)
    cluster.add_argument("--scores", dest="scores")
    cluster.add_argument("--dbNSFP", dest="dbNSFP")
    cluster.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")
    cluster.add_argument("--threshold", dest="threshold", required=False)
    cluster.add_argument("--runs", dest="runs", type=float, default = 1E6)    
    cluster.add_argument("--protein", dest="protein", required=False)
    cluster.add_argument("--inherited", dest="inherited", required=False)
    cluster.add_argument("--pvalues_in", dest="pvalues_in", required=False)
    cluster.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    cluster.add_argument("--histogram_data", dest="dist_file_output", required=False)
    cluster.set_defaults(func=clustering)

    ############################################################################
    # CLI options for entropy_test
    entropy_test = subparsers.add_parser('entropy_test', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    entropy_test.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    entropy_test.add_argument("--scores", dest="scores")
    entropy_test.add_argument("--dbNSFP", dest="dbNSFP")
    entropy_test.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")    
    entropy_test.add_argument("--runs", dest="runs", type=float, default = 1E6)    
    entropy_test.add_argument("--protein", dest="protein", required=False)
    entropy_test.add_argument("--inherited", dest="inherited", required=False)
    entropy_test.add_argument("--pvalues_in", dest="pvalues_in", required=False)
    entropy_test.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    entropy_test.add_argument("--histogram_data", dest="dist_file_output", required=False)    
    entropy_test.set_defaults(func=entropy)

    ############################################################################
    # CLI options for poisson_test
    poisson = subparsers.add_parser('poisson', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    poisson.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    poisson.add_argument("--N", dest="N", type=float)
    poisson.add_argument("--protein", dest="protein", required=False)
    poisson.add_argument("--scores", dest="scores")    
    poisson.add_argument("--dbNSFP", dest="dbNSFP")
    poisson.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")
    poisson.add_argument("--threshold", dest="threshold", required=False)
    # poisson.add_argument("--females", dest="N_females", type=float, required=False)
    # poisson.add_argument("--males", dest="N_males", type=float, required=False)        
    poisson.set_defaults(func=poisson_test)    
    ############################################################################
    # CLI options for clustering
    cluster_multi = subparsers.add_parser('cluster_multi', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster_multi.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    cluster_multi.add_argument("--multimer_dir", dest="multimer_dir", required=True)
    cluster_multi.add_argument("--multimer", dest="multimer", required=True)    
    cluster_multi.add_argument("--scores", dest="scores")
    cluster_multi.add_argument("--chains",  nargs='+', dest="chains")
    cluster_multi.add_argument("--runs", dest="runs", type=float,default = 1E6)    
    cluster_multi.add_argument("--proteins", nargs='+', dest="proteins", required=False)
    cluster_multi.add_argument("--inherited", dest="inherited", required=False)
    cluster_multi.add_argument("--pvalues_in", dest="pvalues_in", required=False)    
    cluster_multi.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    cluster_multi.add_argument("--histogram_data", dest="dist_file_output", required=False)
    cluster_multi.set_defaults(func=clustering_multi)

    ############################################################################
    # CLI options for clustering_coev
    cluster_coev = subparsers.add_parser('cluster_coev', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster_coev.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    cluster_coev.add_argument("--coev_dir", dest="coev_dir", required=True)
    cluster_coev.add_argument("--runs", dest="runs", type=float, default = 1E6)        
    cluster_coev.add_argument("--protein", dest="protein", required=False)
    cluster_coev.add_argument("--pvalues_in", dest="pvalues_in", required=False)    
    cluster_coev.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    cluster_coev.add_argument("--histogram_data", dest="dist_file_output", required=False)
    cluster_coev.set_defaults(func=clustering_coev)
    

    
    ############################################################################
    # CLI options for clustering_1d
    cluster_1d = subparsers.add_parser('cluster_1d', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster_1d.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    cluster_1d.add_argument("--protein", dest="protein", required=False)
    cluster_1d.add_argument("--threshold", dest="threshold", required=False)    
    cluster_1d.add_argument("--runs", dest="runs", type = float, default = 1E6)
    cluster_1d.add_argument("--pvalues_in", dest="pvalues_in", required=False)    
    cluster_1d.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    cluster_1d.add_argument("--scores", dest="scores")
    cluster_1d.add_argument("--histogram_data", dest="dist_file_output", required=False)    
    cluster_1d.set_defaults(func=clustering_1d)
    

    ############################################################################
    # CLI options for identifing transcripts to use
    transcripts = subparsers.add_parser('transcripts', parents=[parent],
        description="Identify transcripts for a gene containing de novo events.")
    transcripts.add_argument("--de-novos", required=True, help="Path to "
        "file listing de novo variants in genes.")
    transcripts.set_defaults(func=find_transcripts)
    group = transcripts.add_mutually_exclusive_group(required=True)
    group.add_argument("--all-transcripts", action="store_true", default=False,
        help="Flag if you want to identify all transcripts with more than "
        "one de novo on it.")
    group.add_argument("--minimise-transcripts", action="store_true",
        default=False, help="Flag if you want to identify the minimal set of "
        "transcripts to contain all de novos.")
    
    ############################################################################
    # CLI options for identifing transcripts to use
    rater = subparsers.add_parser("rates", parents=[parent],
        description="determine mutation rates for genes given transcript IDs.")
    rater.add_argument("--genes", help="Path to file "
        "listing HGNC symbols, with one or more transcript IDs per gene. "
        "The tab-separated input format is gene symbol followed by transcript "
        "ID. Alternative transcripts are listed on separate lines.")
    rater.set_defaults(func=gene_rates)
    
    args = parser.parse_args()
    if 'func' not in args:
        print('Use one of the subcommands: cluster, rates, or transcripts\n')
        parser.print_help()
        sys.exit()
    
    return args

async def runner():
    args = get_options()
    FORMAT = '%(asctime)-15s %(message)s'
    logging.basicConfig(filename=args.log, format=FORMAT, level=logging.INFO)
    
    async with RateLimiter(per_second=15) as ensembl:
        mut_dict = load_mutation_rates(args.rates)
        await args.func(ensembl, mut_dict, args.out, args)
        #with open(args.out, "wt") as output:
        #    await args.func(ensembl, mut_dict, output, args)

def main():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(runner())

if __name__ == '__main__':
    main()
