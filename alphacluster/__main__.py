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
from math import floor, log
from difflib import SequenceMatcher

from alphacluster.rate_limiter import RateLimiter
from alphacluster.load_mutation_rates import load_mutation_rates
from alphacluster.load_de_novos import load_de_novos, load_de_novos_chrom_pos_alt
from alphacluster.load_inherited import load_inherited, load_inherited_controls
from alphacluster.west_weights import generate_expected_buckets, generate_observed_buckets, generate_scores_from_buckets, get_pred_count
from alphacluster.cluster_test import cluster_de_novos, cluster_de_novos_1d, cluster_de_novos_coevol, cluster_de_novos_multi, fishers_method
from alphacluster.load_three_d_locations import load_three_d_locations, load_three_d_multimer, load_three_d_locations_from_pdb
from alphacluster.load_coevol import load_basic_coevol_strength
from alphacluster.load_pvalues import load_pvalues
from alphacluster.load_gene import (load_gene, construct_gene_object,
                                  count_de_novos_per_transcript, minimise_transcripts, get_transcript_ids,  get_de_novos_in_transcript)
from alphacluster.simulate import get_p_value_west

from alphacluster.site_specific_rates import SiteRates
from alphacluster.frameshift_rate import include_frameshift_rates
from alphacluster.log_transform_rates import log_transform
from alphacluster.gencode import Gencode

import operator as op
from functools import reduce

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

def find_transcripts(args, output):
    
    mut_dict = load_mutation_rates(args.rates)
    de_novos = load_de_novos(args.de_novos)
    gencode = load_gencode(de_novos, args.gencode, args.fasta)
        
    output.write("hgnc_symbol\ttranscript_id\tlength\tde_novos\n")
    
    for symbol in sorted(de_novos):
        logging.info(f'checking {symbol}')
        func_events = de_novos[symbol]["missense"] + de_novos[symbol]["nonsense"]
        
        transcripts = gencode[symbol].transcripts
        
        # find the counts per transcript, depending on whether we want to count
        # for all transcripts containing one or more de novos, or to find the
        # minimum set of transcripts to contain the de novos
        try:
            if args.all_transcripts:
                counts = count_de_novos_per_transcript(transcripts, func_events)
            elif args.minimal_transcripts:
                counts = minimise_transcripts(transcripts, func_events)
        except (ValueError, IndexError):
            logging.error(f"error occured with {symbol}")
            continue
        
        # write the transcript details to a file
        for key in counts:
            line = "{}\t{}\t{}\t{}\n".format(symbol, key, counts[key]["len"],
                counts[key]["n"])
            output.write(line)

    
def poisson_test(output, args):
    de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)
    if args.protein in de_novos:
        de_novos = {args.protein : de_novos[args.protein]}
    else:
        return

    for symbol in sorted(de_novos):
        genes = load_gencode([symbol], args.gencode, args.fasta)
        gene = genes[symbol]

        transcript = gene.canonical

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
                with open('position_scores_level.txt', 'a') as the_file:
                    for line in dbNSFP_score_obj.fetch(transcript.get_chrom()[3:],
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
                for line in scores_file.fetch(transcript.get_chrom()[3:],
                                              transcript.get_start()-1,
                                              transcript.get_end()-1):
                    pos = line.split('\t')[1]
                    ref = line.split('\t')[2]
                    alt = line.split('\t')[3]
                    score = line.split('\t')[int(args.score_col)]                

                    if score == '.' :
                        has_no_score = has_no_score + 1
                        continue
                    else:
                        has_score = has_score + 1
                        score = float(score)
                        if int(pos) not in scores:
                            scores[int(pos)] = {}
                        scores[int(pos)][ord(alt)] = float(score)

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
#                scores = {9:{9:10}}
#                print(scores)
#                print(type(scores))
#                print(type(scores[9]))
                missense_events = get_de_novos_in_transcript(transcript, [pos for pos, alt in de_novos[symbol]["missense"]])
                rates = SiteRates(transcript, mut_dict, scores, float(args.threshold))
                try:
                    missense = [pos for pos,alt in de_novos[symbol]["missense"] if pos in missense_events and  scores[int(pos)][ord(alt)] >= float(args.threshold)]
                except KeyError:
                    print("No score exists for this position")
                except:
                    print("Error while loading scores for missense positions")
            else:
                rates = SiteRates(transcript, mut_dict, {})

            missense_events = get_de_novos_in_transcript(transcript, missense)
            weights = rates["missense"]

            miss_codons = [transcript.get_codon_info(i)["codon_number"] for i in missense_events]
            count = len(missense_events)

            prob_aa = weights.get_summed_rate()
            
            pois_rate = poisson.pmf(count,
                                     get_pred_count(prob_aa,
                                                    args.N/2.0,
                                                    args.N/2.0,
                                                    transcript.get_chrom()[3:])) 
            
            
            print("Number of missense variants = " + str(len(rates["missense"])))
            print("Summed rates = " + str(rates["missense"].get_summed_rate()))
            print("Count = " + str(count))
            print("Poisson test = " + str(pois_rate))
                
        #with open(output, "w") as output_file:
        #    output_file.write("gene_id\tpoisson_p_value\n")
        #    output_file.write("{}\t{}\n".format(symbol, pois_rate)) 
        output.write("gene_id\tpoisson_p_value\n")
        output.write("{}\t{}\n".format(symbol, pois_rate)) 

def poisson_multi_test(output, args):
    all_de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)
    de_novos = {}
    count = 0
    pred_count = 0
    for protein in args.proteins:
        print(protein)
        if protein in all_de_novos:
            de_novos = {protein : all_de_novos[protein]}
        else:
            continue
        print(de_novos)    

        for symbol in sorted(de_novos):
            genes = load_gencode([symbol], args.gencode, args.fasta)
            gene = genes[symbol]

            transcript = gene.canonical
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
                    with open('position_scores_level.txt', 'a') as the_file:
                        for line in dbNSFP_score_obj.fetch(transcript.get_chrom()[3:],
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
                    for line in scores_file.fetch(transcript.get_chrom()[3:],
                                                  transcript.get_start()-1,
                                                  transcript.get_end()-1):
                        pos = line.split('\t')[1]
                        ref = line.split('\t')[2]
                        alt = line.split('\t')[3]
                        score = line.split('\t')[int(args.score_col)]                

                        if score == '.' :
                            has_no_score = has_no_score + 1
                            continue
                        else:
                            has_score = has_score + 1
                            score = float(score)
                            if int(pos) not in scores:
                                scores[int(pos)] = {}
                            scores[int(pos)][ord(alt)] = float(score)

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
                    #                scores = {9:{9:10}}
                    #                print(scores)
                    #                print(type(scores))
                    #                print(type(scores[9]))
                    rates = SiteRates(transcript, mut_dict, scores, float(args.threshold))
                    missense = [pos for pos,alt in de_novos[symbol]["missense"] if scores[int(pos)][ord(alt)] >= float(args.threshold) ]
                else:
                    rates = SiteRates(transcript, mut_dict, {})

                missense_events = get_de_novos_in_transcript(transcript, missense)
                weights = rates["missense"]

                miss_codons = [transcript.get_codon_info(i)["codon_number"] for i in missense_events]
                count = count + len(missense_events)
                prob_aa =  weights.get_summed_rate()
                pred_count = pred_count + get_pred_count(prob_aa,
                                           args.N/2.0,
                                           args.N/2.0,
                                           transcript.get_chrom()[3:])
            
                print("Number of missense variants = " + str(len(rates["missense"])))
                print("Summed rates = " + str(rates["missense"].get_summed_rate()))
                print("Count = " + str(len(missense_events)))

        # Total pois rate    
    pois_rate = poisson.pmf(count,pred_count)
    print("Total count = " + str(count))
    print("Total pred = " + str(pred_count))
    print("Poisson test = " + str(pois_rate))
                
    #with open(output, "w") as output_file:
    #    output_file.write("gene_id\tpoisson_p_value\n")
    #    output_file.write("{}\t{}\n".format(symbol, pois_rate)) 
    #    output.write("gene_id\tpoisson_p_value\n")
    #    output.write("{}\t{}\n".format(symbol, pois_rate)) 

        
def genomic_to_residue_position(output, args):
    three_d_locations_list =  load_three_d_locations(args.protein_dir, args.protein)
    three_d_aa = load_three_d_locations_from_pdb(args.protein_dir, args.protein)
    genes = load_gencode([args.protein], args.gencode, args.fasta)
    print("Length is " + str(len(three_d_aa)))

    
    try:
        transcripts = [genes[args.protein].canonical]
    except IndexError as e:
        print(e)
        return None

    for transcript in transcripts:
        #        print(transcript.get_name())
        t_aa = "".join(list(transcript.translate(transcript.get_cds_sequence())))
        p_aa = "".join([row[3] for row in three_d_aa])
        #        print(len(t_aa))
        #        print(len(p_aa))
        #        print("Transcript")
        #        print(t_aa)
        #        print("Uniprot")
        #        print(p_aa)
        try:
            #assert p_aa in t_aa
            print("XO")
            #Find longest common substring
            match = SequenceMatcher(None, p_aa, t_aa).find_longest_match()
            lcs = p_aa[match.a:match.a + match.size]
            print("Longest Common Substring")
            print(lcs)
            
            #find offsets
            #start_offset = t_aa.find(p_aa)
            #end_offset = t_aa[::-1].find(p_aa[::-1])-1
            start_offset = t_aa.find(lcs)
            end_offset = t_aa[::-1].find(lcs[::-1])-1

            print(start_offset)
            print(end_offset)
            if(start_offset != 0):
                print("And then")
                transcript.set_cds_start(transcript.get_cds_start(),start_offset*3)
            if(end_offset != 0):
                transcript.set_cds_end(transcript.get_cds_end(),end_offset*-3)                
        except:
            print("AlphaFold and transcript disagree for transcript " + transcript.get_name())
            output.write("gene_id\tresidue\n")
            output.write(args.protein + "\t" + "disagreement")
            continue

        de_novos = load_de_novos(args.input)
        # if protein specific run, limit the input variants to this protein/gene
        if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
        else:
            return
        print(de_novos)

        missense = [pos for pos,alt in de_novos[args.protein]["missense"]]
        nonsense = [pos for pos,alt in de_novos[args.protein]["nonsense"]]

        missense_events = get_de_novos_in_transcript(transcript, missense)
        nonsense_events = get_de_novos_in_transcript(transcript, nonsense)

        miss_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in missense_events ]
        nons_cds_positions = [ transcript.get_coding_distance(x)['pos'] for x in nonsense_events ]        

        miss_aa_positions = [floor(float(x)/3)+1 for x in miss_cds_positions]
        nons_aa_positions = [floor(float(x)/3)+1 for x in nons_cds_positions]        

        print(sorted(miss_aa_positions))
        print(sorted(nons_aa_positions))
        output.write("gene_id\tresidue\n")
        for aa in miss_aa_positions:
            output.write(args.protein +"\t" + str(aa) + "\n")
        
def clustering(output, args):
    de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)

    #set default iterations
    iterations = 1000000
    
    # if protein specific run, limit the input variants to this protein/gene
    if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
    else:
        return

    for symbol in sorted(de_novos):
        # Handle case where there are less than two de novo missense variants
        if len(de_novos[symbol]["missense"]) < 2:
            output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense", len(de_novos[symbol]["missense"]), "N/A", "N/A"))
            continue

        three_d_aa = load_three_d_locations_from_pdb(args.protein_dir, symbol)

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

        genes = load_gencode([symbol], args.gencode, args.fasta)

        (probs,results) =  cluster_de_novos(symbol,
                                            de_novos[symbol],
                                            three_d_aa,
                                            genes[symbol],
                                            args.scale,
                                            args.threshold,
                                            float(args.degree),
                                            args.dist_file_output,
                                            args.runs,#iterations,
                                            mut_dict,
                                            scores,
                                            args.score_col,
                                            dbNSFP_score_obj,
                                            args.dbNSFP_annotator,
                                            pvalues_in)
        
        if probs is None:
            # If there is a non-return error, don't print a results file
            continue

        output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                   "missense",
                                                   probs["miss_count"],
                                                   probs["miss_dist"],
                                                   probs["miss_prob"]))
        if("miss_prob_comb" in probs):
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                       "missense_comb",
                                                       probs["miss_count"],
                                                       probs["miss_dist_comb"],
                                                       probs["miss_prob_comb"]))
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                       "missense_pois",
                                                       probs["miss_count"],
                                                       probs["miss_dist_pois"],
                                                       probs["miss_prob_pois"]))


def clustering_multi(output, args):

    assert len(args.chains) == len(args.proteins)
    
    # TODO: make load_de_denovos take args.protein as input
    all_de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)
    de_novos = {}
    three_d_locations = []
    proteins_to_chains = {}
    iterations = args.runs

    #Load locations of amino acids from the multimer file
    #three_d_locations +=  load_three_d_multimer(args.multimer_dir, args.multimer)

    three_d_locations = load_three_d_locations_from_pdb(args.multimer_dir, args.multimer)
    print(three_d_locations)
    
    #Create the proteins_to_chains map (e.g. GABRA2 -> A, GABRB1 -> B, ...)
    if "proteins" in args:
        for i in range(len(args.proteins)):
            if args.proteins[i] not in proteins_to_chains:
                proteins_to_chains[args.proteins[i]] = []
            proteins_to_chains[args.proteins[i]] += args.chains[i]
            if args.proteins[i] in all_de_novos:
                de_novos[args.proteins[i]] = all_de_novos[args.proteins[i]]
    else:
        return
    print(proteins_to_chains)
    
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
    
    genes = load_gencode(args.proteins, args.gencode, args.fasta)
    (probs,results) =  cluster_de_novos_multi(args.multimer,
                                                   args.chains,
                                                   args.proteins,
                                                   proteins_to_chains,
                                                   de_novos,
                                                   three_d_locations,
                                                   genes,
                                                   args.scale,
                                                   args.threshold,
                                                   float(args.degree),
                                                   args.dist_file_output,
                                                   args.runs,#iterations,
                                                   mut_dict,
                                                   scores,
                                                   args.score_col,
                                                   pvalues_in)
    #os.rename('current_distribution.txt',  args.multimer + ".missense.multi.dist.txt")
    if probs is None:
        return

        symbol = ""
        for p in args.proteins:
            symbol += p+"_"
        output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                   "missense",
                                                   len(de_novos[symbol]["missense"]),
                                                   probs["miss_dist"],
                                                   probs["miss_prob"]))
        if("miss_prob_comb" in probs):
                    output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                               "missense_comb",
                                                               len(de_novos[symbol]["missense"]),
                                                               probs["miss_dist_comb"],
                                                               probs["miss_prob_comb"]))
                    output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                               "missense_pois",
                                                               len(de_novos[symbol]["missense"]),
                                                               probs["miss_dist_pois"],
                                                               probs["miss_prob_pois"]))
                
def clustering_coev(output, args):
    
    de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)
    
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
        
        coev_relations =  load_basic_coevol_strength(args.coev_dir, symbol)
        genes = load_gencode([symbol], args.gencode, args.fasta)
        probs =  cluster_de_novos_coevol(symbol,
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

            
def clustering_1d(output, args):
    
    de_novos = load_de_novos(args.input)
    mut_dict = load_mutation_rates(args.rates)
    iterations = 1000000

    if "protein" in args and args.protein in de_novos:
            de_novos = {args.protein : de_novos[args.protein]}
    else:
        return

    for symbol in sorted(de_novos):
        if len(de_novos[symbol]["missense"]) < 2:
            output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol, "missense", len(de_novos[symbol]["missense"]), "N/A", "N/A"))
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
            
        genes = load_gencode([symbol], args.gencode, args.fasta)  
        probs =  cluster_de_novos_1d(symbol,
                                     de_novos[symbol],
                                     genes[symbol],#ensembl,
                                     args.scale,
                                     args.threshold,
                                     args.degree,
                                     args.dist_file_output,
                                     args.runs,#iterations,
                                     mut_dict,
                                     scores,
                                     args.score_col,
                                     pvalues_in)
        
        if probs is None:
            continue

        output.write("gene_id\tmutation_category\tevents_n\tdist\tprobability\n")    
        output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                        "missense",
                                                        len(de_novos[symbol]["missense"]),
                                                        probs["miss_dist"],
                                                        probs["miss_prob"]))
        if("miss_dist_comb" in probs):
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                            "missense_comb",
                                                            len(de_novos[symbol]["missense"]),
                                                            probs["miss_dist_comb"],
                                                            probs["miss_prob_comb"]))
            output.write("{}\t{}\t{}\t{}\t{}\n".format(symbol,
                                                            "missense_pois",
                                                            len(de_novos[symbol]["missense"]),
                                                            probs["miss_dist_pois"],
                                                            probs["miss_prob_pois"]))
            
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

def get_mutation_rates(transcripts, mut_dict):
    """ determines mutation rates per functional category for transcripts
    
    Args:
        transcripts: list of transcript IDs for a gene
        mut_dict: dictionary of local sequence context mutation rates
    
    Returns:
        tuple of (rates, merged transcript, and transcript CDS length)
    """
    
    rates = {'missense': 0, 'nonsense': 0, 'splice_lof': 0,
        'splice_region': 0, 'synonymous': 0}
    combined = None
    
    for tx in transcripts:
        
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

def gene_rates(args, output):
    ''' calculate per-consequence mutation rates per gene
    '''
    mut_dict = load_mutation_rates(args.rates)
    transcripts = load_genes(args.genes)
    gencode = load_gencode(transcripts, args.gencode, args.fasta)
    
    header = ['transcript_id', 'chrom', 'length', 'missense_rate', 'nonsense_rate',
        'splice_lof_rate', 'splice_region_rate', 'synonymous_rate']
    output.write('\t'.join(header) + "\n")
    
    for symbol in sorted(transcripts):
        tx_ids = set(transcripts[symbol])
        gene = gencode[symbol]
        txs = [x for x in gene.transcripts if x.get_name().split('.')[0] in tx_ids]
        try:
            rates, tx, length = get_mutation_rates(txs, mut_dict)
            # log transform rates, for consistency with Samocha et al.
            line = "{}\t{}\t{}\t{}".format(symbol, tx.get_chrom(), length, log_transform(rates))
        except (ValueError, KeyError) as error:
            logging.error(f"{symbol}\t{error}\n")
            line = "{}\t{}\tNA\tNA\tNA\tNA\tNA\tNA".format(symbol, tx.get_chrom())
        
        output.write(line + '\n')
    
    output.close()
    include_frameshift_rates(args.out)

def get_options():
    """ get the command line switches
    """
    
    parser = argparse.ArgumentParser(description='alphacluster cli interface')
    
    ############################################################################
    # CLI options in common
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument("--out", default=sys.stdout, help="output filename")
    parent.add_argument("--rates",
        help="optional path to file containing sequence context-based mutation rates.")
    parent.add_argument("--gencode",
        help="optional path to gencode annotations file. If not provided, gene " \
            "coordinates will be obtained via the Ensembl REST API.")
    parent.add_argument("--fasta",
        help="optional path to genome fasta file. If not provided, gene " \
            "coordinates will be obtained via the Ensembl REST API.")
    parent.add_argument("--genome-build", choices=["grch37",
        "GRCh37", "grch38", "GRCh38"], default="grch38", help="Genome build "
        "that the de novo coordinates are based on (GRCh37 or GRCh38")
    parent.add_argument("--log", default=sys.stdout, help="where to write log files")
    
    subparsers = parser.add_subparsers()

    
    
    ############################################################################
    # CLI options for clustering
    cluster = subparsers.add_parser('cluster', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    cluster.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    cluster.add_argument("--protein_dir", dest="protein_dir", required=True)
    cluster.add_argument("--scores", dest="scores")
    cluster.add_argument("--scale", dest="scale", action="store_true")
    cluster.add_argument("--dbNSFP", dest="dbNSFP")
    cluster.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")
    cluster.add_argument("--threshold", dest="threshold", required=False)
    cluster.add_argument("--runs", dest="runs", type=float, default = 1E6)    
    cluster.add_argument("--protein", dest="protein", required=False)
    cluster.add_argument("--inherited", dest="inherited", required=False)
    cluster.add_argument("--pvalues_in", dest="pvalues_in", required=False)
    cluster.add_argument("--degree", dest="degree", default = 0, type = float, required=False)
    cluster.add_argument("--histogram_data", dest="dist_file_output", required=False)
    cluster.add_argument("--score_col", dest="score_col", required=False,default = 5)    
    cluster.set_defaults(func=clustering)
    ############################################################################
    # CLI options for clustering
    gen_to_aa = subparsers.add_parser('gen_to_aa', parents=[parent],
        description="Get animo acids impacted by variants.")
    gen_to_aa.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    gen_to_aa.add_argument("--protein_dir", dest="protein_dir", required=True)
    gen_to_aa.add_argument("--protein", dest="protein", required=False)

    gen_to_aa.set_defaults(func=genomic_to_residue_position)
    
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
    poisson.add_argument("--score_col", dest="score_col", required=False,default = 5)
    poisson.set_defaults(func=poisson_test)    
    ############################################################################
    # CLI options for poisson_multi_test
    poisson_multi = subparsers.add_parser('poisson_multi', parents=[parent],
        description="Tests the proximity of de novo mutations in genes.")
    poisson_multi.add_argument("--in", dest="input", required=True, help="Path to "
        "file listing known mutations in genes. See example file in data folder "
        "for format.")
    poisson_multi.add_argument("--N", dest="N", type=float)
    poisson_multi.add_argument("--proteins", nargs='+', dest="proteins", required=False)
    poisson_multi.add_argument("--scores", dest="scores")    
    poisson_multi.add_argument("--dbNSFP", dest="dbNSFP")
    poisson_multi.add_argument("--dbNSFP_annotator", dest="dbNSFP_annotator")
    poisson_multi.add_argument("--threshold", dest="threshold", required=False)
    # poisson_multi.add_argument("--females", dest="N_females", type=float, required=False)
    # poisson_multi.add_argument("--males", dest="N_males", type=float, required=False)
    poisson_multi.add_argument("--score_col", dest="score_col", required=False,default = 5)
    poisson_multi.set_defaults(func=poisson_multi_test)    
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
    cluster_multi.add_argument("--score_col", dest="score_col", required=False,default = 5)    
    cluster_multi.add_argument("--scale", dest="scale", action="store_true")
    cluster_multi.add_argument("--threshold", dest="threshold", required=False)
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
    cluster_1d.add_argument("--scale", dest="scale", action="store_true")
    cluster_1d.add_argument("--histogram_data", dest="dist_file_output", required=False)
    cluster_1d.add_argument("--score_col", dest="score_col", required=False,default = 5)
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
    # CLI options for getting mutation rates per gene
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

def open_output(path):
    ''' open output, which could be standard out
    '''
    try:
        output = open(path, 'wt')
    except TypeError:
        output = path
    return output

def main():
    args = get_options()
    FORMAT = '%(asctime)-15s %(message)s'
    log = open(args.log, 'at') if isinstance(args.log, str) else args.log
    logging.basicConfig(stream=log, format=FORMAT, level=logging.INFO)
    
    output = open_output(args.out)
    args.func(output, args)


if __name__ == '__main__':
    main()
