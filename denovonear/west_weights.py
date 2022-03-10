import bisect
import gzip
import pysam
import os
from itertools import product, chain

from denovonear.site_specific_rates import SiteRates
from load_gene import best_transcript
from denovonear.gencode import Gencode

def get_annotation_score(gene_scores, annotations_files, annotator, gene, chrom, pos, alt):
    pos = int(pos)
    score = -1

    if gene_scores != {} and gene in gene_scores and annotator in gene_scores[gene]:
        score = float(gene_scores[gene][annotator])        
    elif annotator == "gMVP_rankscore":
        score_file = annotations_files["gMVP_rankscore"]
        score_obj = pysam.TabixFile(score_file)
        try:
            for line in score_obj.fetch(chrom, pos-1, pos):
                values = line.split('\t')
                if values[3] != alt:
                    continue
                score = float(values[14])
        except:
            score = -1
    elif annotator == "MCR_region":
        score_file = annotations_files["MCR_region"]
        score_obj = pysam.TabixFile(score_file)
        try:
            score_obj.fetch("chr"+chrom,pos-1,pos)
            score = 1
        except:
            score = 0
    else:
        dbNSFP_score_file = annotations_files["dbNSFP"]
        dbNSFP_score_obj = pysam.TabixFile(dbNSFP_score_file)
        dbNSFP_header = dbNSFP_score_obj.header[0].split('\t')
            
        if annotator in dbNSFP_header:
            anno_idx = dbNSFP_header.index(annotator)
            alt_idx = dbNSFP_header.index('alt')
            try:
                for line in dbNSFP_score_obj.fetch(chrom, pos-1,pos):
                    values = line.split('\t')
                    #print(values[alt_idx])
                    #print(values[anno_idx])
                    alt_ = values[alt_idx]
                    if alt_ != alt:
                        continue
                    score = float(values[anno_idx])
            except:
                score = -1
    #print(score)
    return score

def get_bucket(gene, chrom, position, alt, gene_scores, bucket_def, annotation_files):
    '''
    bucket_def is of form {"CADD_rankscore":[0,.25,.75,1.00], "gMVP_rankscore":[0,.25,.75,1.00], "sHet":[.1,.9], "MCR": [T,F]}, etc...}
    '''
    bucket = dict.fromkeys(bucket_def.keys())
    scores = {}
    for annotator, ranges in bucket_def.items():
        score = get_annotation_score(gene_scores, annotation_files, annotator, gene, chrom, position, alt)
        #print(gene, chrom, position, alt, annotator, score)
        bucket[annotator] = ranges[bisect.bisect_left(ranges,score)]
        scores[annotator] = score
        #print(ranges[bisect.bisect_left(ranges,score)])
    return tuple(bucket.values()), scores

def get_all_buckets_from_bucket_def(bucket_def):
    all_tuples = ()
    #print(bucket_def)
    prod = product(*list(bucket_def.values()))
    all_buckets = dict.fromkeys([tuple(p) for p in prod], 0)
    #print(all_buckets)
    return all_buckets

def get_exp(prob,
            autosomal,
            x_factor,
            y_factor,
            chrom,
            position,
            alt):
    
    if chrom == "X" or chrom == "chrX":
        exp = prob*autosomal*x_factor
    elif chrom == "Y" or chrom == "chrY":
        exp = prob*autosomal*y_factor
    else:
        exp = prob*autosomal

    return exp


def get_pred_count(prob,
                   N_males,
                   N_females,
                   chrom):

    # Get important factors for scaling probabilities
    autosomal = 2 * (N_males + N_females)
    alpha =  3.4
    male_factor = 2 / (1 + (1 / alpha))
    female_factor = 2 / (1 + alpha)
    male_x_transmissions = N_females
    female_x_transmissions = N_males + N_females
    male_y_transmissions = N_males
    female_y_transmissions = 0
    x_factor = ((male_x_transmissions * male_factor) + (female_x_transmissions * female_factor)) / autosomal
    y_factor = ((male_y_transmissions * male_factor)+ (female_y_transmissions * female_factor)) / autosomal

        
    if chrom == "X" or chrom == "chrX":
        exp = prob*autosomal*x_factor
    elif chrom == "Y" or chrom == "chrY":
        exp = prob*autosomal*y_factor
    else:
        exp = prob*autosomal

    return exp

async def expected_per_bucket(gene_name, genes,mut_dict, N_males, N_females, gene_scores, bucket_def, annotation_files):
    all_buckets_exp = get_all_buckets_from_bucket_def(bucket_def)

    # Get important factors for scaling probabilities
    autosomal = 2 * (N_males + N_females)
    alpha =  3.4
    male_factor = 2 / (1 + (1 / alpha))
    female_factor = 2 / (1 + alpha)
    male_x_transmissions = N_females
    female_x_transmissions = N_males + N_females
    male_y_transmissions = N_males
    female_y_transmissions = 0
    x_factor = ((male_x_transmissions * male_factor) + (female_x_transmissions * female_factor)) / autosomal
    y_factor = ((male_y_transmissions * male_factor)+ (female_y_transmissions * female_factor)) / autosomal

    # Get score for all possible snv
    snv_by_bucket = []
    snv_by_scores = []    
    
    if gene_name == None:
        genes = [ gene for gene in genes ]
    else:
        genes = [genes[gene_name]]
    for gene in genes:
        transcript = gene.canonical
        rates = SiteRates(transcript,mut_dict, {})
        weight = rates["missense"]
        print(gene.symbol)
        for choice in weight:
            position = transcript.get_position_on_chrom(choice["pos"],0)
            alt = choice["alt"]
            prob = choice['prob']
            bucket, scores = get_bucket(gene.symbol, gene.chrom[3:], position, alt, gene_scores, bucket_def, annotation_files)
            exp = get_exp(prob, autosomal, x_factor, y_factor, gene.chrom, position, alt)
            all_buckets_exp[bucket] = all_buckets_exp[bucket]+exp
            snv_by_bucket.append([gene.symbol, gene.chrom[3:], position, alt, bucket])
            snv_by_scores.append([gene.symbol, gene.chrom[3:], position, alt, scores])
    #for gene, chrom, position, alt in all_possible_snv:
    #    bucket = get_bucket(gene, chrom, position, alt, gene_scores, bucket_def, annotation_files)
    #    exp = get_exp(autosomal, x_factor, y_factor, chrom, position, alt)
        #print(exp)
        #print(bucket)
    #    all_buckets_exp[bucket] = all_buckets_exp[bucket] + exp

    return (all_buckets_exp, snv_by_bucket, snv_by_scores)

async def observed_per_bucket(de_novos, gene_scores, bucket_def, annotation_files):
    all_buckets_obs = get_all_buckets_from_bucket_def(bucket_def)
    for gene, dnvs in de_novos.items():
        print(gene)
        for chrom, position, alt in dnvs["missense"]:
            bucket, _ = get_bucket(gene, chrom, position, alt, gene_scores, bucket_def, annotation_files)
            all_buckets_obs[bucket] = all_buckets_obs[bucket] + 1

    return all_buckets_obs

async def get_scores_per_bucket(de_novos, genes, mut_dict, N_males, N_females, gene_scores, bucket_def, annotation_files):
    all_buckets_exp, snv_by_bucket, _ = await expected_per_bucket(None, genes, mut_dict, N_males, N_females, gene_scores, bucket_def, annotation_files)
    all_buckets_obs = await observed_per_bucket(de_novos, gene_scores, bucket_def, annotation_files)

    all_buckets_scores = get_all_buckets_from_bucket_def(bucket_def)

    for bucket, score in all_buckets_scores.items():
        all_buckets_scores[bucket] = all_buckets_obs[bucket]/all_buckets_exp[bucket]

    return (all_buckets_scores, snv_by_bucket)

#def get_score(gene, chrom, position, alt, all_buckets_scores, gene_scores, bucket_def, annotation_files):
#    bucket = get_bucket(gene, chrom, position, alt, gene_scores, bucket_def, annotation_files)
#    score = all_buckets_scores[bucket]
#    return score

def write_bucket_counts(output, bucket_counts):
    #print(bucket_counts)
    with open(output, "w") as output_file:    
        for bucket, count in bucket_counts.items():
            line = '\t'.join([str(bucket), str(count)])
            output_file.write(line + "\n")
            
def write_scores_per_snv(output, snv_by_scores):
    with open(output, "w") as output_file:
        # Header
        line = "gene\tchrom\tpos\talt"
        for annotator in snv_by_scores[0][4].keys():
            line = line + "\t" + annotator
        output_file.write(line + "\n")
        # Values
        for gene, chrom, pos, alt, scores in snv_by_scores:
            line = '\t'.join([gene,str(chrom), str(pos), alt])
            for score in scores.values():
                line = line + "\t" + str(score)
            output_file.write(line + "\n")
    

def write_bucket_per_snv(output, snv_by_bucket):
    with open(output, "w") as output_file:
        output_file.write("gene\tchrom\tpos\talt\tbucket\n")
        for gene, chrom, pos, alt, bucket in snv_by_bucket:
            line = '\t'.join([gene,str(chrom), str(pos), alt,  str(bucket)])
            output_file.write(line + "\n")

def write_all_buckets_scores(output, all_bucket_scores, bucket_def):
    with open(output, "x") as output_file:
        header = '\t'.join([str(x) for x in bucket_def.keys()])
        header = header + "\t score" + "\n"
        output_file.write(header)
        for bucket, score in all_buckets_scores.items():
            line = '\t'.join(bucket) + "\t" + score + "\n"
            output_file.write(line)

    with open(output + ".def", "x") as output_file:
        output_file.write(json.dumps(bucket_def))

def get_gene_scores(bucket_def, annotation_files):
    annotators = bucket_def.keys()
    gene_scores = {}
    dbNSFP_gene_score_file = annotation_files["dbNSFP_gene"]
    with gzip.open(dbNSFP_gene_score_file, 'rt') as f:
        idx = 0
        for line in f:
            if idx == 0:
                dbNSFP_gene_header = line.split('\t')
                gene_idx = dbNSFP_gene_header.index('Gene_name')
                anno_idxs = {}
                for annotator in annotators:
                    if annotator in dbNSFP_gene_header:
                        anno_idxs[annotator] = dbNSFP_gene_header.index(annotator)
            elif len(anno_idxs) > 0:
                values = line.split('\t')
                gene_scores[values[gene_idx]] = {}
                for annotator, anno_idx in anno_idxs.items():
                    score = values[anno_idx]
                    if score == ".":
                        score = 0
                    elif score == "NA":
                        score = -1
                    else:
                        score = float(score)
                    gene_scores[values[gene_idx]][annotator] = score
            else:
                break
            idx=idx+1
    return gene_scores

async def generate_observed_buckets(output_dir,
                                    de_novos,
                                    bucket_def,
                                    annotation_files):
    for key, ranges in bucket_def.items():
        bucket_def[key] = [-1] + bucket_def[key]
        
    gene_scores = get_gene_scores(bucket_def, annotation_files)

    all_buckets_obs = await observed_per_bucket(de_novos,
                                                gene_scores,
                                                bucket_def,
                                                annotation_files)
    write_bucket_counts(output_dir+"/"+"observed.buckets_count.txt", all_buckets_obs)
    
async def generate_expected_buckets(output_score_filename,
                           N_males,
                           N_females,
                           gene_name,
                           genes,
                           mut_dict,
                           bucket_def,
                           annotation_files):
    
    for key, ranges in bucket_def.items():
        bucket_def[key] = [-1] + bucket_def[key]
        
    gene_scores = get_gene_scores(bucket_def, annotation_files)
    all_buckets_exp, snv_by_bucket, snv_by_scores = await expected_per_bucket(gene_name,
                                                               genes,
                                                               mut_dict,
                                                               N_males,
                                                               N_females,
                                                               gene_scores,
                                                               bucket_def,
                                                               annotation_files)
    write_bucket_per_snv(output_score_filename+gene_name+".buckets.txt", snv_by_bucket)
    write_scores_per_snv(output_score_filename+gene_name+".scores.txt", snv_by_scores)
    write_bucket_counts(output_score_filename+gene_name+".buckets_count.txt", all_buckets_exp)

def generate_scores_from_buckets(output_dir,
                                 genes,
                                 obs_bucket_file,
                                 bucket_dir,
                                 bucket_def):
    exp_bucket_counts = {}
    for filename in os.listdir(bucket_dir):
        if filename.endswith(".buckets_count.txt") and filename != "observed.buckets_count.txt":
            with open(bucket_dir + "/" + filename) as f: 
                for line in f:
                    values = line.strip('\n').split('\t')
                    bucket = values[0]
                    count = float(values[1])
                    if bucket not in exp_bucket_counts:
                        exp_bucket_counts[bucket] = 0
                    exp_bucket_counts[bucket] += count 

    obs_bucket_counts = {}#get_all_buckets_from_bucket_def(bucket_def)
    with open(obs_bucket_file) as f:
        for line in f:
            values = line.strip('\n').split('\t')
            bucket = values[0]
            count = float(values[1])
            if bucket not in obs_bucket_counts:
                obs_bucket_counts[bucket] = 0
            obs_bucket_counts[bucket] += count 
        
    all_buckets_OR = {}#get_all_buckets_from_bucket_def(bucket_def)
    all_buckets_PPV = {}#get_all_buckets_from_bucket_def(bucket_def)
    for bucket in exp_bucket_counts.keys():
        if bucket in obs_bucket_counts and obs_bucket_counts[bucket] != 0:
            OR = obs_bucket_counts[bucket]/exp_bucket_counts[bucket]
            all_buckets_OR[bucket] = OR
            all_buckets_PPV[bucket] = (OR-1)/OR        

    write_bucket_counts(output_dir+"/OR_per_bucket.txt", all_buckets_OR)
    write_bucket_counts(output_dir+"/PPV_per_bucket.txt", all_buckets_PPV)    

    snv_by_scores = []        
    with open(output_dir + "/PPV_per_snv.txt", "w+") as o:
        for filename in os.listdir(bucket_dir):
            if filename.endswith(".buckets.txt") and filename != "observed.buckets.txt":
                with open(bucket_dir + "/" + filename) as f:
                    for line in f:
                        values = line.strip('\n').split('\t')
                        bucket = values[4]
                        if bucket != "bucket" and bucket in all_buckets_PPV:
                            o.write(line + '\t' +str(all_buckets_PPV[bucket]) + '\n')
        
    
        

# async def generate_scores(output_score_filename,
#                           de_novos,
#                           N_males,
#                           N_females,
#                           genes,
#                           mut_dict,
#                           bucket_def,
#                           annotation_files):

#     for key, ranges in bucket_def.items():
#         bucket_def[key] = [-1] + bucket_def[key]
        
#     #all_possible_snv = []
#     #score_file = annotation_files["gMVP_rankscore"]
#     #score_obj = pysam.TabixFile(score_file)
#     #for line in score_obj.fetch(1, 1, 69093):
#     #    values = line.split('\t')
#     #    all_possible_snv.append([values[4], values[0], values[1], values[3]])
#     #print(str(values[0])+ " "+ str(values[1] + " "+str(values[2]) + " " + str(values[3])))

#     #print(len(all_possible_snv))

#     gene_scores = get_gene_scores(bucket_def, annotation_files)
#     #print(len(gene_scores))
#     all_buckets_scores, score_by_snv = await get_scores_per_bucket(de_novos,
#                                                                    genes,
#                                                                    mut_dict,
#                                                                    N_males,
#                                                                    N_females,
#                                                                    gene_scores,
#                                                                    bucket_def,
#                                                                    annotation_files)
                                                                   
    
#     write_scores_per_snv(output_score_filename,
#                          all_possible_snv,
#                          all_buckets_scores,
#                          bucket_def,
#                          annotation_files)


'''
Sample run

__main__.py west_weight
--in data/asd.hg38.txt        
--N_males 17750
--N_females 44262
--gMVP_file scores/gMVP.gz
--dbNSFP_file scores/dbNSFP4.2a_grch38.gz
--MCR_file scores/MCR.gz
--bucket_ann gMVP_rankscore CADD_rankscore MCR_region gnomad_pLI
--bucket_cutoffs [.25,.5,.75,.9,1.0] [.25,.5,.75,.9,1.0] [True,False] [.1, .9, 1.0]
'''
