#!/bin/bash

#SBATCH -J SCORES
#SBATCH -o log/SCORES."%j".out
#SBATCH -e log/SCORES."%j".err
#SBATCH --mem=128G
python denovonear/__main__.py west_scores --obs_buckets buckets/observed.buckets_count.txt --buckets_dir buckets/ --out_dir buckets/ --bucket_ann gMVP_rankscore CADD_raw_rankscore MCR_region gnomAD_pLI --bucket_cutoffs 0,.25,.5,.75,.9,1.0 0,.25,.5,.75,.9,1.0 0,1 0,.1,.9,1.0 --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa
