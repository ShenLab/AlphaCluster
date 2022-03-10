#!/bin/bash

#SBATCH -J OBS_BUCKET
#SBATCH -o log/OBS_BUCKET..out
#SBATCH -e log/OBS_BUCKET.err
#SBATCH --mem=64G
python denovonear/__main__.py west_obs_buckets --in data/asd.hg38.txt --gMVP_file scores/gMVP_hg38.2021-02-28.txt.gz --dbNSFP_file scores/dbNSFP4.2a_grch38.gz --dbNSFP_gene_file scores/dbNSFP4.2_gene.gz --MCR_file scores/missenseConstrained.hg38.sorted.bed.gz --bucket_ann gMVP_rankscore CADD_raw_rankscore MCR_region gnomAD_pLI --bucket_cutoffs .25,.5,.75,.9,1.0 .25,.5,.75,.9,1.0 1 .1,.9,1.0 --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --out_dir buckets/
