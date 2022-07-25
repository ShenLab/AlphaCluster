#!/bin/bash

#SBATCH -J BUCKET.$gene
#SBATCH -o log/BUCKET.$gene."%j".out
#SBATCH -e log/BUCKET.$gene."%j".err
#SBATCH --mem=64G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python alphacluster/__main__.py west_exp_buckets --males 17750 --females 4262 --gMVP_file scores/gMVP_hg38.2021-02-28.txt.gz --dbNSFP_file scores/dbNSFP4.2a_grch38.gz --dbNSFP_gene_file scores/dbNSFP4.2_gene.gz --MCR_file scores/missenseConstrained.hg38.sorted.bed.gz --bucket_ann gMVP_rankscore CADD_raw_rankscore MCR_region gnomAD_pLI --bucket_cutoffs 0,.25,.5,.75,.9,1.0 0,.25,.5,.75,.9,1.0 0,1 0,.1,.9,1.0 --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --out buckets/ --gene $gene
