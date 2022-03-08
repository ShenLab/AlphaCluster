#!/bin/bash

#SBATCH -J ASD_POWER_C
#SBATCH -o log/POWER_C."%j".out
#SBATCH -e log/POWER_C."%j".err
#SBATCH --mem=16G
gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
size=$2
cohort_file=$3
dnv_file=$4
disease=$5
tag=$6
echo $gene
echo $disease
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

#select subset of samples and variants
mkdir data/$disease.power_cohort_$tag

grep -e "Affected" $cohort_file | shuf -n $size | awk -F"\t" '{print $4}' | grep -f - $dnv_file | grep -e "missense" >  data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt

#run poisson test
python alphacluster/__main__.py poisson --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --N $size --protein $gene --out data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py poisson --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --N $size --protein $gene --out data/$disease.power_cohort_$tag/$SLURM_JOB_ID.REVEL0.5.poisson_pvalue.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_score --threshold 0.5

python alphacluster/__main__.py poisson --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --N $size --protein $gene --out data/$disease.power_cohort_$tag/$SLURM_JOB_ID.CADD25.poisson_pvalue.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator CADD_phred --threshold 25

mkdir $disease.power_cohort_$tag
mkdir $disease.power_cohort_$tag/$gene

# run tests on these variants

python alphacluster/__main__.py cluster_1d --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.1d.out

python alphacluster/__main__.py cluster_1d --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.1d.gMVP.out

python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.gMVP.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.REVEL.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.CADD.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator CADD_raw_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

# Thresholding runs
python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.REVEL0.5.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.REVEl.threshold0.5.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_score --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.5

python alphacluster/__main__.py cluster --in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.cohort_size_$size.$disease.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$disease.power_cohort_$tag/$SLURM_JOB_ID.CADD25.poisson_pvalue.txt --out $disease.power_cohort_$tag/$gene/$gene.cohort_size_$size.$SLURM_JOB_ID.3d.CADD.threshold25.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator CADD_phred --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25


