#!/bin/bash

#SBATCH -J CLUSTERING
#SBATCH -o log/CLUSTERING."%j".out
#SBATCH -e log/CLUSTERING."%j".err
#SBATCH --mem=16G
gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
size=$2
cohort_file=$3
dnv_file=$4
enrich_pvalues_file=$5
disease=$6
tag=$7
echo $gene
echo $disease
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

#select subset of samples and variants
mkdir data/$disease.results_$tag
mkdir $disease.results_$tag
mkdir $disease.results_$tag/$gene

grep -e "missense" $dnv_file | grep -e "$gene" >  $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt

#run poisson test
python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.CADD15.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --threshold 15

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.CADD20.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --threshold 20

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --threshold 25

# run tests on these variants
#1D
python denovonear/__main__.py cluster_1d --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.1d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa


#1D with CADD
python denovonear/__main__.py cluster_1d --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.1d.CADDscaled.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale

#3D
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa


#3D with CADD score scaling
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.CADDscaled.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale

# Thresholding runs

#3D with CADD threshold 15
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.CADD15.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.CADDthreshold15.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 15

#3D with CADD threshold 20
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.CADD20.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.CADDthreshold20.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 20

#3D with CADD threshold 25
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.CADDthreshold25.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25


#3D with CADD threshold 25 and scale scoring
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.CADDscaled.CADDthreshold25.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25 --scale


