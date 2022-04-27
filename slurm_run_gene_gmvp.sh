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

head -n 1 $dnv_file > $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt
grep -e "missense" $dnv_file | grep -e "$gene" >>  $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt

#run poisson test
python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.5.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --threshold .5 --score_col 14

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.8.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --threshold .8  --score_col 14

python denovonear/__main__.py poisson --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.9.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --threshold .9  --score_col 14

# run tests on these variants
#1D
#python denovonear/__main__.py cluster_1d --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.1d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#1D with DenovoWEST enrichment
#python denovonear/__main__.py cluster_1d --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 10000000 --pvalues_in $enrich_pvalues_file --out $disease.results_$tag/$gene/$size.$gene.$disease.1d.enrich.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#1D with gMVP
python denovonear/__main__.py cluster_1d --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.1d.gMVPscaled.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale  --score_col 14

#3D
#python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#3D with enrichent
#python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $enrich_pvalues_file --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.enrich.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa


#3D with gMVP score scaling
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale --score_col 14

#3D with gMVP score scaling and enrichment
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $enrich_pvalues_file --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.enrich.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale  --score_col 14

# Thresholding runs
#3D with gMVP threshold .5
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.5.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPthreshold.5.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold .5 --score_col 14

#3D with gMVP threshold .8
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.8.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPthreshold.8.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold .8 --score_col 14

#3D with gMVP threshold .9
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.9.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPthreshold.9.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold .9 --score_col 14


#3D with gMVP threshold .8 and scale scoring
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.8.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.8.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold .8 --scale  --score_col 14

#3D with gMVP threshold .5 and scale scoring
python denovonear/__main__.py cluster --in $disease.results_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --pvalues_in $disease.results_$tag/$gene/$size.$gene.$disease.gMVP.5.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.results_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.5.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold .5 --scale  --score_col 14


