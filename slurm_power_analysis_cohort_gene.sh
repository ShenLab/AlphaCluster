#!/bin/bash

#SBATCH -J ASD_POWER_C
#SBATCH -o log/POWER_C."%j".out
#SBATCH -e log/POWER_C."%j".err
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
mkdir data/$disease.power_cohort_$tag
mkdir $disease.power_cohort_$tag
mkdir $disease.power_cohort_$tag/$gene

head -n 1 $dnv_file > $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt

grep -e "Affected" $cohort_file | shuf -n $size | awk -F"\t" '{print $1}' | grep -f - $dnv_file | grep -e "missense" | grep -e "$gene" >>  $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt

#run poisson test
python denovonear/__main__.py poisson --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py poisson --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD15.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --threshold 15

python denovonear/__main__.py poisson --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --threshold 25 --score_col 5

python denovonear/__main__.py poisson --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.7.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz  --threshold .7 --score_col 14


python denovonear/__main__.py poisson --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --N $size --protein $gene --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.85.poisson_pvalue.$SLURM_JOB_ID.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz  --threshold .85 --score_col 14

# run tests on these variants
#1D
#python denovonear/__main__.py cluster_1d --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.1d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#1D with DenovoWEST enrichment
if [[ $size == 31783 || $size == 23425 ]];
then
#       python denovonear/__main__.py cluster_1d --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 5000000 --pvalues_in $enrich_pvalues_file --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.1d.enrich.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa
fi

#1D with CADD
#python denovonear/__main__.py cluster_1d --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.1d.CADD.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale

#1D with gMVP
#python denovonear/__main__.py cluster_1d --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.1d.CADD.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz

#3D
#python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

if [[ $size == 31783 || $size == 23425 ]];
then
#       python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $enrich_pvalues_file --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.enrich.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa
fi

#3D with CADD score scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADDscaled.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale --score_col 5

#3D with gMVP score scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.7.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.$SLURM_JOB_ID.out  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --scale --score_col 14 --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz

# Thresholding runs

#3D with CADD threshold 15
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD15.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold15.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 15

#3D with CADD threshold 20
#python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD20.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold20.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 20

#3D with CADD threshold 25
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold25.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25 --score_col 5

#3D with gMVP threshold 0.7 
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.7.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.7.$SLURM_JOB_ID.out  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.7 --score_col 14 --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz

#3D with gMVP threshold 0.85 
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.85.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.85.$SLURM_JOB_ID.out  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.85 --score_col 14 --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz

#3D with CADD threshold 15 with CADD scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD15.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADDscaled.CADDthreshold15.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 15 --scale --score_col 5

#3D with CADD threshold 15 with CADD scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.CADD25.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.CADDscaled.CADDthreshold25.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25 --scale --score_col 5

#3D with gMVP threshold 0.7 with gMVP scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.7.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.7.$SLURM_JOB_ID.out  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.7 --scale --score_col 14 --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz

#3D with gMVP threshold 0.7 with gMVP scaling
python denovonear/__main__.py cluster --in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 5000000 --pvalues_in $disease.power_cohort_$tag/$gene/$size.$gene.$disease.gMVP.85.poisson_pvalue.$SLURM_JOB_ID.txt --out $disease.power_cohort_$tag/$gene/$size.$gene.$disease.3d.gMVPscaled.gMVPthreshold.85.$SLURM_JOB_ID.out  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.85 --scale --score_col 14 --score /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz 
