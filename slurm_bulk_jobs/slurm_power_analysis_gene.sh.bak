#!/bin/bash

#SBATCH -J ASD_3D
#SBATCH -o log/ASD_3D."%j".out
#SBATCH -e log/ASD_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

count=$2

# select subset of variants

grep -e "$gene" data/asd.released.hg38.txt |grep -e "missense" | shuf -n $count > data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt

mkdir ASD_power
mkdir ASD_power/$gene

# run tests on these variants

python denovonear/__main__.py cluster_1d --in data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_power/$gene/$gene.$SLURM_JOB_ID.1d.out

python denovonear/__main__.py cluster --in data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_power/$gene/$gene.$SLURM_JOB_ID.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_power/$gene/$gene.$SLURM_JOB_ID.3d.gMVP.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_power/$gene/$gene.$SLURM_JOB_ID.3d.REVEL.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/power/$SLURM_JOB_ID.$count.asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_power/$gene/$gene.$SLURM_JOB_ID.3d.CADD.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator CADD_raw_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa


