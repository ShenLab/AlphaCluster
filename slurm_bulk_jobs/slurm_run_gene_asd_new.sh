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

python denovonear/__main__.py cluster_1d --in data/asd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_new/$gene.1d.out

#python denovonear/__main__.py cluster_1d --in data/asd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_new/$gene.1d.gMVP.out

python denovonear/__main__.py cluster --in data/asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_new/$gene.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/asd.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_new/$gene.3d.gMVP.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_new/$gene.3d.REVEL.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#python denovonear/__main__.py cluster --in data/asd.hg38.txt --protein_dir proteins --protein $gene --runs 90000000  --out ASD/$gene.3d.rec.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --males 17750 --females 4262 --pvalues data/asd_enrich_pvalues.txt

mkdir data/$SLURM_JOB_ID
mkdir results
mkdir results/$SLURM_JOB_ID

python denovonear/__main__.py poisson --in data/asd.released.iid.hg3d.txt --N 23425 --protein $gene --out data/$SLURM_JOB_ID/REVEL0.5.poisson_pvalue.txt  --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_score --threshold 0.5

python denovonear/__main__.py cluster --in data/asd.released.iid.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/$SLURM_JOB_ID/REVEL0.5.poisson_pvalue.txt --out results/$SLURM_JOB_ID/$gene.3d.REVEl.threshold0.5.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_score --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 0.5

python denovonear/__main__.py cluster --in data/asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out results/$SLURM_JOB_ID/$gene.3d.REVEL.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/asd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out results/$SLURM_JOB_ID/$gene.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa
