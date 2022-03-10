#!/bin/bash

#SBATCH -J ASD_NDD_3D
#SBATCH -o log/ASD_NDD_3D."%j".out
#SBATCH -e log/ASD_NDD_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python denovonear/__main__.py cluster_1d --in data/asd_ndd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/asd_ndd_released_pvalues.txt --out ASD_NDD_released/$gene.1d.out

python denovonear/__main__.py cluster_1d --in data/asd_ndd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/asd_ndd_released_pvalues.txt --out ASD_NDD_released/$gene.1d.gMVP.out

python denovonear/__main__.py cluster --in data/asd_ndd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_ndd_released_pvalues.txt --out ASD_NDD_released/$gene.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/asd_ndd.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_ndd_released_pvalues.txt --out ASD_NDD_released/$gene.3d.gMVP.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python denovonear/__main__.py cluster --in data/asd_ndd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/asd_released_pvalues.txt --out ASD_NDD_released/$gene.3d.DANN.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator DANN_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#python denovonear/__main__.py cluster --in data/asd_ndd.hg38.txt --protein_dir proteins --protein $gene --runs 90000000  --out ASD_NDD/$gene.3d.rec.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --males 17750 --females 4262 --pvalues data/asd_ndd_enrich_pvalues.txt
