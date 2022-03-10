#!/bin/bash

#SBATCH -J CDH_3D
#SBATCH -o log/CDH_3D."%j".out
#SBATCH -e log/CDH_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

#denovonear cluster --in data/all_aff_updated.txt --protein_dir proteins --protein $gene --out out/$gene.out

python denovonear/__main__.py cluster_1d --in data/cdh.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/cdh_pvalues.txt --out CDH/$gene.1d.out

python denovonear/__main__.py cluster_1d --in data/cdh.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/cdh_pvalues.txt --out CDH/$gene.1d.scores.out

python denovonear/__main__.py cluster --in data/cdh.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/cdh_pvalues.txt --out CDH/$gene.3d.out

python denovonear/__main__.py cluster --in data/cdh.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/cdh_pvalues.txt --out CDH/$gene.3d.scores.out
