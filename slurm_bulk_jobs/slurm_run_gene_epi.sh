#!/bin/bash

#SBATCH -J EPI_3D
#SBATCH -o log/EPI_3D."%j".out
#SBATCH -e log/EPI_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python alphacluster/__main__.py cluster_1d --in data/epi.hg38.txt --protein $gene --pvalues_in data/epi_pvalues.txt --runs 90000000 --out EPI/$gene.1d.out

python alphacluster/__main__.py cluster_1d --in data/epi.hg38.txt --protein $gene --pvalues_in data/epi_pvalues.txt --runs 90000000 --out EPI/$gene.1d.scores.out

python alphacluster/__main__.py cluster --in data/epi.hg38.txt --protein_dir proteins --protein $gene --pvalues_in data/epi_pvalues.txt --runs 90000000 --out EPI/$gene.3d.out

python alphacluster/__main__.py cluster --in data/epi.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --pvalues_in data/epi_pvalues.txt --runs 90000000 --out EPI/$gene.3d.scores.out
