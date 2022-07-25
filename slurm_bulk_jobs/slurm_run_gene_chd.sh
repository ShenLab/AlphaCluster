#!/bin/bash

#SBATCH -J CHD_3D
#SBATCH -o log/CHD_3D."%j".out
#SBATCH -e log/CHD_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

#alphacluster cluster --in data/all_aff_updated.txt --protein_dir proteins --protein $gene --out out/$gene.out

python alphacluster/__main__.py cluster_1d --in data/chd.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/chd_pvalues.txt --out CHD/$gene.1d.out

python alphacluster/__main__.py cluster_1d --in data/chd.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/chd_pvalues.txt --out CHD/$gene.1d.scores.out

python alphacluster/__main__.py cluster --in data/chd.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/chd_pvalues.txt --out CHD/$gene.3d.out

python alphacluster/__main__.py cluster --in data/chd.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/chd_pvalues.txt --out CHD/$gene.3d.scores.out
