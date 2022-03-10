#!/bin/bash

#SBATCH -J COEVOLNEAR
#SBATCH -o log/COEVOLNEAR."%j".out
#SBATCH -e log/COEVOLNEAR."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

denovonear cluster_coev --in data/all_aff_updated.txt --coev_dir coevol --protein $gene --out out_coevol/$gene.coevol.out
