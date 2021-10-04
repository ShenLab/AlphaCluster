#!/bin/bash

#SBATCH -J THREEDNEAR
#SBATCH -o log/THREEDNEAR."%j".out
#SBATCH -e log/THREEDNEAR."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

denovonear cluster --in data/all_aff_updated.txt --protein_dir proteins --protein $gene --out out/$gene.out
