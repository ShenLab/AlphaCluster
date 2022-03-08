#!/bin/bash

#SBATCH -J RENAME_COEV
#SBATCH -o log/RENAME_GEN."%j".out
#SBATCH -e log/RENAME_GEN."%j".err
#SBATCH --mem=16G

file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $transcript_id
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python /scripts/rename_coevol.py 

