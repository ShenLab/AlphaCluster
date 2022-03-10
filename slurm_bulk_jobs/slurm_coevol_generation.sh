#!/bin/bash

#SBATCH -J COEV_GEN
#SBATCH -o log/COEV_GEN."%j".out
#SBATCH -e log/COEV_GEN."%j".err
#SBATCH --mem=16G

transcript_file=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
transcript_id="${transcript_file%.*}"
echo $transcript_id
#echo "${SLURM_ARRAY_TASK_ID}"
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

./scripts/coevol_calc.py -t $transcript_id -o coevol/$transcript_id.coevol.txt

