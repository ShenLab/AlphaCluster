#!/bin/bash

TRANSCRIPT_FILE=$1
LEN=$(wc -l $TRANSCRIPT_FILE | awk '{print $1}')
sbatch --array=1-$LEN slurm_coevol_generation.sh $TRANSCRIPT_FILE
