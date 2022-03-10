#!/bin/bash

FILE=$1
LEN=$(wc -l $TRANSCRIPT_FILE | awk '{print $1}')
sbatch --array=1-$LEN slurm_rename_coevol.sh $FILE
