#!/bin/bash

FILE=$1
LEN=$(wc -l $FILE | awk '{print $1}')
sbatch --array=1-$LEN%20 slurm_rename_coevol.sh $FILE
