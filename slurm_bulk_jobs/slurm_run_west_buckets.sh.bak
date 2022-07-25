#!/bin/bash

GENE_FILE=$1
LEN=$(wc -l $GENE_FILE | awk '{print $1}')
sbatch --array=1-$LEN%200 slurm_run_west_gene_bucket.sh $GENE_FILE 

