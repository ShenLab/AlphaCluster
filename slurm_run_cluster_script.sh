#!/bin/bash

GENE_FILE=$1
LEN=$(wc -l $GENE_FILE | awk '{print $1}')
sbatch --array=1-$LEN%15 slurm_run_gene.sh $GENE_FILE
