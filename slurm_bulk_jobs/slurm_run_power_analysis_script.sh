#!/bin/bash

GENE_FILE=$1
RUNS=$2
LEN=$(wc -l $GENE_FILE | awk '{print $1}')
for RUN in {0..20}
do
for COUNT in 2 3 4 5 6 7 8 9 10
do
    sbatch --array=1-$LEN%40 slurm_power_analysis_gene.sh $GENE_FILE $COUNT
done
done
