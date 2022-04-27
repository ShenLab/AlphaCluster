#!/bin/bash

EPI_GENE_FILE=$1
EPI_LEN=$(wc -l $EPI_GENE_FILE | awk '{print $1}')

sbatch --array=1-$EPI_LEN%30 slurm_run_gene_gmvp.sh $EPI_GENE_FILE 356 data/empty.tx data/epi.hg38.txt data/empty.txt EPI apr4.gMVP
