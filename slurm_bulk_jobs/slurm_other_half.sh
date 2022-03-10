#!/bin/bash

GENE_FILE=$1
LEN=$(wc -l $GENE_FILE | awk '{print $1}')
#sbatch --array=1-$LEN%15 slurm_run_gene_asd.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_ndd.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_chd.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_cdh.sh $GENE_FILE
sbatch --array=1-$LEN%15 slurm_run_gene_asd_ndd.sh $GENE_FILE
