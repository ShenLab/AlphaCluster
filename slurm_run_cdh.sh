#!/bin/bash

CDH_GENE_FILE=$1
CDH_LEN=$(wc -l $CDH_GENE_FILE | awk '{print $1}')

sbatch --array=1-$CDH_LEN%30 slurm_run_gene_gmvp.sh $CDH_GENE_FILE 827 data/empty.tx data/cdh.hg38.txt data/empty.txt CDH apr4.gMVP
