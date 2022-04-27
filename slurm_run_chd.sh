#!/bin/bash

CHD_GENE_FILE=$1
CHD_LEN=$(wc -l $CHD_GENE_FILE | awk '{print $1}')

sbatch --array=1-$CHD_LEN%30 slurm_run_gene_gmvp.sh $CHD_GENE_FILE 3678 data/empty.tx data/chd.hg38.txt data/empty.txt CHD apr4.gMVP
