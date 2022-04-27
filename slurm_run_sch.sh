#!/bin/bash

SCH_GENE_FILE=$1
SCH_LEN=$(wc -l $SCH_GENE_FILE | awk '{print $1}')

sbatch --array=1-$SCH_LEN%30 slurm_run_gene_gmvp.sh $SCH_GENE_FILE 2773 data/empty.tx data/sch.hg38.txt data/empty.txt SCH apr4.gMVP
