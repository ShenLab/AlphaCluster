#!/bin/bash

GENE_FILE=$1
 LEN=$(wc -l $GENE_FILE | awk '{print $1}')
# for RUN in {0..20}
# do
# for COUNT in 2 3 4 5 6 7 8 9 10
# do
#     sbatch --array=1-$LEN%40 slurm_power_analysis_gene.sh $GENE_FILE $COUNT
# done
# done

sbatch --array=1-$LEN%40 slurm_run_gene_asd_coev.sh $GENE_FILE
sbatch --array=1-$LEN%15 slurm_run_gene_ndd_coev.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_asd_ndd.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_chd.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_cdh.sh $GENE_FILE
#sbatch --array=1-$LEN%15 slurm_run_gene_epi.sh $GENE_FILE
