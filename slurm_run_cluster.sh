#!/bin/bash

#ASD_GENE_FILE=$1
NDD_GENE_FILE=$1
#ASD_LEN=$(wc -l $ASD_GENE_FILE | awk '{print $1}')
NDD_LEN=$(wc -l $NDD_GENE_FILE | awk '{print $1}')

#sbatch --array=1-$ASD_LEN%30 slurm_run_gene.sh $ASD_GENE_FILE $SIZE data/asd.released.cohort.2.txt data/asd.released.iid.hg38.txt data/asd_enrich_pvalues.txt ASD mar18

sbatch --array=1-$NDD_LEN%30 slurm_run_gene.sh $NDD_GENE_FILE 31566 data/ndd.released.cohort.txt data/ndd.released.iid.hg38.txt data/ndd_enrich_pvalues.txt NDD mar21
