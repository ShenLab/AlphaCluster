#!/bin/bash

ASD_GENE_FILE=$1
ASD_LEN=$(wc -l $ASD_GENE_FILE | awk '{print $1}')

#sbatch --array=1-$ASD_LEN%30 slurm_run_gene_gmvp.sh $ASD_GENE_FILE 23425 data/asd.released.cohort.txt data/asd.released.iid.hg38.txt data/asd_enrich_pvalues.txt ASD apr8.gMVP

#sbatch --array=1-$ASD_LEN%100 slurm_run_gene.sh $ASD_GENE_FILE 23425 data/asd.released.cohort.txt data/asd.released.iid.hg38.txt data/asd_enrich_pvalues.txt ASD apr25

sbatch --array=1-$ASD_LEN%100 slurm_run_gene.sh $ASD_GENE_FILE 23425 data/asd.released.cohort.txt data/asd.clcn4.iid.hg38.txt data/asd_enrich_pvalues.txt ASD apr27
