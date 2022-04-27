#!/bin/bash

ASD_GENE_FILE=$1
NDD_GENE_FILE=$2
ASD_LEN=$(wc -l $ASD_GENE_FILE | awk '{print $1}')
NDD_LEN=$(wc -l $NDD_GENE_FILE | awk '{print $1}')
for RUN in {0..100}
do
    for SIZE in 100 500 1000 5000 10000 15000 20000 23425
    do
    	sbatch --array=1-$ASD_LEN%30 slurm_power_analysis_cohort_gene.sh $ASD_GENE_FILE $SIZE data/asd.released.cohort.2.txt data/asd.released.iid.hg38.txt data/asd_enrich_pvalues.txt ASD apr25
    done

    for SIZE in 100 500 1000 5000 10000 15000 20000 25000 30000 31783
    do
	sbatch --array=1-$NDD_LEN%30 slurm_power_analysis_cohort_gene.sh $NDD_GENE_FILE $SIZE data/ndd.released.cohort.txt data/ndd.released.iid.hg38.txt data/ndd_enrich_pvalues.txt NDD apr25
    done
done
