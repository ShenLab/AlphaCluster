#!/bin/bash

GENE_FILE=$1
RUNS=$2
LEN=$(wc -l $GENE_FILE | awk '{print $1}')
for RUN in {0..100}
do
    for SIZE in 100 500 1000 5000 10000 15000 20000 25000 30000
    do
    	sbatch --array=1-$LEN%30 slurm_power_analysis_cohort_gene.sh $GENE_FILE $SIZE data/asd.released.cohort.txt data/asd.released.iid.hg38.txt ASD feb23
	#sbatch --array=1-$LEN%30 slurm_power_analysis_cohort_gene.sh $GENE_FILE $SIZE data/ndd.released.cohort.txt data/ndd.released.iid.hg38.txt NDD feb23
    done
done
