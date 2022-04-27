#!/bin/bash

ASD_GENE_FILE=$1
NDD_GENE_FILE=$2
ASD_LEN=$(wc -l $ASD_GENE_FILE | awk '{print $1}')
NDD_LEN=$(wc -l $NDD_GENE_FILE | awk '{print $1}')
for RUN in {0..100}
do
    for SIZE in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
    do
    	sbatch --array=1-$ASD_LEN%30 slurm_power_analysis_dnvs_gene.sh $ASD_GENE_FILE $SIZE data/asd.released.iid.hg38.txt data/asd_enrich_pvalues.txt ASD mar18
    done

    for SIZE in  2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
    do
	sbatch --array=1-$NDD_LEN%30 slurm_power_analysis_dnvs_gene.sh $NDD_GENE_FILE $SIZE data/ndd.released.iid.hg38.txt data/ndd_enrich_pvalues.txt NDD mar18
    done
done
