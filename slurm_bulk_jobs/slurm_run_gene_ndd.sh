#!/bin/bash

#SBATCH -J NDD_3D
#SBATCH -o log/NDD_3D."%j".out
#SBATCH -e log/NDD_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python alphacluster/__main__.py cluster_1d --in data/ndd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/ndd_released_pvalues.txt --out NDD_known/$gene.1d.out

#python alphacluster/__main__.py cluster_1d --in data/ndd.released.hg38.txt --protein $gene --runs 90000000 --pvalues_in data/ndd_released_pvalues.txt --out NDD_released/$gene.1d.gMVP.out

python alphacluster/__main__.py cluster --in data/ndd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/ndd_released_pvalues.txt --out NDD_known/$gene.3d.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py cluster --in data/ndd.released.hg38.txt --scores /share/terra/rsrc/hg38/gMVP/gMVP_hg38.2021-02-28.txt.gz --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/ndd_released_pvalues.txt --out NDD_known/$gene.3d.gMVP.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

python alphacluster/__main__.py cluster --in data/ndd.released.hg38.txt --protein_dir proteins --protein $gene --runs 90000000 --pvalues_in data/ndd_released_pvalues.txt --out NDD_known/$gene.3d.REVEL.out --dbNSFP scores/dbNSFP4.2a_grch38.gz --dbNSFP_annotator REVEL_rankscore --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#python alphacluster/__main__.py cluster --in data/ndd.hg38.txt --protein_dir proteins --protein $gene --runs 90000000  --out NDD/$gene.3d.rec.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --males 17750 --females 4262 --pvalues data/ndd_enrich_pvalues.txt
