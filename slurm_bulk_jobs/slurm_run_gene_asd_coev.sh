#!/bin/bash

#SBATCH -J ASD_3D
#SBATCH -o log/ASD_3D."%j".out
#SBATCH -e log/ASD_3D."%j".err
#SBATCH --mem=16G

gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
echo $gene
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

python alphacluster/__main__.py cluster_coev --in data/asd.released.hg38.txt --coev_dir coevol --protein $gene  --out ASD_coev/$gene.coev.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --run 900000
