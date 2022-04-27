#!/bin/bash

#SBATCH -J ASD_POWER_C
#SBATCH -o log/POWER_C."%j".out
#SBATCH -e log/POWER_C."%j".err
#SBATCH --mem=16G
gene=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $1)
size=$2
dnv_file=$3
enrich_pvalues_file=$4
disease=$5
tag=$6
echo $gene
echo $disease
echo $SLURM_JOB_ID
echo $gene >&2
echo "${SLURM_ARRAY_TASK_ID}" >&2
echo $SLURM_JOB_ID >&2

#select subset of samples and variants
mkdir data/$disease.power_dnvs_$tag
mkdir $disease.power_dnvs_$tag
mkdir $disease.power_dnvs_$tag/$gene

grep -e $gene $dnv_file | grep -e "missense" | shuf -n $size >  $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt


# run tests on these variants
#1D
python denovonear/__main__.py cluster_1d --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 10000000  --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.1d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#1D with DenovoWEST enrichment
if [[ $size == 31783 || $size == 23425 ]];
then
       python denovonear/__main__.py cluster_1d --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein $gene --runs 10000000 --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.1d.enrich.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa
fi

#1D with CADD
python denovonear/__main__.py cluster_1d --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --protein $gene --runs 10000000  --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.1d.CADD.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

#3D
python denovonear/__main__.py cluster --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.3d.$SLURM_JOB_ID.out --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa


#3D with CADD score scaling
python denovonear/__main__.py cluster --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000 --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.3d.CADD.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa

# Thresholding runs

#3D with CADD threshold 15
python denovonear/__main__.py cluster --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000  --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold15.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 15

#3D with CADD threshold 20
python denovonear/__main__.py cluster --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000  --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold20.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 20

#3D with CADD threshold 25
python denovonear/__main__.py cluster --in $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.released.hg38.$SLURM_JOB_ID.txt --protein_dir proteins --protein $gene --runs 10000000  --out $disease.power_dnvs_$tag/$gene/$size.$gene.$disease.3d.CADD.threshold25.$SLURM_JOB_ID.out --scores  /share/terra/rsrc/hg38/cadd/cadd_1.6/whole_genome_SNVs.tsv.gz --gencode data/gencode.v38.annotation.gtf.gz --fasta data/genome.hg38rg.fa --threshold 25


