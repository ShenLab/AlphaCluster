![travis](https://travis-ci.org/jeremymcrae/denovonear.svg?branch=master)

<p align="center">
  <img src="img/img1.gif">
</p>

### 3dnear

This codebase determines the statistical significance of clustering of
variants on coding sequence of a gene. It is designed to be used as part
of a disease cohort study, in which variants are taken from the affected
individuals.

This codebase was built on top of Jeremy McRae's `denovonear`, which tests
clustering for de novo single-nucleotide variants, where clustering is
measured by 1 dimmensional level (the distance metric is base-pair
seperation). Our code greatly extends the initial test by:
- incorporating 3-dimmensional clustering analysis
- incorporating co-evolutionary clustering analysis
- extending to include inherited variants

### Install
```sh
pip install denovonear
```


### Usage
The core usage options are `cluster`, `cluster-1d` and (pending) `cluster-coev`.
The user can run these commands to perform the cluster analysis in 3-dimmensional
space, 1-dimmensional space, or "co-evolution" space.

#### 3D space
The 3-dimmensional analysis of the package is acheived in the following command:
```sh
denovonear cluster \
   --in data/example_de_novos.txt \
   --protein_dir protein
   --protein CDH8
   --out output.txt
```

This will run the 3d clustering analysis for the given protein (CHD8) with
xyz coordinates for each amino acid in the `protein_dir` folder. The system
assumes that there will be an `.xyz` file named `CHD8.xyz` in the `protein_dir`
folder.

#### 1D space
The 1 dimmensional analysis is acheived via the following command:
```sh
denovonear cluster-1d \
   --in data/example_de_novos.txt \
   --protein CDH8
   --out output.txt
```

#### Coevolutionary space
The coevolutionary analysis is acheived via the following command:
```sh
denovonear cluster-coevol \
   --in data/example_de_novos.txt \
   --coevol_dir coevol
   --protein CDH8
   --out output.txt
```

#### Output

The tab-separated output file will contain one row per gene/transcript, with
each line containing a transcript ID or gene symbol, a log10 transformed
missense mutation rate, a log10 transformed nonsense mutation rate, and a log10
transformed synonymous mutation rate.

### Helper scripts

#### PDB to XYZ file
