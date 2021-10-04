# 1) Download AlphaFold2 files
wget https://ftp.ebi.ac.uk/pub/databases/alphafold/UP000005640_9606_HUMAN.tar

# 2) Extract PDB.GZ files to protein folder
mkdir ../proteins
tar -xf UP000005640_9606_HUMAN.tar -C ../proteins
rm ../proteins/*.cif.gz

# 3) Strip down PDG.GZ files to be of form "uniport.pdb.gz" 
cd ../proteins
rename 's/AF-//' *.pdb.gz
rename 's/-F*-model_v*//' *.pdb.gz
cd ../scripts

# 4) Rename PDB.GZ files from "uniport.pdb.gz" to "hgnc_symbol.pdb.gz
python uniport_to_gene.py

# 5) Generate xyz files from pdb.gz files
generate_xyz_files.sh
