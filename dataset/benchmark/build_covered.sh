
test

for n in Cancer ASD CHD Control PTEN TP53
do
    python ./build_covered_var_set.py --input ./data/GRCh38_$n.csv
done

python ./build_covered_var_set.py --input ./train/ClinVar/GRCh38_ClinVar.csv
python ./build_covered_var_set.py --input ./train/UniProt/GRCh38_UniProt.csv
python ./build_covered_var_set.py --input ./train/HGMD/GRCh38_HGMD.csv
python ./build_covered_var_set.py --input ./train/DiscovEHR/GRCh38_DiscovEHR.csv
