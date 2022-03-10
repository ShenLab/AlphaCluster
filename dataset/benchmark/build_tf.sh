
script=./build_tf_dataset.py
output_dir=./tf/v1
mkdir -p ${output_dir}

build() {
	python ${script} --input ./data/GRCh38_${1}.covered.csv --output ${output_dir}/${1} \
		--cpu ${2} --af 1e-3
}


#build BRCA1 1
build TP53 1
build PTEN 1
build CHD 1
build ASD 1
build Control 1

for i in `seq 0 4`
do
	python ${script} --input ./train/list/r${i}.csv --output ${output_dir}/train_r${i} --cpu 10
done
