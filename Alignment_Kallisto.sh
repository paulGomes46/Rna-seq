#!/bin/bash

dataset="scl6"
index="sinensis_v2.idx"

# Initialize the array
read_prefix=()
index_dir="/home/paul/rna_seq/genome_index"
output_dir="/home/paul/rna_seq/data/${dataset}/alignment"
trimmed_sequences="/home/paul/rna_seq/data/${dataset}/trimmed_sequences"

# Retrieve the prefix of all the files in the specified folder
for file in "${trimmed_sequences}"/*.fastq.gz; do
# Extract the first 6 characters as the sample name
    sample_name=$(basename "$file" | cut -c1-11)
    # Add to the array if it's not already present
    if [[ ! ${read_prefix[@]} =~ $sample_name ]]
    then
        read_prefix+=($sample_name)
        echo $sample_name
    fi
done


# Kallisto 
for element in ${read_prefix[@]}
do
        echo "start : ${element}"
        kallisto quant \
        -b 100 \
        -i ${index_dir}/${index} \
        -o ${output_dir}/${element} \
        ${trimmed_sequences}/${element}_1* ${trimmed_sequences}/${element}_2*
done


