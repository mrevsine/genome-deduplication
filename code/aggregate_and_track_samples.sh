#!/usr/bin/env bash

# Produces a corresponding txt file for each input bed
# e.g. for the above example, produces human.dev.txt and human.train.txt


# module load bedtools  # Comment out if not on HPC system

#in_bed=tests/out/dummy/dummy_1.samples.bed
#in_fa=tests/test-data/dummy_1.fa
#out_file=foo.txt

#in_bed=tests/out/small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.bed
#in_fa=/scratch4/mschatz1/mrevsin1/ncbi_genomes/genomes/GCA_001775245.1_ASM177524v1_genomic.fna.gz
#out_file=tests/out_small_input/GCA_001775245.1_ASM177524v1_genomic.samples.dev.txt

in_dir=$1
random_seed=$2
no_txt=${3:-false}

train_bed=${in_dir}/all_train.bed
dev_bed=${in_dir}/all_dev.bed
train_txt=${in_dir}/all_train.txt
dev_txt=${in_dir}/all_dev.txt
touch $train_bed
touch $dev_bed

train_files=$(ls ${in_dir}/*.samples.train.bed)
dev_files=$(ls ${in_dir}/*.samples.dev.bed)

#For each file, getfasta, then pipe that txt to a bedfile with [chrom, start, end, sample_id, sequence]
for file in $train_files; do
    file_basename=$(basename "$file" .samples.train.bed)
    fasta_file=$(grep $file_basename ${in_dir}/basename_fasta_match.txt | cut -f2)
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab > $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $train_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $train_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab > $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $train_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $train_bed
    fi
done

for file in $dev_files; do
    file_basename=$(basename "$file" .samples.dev.bed)
    fasta_file=$(grep $file_basename ${in_dir}/basename_fasta_match.txt | cut -f2)
    temp_bed=$(mktemp)
    if [[ "$fasta_file" == *.gz ]]; then
        temp_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
        gunzip -c "$fasta_file" > "$temp_fasta"
        bedtools getfasta -fi "$temp_fasta" -bed "$file" -tab >> $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $dev_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $dev_bed
        rm "$temp_fasta"
    else
        bedtools getfasta -fi "$fasta_file" -bed "$file" -tab >> $temp_bed
        if [ "$no_txt" = false ]; then
            cut -f2 $temp_bed | tr '[:lower:]' '[:upper:]' >> $dev_txt
        fi
        awk -v basename="$file_basename" '{print $0 "\t" basename}' $file >> $dev_bed
    fi
done

#Now shuffle the train and dev bed files
shuf --random-source="$train_bed" -o ${train_bed}.shuf $train_bed
if [ "$no_txt" = false ]; then
    shuf --random-source="$train_bed" -o ${train_txt}.shuf $train_txt
fi
shuf --random-source="$dev_bed" -o ${dev_bed}.shuf $dev_bed
if [ "$no_txt" = false ]; then
    shuf --random-source="$dev_bed" -o ${dev_txt}.shuf $dev_txt
fi
mv ${train_bed}.shuf $train_bed
mv ${dev_bed}.shuf $dev_bed
if [ "$no_txt" = false ]; then
    mv ${train_txt}.shuf $train_txt
    mv ${dev_txt}.shuf $dev_txt
fi
