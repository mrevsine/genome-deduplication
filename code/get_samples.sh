#!/usr/bin/env bash


# Input arg
data_folder=$1

# Global args
dev_pct=0.1
max_dev_samples=100000
shuf_seed=123


while getopts ":i:p:m:s:" opt; do
    case $opt in
        i)
            data_folder="$OPTARG"
            ;;
        p)
            dev_pct="$OPTARG"
            ;;
        m)
            max_dev_samples="$OPTARG"
            ;;
        s)
            shuf_seed="$OPTARG"
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            echo "Usage: $0 -i input_folder -p dev_pct -m max_dev_samples -s seed"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            echo "Usage: $0 -i input_folder -p dev_pct -m max_dev_samples -s seed"
            exit 1
            ;;
    esac
done

# Check if all required arguments are provided and not empty
if [ -z "$data_folder" ] || [ -z "$dev_pct" ] || [ -z "$max_dev_samples" ] || [ -z "$shuf_seed" ]; then
    echo "Error: All arguments (-i, -p, -m, -s) are required and cannot be empty."
    echo "Usage: $0 -i input_folder -p output_folder -m max_dev_samples -s seed"
    exit 1
fi


#Check if paritioned files exist already and check for confirmation to overwrite
if ls ${data_folder}/*.samples.dev.bed 1> /dev/null 2>&1 || ls ${data_folder}/*.samples.train.bed 1> /dev/null 2>&1; then
    read -p "Partitioned sample files already exist. Do you want to overwrite them? (y/n) " choice
    case "$choice" in
        y|Y ) echo "Overwriting existing partitioned sample files..." && rm ${data_folder}/*.samples.dev.bed ${data_folder}/*.samples.train.bed;;
        n|N ) echo "Exiting script."; exit 1;;
        * ) echo "Invalid choice. Exiting script."; exit 1;;
    esac
fi

#Check if final output files exist already and check for confirmation to overwrite
if [ -f "${data_folder}/all_train.bed" ] || [ -f "${data_folder}/all_dev.bed" ]; then
    read -p "Combined and shuffled files already exist. Do you want to overwrite them? (y/n) " choice
    case "$choice" in
        y|Y ) echo "Overwriting existing files..." && rm "${data_folder}/all_train.bed" "${data_folder}/all_dev.bed";;
        n|N ) echo "Exiting script."; exit 1;;
        * ) echo "Invalid choice. Exiting script."; exit 1;;
    esac
fi

## Create train/dev txt files for each sample in the input directory

# Get all samples.bed files in the input directory
sample_beds=$(ls $data_folder/*.samples.bed)

for sample_bed in $sample_beds; do
	
	# Get all needed file names
	basename=$(basename $sample_bed | rev | cut -d'.' -f3- | rev)
	sample_fasta=$(grep "^${basename}" ${data_folder}/basename_fasta_match.txt | cut -f2) # query corresponding fasta for this sample from basename_fasta_match.txt
	dev_sample_bed=${data_folder}/${basename}.samples.dev.bed
	train_sample_bed=${data_folder}/${basename}.samples.train.bed

	# Get training/dev beds/txts for each sample bed
	echo "Partitioning samples"
	./code/partition_samples.sh $sample_bed $dev_pct $max_dev_samples $shuf_seed
done
echo "Finished partitioning samples into train/dev sets. Aggregating and retrieving sequences now"
./code/aggregate_and_track_samples.sh $data_folder $shuf_seed 
