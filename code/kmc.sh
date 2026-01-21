#!/usr/bin/env bash

# Input args
in_bed=$1
kmer_size=$2
threads=$3
mem=$4

# Extrapolate results directory from input bed file
results_dir=$(dirname $in_bed)

# Make sure basename_fasta_match exists
fasta_match_file=$results_dir/basename_fasta_match.txt
if [[ ! -f $fasta_match_file ]]; then
        echo "ERROR: could not find file matching species name to fasta file in the input directory"
        exit 1
fi

# Get corresponding fasta for this .samples.bed file
sample=$(basename "$in_bed" ".samples.bed")
sample_fasta=$(grep -e "^$sample[[:space:]]" "$fasta_match_file" | head -1 | cut -f2)
if [[ -z "$sample_fasta" ]]; then
	echo "ERROR: could not find fasta file for sample $sample"
	exit 1
fi

# Extract sequences from input bed and corresponding fasta
sequences_fasta=$(mktemp -t sample_fasta.XXXXXX).fa
if [[ "$sample_fasta" == *.gz ]]; then
	#echo "Processing sample kmers for $sample"
	unzipped_sample_fasta=$(mktemp -t temp_fasta.XXXXXX).fa
	gunzip -c "$sample_fasta" > "$unzipped_sample_fasta"
	bedtools getfasta -fi "$unzipped_sample_fasta" -bed "$in_bed" >> $sequences_fasta
	rm "$unzipped_sample_fasta"
else
	#echo "Processing sample kmers for $sample"
	bedtools getfasta -fi "$sample_fasta" -bed "$in_bed" >> $sequences_fasta
fi

# Run KMC on sequences
working_dir=$results_dir/kmc
if [[ ! -d "$working_dir" ]]; then
	mkdir -p $working_dir
fi
kmc -k${kmer_size} -cs2 -fa -b -ci1 -t${threads} -m${mem} $sequences_fasta $working_dir/$sample $working_dir
kmc_dump $working_dir/$sample $working_dir/${sample}.kmers.txt

# Final clean up
rm $sequences_fasta

