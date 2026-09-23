#!/usr/bin/env bash

input_fasta=$1
files_dir=$2
awk 'BEGIN {srand()}
/^>/ {
    if (seq != "") {
        printf "%020d\t%s\t%s\n", int(rand()*10^15), header, seq
    }
    header = $0
    seq = ""
    next
}
{
    seq = seq $0
}
END {
    if (seq != "") {
        printf "%020d\t%s\t%s\n", int(rand()*10^15), header, seq
    }
}' $input_fasta > $files_dir/all_sequences_labeled.txt

sort --parallel=8 --buffer-size=80% -k1,1n $files_dir/all_sequences_labeled.txt | \
cut -f2,3 | \
tr '\t' '\n' > $files_dir/shuffled.fasta
