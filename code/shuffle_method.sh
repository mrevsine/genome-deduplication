#!/usr/bin/env bash

ignored=false
masks=false
add_originally_ignored=false

# dedup defaults
evaluation_method="per_kmer"
agnostic_retain_pct="0.0"
kmer="32"
sample_len="1000"
min_sample_len=""
outdir=""
seen_kmers=""
per_kmer_retain_pct="0.0"
ambiguous_base_threshold="0.0"
duplicate_base_threshold="0.0"
save_every="0"
no_overlap=false
no_save_kmers_at_end=false
print_config=false
no_write_config=false
seed="123"

usage() {
    echo "Usage: $0 [wrapper options] [dedup options] file_list indir"
    echo ""
    echo "Wrapper options:"
    echo "  -i                      generate final ignored bed files"
    echo "  -z                      generate final masked bed files"
    echo "  -x                      add originally ignored files"
    echo ""
    echo "Dedup options:"
    echo "  -e <str>                Evaluation method: per_kmer | per_sample_agnostic | per_sample_threshold"
    echo "  -a <float>              Likelihood to retain a sample with any duplication/ambiguity in per_sample_agnostic mode [0,1]"
    echo "  -k <int>                K-mer size [1,32]"
    echo "  -l <int>                Sample length"
    echo "  -m <int>                Minimum sample length"
    echo "  -o <path>               Output directory"
    echo "  -p <path>               Pickle file containing seen kmers"
    echo "  -r <float>              Likelihood to allow duplicate kmer [0,1] in per_kmer mode"
    echo "  -b <float>              Allowed ambiguous base threshold in per_sample_threshold mode [0,1]"
    echo "  -d <float>              Allowed duplicate base threshold in per_sample_threshold mode [0,1]"
    echo "  -s <int>                Save seen kmers every n samples"
    echo "  --no-overlap            Keep neighboring samples discrete"
    echo "  --no-save-kmers-at-end  Don't save kmers at program end"
    echo "  --print-config          Print run arguments as JSON to stdout"
    echo "  --no-write-config       Do not write run arguments to config.json"
    echo "  --seed <int>            Random seed"
    echo "  -h, --help              Show this help message"
}

# handle long options first
args=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-overlap)
            no_overlap=true
            shift
            ;;
        --no-save-kmers-at-end)
            no_save_kmers_at_end=true
            shift
            ;;
        --print-config)
            print_config=true
            shift
            ;;
        --no-write-config)
            no_write_config=true
            shift
            ;;
        --help)
            usage
            exit 0
            ;;
        --*)
            echo "Invalid option: $1" >&2
            usage >&2
            exit 1
            ;;
        *)
            args+=("$1")
            shift
            ;;
    esac
done

set -- "${args[@]}"

while getopts ":izxe:a:k:l:m:o:p:r:b:d:s:h" opt; do
    case $opt in
        i)
            ignored=true
            ;;
        z)
            masks=true
            ;;
        x)
            add_originally_ignored=true
            ;;
        e)
            evaluation_method="$OPTARG"
            ;;
        a)
            agnostic_retain_pct="$OPTARG"
            ;;
        k)
            kmer="$OPTARG"
            ;;
        l)
            sample_len="$OPTARG"
            ;;
        m)
            min_sample_len="$OPTARG"
            ;;
        o)
            outdir="$OPTARG"
            ;;
        p)
            seen_kmers="$OPTARG"
            ;;
        r)
            per_kmer_retain_pct="$OPTARG"
            ;;
        b)
            ambiguous_base_threshold="$OPTARG"
            ;;
        d)
            duplicate_base_threshold="$OPTARG"
            ;;
        s)
            save_every="$OPTARG"
            ;;
        h)
            usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage >&2
            exit 1
            ;;
    esac
done

shift $((OPTIND - 1))

file_list="$1"
indir="$2"

if [ -z "$file_list" ] || [ -z "$indir" ]; then
    usage >&2
    exit 1
fi

if [ -z "$min_sample_len" ]; then
    min_sample_len="$sample_len"
fi

outdir="$indir/final_files"
combined_file_name="combined_deduped_samples.fasta"
CURRENT_DEDUP_COMMAND="dedup4"
ROOTDIR="$(pwd)"


if [ -z "$file_list" ] || [ -z "$indir" ]; then
    usage >&2
    exit 1
fi


# Check if files have already been deduplicated
while IFS= read -r file; do
    base_name="$(basename "$file")"
    base_name="${base_name%.*}"

    deduped_subdir="$indir/$base_name"
    deduped_file="$deduped_subdir/$base_name.samples.bed"
    echo $deduped_file
    if [ ! -f "$deduped_file" ]; then
        echo "Deduplicated file for $base_name not found. Please run deduplication step first."
        exit 1
    fi
done < "$file_list"

# Merge deduplicated bed files for each sample
while IFS= read -r file; do
    base_name="$(basename "$file")"
    base_name="${base_name%.*}"

    deduped_subdir="$indir/$base_name"
    fasta_match="$deduped_subdir/basename_fasta_match.txt"
    fasta_file="$(cut -f2 "$fasta_match")"
    deduped_file="$deduped_subdir/$base_name.samples.bed"
    echo "Merging deduplicated bed file for $base_name"
    if [ $add_originally_ignored = true ]; then
        ignored_file="$deduped_subdir/$base_name.ignored.bed"
        cat "$deduped_file" > "$deduped_subdir/$base_name.samples_with_ignored.bed"
        cat "$ignored_file" >> "$deduped_subdir/$base_name.samples_with_ignored.bed"
        bedtools sort -i "$deduped_subdir/$base_name.samples_with_ignored.bed" > "$deduped_subdir/$base_name.sorted_samples_with_ignored.bed"
        bedtools merge -i "$deduped_subdir/$base_name.sorted_samples_with_ignored.bed" > "$indir/$base_name.merged.bed"
        bedtools getfasta -fi "$fasta_file" -bed "$indir/$base_name.merged.bed" -fo "$indir/$base_name.merged.fasta"
        rm "$deduped_subdir/$base_name.samples_with_ignored.bed" "$deduped_subdir/$base_name.sorted_samples_with_ignored.bed"
    else
        bedtools merge -i "$deduped_file" > "$indir/$base_name.merged.bed"
        bedtools getfasta -fi "$fasta_file" -bed "$indir/$base_name.merged.bed" -fo "$indir/$base_name.merged.fasta"
    fi
    
done < "$file_list"

# Combine all merged fasta files into a single file
: > "$indir/$combined_file_name"
for file in "$indir"/*.merged.fasta; do
    [ -e "$file" ] || continue
    sed "s/^>/>$(basename "$file")::/" "$file" >> "$indir/$combined_file_name"
done

# Shuffle the combined file
sh "$ROOTDIR/code/shuffle_awk.sh" "$indir/$combined_file_name" "$indir"
rm -f "$indir/all_sequences_labeled.txt"





# Deduplicate the shuffled file
dedup_args=(
    -e "$evaluation_method"
    -a "$agnostic_retain_pct"
    -k "$kmer"
    -l "$sample_len"
    -m "$min_sample_len"
    -o "$outdir"
    -r "$per_kmer_retain_pct"
    -b "$ambiguous_base_threshold"
    -d "$duplicate_base_threshold"
    -s "$save_every"
    --seed "$seed"
)

if [ -n "$seen_kmers" ]; then
    dedup_args+=( -p "$seen_kmers" )
fi

if [ "$no_overlap" = true ]; then
    dedup_args+=( --no-overlap )
fi

if [ "$no_save_kmers_at_end" = true ]; then
    dedup_args+=( --no-save-kmers-at-end )
fi

if [ "$print_config" = true ]; then
    dedup_args+=( --print-config )
fi

if [ "$no_write_config" = true ]; then
    dedup_args+=( --no-write-config )
fi

"$ROOTDIR/code/$CURRENT_DEDUP_COMMAND" "${dedup_args[@]}" "$indir/shuffled.fasta"

# Re-sort to individual files based on sample name
while IFS= read -r line; do
    file_part=${line%%::*}
    rest=${line#*::}

    # remove .merged.fasta from the filename
    file="${file_part%.merged.fasta}"

    contig=${rest%%:*}
    rest=${rest#*:}

    realstart=${rest%%-*}
    rest=${rest#*-}

    realend=${rest%%$'\t'*}
    rest=${rest#*$'\t'}
    sequencestart=${rest%%$'\t'*}
    sequenceend=${rest#*$'\t'}

    true_start=$((realstart + sequencestart))
    true_end=$((realstart + sequenceend))

    outfile="$outdir/$file.final.samples.bed"
    printf "%s\t%s\t%s\n" "$contig" "$true_start" "$true_end" >> "$outfile"
done < "$outdir/shuffled.samples.bed"


if [ "$ignored" = true ]; then
    echo "Generating final ignored bed files..."
    while IFS= read -r line; do
        file_part=${line%%::*}
        rest=${line#*::}

        # remove .merged.fasta from the filename
        file="${file_part%.merged.fasta}"

        contig=${rest%%:*}
        rest=${rest#*:}

        realstart=${rest%%-*}
        rest=${rest#*-}

        realend=${rest%%$'\t'*}
        rest=${rest#*$'\t'}
        sequencestart=${rest%%$'\t'*}
        sequenceend=${rest#*$'\t'}

        true_start=$((realstart + sequencestart))
        true_end=$((realstart + sequenceend))

        outfile="$outdir/$file.final.ignored.bed"
        printf "%s\t%s\t%s\n" "$contig" "$true_start" "$true_end" >> "$outfile"
    done < "$outdir/shuffled.ignored.bed"
fi

if [ "$masks" = true ]; then
    echo "Generating final masked bed files..."
    while IFS= read -r line; do
        file_part=${line%%::*}
        rest=${line#*::}

        # remove .merged.fasta from the filename
        file="${file_part%.merged.fasta}"

        contig=${rest%%:*}
        rest=${rest#*:}

        realstart=${rest%%-*}
        rest=${rest#*-}

        realend=${rest%%$'\t'*}
        rest=${rest#*$'\t'}
        sequencestart=${rest%%$'\t'*}
        sequenceend=${rest#*$'\t'}

        true_start=$((realstart + sequencestart))
        true_end=$((realstart + sequenceend))

        outfile="$outdir/$file.final.masks.bed"
        printf "%s\t%s\t%s\n" "$contig" "$true_start" "$true_end" >> "$outfile"
    done < "$outdir/shuffled.masks.bed"
fi