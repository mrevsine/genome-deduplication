#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: individual mode without --config executes all three GCA inputs.
# Expected: no wrapper output and three exact recorded per-file invocations.
output="$(run_synopsis individual "$input_list")"
assert_md5_text "$output" '' 'individual wrapper output'
expected_output=''
while IFS= read -r fasta; do
	base="${fasta##*/}"
	base="${base%.gz}"
	base="${base%.*}"
	line="-e per_kmer -o ./$base $fasta"
	if [ -n "$expected_output" ]; then
		expected_output="$expected_output
$line"
	else
		expected_output="$line"
	fi
done < "$input_list"
assert_md5_text "$(cat "$dedup_log")" "$expected_output" 'individual deduplication arguments'

# Test: -y without --config is rejected.
# Expected: nonzero exit status and an explicit --config requirement.
assert_fails_with '-y requires --config' individual -y "$input_list"
# Test: the removed standalone batch mode is rejected.
# Expected: nonzero exit status identifying batch as an unrecognized argument.
assert_fails_with 'unrecognized argument: batch' batch "$input_list"

printf '%s\n' 'PASS: individual mode'

printf '%s\n' '--- individual with --config: empty Slurm config ---'

# Test: individual --config with no Slurm variables.
# Expected: only the computed array directive, input-selection code, and dedup command.
command_args=(
	-e per_kmer
	-d 0.25
	-k 31
	-l 100
	-m 50
	-p /tmp/test-seen
	-s 10
	-v 5
	-r
	--save_kmers_at_end
	--write_ambiguous_beds
	--write_ignored_beds
	--write_masked_beds
	--seed 42
)

write_expected_script() {
	local expected_script="$1"
	shift
	local has_slurm_lines=0
	[ "$#" -gt 0 ] && has_slurm_lines=1

	{
		printf '%s\n' '#!/bin/bash'
		if [ "$has_slurm_lines" -eq 1 ]; then
			for line in "$@"; do
				printf '%s\n' "$line"
			done
		fi
		if [ "$has_slurm_lines" -eq 1 ]; then
			printf '%s\n' '#SBATCH --array=0-2%2'
		else
			printf '%s\n' '#SBATCH --array=0-2'
		fi
		if [ "$has_slurm_lines" -eq 1 ]; then
			printf '\n'
			printf '%s\n' 'echo setup-from-test-config'
		else
			printf '\n'
		fi
		printf '\n'
		printf '%s\n' '## Using SLURM_ARRAY_TASK_ID to select the input file from the list'
		printf 'FASTA_LIST="%s"\n' "$normalized_input_list"
		printf '%s\n' 'FASTA_PATH=$(awk -v index="$SLURM_ARRAY_TASK_ID" '\''NF && $0 !~ /^[[:space:]]*#/ { if (count++ == index) { print; exit } }'\'' "$FASTA_LIST")'
		printf '%s\n' 'BASENAME=$(basename "$FASTA_PATH")'
		printf '%s\n' 'BASENAME="${BASENAME%.gz}"'
		printf '%s\n' 'BASENAME="${BASENAME%.*}"'
		printf '\n'
		printf '"%s"' "$DEDUP_BIN"
		for arg in "${command_args[@]}"; do
			printf ' "%s"' "$arg"
		done
		printf ' -o "/tmp/test-output/${BASENAME}" "$FASTA_PATH"\n'
	} > "$expected_script"
}

normalized_input_list="$(cd "$(dirname "$input_list")" && pwd)/$(basename "$input_list")"

# Test: empty Slurm config generates the expected batch script.
# Expected: exact MD5 match for the complete generated script.
empty_output="$(run_synopsis individual \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds \
	--write_ignored_beds --write_masked_beds --seed 42 \
	--config "$empty_config_file" "$input_list")"
assert_md5_text "$empty_output" 'Generated: synopsis_batch_kmer.sh' 'empty-config output'
empty_expected="$test_dir/expected_empty_batch.sh"
write_expected_script "$empty_expected"
assert_md5_file "$generated_script" "$empty_expected" 'empty-config generated script'

printf '%s\n' '--- individual with --config: all Slurm variables and setup ---'

# Test: individual --config with every Slurm variable and setup code.
# Expected: all #SBATCH directives, setup code before input selection, and the exact dedup command.
all_output="$(run_synopsis individual \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds \
	--write_ignored_beds --write_masked_beds --seed 42 \
	--config "$all_config_file" "$input_list")"
assert_md5_text "$all_output" 'Generated: synopsis_batch_kmer.sh' 'all-config output'
all_expected="$test_dir/expected_all_batch.sh"
write_expected_script "$all_expected" \
	'#SBATCH --job-name="synopsis-test"' \
	'#SBATCH --account="account-test"' \
	'#SBATCH --partition="partition-test"' \
	'#SBATCH --time="00:05:00"' \
	'#SBATCH --mem="80"' \
	'#SBATCH --cpus-per-task="2"' \
	'#SBATCH --ntasks="1"' \
	'#SBATCH --nodes="1"' \
	'#SBATCH --gres="gpu:1"' \
	'#SBATCH --output="slurm-%A_%a.out"' \
	'#SBATCH --error="slurm-%A_%a.err"' \
	'#SBATCH --mail-user="test@example.com"' \
	'#SBATCH --mail-type="END,FAIL"'
assert_md5_file "$generated_script" "$all_expected" 'all-config generated script'

# Test: -y submits the generated script for the full config.
# Expected: sbatch receives exactly synopsis_batch_kmer.sh.
export SYNOPSIS_SUBMISSION_LOG="$submission_log"
PATH="$test_dir/bin:$PATH" \
	run_synopsis individual --config "$all_config_file" -y "$input_list" >/dev/null
assert_md5_text "$(cat "$submission_log")" 'synopsis_batch_kmer.sh' 'sbatch argument'

printf '%s\n' 'PASS: individual mode with --config'
