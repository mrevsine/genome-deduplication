#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: dedup6 halts when an existing output directory receives a non-affirmative answer.
# Expected: exit status 1, overwrite prompt, and unchanged existing contents.
overwrite_dir="$test_dir/overwrite-output"
marker_file="$overwrite_dir/sentinel.txt"
mkdir -p "$overwrite_dir"
printf '%s\n' 'preserve this file' > "$marker_file"

set +e
prompt_output="$(printf '%s\n' 'n' | "$ROOTDIR/code/dedup6" -e per_kmer \
    -o "$overwrite_dir" "$ROOTDIR/tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz" 2>&1)"
exit_code=$?
set -e

assert_equal "$exit_code" 1 'overwrite rejection exit status'
assert_contains "$prompt_output" 'Output directory already exists. Continue and potentially overwrite files? (y/n):'
assert_equal "$(cat "$marker_file")" 'preserve this file' 'existing output preservation'

printf '%s\n' 'PASS: overwrite halt'
