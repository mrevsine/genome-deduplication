#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: direct mode executes the deduplication binary with the input path.
# Expected: no wrapper output and one exact recorded deduplication invocation.
output="$(run_synopsis "$ROOTDIR/tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz")"
assert_md5_text "$output" '' 'direct wrapper output'
expected="-e per_kmer $ROOTDIR/tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz"
assert_md5_text "$(cat "$dedup_log")" "$expected" 'direct deduplication arguments'

printf '%s\n' 'PASS: direct mode'
