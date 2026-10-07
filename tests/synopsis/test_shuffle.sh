#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: shuffle mode executes with sample evaluation, all shared options, and -i/-z/-x.
# Expected: no wrapper output and one exact recorded shuffle invocation.
output="$(run_synopsis shuffle sample \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -r -i -z -x "$input_list" "$test_dir")"
assert_md5_text "$output" '' 'shuffle wrapper output'
expected="-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen -s 10 -r -i -z -x -e per_sample_threshold $input_list $test_dir"
assert_md5_text "$(cat "$shuffle_log")" "$expected" 'shuffle arguments'

printf '%s\n' 'PASS: shuffle mode'
