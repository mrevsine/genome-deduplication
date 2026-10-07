#!/usr/bin/env bash

set -eu

test_dir="$(cd "$(dirname "$0")" && pwd)"
# Test: every subcommand test file completes successfully.
# Expected: each file reports PASS and the runner exits successfully.
for test_file in "$test_dir"/test_*.sh; do
    "$test_file"
done

printf '%s\n' 'PASS: all synopsis subcommand tests'
