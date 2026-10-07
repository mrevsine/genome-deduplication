#!/usr/bin/env bash

set -eu

ROOTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
SYNOPSIS="$ROOTDIR/code/synopsis.sh"
TARGET_SHELL="${SYNOPSIS_SHELL:-zsh}"

if ! command -v "$TARGET_SHELL" >/dev/null 2>&1; then
    echo "SKIP: $TARGET_SHELL is not installed" >&2
    exit 0
fi

# Create the isolated fixture directory and command stubs used by each test.
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/synopsis-tests.XXXXXX")"
trap 'rm -rf "$test_dir"' EXIT

input_list="$test_dir/gca_inputs.txt"
empty_config_file="$ROOTDIR/tests/synopsis/fixtures/slurm/slurm_empty.conf"
all_config_file="$ROOTDIR/tests/synopsis/fixtures/slurm/slurm_all.conf"
generated_script="$test_dir/synopsis_batch_kmer.sh"
submission_log="$test_dir/submission.log"
dedup_log="$test_dir/dedup.log"
shuffle_log="$test_dir/shuffle.log"
mkdir "$test_dir/bin"

while IFS= read -r filename; do
    printf '%s/%s\n' "$ROOTDIR" "$filename"
done < "$ROOTDIR/tests/test-data/small_input.txt" > "$input_list"

cat > "$test_dir/bin/fake-dedup" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SYNOPSIS_DEDUP_LOG"
EOF
chmod +x "$test_dir/bin/fake-dedup"

cat > "$test_dir/bin/fake-shuffle" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$*" >> "$SYNOPSIS_SHUFFLE_LOG"
EOF
chmod +x "$test_dir/bin/fake-shuffle"

cat > "$test_dir/bin/sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$1" > "$SYNOPSIS_SUBMISSION_LOG"
EOF
chmod +x "$test_dir/bin/sbatch"

export DEDUP_BIN="$test_dir/bin/fake-dedup"
export SYNOPSIS_DEDUP_LOG="$dedup_log"
export SHUFFLE_SCRIPT="$test_dir/bin/fake-shuffle"
export SYNOPSIS_SHUFFLE_LOG="$shuffle_log"

# Run synopsis from the fixture directory and preserve its exit status/output.
run_synopsis() {
    (cd "$test_dir" && "$TARGET_SHELL" "$SYNOPSIS" "$@")
}

# Abort the current test with a formatted failure message.
fail() {
    echo "FAIL: $*" >&2
    exit 1
}

# Assert that a string contains the expected substring.
assert_contains() {
    local text="$1" expected="$2"
    case "$text" in
        *"$expected"*) ;;
        *) fail "expected output to contain: $expected" ;;
    esac
}

# Assert that two strings are exactly equal and label the comparison.
assert_equal() {
    local actual="$1" expected="$2" label="$3"
    if [ "$actual" != "$expected" ]; then
        printf 'FAIL: %s\nexpected:\n%s\nactual:\n%s\n' "$label" "$expected" "$actual" >&2
        exit 1
    fi
}

# Return the MD5 digest of text using the platform's available MD5 utility.
md5_text() {
    if command -v md5 >/dev/null 2>&1; then
        printf '%s' "$1" | md5 -q
    else
        printf '%s' "$1" | md5sum | awk '{print $1}'
    fi
}

# Assert that two text values have identical MD5 digests.
assert_md5_text() {
    local actual="$1" expected="$2" label="$3"
    assert_equal "$(md5_text "$actual")" "$(md5_text "$expected")" "$label"
}

# Assert that two files have identical MD5 digests.
assert_md5_file() {
    local actual_file="$1" expected_file="$2" label="$3"
    local actual_hash expected_hash
    actual_hash="$(md5 -q "$actual_file" 2>/dev/null || md5sum "$actual_file" | awk '{print $1}')"
    expected_hash="$(md5 -q "$expected_file" 2>/dev/null || md5sum "$expected_file" | awk '{print $1}')"
    assert_equal "$actual_hash" "$expected_hash" "$label"
}

# Assert that a command fails and emits the expected diagnostic substring.
assert_fails_with() {
    local expected="$1"
    shift
    local output
    if output="$(run_synopsis "$@" 2>&1)"; then
        fail "expected command to fail: $*"
    fi
    assert_contains "$output" "$expected"
}

# Assert that a string does not contain an unexpected substring.
assert_not_contains() {
    local text="$1" unexpected="$2"
    case "$text" in
        *"$unexpected"*) fail "expected output not to contain: $unexpected" ;;
        *) ;;
    esac
}
