#!/usr/bin/env bash

. "$(dirname "$0")/common.sh"

# Test: direct dedup6 output is reproducible for every real GCA fixture.
# Expected: each run produces the same samples BED, basename map, and config as the checked-in fixtures.
normalize_config() {
    local config_file="$1"
    sed -E \
        -e 's#"output_dir": "[^"]+"#"output_dir": "<OUTPUT_DIR>"#' \
        -e 's#"seen_kmers": "[^"]+"#"seen_kmers": "<SEEN_KMERS>"#' \
        "$config_file"
}

compare_output_tree() {
    local expected_root="$1" actual_root="$2" label="$3"
    local expected_dir
    for expected_dir in "$expected_root"/*; do
        fixture_name="$(basename "$expected_dir")"
        [ "$fixture_name" = seen ] && continue
        actual_dir="$actual_root/$fixture_name"

        for expected_file in "$expected_dir"/*; do
            filename="$(basename "$expected_file")"
            actual_file="$actual_dir/$filename"
            [ -f "$actual_file" ] || fail "missing $label output: $actual_file"
            if [ "$filename" = config.json ]; then
                assert_md5_text "$(normalize_config "$actual_file")" "$(normalize_config "$expected_file")" "normalized $label config for $fixture_name"
            else
                assert_md5_file "$actual_file" "$expected_file" "$label output $fixture_name/$filename"
            fi
        done
    done
}

compare_flat_output() {
    local expected_dir="$1" actual_dir="$2" label="$3"
    local expected_file filename actual_file
    for expected_file in "$expected_dir"/*; do
        filename="$(basename "$expected_file")"
        actual_file="$actual_dir/$filename"
        [ -f "$actual_file" ] || fail "missing $label output: $actual_file"
        if [ "$filename" = config.json ]; then
            assert_md5_text "$(normalize_config "$actual_file")" "$(normalize_config "$expected_file")" "normalized $label config"
        else
            assert_md5_file "$actual_file" "$expected_file" "$label output $filename"
        fi
    done
}

assert_parameterized_config() {
    local config_file="$1" output_dir="$2" seen_kmers="$3" evaluation_method="$4" label="$5"
    local result
    result="$(jq -e \
        --arg output_dir "$output_dir" \
        --arg seen_kmers "$seen_kmers" \
        --arg evaluation_method "$evaluation_method" \
        '(
            .dedup_param == 0.25 and
            .evaluation_method == $evaluation_method and
            .kmer == 31 and
            .sample_len == 100 and
            .min_sample_len == 50 and
            .output_dir == $output_dir and
            .seen_kmers == $seen_kmers and
            .retain_info == true and
            .save_every == 10 and
            .overlap == 5 and
            .save_kmers_at_end == true and
            .write_ambiguous_beds == true and
            .write_ignored_beds == true and
            .write_masked_beds == true and
            .seed == 42
        )' "$config_file")" || fail "invalid parameterized config: $config_file"
    assert_equal "$result" true "$label config values"
}

run_default_outputs() {
    local actual_root="$test_dir/default_outputs"
    local input_file fixture_name actual_dir
    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        actual_dir="$actual_root/$fixture_name"
        (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" -e per_kmer -o "$actual_dir" "tests/test-data/$fixture_name.fna.gz" >/dev/null)
    done
    compare_output_tree "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/default" "$actual_root" default
}

run_parameterized_outputs() {
    local actual_root="$test_dir/parameterized_outputs"
    local seen_root="$test_dir/parameterized_seen"
    local input_file fixture_name actual_dir
    mkdir -p "$seen_root"
    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        actual_dir="$actual_root/$fixture_name"
        (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
            -e per_kmer -d 0.25 -k 31 -l 100 -m 50 -o "$actual_dir" \
            -p "$seen_root" -s 10 -v 5 -r --save_kmers_at_end \
            --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 \
            "tests/test-data/$fixture_name.fna.gz" >/dev/null)
    done
    compare_output_tree "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/parameterized" "$actual_root" parameterized
    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        assert_parameterized_config "$actual_root/$fixture_name/config.json" "$actual_root/$fixture_name" "$seen_root" per_kmer "per-kmer $fixture_name"
    done
}

run_list_outputs() {
    local actual_root="$test_dir/list_outputs"
    local actual_default="$actual_root/default"
    local actual_parameterized="$actual_root/parameterized"
    local seen_root="$test_dir/list_seen"
    mkdir -p "$seen_root"

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" -e per_kmer -o "$actual_default" tests/test-data/small_input.txt >/dev/null)
    compare_flat_output "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/list_default" "$actual_default" 'list default'

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
        -e per_kmer -d 0.25 -k 31 -l 100 -m 50 -o "$actual_parameterized" \
        -p "$seen_root" -s 10 -v 5 -r --save_kmers_at_end \
        --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 \
        tests/test-data/small_input.txt >/dev/null)
    compare_flat_output "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/list_parameterized" "$actual_parameterized" 'list parameterized'
    assert_parameterized_config "$actual_parameterized/config.json" "$actual_parameterized" "$seen_root" per_kmer 'list per-kmer'
}

run_sample_outputs() {
    local actual_root="$test_dir/sample_outputs"
    local parameterized_root="$test_dir/sample_parameterized_outputs"
    local resume_root="$test_dir/sample_resume"
    local input_file fixture_name actual_dir

    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        actual_dir="$actual_root/$fixture_name"
        (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" -e per_sample_threshold -o "$actual_dir" "tests/test-data/$fixture_name.fna.gz" >/dev/null)
    done
    compare_output_tree "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/sample_default" "$actual_root" 'sample default'

    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        local resume_dir="$resume_root/$fixture_name"
        (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
            -e per_sample_threshold --save_kmers_at_end -o "$resume_dir" \
            "tests/test-data/$fixture_name.fna.gz" >/dev/null)
        actual_dir="$parameterized_root/$fixture_name"
        (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
            -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o "$actual_dir" \
            -p "$resume_dir/final.kmers.bin" -s 10 -v 5 -r --save_kmers_at_end \
            --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 \
            "tests/test-data/$fixture_name.fna.gz" >/dev/null)
    done
    compare_output_tree "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/sample_parameterized" "$parameterized_root" 'sample parameterized'
    for input_file in "$ROOTDIR"/tests/test-data/GCA_*.fna.gz; do
        fixture_name="$(basename "$input_file" .fna.gz)"
        assert_parameterized_config "$parameterized_root/$fixture_name/config.json" "$parameterized_root/$fixture_name" "$resume_root/$fixture_name/final.kmers.bin" per_sample_threshold "per-sample $fixture_name"
    done
}

run_list_sample_outputs() {
    local actual_root="$test_dir/list_sample_outputs"
    local actual_default="$actual_root/default"
    local actual_parameterized="$actual_root/parameterized"
    local resume_dir="$test_dir/list_sample_resume"

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" -e per_sample_threshold -o "$actual_default" tests/test-data/small_input.txt >/dev/null)
    compare_flat_output "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/list_sample_default" "$actual_default" 'list sample default'

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" -e per_sample_threshold --save_kmers_at_end -o "$resume_dir" tests/test-data/small_input.txt >/dev/null)
    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
        -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o "$actual_parameterized" \
        -p "$resume_dir/final.kmers.bin" -s 10 -v 5 -r --save_kmers_at_end \
        --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 \
        tests/test-data/small_input.txt >/dev/null)
    compare_flat_output "$ROOTDIR/tests/synopsis/fixtures/expected_outputs/list_sample_parameterized" "$actual_parameterized" 'list sample parameterized'
    assert_parameterized_config "$actual_parameterized/config.json" "$actual_parameterized" "$resume_dir/final.kmers.bin" per_sample_threshold 'list per-sample'
}

run_kmer_resume() {
    local source_dir="$test_dir/resume_source"
    local resumed_dir="$test_dir/resume_output"
    local input_file='tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz'

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
        -e per_kmer --save_kmers_at_end -o "$source_dir" "$input_file" >/dev/null)
    [ -s "$source_dir/final.kmers.bin" ] || fail 'resume checkpoint was not created'

    (cd "$ROOTDIR" && "$ROOTDIR/code/dedup6" \
        -e per_kmer -p "$source_dir/final.kmers.bin" -o "$resumed_dir" "$input_file" >/dev/null)
    [ -f "$resumed_dir/GCA_001775245.1_ASM177524v1_genomic.samples.bed" ] || fail 'resumed samples output was not created'
    [ -f "$resumed_dir/basename_fasta_match.txt" ] || fail 'resumed basename map was not created'
    [ -f "$resumed_dir/config.json" ] || fail 'resumed config was not created'
}

# Test: default dedup6 output is reproducible for every real GCA fixture.
# Expected: each run matches the checked-in three-file default fixture tree.
run_default_outputs

# Test: every shared dedup6 parameter and output-writing flag is reproducible.
# Expected: each run matches the checked-in seven-file parameterized fixture tree.
run_parameterized_outputs

# Test: dedup6 accepts the complete files.txt list directly.
# Expected: combined default and full-parameter output trees match their fixtures.
run_list_outputs

# Test: real per-sample output is reproducible for individual and list inputs.
# Expected: default and full-parameter output trees match their MD5 fixtures.
run_sample_outputs
run_list_sample_outputs

# Test: k-mer resume reloads the checkpoint written by --save_kmers_at_end.
# Expected: final.kmers.bin exists and a second run using -p completes with output files.
run_kmer_resume

printf '%s\n' 'PASS: direct dedup6 output fixtures'
