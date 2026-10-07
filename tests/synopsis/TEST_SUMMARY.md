# Synopsis Test Summary

## Run the Suite

From the repository root:

```sh
tests/synopsis/run.sh
```

The runner executes every `test_*.sh` file in `tests/synopsis/`.

## Real Deduplication Output

`test_outputs.sh` runs `code/dedup6` directly, without `synopsis.sh`, for each GCA input:

```text
tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz
tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz
```

It runs two parameter sets:

1. Default `per_kmer` execution.
2. Full parameterized execution using `-d`, `-k`, `-l`, `-m`, `-o`, `-p`, `-s`, `-v`, `-r`, `--save_kmers_at_end`, all three BED-writing flags, and `--seed`.

Each parameter set is run in two input modes:

- Once per GCA FASTA file.
- Once with `tests/test-data/small_input.txt` passed directly as the complete file list.

The list-input runs produce combined output directories and verify every generated file, including one samples BED per listed GCA file.

Per-sample evaluation is also run directly with `-e per_sample_threshold` for the same per-file and full-list cases. These outputs are stored under `sample_default`, `sample_parameterized`, `list_sample_default`, and `list_sample_parameterized`.

Commands used to generate and verify the per-file default fixtures:

```sh
./code/dedup6 -e per_kmer \
	-o tests/synopsis/fixtures/expected_outputs/default/GCA_001775245.1_ASM177524v1_genomic \
	tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
```

The same command is run for these two additional inputs, with the corresponding literal output directory:

```sh
./code/dedup6 -e per_kmer \
	-o tests/synopsis/fixtures/expected_outputs/default/GCA_001786795.1_ASM178679v1_genomic \
	tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz

./code/dedup6 -e per_kmer \
	-o tests/synopsis/fixtures/expected_outputs/default/GCA_001820145.1_ASM182014v1_genomic \
	tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz
```

The expected output is every file in each literal `default/` directory, compared by MD5.

Commands used to generate and verify the per-file parameterized fixtures:

```sh
./code/dedup6 -e per_kmer -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/parameterized/GCA_001775245.1_ASM177524v1_genomic -p /tmp/synopsis-kmers-775245 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz

./code/dedup6 -e per_kmer -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/parameterized/GCA_001786795.1_ASM178679v1_genomic -p /tmp/synopsis-kmers-1786795 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz

./code/dedup6 -e per_kmer -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/parameterized/GCA_001820145.1_ASM182014v1_genomic -p /tmp/synopsis-kmers-1820145 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz
```

Expected output: all BED files, `basename_fasta_match.txt`, `config.json`, and `final.kmers.bin` in the matching parameterized fixture directory, each compared by MD5.

Per-sample commands for the same three GCA inputs:

```sh
./code/dedup6 -e per_sample_threshold -o tests/synopsis/fixtures/expected_outputs/sample_default/GCA_001775245.1_ASM177524v1_genomic tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz

./code/dedup6 -e per_sample_threshold -o tests/synopsis/fixtures/expected_outputs/sample_default/GCA_001786795.1_ASM178679v1_genomic tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz

./code/dedup6 -e per_sample_threshold -o tests/synopsis/fixtures/expected_outputs/sample_default/GCA_001820145.1_ASM182014v1_genomic tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz

./code/dedup6 -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/sample_parameterized/GCA_001775245.1_ASM177524v1_genomic -p /tmp/synopsis-sample-kmers-775245 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz

./code/dedup6 -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/sample_parameterized/GCA_001786795.1_ASM178679v1_genomic -p /tmp/synopsis-sample-kmers-1786795 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz

./code/dedup6 -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/sample_parameterized/GCA_001820145.1_ASM182014v1_genomic -p /tmp/synopsis-sample-kmers-1820145 -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz
```

Per-sample commands for the complete file list:

```sh
./code/dedup6 -e per_sample_threshold -o tests/synopsis/fixtures/expected_outputs/list_sample_default tests/test-data/small_input.txt

./code/dedup6 -e per_sample_threshold -d 0.25 -k 31 -l 100 -m 50 -o tests/synopsis/fixtures/expected_outputs/list_sample_parameterized -p /tmp/synopsis-list-sample-kmers -s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 tests/test-data/small_input.txt
```

Expected output: every file in the corresponding sample fixture directory is compared by MD5.

Commands used for the complete `files.txt` list:

```sh
./code/dedup6 -e per_kmer \
	-o tests/synopsis/fixtures/expected_outputs/list_default \
	tests/test-data/small_input.txt

./code/dedup6 -e per_kmer -d 0.25 -k 31 -l 100 -m 50 \
	-o tests/synopsis/fixtures/expected_outputs/list_parameterized \
	-p /tmp/synopsis-list-kmers \
	-s 10 -v 5 -r --save_kmers_at_end \
	--write_ambiguous_beds --write_ignored_beds --write_masked_beds --seed 42 \
	tests/test-data/small_input.txt
```

Expected output: the five default list files, or the full parameterized list output set, with every file compared by MD5.

K-mer resume is tested with two independent commands:

```sh
rm -rf /tmp/synopsis-resume-source /tmp/synopsis-resume-output

./code/dedup6 -e per_kmer --save_kmers_at_end \
	-o /tmp/synopsis-resume-source \
	tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz

./code/dedup6 -e per_kmer \
	-p /tmp/synopsis-resume-source/final.kmers.bin \
	-o /tmp/synopsis-resume-output \
	tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
```

Expected output: the first command creates a non-empty `final.kmers.bin`; the second command successfully creates the samples BED, basename map, and config files.

Expected outputs are stored in:

```text
tests/synopsis/fixtures/expected_outputs/default/
tests/synopsis/fixtures/expected_outputs/parameterized/
tests/synopsis/fixtures/expected_outputs/list_default/
tests/synopsis/fixtures/expected_outputs/list_parameterized/
tests/synopsis/fixtures/expected_outputs/sample_default/
tests/synopsis/fixtures/expected_outputs/sample_parameterized/
tests/synopsis/fixtures/expected_outputs/list_sample_default/
tests/synopsis/fixtures/expected_outputs/list_sample_parameterized/
```

Each GCA directory contains the expected BED files, `basename_fasta_match.txt`, `config.json`, and, for the parameterized run, `final.kmers.bin`. Every file is compared by MD5. The path-dependent `output_dir` and `seen_kmers` fields in `config.json` are normalized before the file comparison. Parameterized configs are also parsed with `jq` and checked field-by-field for the requested evaluation method, numeric values, paths, booleans, output flags, and seed.

## Wrapper Tests

### `test_direct.sh`

Runs direct mode through `synopsis.sh` with a real GCA input.

Command under test:

```sh
code/synopsis.sh tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
```

Checks:

- Expected stdout: empty.
- Expected dedup6 arguments: `-e per_kmer tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz`.

### `test_individual.sh`

Contains both individual-mode paths because `--config` changes individual mode into batch generation.

Checks without `--config`:

Command under test:

```sh
code/synopsis.sh individual tests/test-data/small_input.txt
```

- All three GCA files are executed.
- Each deduplication invocation has the expected input and per-file output directory.
- Expected stdout: empty.
- Expected dedup6 arguments:

```text
-e per_kmer -o ./GCA_001775245.1_ASM177524v1_genomic /Users/alexanderostrovsky/Desktop/genome-deduplication/tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
-e per_kmer -o ./GCA_001786795.1_ASM178679v1_genomic /Users/alexanderostrovsky/Desktop/genome-deduplication/tests/test-data/GCA_001786795.1_ASM178679v1_genomic.fna.gz
-e per_kmer -o ./GCA_001820145.1_ASM182014v1_genomic /Users/alexanderostrovsky/Desktop/genome-deduplication/tests/test-data/GCA_001820145.1_ASM182014v1_genomic.fna.gz
```

Invalid-command checks:

```sh
code/synopsis.sh individual -y tests/test-data/small_input.txt
code/synopsis.sh batch tests/test-data/small_input.txt
```

Expected output: nonzero exit status; the first reports `-y requires --config`, and the second reports `batch` as unrecognized.

Checks with `--config`:

Empty-config command:

```sh
code/synopsis.sh individual \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds \
	--write_ignored_beds --write_masked_beds --seed 42 \
	--config tests/synopsis/fixtures/slurm/slurm_empty.conf tests/test-data/small_input.txt
```

Full-config command:

```sh
code/synopsis.sh individual \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -v 5 -r --save_kmers_at_end --write_ambiguous_beds \
	--write_ignored_beds --write_masked_beds --seed 42 \
	--config tests/synopsis/fixtures/slurm/slurm_all.conf tests/test-data/small_input.txt
```

- An empty Slurm config produces only the computed array directive and command setup.
- A full Slurm config emits every supported `#SBATCH` variable.
- Setup code from the config appears before the generated command.
- Every shared deduplication parameter appears in the generated command in the expected order.
- Expected output: `Generated: synopsis_batch_kmer.sh` and an MD5 match for the complete generated script.

Submission command:

```sh
code/synopsis.sh individual --config tests/synopsis/fixtures/slurm/slurm_all.conf \
	-y tests/test-data/small_input.txt
```

Expected scheduler command: `sbatch synopsis_batch_kmer.sh`.

Generated scripts are compared byte-for-byte through MD5 against the expected script constructed by the test.

### `test_shuffle.sh`

Runs shuffle mode with the `sample` evaluation method and all shuffle flags:

```sh
code/synopsis.sh shuffle sample \
	-d 0.25 -k 31 -l 100 -m 50 -o /tmp/test-output -p /tmp/test-seen \
	-s 10 -r -i -z -x tests/test-data/small_input.txt /tmp/test-input-dir
```

```text
-i -z -x
```

Expected output: empty wrapper stdout and exact fake-shuffle arguments ending in `-e per_sample_threshold`, the list path, and the input directory.

### `test_overwrite.sh`

Creates an existing output directory containing a sentinel file, then runs:

```sh
printf '%s\n' 'n' | ./code/dedup6 -e per_kmer \
	-o /tmp/synopsis-overwrite-output \
	tests/test-data/GCA_001775245.1_ASM177524v1_genomic.fna.gz
```

Expected output: `dedup6` prints its overwrite prompt, exits with status 1, and leaves the existing sentinel file unchanged.

### `run.sh`

Runs every subcommand test file and requires each test process to complete successfully.

Command:

```sh
tests/synopsis/run.sh
```

Expected output: every test reports `PASS` and the command exits with status 0.

## Comparison Locations

When a test runs, temporary actual outputs and command logs are created under a system temporary directory such as:

```text
/var/folders/.../synopsis-tests.XXXXXX/
```

Those files are removed when the test exits. For persistent expected results, inspect:

- [default expected outputs](fixtures/expected_outputs/default)
- [parameterized expected outputs](fixtures/expected_outputs/parameterized)
- [Slurm fixtures](fixtures/slurm)
- [real-output test](test_outputs.sh)
- [individual and batch-wrapper test](test_individual.sh)
