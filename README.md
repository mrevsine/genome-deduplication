# Synopsis

`synopsis` is the command-line entry point for running FASTA deduplication. It supports running a single dataset, running many files independently in sequence, running many files independently in parallel on a Slurm cluster, or merging/shuffling previously deduplicated files together.

## Requirements

- Bash 4+ or Zsh
- bedtools
- C++17 compiler and zlib

## Usage

```
synopsis [kmer|sample-threshold] [options] <input file|list>
synopsis individual [kmer|sample-threshold] [options] <input file|list>
synopsis batch [kmer|sample-threshold] --config <config_file> [-y] [options] <input file|list>
synopsis shuffle [kmer|sample-threshold] [options] <file_list> <indir>
```

### Modes

| Mode | Description |
|---|---|
| **direct** (default) | Runs deduplication once on the given input. If the input is a `.txt` list of fasta paths, it's treated as a single combined dataset. |
| **individual** | Runs deduplication independently on each file in a list, one after another. Each file gets its own output subdirectory. |
| **batch** | Generates a Slurm array job that runs deduplication independently on each file in a list, in parallel on a cluster. Requires `--config`. Can optionally submit the job immediately with `-y`. |
| **shuffle** | Merges adjacent samples from a directory of previously, independently deduplicated files, shuffles them together, re-deduplicates, and generates new samples. |

### Choosing a mode and evaluation method

The first one or two arguments (before any flags) are positional:

1. If the first argument is a mode word (`individual`, `batch`, `shuffle`), it's used as the mode.
2. The next argument, if it matches a valid evaluation method (`kmer` or `sample-threshold`), is used as the evaluation method.
3. If no evaluation method is given, it defaults to `kmer`.
4. Whatever comes next must either be a flag (starts with `-`) or an existing input file with a recognized extension — otherwise `synopsis` exits with an error (this catches typos in a misspelled mode or eval-method word).

`sample-agnostic` is no longer supported as an evaluation method.

Examples:

```sh
synopsis myfile.fasta                      # direct mode, kmer
synopsis sample-threshold myfile.fasta     # direct mode, sample-threshold
synopsis individual files.txt              # individual mode, kmer
synopsis individual sample-threshold files.txt
synopsis batch files.txt --config cfg.conf # batch mode, kmer
synopsis shuffle sample-threshold list.txt indir/
```

## Input detection

The final positional argument is inspected by extension:

| Extension | Treated as |
|---|---|
| `.txt` | A list of fasta file paths, one per line |
| `.fa`, `.fna`, `.fasta`, `.fa.gz`, `.fna.gz`, `.fasta.gz` | A single fasta file |
| anything else (including `.fq`/`.fastq`) | Not supported |

When a `.txt` list is given to `individual` or `batch` mode:

- Blank lines are skipped.
- Lines beginning with `#` are treated as comments and skipped.
- Every other line is assumed to be a **full, absolute path** to a fasta file. Relative paths are not resolved against anything — use absolute paths.

## Output directories

`-o <path>` sets the base output directory and is shared across all modes.

- In **direct** mode, `-o` is used exactly as given.
- In **individual** and **batch** modes, each fasta file gets its own subdirectory under `-o`, named after that file's basename with its extension(s) stripped (including both parts of a `.gz`-compressed double extension, e.g. `sample1.fasta.gz` → `sample1`):

  ```
  <outdir>/<basename_without_extension>/
  ```

## Dedup options (shared across all modes)

| Flag | Description |
|---|---|
| `-d <float>` | dedup_param |
| `-k <int>` | kmer |
| `-l <int>` | sample_len |
| `-m <int>` | min_sample_len |
| `-o <path>` | output_dir |
| `-p <path>` | seen_kmers |
| `-s <int>` | save_every |
| `-v <int>` | overlap |
| `-r` | retain_info |
| `--save_kmers_at_end` | |
| `--write_ambiguous_beds` | |
| `--write_ignored_beds` | |
| `--write_masked_beds` | |
| `--seed <int>` | |
| `-h`, `--help` | Show usage |

## Shuffle mode options

```
synopsis shuffle [kmer|sample-threshold] [options] <file_list> <indir>
```

| Flag | Description |
|---|---|
| `-i` | Generate final ignored bed files |
| `-z` | Generate final masked bed files |
| `-x` | Add originally ignored regions before merging |

`<file_list>` and `<indir>` specify the previously, independently deduplicated files to merge, shuffle, and re-deduplicate. The evaluation method chosen (default `kmer`) governs the *inner* re-deduplication step. Nested modes are not allowed — you cannot run `synopsis shuffle individual` or `synopsis shuffle batch`.

## Batch mode

```
synopsis batch [kmer|sample-threshold] --config <config_file> [-y] [options] <input file|list>
```

Batch mode resolves the input into a list of fasta files and generates a Slurm script with one array task per file. Each array task deduplicates a single assigned fasta file, writing to its own output subdirectory as described above.

### `--config <file>` (required)

Points to a Slurm configuration file (see format below) describing `#SBATCH` settings and any environment setup needed on the cluster (module loads, conda/venv activation, etc).

### `-y`

If given, after generating the Slurm script, `synopsis` additionally submits it immediately via `sbatch`. Without `-y`, the script is only written to disk — nothing is submitted.

### Generated script

The output script is written to the **caller's current working directory** and named:

```
synopsis_batch_<eval_method>.sh
```

e.g. `synopsis_batch_kmer.sh` or `synopsis_batch_sample-threshold.sh`.

The `--array` range is computed automatically based on the number of resolved input files. Set `max_concurrent_jobs` to a positive integer to limit how many array tasks run concurrently; it is emitted as the Slurm `%<number>` array limit.

### Slurm config file format

The config file has two sections, each introduced by a literal marker line, with an optional comment line immediately below purely for visual separation:

```
###SLURM VARIABLES###
#_____________________________
job-name=""
account=""
partition=""
time=""
mem=""
cpus-per-task=""
ntasks=""
nodes=""
gres=""
output=""
error=""
mail-user=""
mail-type=""
max_concurrent_jobs=""

###LOCAL SETUP (ENVIRONMENTS, MODULES, ETC)###
#_____________________________
module load anaconda3
source activate my_env
```

**Section 1 — `###SLURM VARIABLES###`:**

- Each line must be of the form `key="value"`.
- `key` must be one of the recognized Slurm long-option names listed above (any other key is an error).
- Any key left as `key=""` (empty value) is skipped — no corresponding `#SBATCH` line is generated.
- Any key with a non-empty value generates a `#SBATCH --<key>="<value>"` line in the output script.
- `max_concurrent_jobs` is a special key: it must be a positive integer and is appended to the computed array range as `%<number>` rather than emitted as its own `#SBATCH` option.
- `array` is **not** a valid key here; it is always computed and injected automatically.

**Section 2 — `###LOCAL SETUP (ENVIRONMENTS, MODULES, ETC)###`:**

- Everything after this marker is copied verbatim into the generated script, and runs once per array task, before deduplication runs on that task's assigned file. You do not need to write any array-handling logic here (no need to reference `$SLURM_ARRAY_TASK_ID`, look up the input file, or build output paths) — `synopsis` generates and appends all of that automatically, immediately after this section. Use this section only for one-time-per-task setup, e.g. module loads, environment activation, exports.
