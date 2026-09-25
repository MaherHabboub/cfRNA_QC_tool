# Build a Kraken2 database

This standalone workflow downloads references, builds a full Kraken2 database,
and validates it in the same configured storage directory. Run it once before
using the HPC tool's Kraken module. It is not submitted automatically by the QC
pipeline.

## Software, reference contents, and original snapshot

The supplied workflow uses `env/software/doduo` and
`Kraken2/2.1.3-gompi-2023a`. Software loads are in `environment.sh`; adapt those
module names for your cluster. If `dustmasker` is missing, add the appropriate
BLAST+ module to `KRAKEN_EXTRA_MODULES` there. The environment check verifies
commands and the Kraken helper scripts used for library preparation.

The default libraries are **archaea, viral, fungi, human, and bacteria**.
Bacteria use a frozen NCBI RefSeq assembly manifest selecting Complete Genome
and Chromosome assemblies, downloaded over restartable HTTPS. Taxonomy and the
other libraries retain Kraken2's FTP download mode. Low-complexity masking is
enabled. The build uses k-mer length **35**, minimizer length **31**, and **7**
minimizer spaces, without a compact-hash size cap or fast-build mode.

The original database label was
`cfrna_k2_bacteria_archaea_viral_human_fungi_20260828_full`: **28 August 2026** is
its snapshot label, not a verified build-completion date. Historical completion
dates are unavailable in the supplied scripts. A new download uses the reference
data available when it runs; these scripts do not recreate that historical
snapshot from its label alone. Completed bacterial manifests are reused on
restart rather than refreshed.

## Requirements and configuration

Use a Linux HPC system with Bash, Slurm, environment modules, GNU coreutils/findutils,
`awk`, `gzip`, `wget`, `rsync`, and GNU `/usr/bin/time`. Compute nodes need access
to the reference download servers. Check available storage first: this is an
uncapped database, and references, compressed downloads, intermediate files, and
the final index all consume space. Actual storage and memory requirements vary
as reference collections grow.

Edit `config.sh`:

```bash
KRAKEN_ROOT="/absolute/path/to/kraken_storage"
KRAKEN_DB_NAME="cfrna_k2_bacteria_archaea_viral_human_fungi_full"
KRAKEN_SHARED_GROUP=""  # optional group you belong to
```

The example root is intentionally rejected until edited. The workflow derives:

```text
<KRAKEN_ROOT>/database/<KRAKEN_DB_NAME>/   references, index, validation, provenance
<KRAKEN_ROOT>/logs/<KRAKEN_DB_NAME>/       Slurm output and submitted_jobs.tsv
```

Keep the same config when resuming. Use a new database name for a fresh reference
snapshot or different build settings. Do not run two workflows concurrently
against the same database directory. Config files are sourced as shell scripts;
use your own trusted config.

An empty shared group preserves normal ownership. When set, group ownership and
setgid directories allow shared read access. Data files receive mode `0644`;
directories receive `0755`, or `2755` with a shared group. Permissions are applied
only to this database and its log directory.

Resources remain in the batch scripts, not the config:

| Stage | CPUs | Memory | Time limit |
|---|---:|---:|---:|
| Download | 16 | 32G | 3 days |
| Build | 32 | 200G | 3 days |
| Validation | 4 | 200G | 1 day |

The bacterial downloader limits concurrent transfers to eight even when more
CPUs are allocated. Adjust resource directives for your cluster if needed.

## Submit the workflow

From the project root:

```bash
bash Kraken_DB/submit_all.sh
# Or use a separate config file:
bash Kraken_DB/submit_all.sh "/absolute/path/to/my settings.sh"
```

Setup and environment checks run first. The submitter then schedules download,
build after successful download, and validation after successful build. Failed
prerequisites cancel dependent jobs using `--kill-on-invalid-dep=yes`. Job IDs
are printed and appended to `submitted_jobs.tsv`. Logs are named
`02_download_database_<job-id>.out`/`.err`, with corresponding names for stages
03 and 04. Monitor jobs with `squeue -u "$USER"` and inspect those logs.

The submitter passes absolute config and bundle paths to every batch job,
including jobs executed from Slurm's spool directory. Keep the scripts and config
available, and do not edit them while jobs are queued or running.

## Individual stages and restarting

Run setup/checks before submitting individual stages. For example:

```bash
BUNDLE="$(cd Kraken_DB && pwd)"
CONFIG="$BUNDLE/config.sh"
bash "$BUNDLE/00_setup_directories.sh" "$CONFIG"
bash "$BUNDLE/01_check_environment.sh" "$CONFIG"

# Obtain the derived log path from the same config.
source "$BUNDLE/environment.sh"
load_kraken_config "$CONFIG"
sbatch --output="$KRAKEN_LOG_DIR/03_build_database_%j.out" \
       --error="$KRAKEN_LOG_DIR/03_build_database_%j.err" \
       "$BUNDLE/03_build_database.sbatch" "$CONFIG" "$BUNDLE"
```

Use `02_download_database.sbatch` or `04_validate_database.sbatch` and matching
log names for the other stages. Wait for the preceding stage to succeed before
submitting the next manually. The second argument, `"$BUNDLE"`, is required for
direct `sbatch` submissions so the spooled script can locate shared code.

After a failure, inspect the logs, resolve the cause, ensure no jobs for this
database are still running, and rerun `submit_all.sh` with the same config.
Downloads resume incomplete files and reuse completed libraries; a completed
build with all three nonempty index files is skipped. Validation runs again.
Do not manually add completion markers to bypass an unsuccessful stage.

## Validation, dates, and output

Validation requires nonempty `hash.k2d`, `opts.k2d`, and `taxo.k2d`, and checks
`kraken2-inspect` for each configured library's major taxon. This checks database
structure and expected contents; it is not a classification accuracy benchmark.

The database directory contains:

- `inspect.txt` and `expected_taxa.inspect.tsv`: inspected database contents.
- `checksums.sha256`: SHA-256 checksums for the three runtime index files.
- `reference_metadata_checksums.sha256`: checksums for reference metadata,
  including the frozen bacterial HTTPS manifest and taxonomy files.
- `provenance/`: UTC stage start/completion timestamps, actual Kraken2 version
  output, loaded modules, and settings for each stage.
- `database_manifest.tsv`: database name, libraries, build parameters, software
  version, stage dates, and latest successful validation time.
- `.database_validation_complete`: published only after successful validation;
  a failed revalidation removes any previous success marker.

Completed-stage timestamps are preserved when resuming. Validation records both
its first successful completion and its most recent successful run. Missing
historical build timestamps are reported as unknown, not inferred from the name.
Checksum paths are relative to the database directory. To check them there:

```bash
cd "/absolute/path/to/kraken_storage/database/cfrna_k2_bacteria_archaea_viral_human_fungi_full"
sha256sum -c checksums.sha256
sha256sum -c reference_metadata_checksums.sha256
```

## Use with the HPC QC tool

After successful validation, set the HPC QC config's `KRAKEN_DB` to the directory
containing `hash.k2d`, `opts.k2d`, `taxo.k2d`, and the validation marker. Its
`KRAKEN_CONFIG` must provide `load_kraken_modules`; this bundle's `environment.sh`
provides that function and can be used as the module loader:

```bash
KRAKEN_DB="/absolute/path/to/kraken_storage/database/cfrna_k2_bacteria_archaea_viral_human_fungi_full"
KRAKEN_CONFIG="/absolute/path/to/cfRNA_QC_tool/Kraken_DB/environment.sh"
```

The existing HPC module and its defaults are unchanged. The database-building
config is separate from the HPC QC config.

## Local verification

```bash
python3 Kraken_DB/test_workflow.py
for script in Kraken_DB/*.sh Kraken_DB/*.sbatch; do bash -n "$script" || break; done
```

Tests use isolated temporary copies, executable mocks, and tiny synthetic
references. They do not download real references, run Kraken2, or submit real
Slurm jobs. Perform a real HPC smoke test to confirm cluster modules, network
access, scheduler behavior, and resource requirements.
