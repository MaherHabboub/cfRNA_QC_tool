# cfRNA QC Tool

This repository contains a modular quality control workflow for cell-free RNA
sequencing (cfRNA-seq) data. It was developed as part of a master's thesis
project and is intended as a research workflow for comprehensive QC assessment,
not as fully packaged production software.

cfRNA-seq data can show problems that are not visible from one metric alone.
This workflow therefore combines sequencing-level, alignment-level,
quantification-level, per-sample, and cohort-level checks. The goal is to help
researchers inspect standard RNA-seq QC metrics together with cfRNA-specific
patterns such as insert size distributions, gene body coverage, exon-intron
drop-off, splice junction support, biotype composition, platelet signals, and
possible DNA contamination patterns.

![cfRNA QC workflow overview](cfRNA_QC_workflow.jpg)

## Repository Structure

```text
cfRNA_QC_tool/
├── HPC_tool/
│   ├── submit_HPC_QC.sh
│   ├── requirements.yml
│   ├── config/
│   │   ├── cohortA_10samples_config.sh
│   │   ├── cohortA_10samples_samplesheet.tsv
│   │   └── fastq_only_samplesheet.example.tsv
│   ├── modules/
│   └── wrappers/                 # generated at submission time
├── Local_tool/
│   ├── requirements.yml
│   ├── resources/
│   └── scripts/
├── example_reports/
└── README.md
```

The workflow is split into two main parts:

- `HPC_tool/`: Bash/Slurm workflow intended for the UGent HPC cluster.
- `Local_tool/`: downstream local R scripts for count-based and cohort-level QC.

## Installation

The repository includes separate Conda environment files for the HPC and local
parts of the workflow.

```bash
conda env create -f HPC_tool/requirements.yml
conda env create -f Local_tool/requirements.yml
```

On the UGent HPC, many command-line tools may already be provided through
environment modules. The HPC scripts currently use `module load` statements for
tools such as STAR, HTSeq, FastQC, MultiQC, Picard, RSeQC, BEDTools, SAMtools,
Kent/UCSC tools, and Anaconda. Software loads and versions are specified in
each shell module's **Software environment** section:

- `Star_Alignment.sh` loads `env/software/doduo` and `STAR/2.7.11b-GCC-13.2.0`.
- `HTSeq_Counts.sh` loads `env/software/doduo`, `HTSeq/2.0.7-foss-2023a`, and
  `SAMtools/1.19.2-GCC-13.2.0`.

To adapt to another cluster, edit those sections to match its available
software modules. The Conda file documents the required software stack;
using Conda instead also requires adapting the shell modules' environment
setup. Installing the Conda environment alone does not replace `module load`.

## Required Inputs

The HPC workflow can align FASTQs with STAR and quantify each sample with
HTSeq. Both are enabled by default and independently optional. Set
`STAR_ENABLED="no"` to use existing alignment files, and `HTSEQ_ENABLED="no"`
to skip counting or use a separate quantifier.

Required HPC inputs:

- FASTQ files for each sample when STAR or FastQC is enabled (R1 for SE; R1/R2 for PE).
- An existing STAR index directory when STAR is enabled; index building is external.
- Existing BAM files when STAR is disabled and BAM-based modules are enabled.
- Existing STAR `Log.final.out` files when STAR is disabled and mapping/Kraken is enabled.
- Existing STAR unmapped mate FASTQs beside the final log when STAR is disabled and Kraken is enabled.
- Optional existing STAR `SJ.out.tab` files for the STAR junction-table summary.
- A reference GTF annotation file.
- A reference exon BED file when RSeQC strandedness inference is enabled.
- A tab-separated samplesheet with one row per sample.
- A Bash config file pointing to all required paths.

Required local inputs:

- A gene count matrix, usually HTSeq-style counts.
- Optional sample metadata.
- Cached gene annotation / biotype table.
- Gene set files for platelet, cell-type, and sex-marker QC.
- Aggregated HPC QC summary output if running HPC metric versus PCA correlation.

## HPC Workflow

The HPC workflow is controlled by:

```text
HPC_tool/submit_HPC_QC.sh
```

This script reads a user-provided config file and samplesheet, generates Slurm
wrapper scripts, submits jobs, records submitted job IDs, and sets dependencies
between workflow steps.

The HPC workflow includes:

- Optional STAR alignment, one Slurm job per sample.
- Optional HTSeq counting, one Slurm job per sample on its full genomic BAM.
- Optional BAM downsampling for duplication and gene-body coverage.
- GTF to BED12 annotation conversion.
- Drop-off bin creation.
- FastQC on raw FASTQ files.
- STAR mapping statistics parsing.
- Picard duplication metrics.
- Insert size distribution analysis.
- RSeQC gene body coverage.
- RSeQC read distribution.
- STAR splice junction summary.
- RSeQC strandedness inference.
- Exon-intron drop-off analysis.
- Kraken2 microbial screening of STAR-unmapped reads.
- MultiQC report generation.
- Final aggregation of QC metrics and complete HTSeq count matrices.

Every enabled analysis module runs as a separate job per sample, including
downsampling, mapping-statistics parsing, and splice-junction calculations.
GTF-to-BED12 and drop-off-bin creation each run once and independently; both
read the supplied GTF. A sample can advance without waiting for another
sample's alignment or downsampling.

Each job waits only for the jobs that produce its required files:

| Job | Required preceding jobs |
|---|---|
| FastQC, STAR, GTF-to-BED12, drop-off-bin creation | None |
| HTSeq, mapping statistics, splice calculations, strandedness, Kraken | That sample's STAR job, if enabled |
| Genomic insert-size analysis | That sample's STAR job, if enabled |
| Transcriptome insert-size analysis | None; its transcriptome BAM is supplied separately |
| Downsampling | That sample's STAR job, if enabled |
| Duplication | That sample's downsampling job when enabled; otherwise its STAR job, if enabled |
| Gene-body coverage | GTF-to-BED12 plus that sample's downsampling job, or its STAR job when downsampling is disabled |
| Read distribution | GTF-to-BED12 plus that sample's STAR job, if enabled |
| Drop-off analysis | Drop-off-bin creation plus that sample's STAR job, if enabled |

With `STAR_ENABLED="no"`, supplied alignment files satisfy alignment
prerequisites without a scheduler wait. Providing BAM paths does not
automatically disable STAR. Jobs with no dependencies are immediately eligible
for scheduling; actual start times depend on available Slurm resources.
Only duplication and gene-body coverage use the BAM selected by downsampling;
HTSeq and other genomic-BAM consumers use the full supplied or generated BAM.

MultiQC waits with `afterany` for **every** submitted preparation and analysis
job, including both annotation jobs, to finish, fail, or be cancelled. Inside
that final job, `Make_multiqc_custom_content.sh` first combines current splice
results and creates cohort tables/plots and custom content; MultiQC then builds
the report. Aggregation runs after MultiQC finishes, combining downsampling
manifests and HTSeq counts and creating the transfer bundle.
Input-producing dependencies use `afterok`; impossible dependencies are
cancelled with `--kill-on-invalid-dep=yes`. A failed sample therefore does not
block other samples or prevent final reporting of available results.

STAR and HTSeq can each be disabled independently, as can the ten QC modules:
FastQC, mapping, duplication, insert size, gene-body coverage, read distribution,
splice junctions, strandedness, drop-off, and Kraken2 microbial screening.
They default to `yes`; set the matching `*_ENABLED` variable to `no` to omit
its jobs. GTF-to-BED12, drop-off bin
creation, MultiQC (including custom content), and aggregation always run.
MultiQC and aggregation honor the same switches, so old files from a disabled
module in a reused `OUTDIR` are not included in the current report or summary.

### Alignment-input compatibility

The workflow distinguishes STAR-specific files from generic alignment files.
`Map.sh` requires STAR `Log.final.out`; the STAR junction-table component of
`Splice_Junction.sh` requires `SJ.out.tab`; `Insert_Size_Distribution_Transcriptome.sh`
requires a STAR transcriptome-coordinate BAM; and `Kraken.sh` requires STAR
unmapped FASTQ output (plus the STAR final log). Those modules cannot be
substituted with outputs from another aligner without adding an aligner-specific
adapter.

HTSeq and the BAM-driven portions of downsampling, duplication, genomic insert-size,
gene-body coverage, read distribution, strandedness, drop-off, and the
BAM-derived spliced-read fraction are aligner-agnostic. They accept valid BAMs
from common aligners such as STAR, HISAT2, and Bowtie2. The workflow validates
the BAM with `samtools quickcheck` and requires normal `@SQ` reference records.
For HTSeq, duplication, gene-body coverage, read distribution, and drop-off, an input
not declared coordinate-sorted is sorted temporarily in the module output
directory; the original BAM is not changed.

The BAM reference sequence names and genome build must still match the GTF,
BED12, exon BED, and drop-off bins. For example, a BAM using `1`/`2` chromosome
names is not compatible with annotations using `chr1`/`chr2` until one side is
made consistent. This is a reference-compatibility requirement, not a
STAR-specific one.

### HPC Module Reference

Every module receives the config file as its first argument. Modules that work
per sample read the required FASTQ, BAM, STAR-log, or splice-junction path from
the corresponding row in `SAMPLESHEET`; their outputs are written below
`OUTDIR`. The submitter supplies a sample ID for per-sample jobs, but each of
those modules can also be run manually without a sample ID to process all rows.
The selectable analysis modules use the switches documented below; the
annotation, reporting, and aggregation modules remain mandatory.

The shell modules follow the same structure: config loading, required inputs,
paths, software environment, and sample processing. STAR/HTSeq can write to
custom locations through `STAR_OUTDIR` and `HTSEQ_OUTDIR`.

| Module | Required file inputs | Main outputs | Function |
|---|---|---|---|
| `Star_Alignment.sh` | FASTQ R1 (plus R2 for PE), `STAR_INDEX`, `GTF` | `STAR_OUTDIR/<sample>/<sample>.` followed by `Aligned.sortedByCoord.out.bam`, `Log.final.out`, `SJ.out.tab`, and unmapped mates | Runs the supplied two-pass STAR alignment settings separately for each sample; accepts gzip or plain FASTQs. Does not generate transcriptome BAMs. |
| `HTSeq_Counts.sh` | Full genomic BAM, `GTF` | `HTSEQ_OUTDIR/<sample>/<sample>_htseq_counts.txt` and `.complete` run record | Counts each sample separately with position order, nonunique=none, and configurable strandedness. Aggregation combines the count files. |
| `Downsample.sh` | Valid BAM with `@SQ` header records | `downsampled_bams/manifests/<run-id>/<sample>.tsv`; downsampled BAMs only for samples above the target | Runs independently per sample using deterministic `samtools view -s` sampling. Samples at or below the target retain their original BAM path. Duplication and gene-body coverage read only their sample's atomic manifest; aggregation creates the cohort manifest later. |
| `GTF_to_BED12.sh` | Config `GTF` | `annotation/<gtf-prefix>.genePred`, `.bed12.bed`, and `annotation/BED12.path.txt` | Converts the reference GTF to validated BED12 annotation for RSeQC. The path-record file tells downstream modules the exact BED12 filename generated. |
| `Make_Dropoff_Bins.sh` | Config `GTF` | `annotation/exon_intron_bins/exon_intron_bins.bed` and the raw transcript-level bin BED | Builds deduplicated 50 bp exon- and intron-side bins around exon–intron boundaries for drop-off QC. |
| `Fastqc.sh` | Samplesheet `fastq_r1`; also `fastq_r2` for `PE` samples | `fastqc/raw/<sample>/` FastQC HTML/ZIP reports and `<sample>.fastqc_parsed_metrics.tsv` | Runs FastQC and extracts last-10-base quality, minimum positional quality, GC peak, and maximum adapter content. |
| `Map.sh` | Samplesheet `star_log` | `mapping/<sample>/<sample>.Log.final.out` and `.mapping_summary.tsv` | Copies the STAR final log and extracts mapping, multimapping, and unmapped-read metrics. |
| `Duplication.sh` | Valid BAM with `@SQ` header records | `duplication/<sample>/<sample>.markdup.metrics.txt` and `.duplication_summary.tsv` | Runs Picard MarkDuplicates and records the alignment-based duplicate fraction. A non-coordinate-sorted input is sorted temporarily. |
| `Insert_Size_Distribution_Genomic.sh` | Samplesheet `bam` for `PE` samples without `transcriptome_bam` | `insert_size_distribution/<sample>/` histogram TSV, summary TSV, and histogram PNG | Uses genomic-coordinate paired-end spans to derive insert sizes and summarize cfRNA-relevant size windows and 167 bp peak enrichment. |
| `Insert_Size_Distribution_Transcriptome.sh` | Samplesheet `transcriptome_bam` for `PE` samples | The same `insert_size_distribution/<sample>/` layout, plus classification and ambiguous-pair-example TSVs | Uses STAR transcript-coordinate alignments. It collapses placements by original read pair, accepts one distinct absolute TLEN, and excludes pairs with conflicting transcript-placement lengths. |
| `Genebody.sh` | Valid BAM; generated `annotation/BED12.path.txt` and referenced BED12 | `gene_body_coverage/<sample>.geneBodyCoverage.txt` and RSeQC companion outputs | Builds a BAM index and runs RSeQC gene-body coverage to assess 5′–3′ coverage bias. A non-coordinate-sorted input is sorted temporarily. |
| `Read_Distribution.sh` | Valid BAM; generated `annotation/BED12.path.txt` and referenced BED12 | `read_distribution/<sample>/<sample>.read_distribution.txt` | Runs RSeQC feature distribution to quantify reads in CDS exons, UTR exons, introns, and intergenic regions. A non-coordinate-sorted input is sorted temporarily. |
| `Splice_Junction.sh` | Valid BAM and `condition`; optional STAR `sj_tab` and `star_log` | Per-sample optional STAR junction summary and BAM-derived spliced-read-fraction TSVs in `splice_junctions/<sample>/` | Calculates the fraction of primary, mapped, non-duplicate, QC-passing BAM alignments with MAPQ ≥30 that cross a splice junction. Cohort summaries and plots are created during custom MultiQC preparation. |
| `Strandedness.sh` | Valid BAM; config `EXON_BED` | `strandedness/<sample>/<sample>_RSeQC_output_all.txt` and `_RSeQC_output.txt` | Runs RSeQC library-orientation inference and writes a compact strandedness result. |
| `Dropoff.sh` | Valid BAM; generated `annotation/exon_intron_bins/exon_intron_bins.bed` | `dropoff/<sample>/` bin-coverage TSV, normalized drop-off profile TSV, and PNG | Counts split-read coverage across exon–intron boundary bins and visualizes normalized exon-to-intron drop-off. A non-coordinate-sorted input is sorted temporarily for `bedtools coverage -sorted`. |
| `Kraken.sh` | Samplesheet `star_log`; STAR `Unmapped.out.mate1` (and `Unmapped.out.mate2` for `PE`) | `kraken/results/<sample>/` Kraken report, compressed per-fragment calls, microbial summary, and taxon TSVs | Classifies STAR-unmapped reads against the configured Kraken2 database. Single-end samples use mate 1; paired-end samples use mate 1 and mate 2 with Kraken2 paired mode after matching-record validation. Cohort plots and the MultiQC section are created by `Make_multiqc_custom_content.sh`. |
| `Make_multiqc_custom_content.sh` | Existing QC TSV/PNG outputs under `OUTDIR` | Splice cohort tables and plots; `multiqc/custom_content/` MultiQC YAML tables, compact drop-off TSV, and combined PNGs | Assembles current-run splice results and converts pipeline-specific metrics and plots into MultiQC custom-content files. It is called automatically by `Multiqc.sh`. |
| `Multiqc.sh` | Existing QC outputs under `OUTDIR`; `Make_multiqc_custom_content.sh` | `multiqc/hpc_qc_multiqc_report.html`, report-data directory, staged input, and config files | Stages standard and custom outputs, then generates the combined MultiQC report. |
| `Aggregate.sh` | Resolved `SAMPLESHEET`; available QC outputs; optional MultiQC report; current-run HTSeq counts when enabled | `summary/hpc_qc_summary.tsv`, `summary/hpc_qc_transfer_bundle.zip`, combined downsampling manifest when enabled, and complete HTSeq matrices | Combines QC metrics and successful current-run downsampling manifests and packages the MultiQC report/data. With HTSeq enabled, combines all samples' counts or includes a missing/invalid-count notice. With HTSeq disabled, produces the QC bundle without counts. |

### Example HPC Config

An example config is provided at:

```text
HPC_tool/config/cohortA_10samples_config.sh
```

The config is a Bash file. A minimal version looks like this:

```bash
#!/bin/bash

# Optional Slurm-cluster setup. Leave empty if the cluster is selected before
# submitting the workflow.
CLUSTER_MODULE="cluster/doduo"
CLUSTER_ENV_MODULE="env/software/doduo"

# Independent optional alignment/counting modules. Existing configs also
# default to yes: explicitly disable STAR when reusing existing alignments.
STAR_ENABLED="yes"
HTSEQ_ENABLED="yes"
STAR_INDEX="/path/to/STAR_index"
HTSEQ_STRANDED="yes"  # yes, no, or reverse; no automatic inference

# CPU, memory, and wall-time allocations are set in submit_HPC_QC.sh.
# Software loads and versions are set inside the corresponding shell modules.

# Enabled by default. Only duplication and gene-body coverage use the selected
# BAM; other BAM-based modules use the full supplied or STAR-generated BAM.
DOWNSAMPLE_ENABLED="yes"
DOWNSAMPLE_TARGET_ALIGNMENTS=1000000
DOWNSAMPLE_SEED=42

# Core analysis modules. All default to yes; set any one to no to skip it.
FASTQC_ENABLED="yes"
MAPPING_ENABLED="yes"
DUPLICATION_ENABLED="yes"
INSERT_SIZE_ENABLED="yes"
GENEBODY_ENABLED="yes"
READ_DISTRIBUTION_ENABLED="yes"
SPLICE_JUNCTION_ENABLED="yes"
STRANDEDNESS_ENABLED="yes"
DROPOFF_ENABLED="yes"
KRAKEN_ENABLED="yes"

# Tested Kraken2 installation and cfRNA database. Override on another cluster.
KRAKEN_CONFIG="/path/to/kraken/config.sh"
KRAKEN_DB="/path/to/kraken2_database"
KRAKEN_TOP_TAXA=20
KRAKEN_TOP_GENERA=15

SAMPLESHEET="/path/to/samplesheet.tsv"

GTF="/path/to/reference.gtf"
EXON_BED="/path/to/exons_sorted_merged.bed"

OUTDIR="/path/to/qc_output"
STAR_OUTDIR="${OUTDIR}/star"
HTSEQ_OUTDIR="${OUTDIR}/htseq"
```

The workflow generates the BED12 annotation from `GTF` during its first step.
Downstream modules read the generated path automatically, so `BED12` must not
be set in the config.

STAR and HTSeq reuse `GTF`; it must match the alignment reference/index.
Settings are organized by purpose:

| Setting | Where to edit |
|---|---|
| Input/output paths, analysis parameters, and enable/disable switches | Config file |
| CPUs, memory, and wall time for Slurm jobs | `HPC_tool/submit_HPC_QC.sh` |
| Software module loads and versions | **Software environment** section inside each `HPC_tool/modules/*.sh` script |

`CLUSTER_MODULE` and `CLUSTER_ENV_MODULE` remain optional config settings for
selecting the cluster on the submission host; they do not select STAR/HTSeq
software versions inside jobs.

STAR defaults to 8 CPUs/60G/36 hours and HTSeq to 1 CPU/60G/36 hours. FastQC
and downsampling thread allocations also live in the submitter. These modules read
`SLURM_CPUS_PER_TASK`, keeping tool threads consistent with the allocated CPUs;
HTSeq uses this for temporary sorting and retains single-process counting.
Standalone calls without Slurm use the modules' default thread counts.
Old resource entries in a config no longer override the submitter's allocations.
Samplesheet validation is embedded in the submission script, and standalone
STAR/HTSeq parsing is embedded in those modules; no separate helper directory
is needed. The submission host and standalone STAR/HTSeq calls need `python3`
for this parsing.

With `DOWNSAMPLE_ENABLED="yes"`, each `Downsample.sh` job writes
`<OUTDIR>/downsampled_bams/manifests/<run-id>/<sample>.tsv` before that sample's
duplication and gene-body coverage run. No sample waits for a cohort manifest.
BAMs with no more than the target number of alignments
are not copied or sampled; the manifest records their original BAM as selected.
For larger BAMs, `samtools view -s` retains an approximately target-sized,
deterministic subset using `DOWNSAMPLE_SEED`. Set `DOWNSAMPLE_ENABLED="no"` to
skip the job and run duplication and gene-body coverage on the full BAMs.

At the end, aggregation merges successful current-run records in samplesheet
order into `<OUTDIR>/downsampled_bams/downsampling_manifest.tsv`, retaining
its existing columns. Missing records are reported and omitted; records from
other runs are not reused. An all-failed run produces a header-only manifest.
Mapping and splice jobs also record completion against the run ID. Their
reporting readers exclude stale/failed sample outputs. Splice cohort plots
retain the existing calculations and omit zero-qualifying-read samples; with
no usable samples, reporting warns and continues without a plot.

All `*_ENABLED` values must be exactly `yes` or `no`; unspecified values default
to `yes` for existing configs. A disabled module is neither submitted nor read
by MultiQC or aggregation where applicable. The aggregate TSV keeps its existing columns and
writes `NA` for metrics from disabled modules.

Kraken2 reads its database and module-loader paths from `KRAKEN_DB` and
`KRAKEN_CONFIG`. The example values point to the tested HPC installation; users
on another system must provide paths to their own Kraken2 installation and
database. `KRAKEN_TOP_TAXA` controls the number of per-sample genus/species rows
retained, while `KRAKEN_TOP_GENERA` controls the cohort visualisation. For `PE`
samples, STAR must have produced both `Unmapped.out.mate1` and
`Unmapped.out.mate2` in the directory containing that sample's `star_log`.

#### Building a Kraken2 database

If you need a database, the standalone [Kraken_DB workflow](Kraken_DB/README.md)
downloads references, builds Kraken2's index, and validates it in one configured
storage location. Its README covers setup, submission, restarting, software
versions, build dates, and connecting the validated database to the HPC tool.

Use absolute paths where possible, especially on HPC systems. Relative sample
file paths resolve from the samplesheet directory; relative config input/output
paths resolve from the directory where the submitter is invoked.

### Example Samplesheet

The samplesheet is tab-separated and contains one row per sample. Sample IDs
must start with a letter or digit and contain only letters, digits, `_`, `.`,
or `-`; IDs must be unique and cannot be `ALL`. Layout is `PE` or `SE`, and
condition is required.

With STAR enabled, only FASTQs and sample metadata are needed. For example
(`HPC_tool/config/fastq_only_samplesheet.example.tsv`):

```text
sample_id	fastq_r1	fastq_r2	layout	condition
sample_01	/path/S1_R1.fastq.gz	/path/S1_R2.fastq.gz	PE	CONTROL
sample_02	/path/S2.fastq.gz	NA	SE	CASE
```

The `fastq_r2` column can be omitted for an entirely single-end cohort. Paired
mates must both be gzip-compressed (`.gz`) or both plain FASTQ. Existing
eight-/nine-column sheets remain supported:

```text
sample_id	fastq_r1	fastq_r2	bam	star_log	sj_tab	layout	condition	transcriptome_bam
sample_01	/path/S1_R1.fastq.gz	/path/S1_R2.fastq.gz	/path/S1.bam	/path/S1.Log.final.out	/path/S1.SJ.out.tab	PE	CONTROL	/path/S1.Aligned.toTranscriptome.out.bam
sample_02	/path/S2_R1.fastq.gz	/path/S2_R2.fastq.gz	/path/S2.bam	/path/S2.Log.final.out	/path/S2.SJ.out.tab	PE	CASE	NA
```

For single-end data, set `layout=SE` and use `NA`, `.`, or an empty field for
the unused mate. The submitter normalizes missing fields to `NA`. With STAR
enabled, supplied genomic BAM/log/junction paths are replaced in the runtime
sheet by generated STAR paths; the original input sheet and files are unchanged.

For a non-STAR alignment workflow, set `STAR_ENABLED="no"` and provide the BAM in
`bam`, with `MAPPING_ENABLED="no"` and `KRAKEN_ENABLED="no"`. The `star_log`
and `sj_tab` fields may then be `NA`; `Splice_Junction.sh` will still calculate
its BAM-derived spliced-read fraction, but will omit the optional STAR
junction-table summary. Do not supply a non-STAR BAM in `transcriptome_bam`:
the transcriptome insert-size method is explicitly a STAR transcriptome-BAM
method, so use the genomic insert-size method instead.

With STAR and FastQC both disabled, FASTQs are not required. HTSeq can remain
enabled on an externally generated BAM, independently of the aligner. Set
`HTSEQ_ENABLED="no"` to use another quantifier; importing external count
matrices into the transfer bundle is not part of this workflow.

`transcriptome_bam` is an optional ninth column used only for Insert Size
Distribution. Existing eight-column samplesheets remain valid and use the
genomic-coordinate method. For a non-empty value other than `NA` or `.`, the
submitter uses the transcriptome method (`24 h`, `40G`, `4 CPUs`) for that
sample; an invalid supplied path fails rather than silently falling back. When
the field is empty, `NA`, or `.`, it uses the genomic method (`2 h`, `20G`,
`2 CPUs`). Insert-size analysis applies to PE samples. The `bam` field is the
genomic BAM used by the other BAM-based modules. When STAR is enabled, that
path is generated automatically; the integrated STAR command does not create
a transcriptome BAM, so supply one separately to use the transcriptome method.

### Run the HPC Workflow

From the repository root:

```bash
bash HPC_tool/submit_HPC_QC.sh HPC_tool/config/cohortA_10samples_config.sh
```

The script generates a fresh Slurm wrapper for every submitted job and writes
them to:

```text
HPC_tool/wrappers/run.<id>/
```

Wrappers contain the absolute tool and runtime-config paths calculated at
submission time. Each submission saves a normalized nine-column samplesheet
and config snapshot under `<OUTDIR>/logs/run.<id>/`. Every job uses these
resolved inputs. Wrappers are generated artifacts, intentionally untracked,
and should not be reused after moving the repository or output directory.
Use a separate output directory for concurrently running workflows.

For a manual rerun after submission, use the printed runtime config so that
generated BAM paths and the HTSeq run ID remain consistent:

```bash
bash HPC_tool/modules/HTSeq_Counts.sh "/path/to/qc_output/logs/run.<id>/config.sh" sample_01
bash HPC_tool/modules/Aggregate.sh "/path/to/qc_output/logs/run.<id>/config.sh"
```

Standalone STAR and HTSeq calls also accept the original config and optional
sample ID; omitting the ID processes samples sequentially and separately.
Downsampling, mapping, and splice-junction modules support the same interface.
Standalone calls without a run ID use `manual` for their completion records
and downsampling manifest directory. After standalone sample processing, run
MultiQC and then aggregation to create cohort outputs; sample jobs do not
create cohort plots or merge shared manifest files.
Standalone HTSeq completion records use the run ID `manual`. For aggregation
of FASTQ-only sheets, use the submitter's resolved runtime config.

It writes submission logs to:

```text
<OUTDIR>/logs/submitted_jobs.tsv
```

## Local Workflow

The local workflow contains R scripts for count-based and cohort-level QC
interpretation. These scripts are mainly in:

```text
Local_tool/scripts/
```

The local analysis includes:

- DESeq2-based count normalization.
- PCA on normalized expression matrices.
- PCA scree and scatter plots.
- Gene biotype composition analysis.
- Platelet contamination scoring.
- Other cell-type contamination scoring, including erythroid and endothelial
  gene sets.
- Sex inference using XIST and Y-chromosome marker genes.
- Correlation analysis between HPC-derived QC metrics and principal components.
- HTML report generation.

Activate the local environment before running the scripts:

```bash
conda activate cfrna-qc-local
```

### Run Individual Local Scripts

Example normalization step:

```bash
Rscript Local_tool/scripts/01_normalize_log.R \
  --counts /path/to/count_matrix.tsv \
  --out /path/to/local_qc_results/norm_log
```

Example PCA step:

```bash
Rscript Local_tool/scripts/03_pca_top_variable_genes.R \
  --expr /path/to/local_qc_results/norm_log/counts_deseq2_log2norm_plus1.tsv \
  --out /path/to/local_qc_results/pca \
  --top_n 1000
```

Example HTML report generation:

```bash
Rscript Local_tool/scripts/08_generate_html_report.R \
  --results /path/to/local_qc_results \
  --template Local_tool/scripts/08_report_template.Rmd \
  --out /path/to/local_qc_results/report \
  --title "Local cfRNA QC Report"
```

### Configure the Local Workflow

The repository includes a master runner:

```text
Local_tool/scripts/run_all_local_qc.R
```

This script runs the existing local QC scripts in order, records pass/fail or
skipped status, and renders the HTML report. It resolves bundled scripts and
resources from its own repository location, so it can be run from any working
directory.

Create a personal configuration by copying:

```text
Local_tool/config/local_qc_config.example.R
```

For example:

```bash
cp Local_tool/config/local_qc_config.example.R \
  /path/to/local_qc_config.R
```

Set `COUNTS` and `OUT_ROOT`. `METADATA` and `HPC_SUMMARY` are optional. The
config exposes `BIOTYPE_ENABLED`, `PCA_SSGSEA_ENABLED`, and
`SEX_INFERENCE_ENABLED`; all default to `TRUE`. Normalization always runs.
`PCA_SSGSEA_ENABLED=FALSE` skips both PCA and ssGSEA. HPC metric correlation
runs only when `HPC_SUMMARY` is non-empty and PCA is enabled. The generated HTML
report omits sections for disabled modules.

### Launch from a Terminal

```bash
Rscript Local_tool/scripts/run_all_local_qc.R \
  --config /path/to/local_qc_config.R
```

### Launch from RStudio

For users who prefer not to use a terminal, open this launcher in RStudio:

```text
Local_tool/launch_local_qc_from_rstudio.R
```

At the top of the file, set:

- `local_tool_dir` to the `Local_tool` directory in your cloned repository.
- `config_path` to your edited personal local-QC configuration file.

Then click **Source** in RStudio. The launcher invokes the same master runner
as the terminal command, so it produces the same run logs, module manifest, and
HTML report. It uses the R installation currently open in RStudio.

## Outputs

The HPC workflow writes module-specific output directories under `OUTDIR`.
Typical outputs include:

- Per-sample STAR genomic BAMs, logs, splice-junction tables, and unmapped mates.
- Per-sample HTSeq counts and complete cohort count matrices when enabled.
- FastQC reports and parsed FastQC metrics.
- A downsampling manifest and, when needed, selected downsampled BAMs.
- Mapping summaries parsed from STAR logs.
- Picard duplication metrics.
- Insert size distribution summaries and plots.
- Gene body coverage outputs.
- Read distribution summaries.
- Splice junction summaries.
- Strandedness inference outputs.
- Exon-intron drop-off profiles and plots.
- A MultiQC HTML report.
- A final aggregated HPC QC summary TSV.
- A transfer bundle.

With HTSeq enabled and valid counts for every sample, the transfer ZIP includes
`counts/htseq_counts_combined_all.tsv` (including HTSeq `__` summary rows),
`counts/htseq_counts_combined_genes_only.tsv`, and `counts/htseq_counts_status.tsv`.
Matrices are also saved directly under `HTSEQ_OUTDIR`. They use `gene_id` as
the first column and preserve samplesheet sample order; the genes-only matrix
is ready to use as the local workflow's `COUNTS` input.

If any sample lacks valid counts completed in the current run, aggregation
still creates the QC summary and bundle. It includes the status file listing
affected samples and omits both matrices; previous combined matrices are
removed from `HTSEQ_OUTDIR`. No partial matrix is created. With
`HTSEQ_ENABLED="no"`, aggregation requires no count files and excludes all
HTSeq outputs, including stale counts, from the new bundle.

The local workflow writes results under the selected `out_root`. Typical outputs
include:

- Normalized expression matrices.
- PCA scores, loadings, variance tables, and plots.
- Biotype composition tables and plots.
- ssGSEA contamination score tables and plots.
- Sex-marker expression summaries and plots.
- HPC metric versus principal component correlation tables and scatter plots.
- A local HTML QC report.
- Run logs and step status tables when using the master runner.

Example generated reports are included in:

```text
example_reports/
```
