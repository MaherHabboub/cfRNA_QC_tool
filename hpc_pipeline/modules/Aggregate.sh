#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Aggregate.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

SUMMARY_DIR="${OUTDIR}/summary"
mkdir -p "$SUMMARY_DIR"

SUMMARY_TSV="${SUMMARY_DIR}/hpc_qc_summary.tsv"
ZIP_OUT="${SUMMARY_DIR}/hpc_qc_transfer_bundle.zip"

MULTIQC_REPORT="${OUTDIR}/multiqc/hpc_qc_multiqc_report.html"
MULTIQC_DATA_DIR="${OUTDIR}/multiqc/hpc_qc_multiqc_report_data"

echo "Aggregating HPC QC outputs"
echo "SAMPLESHEET: $SAMPLESHEET"
echo "OUTDIR: $OUTDIR"
echo "SUMMARY_TSV: $SUMMARY_TSV"
echo "ZIP_OUT: $ZIP_OUT"

module purge
module load Anaconda3/2024.06-1

python - "$SAMPLESHEET" "$OUTDIR" "$SUMMARY_TSV" "$ZIP_OUT" "$MULTIQC_REPORT" "$MULTIQC_DATA_DIR" <<'PY'
import sys
import os
import glob
import zipfile
from pathlib import Path

import pandas as pd
import numpy as np

samplesheet = Path(sys.argv[1])
outdir = Path(sys.argv[2])
summary_tsv = Path(sys.argv[3])
zip_out = Path(sys.argv[4])
multiqc_report = Path(sys.argv[5])
multiqc_data_dir = Path(sys.argv[6])

summary_tsv.parent.mkdir(parents=True, exist_ok=True)

# -----------------------------
# Helpers
# -----------------------------

def parse_num(x):
    """Convert numbers and percentage strings to float when possible."""
    if x is None:
        return np.nan
    if pd.isna(x):
        return np.nan

    s = str(x).strip()

    if s == "" or s.upper() == "NA":
        return np.nan

    s = s.replace(",", "")

    if s.endswith("%"):
        s = s[:-1].strip()

    try:
        return float(s)
    except Exception:
        return np.nan


def clean_sample_from_path(path):
    """Infer sample from filenames like SAMPLE.metric.tsv."""
    name = Path(path).name

    suffixes = [
        ".mapping_summary.tsv",
        ".duplication_summary.tsv",
        ".fragment_size_summary.tsv",
        ".splice_junction_summary.tsv",
        ".fastqc_parsed_metrics.tsv",
        ".read_distribution.txt",
    ]

    for suf in suffixes:
        if name.endswith(suf):
            return name[:-len(suf)]

    return Path(path).stem


def read_tsv_if_exists(path):
    path = Path(path)

    if not path.exists():
        return pd.DataFrame()

    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame()


def add_file_to_zip(zf, file_path, arcname=None):
    file_path = Path(file_path)

    if not file_path.exists():
        return

    if arcname is None:
        arcname = file_path.name

    zf.write(file_path, arcname)


def add_dir_to_zip(zf, dir_path, arc_prefix):
    dir_path = Path(dir_path)

    if not dir_path.exists():
        return

    for root, _, files in os.walk(dir_path):
        root = Path(root)

        for f in files:
            p = root / f
            rel = p.relative_to(dir_path)
            zf.write(p, str(Path(arc_prefix) / rel))


# -----------------------------
# Load samplesheet
# -----------------------------

samples = pd.read_csv(samplesheet, sep="\t", dtype=str).fillna("NA")

required_cols = ["sample_id", "layout", "condition"]

for col in required_cols:
    if col not in samples.columns:
        raise ValueError(f"Samplesheet missing required column: {col}")

summary = samples[["sample_id", "condition", "layout"]].copy()
summary = summary.rename(columns={"sample_id": "sample"})

# -----------------------------
# Initialize selected metric columns
# -----------------------------

metric_cols = [
    "fastqc_mean_q_last10bp",
    "fastqc_max_adapter_percent",
    "fastqc_gc_peak_percent",

    "uniquely_mapped_pct",
    "mapped_pct",
    "unmapped_pct",
    "unmapped_too_short_pct",
    "multi_mapped_pct",

    "percent_duplication",

    "fraction_20_120",
    "fraction_150_180",
    "fraction_300_1000",
    "peak167_enrichment",

    "read_dist_total_exonic_pct",
    "read_dist_cds_exons_pct",
    "read_dist_utr_exons_pct",
    "read_dist_intronic_pct",
    "read_dist_intergenic_pct",

    "total_junctions",
    "fraction_annotated",
    "fraction_novel",
    "sum_unique_support",
    "splice_junctions_per_million_reads",
]

for col in metric_cols:
    summary[col] = np.nan

summary = summary.set_index("sample", drop=False)

# -----------------------------
# FastQC custom parsed metrics
# -----------------------------

fastqc_files = glob.glob(
    str(outdir / "fastqc" / "raw" / "**" / "*.fastqc_parsed_metrics.tsv"),
    recursive=True
)

fq_rows = []

for f in fastqc_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    fq_rows.append(df)

if fq_rows:
    fq = pd.concat(fq_rows, ignore_index=True)

    for sample, sdf in fq.groupby("sample"):
        if sample not in summary.index:
            continue

        # For paired-end data, average R1/R2 to keep one clean sample-level metric.
        # For single-end data, this is simply the R1 value.
        summary.loc[sample, "fastqc_mean_q_last10bp"] = pd.to_numeric(
            sdf["mean_q_last10bp"],
            errors="coerce"
        ).mean()

        summary.loc[sample, "fastqc_max_adapter_percent"] = pd.to_numeric(
            sdf["max_adapter_percent"],
            errors="coerce"
        ).mean()

        summary.loc[sample, "fastqc_gc_peak_percent"] = pd.to_numeric(
            sdf["gc_peak_percent"],
            errors="coerce"
        ).mean()

# -----------------------------
# Mapping summary
# -----------------------------

mapping_files = glob.glob(
    str(outdir / "mapping" / "**" / "*.mapping_summary.tsv"),
    recursive=True
)

for f in mapping_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    for _, row in df.iterrows():
        sample = str(row.get("sample", clean_sample_from_path(f)))

        if sample not in summary.index:
            continue

        summary.loc[sample, "uniquely_mapped_pct"] = parse_num(row.get("uniquely_mapped_pct"))
        summary.loc[sample, "mapped_pct"] = parse_num(row.get("mapped_pct"))
        summary.loc[sample, "unmapped_pct"] = parse_num(row.get("unmapped_pct"))
        summary.loc[sample, "unmapped_too_short_pct"] = parse_num(row.get("unmapped_too_short_pct"))
        summary.loc[sample, "multi_mapped_pct"] = parse_num(row.get("multi_mapped_pct"))

# Keep input_reads internally for splice junction normalization
input_reads_by_sample = {}

for f in mapping_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    for _, row in df.iterrows():
        sample = str(row.get("sample", clean_sample_from_path(f)))
        input_reads_by_sample[sample] = parse_num(row.get("input_reads"))

# -----------------------------
# Duplication summary
# -----------------------------

dup_files = glob.glob(
    str(outdir / "duplication" / "**" / "*.duplication_summary.tsv"),
    recursive=True
)

for f in dup_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    for _, row in df.iterrows():
        sample = str(row.get("sample", clean_sample_from_path(f)))

        if sample not in summary.index:
            continue

        val = parse_num(row.get("percent_duplication"))

        # Picard PERCENT_DUPLICATION is often a fraction, e.g. 0.61.
        # Convert to percent units for consistency with mapping percentages.
        if not pd.isna(val) and val <= 1:
            val = val * 100.0

        summary.loc[sample, "percent_duplication"] = val

# -----------------------------
# Fragment size summary
# -----------------------------

frag_files = glob.glob(
    str(outdir / "fragment_size" / "**" / "*.fragment_size_summary.tsv"),
    recursive=True
)

for f in frag_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    for _, row in df.iterrows():
        sample = str(row.get("sample", clean_sample_from_path(f)))

        if sample not in summary.index:
            continue

        for col in [
            "fraction_20_120",
            "fraction_150_180",
            "fraction_300_1000",
            "peak167_enrichment",
        ]:
            summary.loc[sample, col] = parse_num(row.get(col))

# -----------------------------
# Read distribution parsing
# -----------------------------

def parse_read_distribution_file(path):
    """
    Parses RSeQC read_distribution.py output.
    Returns selected percentages based on Total Tags when available.
    """
    result = {
        "read_dist_total_exonic_pct": np.nan,
        "read_dist_cds_exons_pct": np.nan,
        "read_dist_utr_exons_pct": np.nan,
        "read_dist_intronic_pct": np.nan,
        "read_dist_intergenic_pct": np.nan,
    }

    if not Path(path).exists():
        return result

    total_tags = np.nan
    rows = {}

    with open(path, "r", errors="replace") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith("Total Tags"):
                parts = line.replace(":", " ").split()
                nums = [parse_num(x) for x in parts]
                nums = [x for x in nums if not pd.isna(x)]

                if nums:
                    total_tags = nums[-1]

                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            group = parts[0]

            # RSeQC table: Group Total_bases Tag_count Tags/Kb
            tag_count = parse_num(parts[2])

            if not pd.isna(tag_count):
                rows[group] = tag_count

    if pd.isna(total_tags) or total_tags == 0:
        # Fallback: use sum of known rows.
        total_tags = sum(rows.values()) if rows else np.nan

    if pd.isna(total_tags) or total_tags == 0:
        return result

    cds = rows.get("CDS_Exons", 0.0)

    utr = 0.0
    for key in ["5UTR_Exons", "3UTR_Exons", "5'UTR_Exons", "3'UTR_Exons"]:
        utr += rows.get(key, 0.0)

    introns = rows.get("Introns", np.nan)
    intergenic = rows.get("Intergenic", np.nan)

    total_exonic = cds + utr

    result["read_dist_cds_exons_pct"] = cds / total_tags * 100.0
    result["read_dist_utr_exons_pct"] = utr / total_tags * 100.0
    result["read_dist_total_exonic_pct"] = total_exonic / total_tags * 100.0

    if not pd.isna(introns):
        result["read_dist_intronic_pct"] = introns / total_tags * 100.0

    if not pd.isna(intergenic):
        result["read_dist_intergenic_pct"] = intergenic / total_tags * 100.0

    return result


rd_files = glob.glob(
    str(outdir / "read_distribution" / "**" / "*.read_distribution.txt"),
    recursive=True
)

for f in rd_files:
    sample = clean_sample_from_path(f)

    if sample not in summary.index:
        continue

    vals = parse_read_distribution_file(f)

    for k, v in vals.items():
        summary.loc[sample, k] = v

# -----------------------------
# Splice junction summary
# -----------------------------

sj_files = glob.glob(
    str(outdir / "splice_junctions" / "**" / "*.splice_junction_summary.tsv"),
    recursive=True
)

for f in sj_files:
    df = read_tsv_if_exists(f)

    if df.empty:
        continue

    for _, row in df.iterrows():
        sample = str(row.get("sample", clean_sample_from_path(f)))

        if sample not in summary.index:
            continue

        total_junctions = parse_num(row.get("total_junctions"))
        sum_unique_support = parse_num(row.get("sum_unique_support"))

        summary.loc[sample, "total_junctions"] = total_junctions
        summary.loc[sample, "fraction_annotated"] = parse_num(row.get("fraction_annotated"))
        summary.loc[sample, "fraction_novel"] = parse_num(row.get("fraction_novel"))
        summary.loc[sample, "sum_unique_support"] = sum_unique_support

        input_reads = input_reads_by_sample.get(sample, np.nan)

        if not pd.isna(total_junctions) and not pd.isna(input_reads) and input_reads > 0:
            summary.loc[sample, "splice_junctions_per_million_reads"] = (
                total_junctions / input_reads * 1_000_000.0
            )

# -----------------------------
# Write main summary
# -----------------------------

summary = summary.reset_index(drop=True)

final_cols = ["sample", "condition", "layout"] + metric_cols
summary = summary[final_cols]

summary.to_csv(summary_tsv, sep="\t", index=False, na_rep="NA")

print(f"Wrote summary table: {summary_tsv}")

# -----------------------------
# Create ZIP bundle
# -----------------------------

if zip_out.exists():
    zip_out.unlink()

with zipfile.ZipFile(zip_out, "w", compression=zipfile.ZIP_DEFLATED) as zf:
    add_file_to_zip(zf, summary_tsv, "hpc_qc_summary.tsv")

    if multiqc_report.exists():
        add_file_to_zip(zf, multiqc_report, "multiqc/hpc_qc_multiqc_report.html")
    else:
        print(f"WARNING: MultiQC report not found: {multiqc_report}")

    if multiqc_data_dir.exists():
        add_dir_to_zip(zf, multiqc_data_dir, "multiqc/hpc_qc_multiqc_report_data")

print(f"Wrote ZIP bundle: {zip_out}")
PY

echo
echo "Aggregation complete."
echo "Main summary:"
echo "$SUMMARY_TSV"
echo
echo "Download bundle:"
echo "$ZIP_OUT"