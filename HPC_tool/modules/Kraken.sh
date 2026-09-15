#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"
TARGET_SAMPLE="${2:-ALL}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash Kraken.sh path/to/config.sh [sample_id]" >&2
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG" >&2
    exit 1
fi

source "$CONFIG"

# -----------------------------
# Required config variables
# -----------------------------
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"
: "${OUTDIR:?ERROR: OUTDIR not set in config}"

# -----------------------------
# Kraken2 settings
# -----------------------------
# These defaults retain the tested HPC installation and database. Override any
# of them in the config when running on another cluster or with another index.
KRAKEN_CONFIG="${KRAKEN_CONFIG:-/scratch/gent/vo/000/gvo00027/tools/kraken/config.sh}"
KRAKEN_DB="${KRAKEN_DB:-/scratch/gent/vo/000/gvo00027/tools/kraken/database/cfrna_k2_bacteria_archaea_viral_human_fungi_20260828_full}"
KRAKEN_RESULTS_DIR="${KRAKEN_RESULTS_DIR:-${OUTDIR}/kraken/results}"
KRAKEN_PYTHON="${KRAKEN_PYTHON:-/usr/bin/python3}"
KRAKEN_CONFIDENCE="${KRAKEN_CONFIDENCE:-0.0}"
KRAKEN_MINIMUM_HIT_GROUPS="${KRAKEN_MINIMUM_HIT_GROUPS:-2}"
KRAKEN_TOP_TAXA="${KRAKEN_TOP_TAXA:-20}"
KRAKEN_THREADS="${SLURM_CPUS_PER_TASK:-${KRAKEN_THREADS:-4}}"

if [[ ! -s "$SAMPLESHEET" ]]; then
    echo "ERROR: Samplesheet is missing or empty: $SAMPLESHEET" >&2
    exit 1
fi

if [[ ! -r "$KRAKEN_CONFIG" ]]; then
    echo "ERROR: Kraken configuration is unreadable: $KRAKEN_CONFIG" >&2
    exit 1
fi

if ! [[ "$KRAKEN_THREADS" =~ ^[1-9][0-9]*$ ]] \
    || ! [[ "$KRAKEN_TOP_TAXA" =~ ^[1-9][0-9]*$ ]] \
    || ! [[ "$KRAKEN_MINIMUM_HIT_GROUPS" =~ ^[0-9]+$ ]]; then
    echo "ERROR: Kraken thread, top-taxa, and minimum-hit-group settings must be valid integers." >&2
    exit 1
fi

source "$KRAKEN_CONFIG"
load_kraken_modules

for required_command in kraken2 "$KRAKEN_PYTHON" gzip; do
    if ! command -v "$required_command" >/dev/null 2>&1; then
        echo "ERROR: Required command is unavailable: $required_command" >&2
        exit 1
    fi
done

for required_path in \
    "$KRAKEN_DB/hash.k2d" \
    "$KRAKEN_DB/opts.k2d" \
    "$KRAKEN_DB/taxo.k2d" \
    "$KRAKEN_DB/.database_validation_complete"; do
    if [[ ! -r "$required_path" ]]; then
        echo "ERROR: Required Kraken database file is missing or unreadable: $required_path" >&2
        exit 1
    fi
done

mkdir -p "$KRAKEN_RESULTS_DIR"

process_sample() {
    local sample="$1"
    local star_log="$2"
    local layout="$3"
    local condition="$4"
    local map_dir unmapped_r1 unmapped_r2 sample_dir work_dir
    local total_input_fragments r1_lines r2_lines unmapped_fragments
    local input_source job_label

    if [[ ! -r "$star_log" ]]; then
        echo "ERROR: STAR log is unreadable for $sample: $star_log" >&2
        return 1
    fi

    map_dir="$(dirname "$star_log")"
    unmapped_r1="${map_dir}/${sample}.Unmapped.out.mate1"
    unmapped_r2="${map_dir}/${sample}.Unmapped.out.mate2"

    if [[ ! -s "$unmapped_r1" ]]; then
        echo "ERROR: STAR unmapped mate1 file is missing or empty for $sample: $unmapped_r1" >&2
        return 1
    fi

    case "$layout" in
        SE)
            input_source="STAR_Unmapped.out.mate1"
            ;;
        PE)
            if [[ ! -s "$unmapped_r2" ]]; then
                echo "ERROR: STAR unmapped mate2 file is missing or empty for $sample: $unmapped_r2" >&2
                return 1
            fi
            input_source="STAR_Unmapped.out.mate1,STAR_Unmapped.out.mate2"
            ;;
        *)
            echo "ERROR: $sample has unsupported layout '$layout'; expected SE or PE." >&2
            return 1
            ;;
    esac

    total_input_fragments="$(awk -F '|' '
        /Number of input reads/ {
            value = $2
            gsub(/[[:space:],]/, "", value)
            print value
            exit
        }
    ' "$star_log")"

    if [[ -z "$total_input_fragments" || ! "$total_input_fragments" =~ ^[0-9]+$ ]]; then
        echo "ERROR: Could not parse Number of input reads from $star_log" >&2
        return 1
    fi

    r1_lines="$(wc -l < "$unmapped_r1")"
    if (( r1_lines % 4 != 0 )); then
        echo "ERROR: mate1 FASTQ line count is not divisible by four for $sample." >&2
        return 1
    fi
    unmapped_fragments=$((r1_lines / 4))

    if [[ "$layout" == "PE" ]]; then
        r2_lines="$(wc -l < "$unmapped_r2")"
        if (( r2_lines % 4 != 0 )); then
            echo "ERROR: mate2 FASTQ line count is not divisible by four for $sample." >&2
            return 1
        fi
        if (( r2_lines / 4 != unmapped_fragments )); then
            echo "ERROR: STAR unmapped mate files have different record counts for $sample." >&2
            return 1
        fi
    fi

    if (( unmapped_fragments == 0 )); then
        echo "ERROR: STAR unmapped FASTQ contains no usable fragments for $sample." >&2
        return 1
    fi

    sample_dir="${KRAKEN_RESULTS_DIR}/${sample}"
    job_label="${SLURM_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-0}"
    work_dir="${sample_dir}/work_${job_label}"
    mkdir -p "$sample_dir" "$work_dir"

    if [[ -f "${sample_dir}/.kraken_analysis_complete" && -s "${sample_dir}/${sample}.microbial_summary.tsv" ]]; then
        echo "Kraken analysis already complete: $sample"
        return 0
    fi

    echo "============================================================"
    echo "Kraken2 microbial screening of STAR-unmapped reads"
    echo "Sample: $sample"
    echo "Condition: $condition"
    echo "Layout: $layout"
    echo "STAR log: $star_log"
    echo "STAR unmapped reads: $input_source"
    echo "STAR total input fragments: $total_input_fragments"
    echo "Kraken input fragments: $unmapped_fragments"
    echo "Database: $KRAKEN_DB"

    {
        printf 'metric\tvalue\n'
        printf 'sample\t%s\n' "$sample"
        printf 'condition\t%s\n' "$condition"
        printf 'layout\t%s\n' "$layout"
        printf 'input_source\t%s\n' "$input_source"
        printf 'star_total_input_fragments\t%s\n' "$total_input_fragments"
        printf 'kraken_input_fragments\t%s\n' "$unmapped_fragments"
        printf 'star_unmapped_mate1\t%s\n' "$unmapped_r1"
        if [[ "$layout" == "PE" ]]; then
            printf 'star_unmapped_mate2\t%s\n' "$unmapped_r2"
        fi
    } > "${sample_dir}/${sample}.extraction_stats.tsv"

    local report raw_output output_gz summary all_taxa top_genera top_species
    report="${sample_dir}/${sample}.kraken.report.tsv"
    raw_output="${work_dir}/${sample}.kraken.output.tsv"
    output_gz="${sample_dir}/${sample}.kraken.output.tsv.gz"
    summary="${sample_dir}/${sample}.microbial_summary.tsv"
    all_taxa="${sample_dir}/${sample}.all_taxa.tsv"
    top_genera="${sample_dir}/${sample}.top_genera.tsv"
    top_species="${sample_dir}/${sample}.top_species.tsv"

    if [[ "$layout" == "PE" ]]; then
        /usr/bin/time -v kraken2 \
            --db "$KRAKEN_DB" \
            --threads "$KRAKEN_THREADS" \
            --confidence "$KRAKEN_CONFIDENCE" \
            --minimum-hit-groups "$KRAKEN_MINIMUM_HIT_GROUPS" \
            --use-names \
            --report-minimizer-data \
            --paired \
            --report "$report" \
            --output "$raw_output" \
            "$unmapped_r1" "$unmapped_r2"
    else
        /usr/bin/time -v kraken2 \
            --db "$KRAKEN_DB" \
            --threads "$KRAKEN_THREADS" \
            --confidence "$KRAKEN_CONFIDENCE" \
            --minimum-hit-groups "$KRAKEN_MINIMUM_HIT_GROUPS" \
            --use-names \
            --report-minimizer-data \
            --report "$report" \
            --output "$raw_output" \
            "$unmapped_r1"
    fi

    [[ -s "$report" ]] || { echo "ERROR: Kraken report is missing or empty for $sample." >&2; return 1; }
    [[ -s "$raw_output" ]] || { echo "ERROR: Kraken per-fragment output is missing or empty for $sample." >&2; return 1; }
    gzip -c "$raw_output" > "$output_gz"
    rm -f "$raw_output"

    "$KRAKEN_PYTHON" - \
        "$report" "$sample" "$layout" "$condition" "$total_input_fragments" "$unmapped_fragments" \
        "$summary" "$all_taxa" "$top_genera" "$top_species" "$KRAKEN_TOP_TAXA" <<'PY'
import csv
import sys

MICROBIAL_ROOTS = {2: "Bacteria", 2157: "Archaea", 10239: "Viruses", 4751: "Fungi"}
HUMAN_TAXID = 9606

def percentage(numerator, denominator):
    return 0.0 if denominator == 0 else 100.0 * numerator / denominator

def write_tsv(path, fieldnames, rows):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

(report_path, sample, layout, condition, total_fragments, unmapped_fragments,
 summary_path, all_taxa_path, top_genera_path, top_species_path, top_n) = sys.argv[1:]
total_fragments, unmapped_fragments, top_n = map(int, (total_fragments, unmapped_fragments, top_n))
rows, stack = [], []
with open(report_path, encoding="utf-8") as handle:
    for line_number, line in enumerate(handle, 1):
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 8:
            raise ValueError(f"Expected an 8-column minimizer report, but line {line_number} contains {len(fields)} columns")
        raw_name = fields[7]
        indentation = len(raw_name) - len(raw_name.lstrip(" "))
        while stack and stack[-1]["indentation"] >= indentation:
            stack.pop()
        taxid = int(fields[6])
        microbial_domain = MICROBIAL_ROOTS[taxid] if taxid in MICROBIAL_ROOTS else (stack[-1]["microbial_domain"] if stack else "")
        row = {"report_percent": float(fields[0]), "clade_fragments": int(fields[1]), "direct_fragments": int(fields[2]), "minimizers": int(fields[3]), "distinct_minimizers": int(fields[4]), "rank": fields[5], "taxid": taxid, "name": raw_name.strip(), "indentation": indentation, "microbial_domain": microbial_domain}
        rows.append(row)
        stack.append(row)

rows_by_taxid = {row["taxid"]: row for row in rows}
unclassified = rows_by_taxid.get(0, {"clade_fragments": 0})["clade_fragments"]
classified = rows_by_taxid.get(1, {"clade_fragments": unmapped_fragments - unclassified})["clade_fragments"]
if classified + unclassified != unmapped_fragments:
    raise ValueError(f"Kraken report does not reconcile with the extracted input: classified={classified}, unclassified={unclassified}, input={unmapped_fragments}")
domain_counts = {domain: rows_by_taxid.get(taxid, {"clade_fragments": 0})["clade_fragments"] for taxid, domain in MICROBIAL_ROOTS.items()}
bacteria, archaea, viruses, fungi = (domain_counts[key] for key in ("Bacteria", "Archaea", "Viruses", "Fungi"))
microbial = bacteria + archaea + viruses + fungi
residual_human = rows_by_taxid.get(HUMAN_TAXID, {"clade_fragments": 0})["clade_fragments"]
other_classified = classified - microbial - residual_human
if other_classified < 0:
    raise ValueError(f"High-level Kraken counts overlap unexpectedly: classified={classified}, microbial={microbial}, human={residual_human}")
microbial_rows = [row for row in rows if row["microbial_domain"] in MICROBIAL_ROOTS.values()]
key = lambda row: (row["clade_fragments"], row["distinct_minimizers"])
genera = sorted([row for row in microbial_rows if row["rank"] == "G"], key=key, reverse=True)[:top_n]
species = sorted([row for row in microbial_rows if row["rank"] == "S"], key=key, reverse=True)[:top_n]
top_genus, top_species = (genera[0] if genera else None), (species[0] if species else None)
summary_fields = ["sample", "condition", "layout", "total_input_fragments", "star_unmapped_fragments", "both_unmapped_fragments", "unmapped_percent_total", "kraken_classified_fragments", "classified_percent_unmapped", "kraken_unclassified_fragments", "unclassified_percent_unmapped", "residual_human_fragments", "residual_human_percent_total", "residual_human_percent_unmapped", "microbial_fragments", "microbial_percent_total", "microbial_percent_unmapped", "bacterial_fragments", "bacterial_percent_total", "bacterial_percent_microbial", "archaeal_fragments", "archaeal_percent_total", "archaeal_percent_microbial", "viral_fragments", "viral_percent_total", "viral_percent_microbial", "fungal_fragments", "fungal_percent_total", "fungal_percent_microbial", "other_classified_fragments", "top_genus", "top_genus_fragments", "top_genus_percent_microbial", "top_genus_distinct_minimizers", "top_species", "top_species_fragments", "top_species_percent_microbial", "top_species_distinct_minimizers"]
summary_row = {"sample": sample, "condition": condition, "layout": layout, "total_input_fragments": total_fragments, "star_unmapped_fragments": unmapped_fragments, "both_unmapped_fragments": unmapped_fragments, "unmapped_percent_total": f"{percentage(unmapped_fragments, total_fragments):.6f}", "kraken_classified_fragments": classified, "classified_percent_unmapped": f"{percentage(classified, unmapped_fragments):.6f}", "kraken_unclassified_fragments": unclassified, "unclassified_percent_unmapped": f"{percentage(unclassified, unmapped_fragments):.6f}", "residual_human_fragments": residual_human, "residual_human_percent_total": f"{percentage(residual_human, total_fragments):.6f}", "residual_human_percent_unmapped": f"{percentage(residual_human, unmapped_fragments):.6f}", "microbial_fragments": microbial, "microbial_percent_total": f"{percentage(microbial, total_fragments):.6f}", "microbial_percent_unmapped": f"{percentage(microbial, unmapped_fragments):.6f}", "bacterial_fragments": bacteria, "bacterial_percent_total": f"{percentage(bacteria, total_fragments):.6f}", "bacterial_percent_microbial": f"{percentage(bacteria, microbial):.6f}", "archaeal_fragments": archaea, "archaeal_percent_total": f"{percentage(archaea, total_fragments):.6f}", "archaeal_percent_microbial": f"{percentage(archaea, microbial):.6f}", "viral_fragments": viruses, "viral_percent_total": f"{percentage(viruses, total_fragments):.6f}", "viral_percent_microbial": f"{percentage(viruses, microbial):.6f}", "fungal_fragments": fungi, "fungal_percent_total": f"{percentage(fungi, total_fragments):.6f}", "fungal_percent_microbial": f"{percentage(fungi, microbial):.6f}", "other_classified_fragments": other_classified, "top_genus": top_genus["name"] if top_genus else "NA", "top_genus_fragments": top_genus["clade_fragments"] if top_genus else 0, "top_genus_percent_microbial": f"{percentage(top_genus['clade_fragments'], microbial):.6f}" if top_genus else "0.000000", "top_genus_distinct_minimizers": top_genus["distinct_minimizers"] if top_genus else 0, "top_species": top_species["name"] if top_species else "NA", "top_species_fragments": top_species["clade_fragments"] if top_species else 0, "top_species_percent_microbial": f"{percentage(top_species['clade_fragments'], microbial):.6f}" if top_species else "0.000000", "top_species_distinct_minimizers": top_species["distinct_minimizers"] if top_species else 0}
write_tsv(summary_path, summary_fields, [summary_row])
all_taxa_fields = ["sample", "report_percent_unmapped", "clade_fragments", "direct_fragments", "minimizers", "distinct_minimizers", "rank", "taxid", "name", "microbial_domain"]
write_tsv(all_taxa_path, all_taxa_fields, [{"sample": sample, "report_percent_unmapped": f"{row['report_percent']:.6f}", **row} for row in rows])
top_taxa_fields = ["sample", "rank", "microbial_domain", "taxid", "name", "clade_fragments", "direct_fragments", "percent_microbial", "percent_total_input", "minimizers", "distinct_minimizers"]
def format_top_taxa(taxa):
    return [{"sample": sample, "rank": row["rank"], "microbial_domain": row["microbial_domain"], "taxid": row["taxid"], "name": row["name"], "clade_fragments": row["clade_fragments"], "direct_fragments": row["direct_fragments"], "percent_microbial": f"{percentage(row['clade_fragments'], microbial):.6f}", "percent_total_input": f"{percentage(row['clade_fragments'], total_fragments):.6f}", "minimizers": row["minimizers"], "distinct_minimizers": row["distinct_minimizers"]} for row in taxa]
write_tsv(top_genera_path, top_taxa_fields, format_top_taxa(genera))
write_tsv(top_species_path, top_taxa_fields, format_top_taxa(species))
PY

    {
        printf 'field\tvalue\n'
        printf 'sample\t%s\n' "$sample"
        printf 'condition\t%s\n' "$condition"
        printf 'layout\t%s\n' "$layout"
        printf 'input_source\t%s\n' "$input_source"
        printf 'kraken_input_fragments\t%s\n' "$unmapped_fragments"
        printf 'completed_at\t%s\n' "$(date --iso-8601=seconds)"
        printf 'hostname\t%s\n' "$(hostname)"
        printf 'slurm_job_id\t%s\n' "${SLURM_JOB_ID:-manual}"
        printf 'kraken_version\t%s\n' "$(kraken2 --version | head -n 1)"
        printf 'database\t%s\n' "$KRAKEN_DB"
        printf 'threads\t%s\n' "$KRAKEN_THREADS"
        printf 'confidence\t%s\n' "$KRAKEN_CONFIDENCE"
        printf 'minimum_hit_groups\t%s\n' "$KRAKEN_MINIMUM_HIT_GROUPS"
    } > "${sample_dir}/${sample}.run_manifest.tsv"

    rm -rf "$work_dir"
    touch "${sample_dir}/.kraken_analysis_complete"
    echo "Kraken2 analysis completed successfully: $sample"
}

found_target="no"
while IFS=$'\t' read -r SAMPLE FASTQ1 FASTQ2 BAM STARLOG SJTAB LAYOUT CONDITION TRANSCRIPTOME_BAM; do
    [[ -z "${SAMPLE:-}" ]] && continue
    if [[ "$TARGET_SAMPLE" != "ALL" && "$SAMPLE" != "$TARGET_SAMPLE" ]]; then
        continue
    fi
    found_target="yes"
    process_sample "$SAMPLE" "$STARLOG" "$LAYOUT" "$CONDITION"
done < <(tail -n +2 "$SAMPLESHEET")

if [[ "$found_target" != "yes" ]]; then
    echo "ERROR: Sample '$TARGET_SAMPLE' was not found in $SAMPLESHEET" >&2
    exit 1
fi
