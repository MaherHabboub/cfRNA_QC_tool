#!/bin/bash

set -euo pipefail

CONFIG="${1:-}"

if [[ -z "$CONFIG" ]]; then
    echo "Usage: bash QC_make_multiqc_custom_content.sh path/to/config.sh"
    exit 1
fi

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

KRAKEN_TOP_GENERA="${KRAKEN_TOP_GENERA:-15}"
if ! [[ "$KRAKEN_TOP_GENERA" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: KRAKEN_TOP_GENERA must be a positive integer: $KRAKEN_TOP_GENERA" >&2
    exit 1
fi

MODULE_SWITCHES=(
    FASTQC_ENABLED
    MAPPING_ENABLED
    DUPLICATION_ENABLED
    INSERT_SIZE_ENABLED
    GENEBODY_ENABLED
    READ_DISTRIBUTION_ENABLED
    SPLICE_JUNCTION_ENABLED
    STRANDEDNESS_ENABLED
    DROPOFF_ENABLED
    KRAKEN_ENABLED
)

for switch_name in "${MODULE_SWITCHES[@]}"; do
    switch_value="${!switch_name:-yes}"
    case "$switch_value" in
        yes|no) ;;
        *)
            echo "ERROR: $switch_name must be 'yes' or 'no': $switch_value" >&2
            exit 1
            ;;
    esac
    printf -v "$switch_name" '%s' "$switch_value"
done

# -----------------------------
# Required config variables
# -----------------------------
: "${OUTDIR:?ERROR: OUTDIR not set in config}"
: "${SAMPLESHEET:?ERROR: SAMPLESHEET not set in config}"

# -----------------------------
# Paths
# -----------------------------
CUSTOM_MQC_DIR="${OUTDIR}/multiqc/custom_content"

# -----------------------------
# Software environment
# -----------------------------
module purge
module load Anaconda3/2024.06-1

mkdir -p "$CUSTOM_MQC_DIR"
find "$CUSTOM_MQC_DIR" -type f -delete 2>/dev/null || true

echo "Creating MultiQC custom content..."
echo "Custom content folder: $CUSTOM_MQC_DIR"

python - "$OUTDIR" "$CUSTOM_MQC_DIR" "$SAMPLESHEET" \
    "$FASTQC_ENABLED" "$MAPPING_ENABLED" "$INSERT_SIZE_ENABLED" \
    "$SPLICE_JUNCTION_ENABLED" "$DROPOFF_ENABLED" "$KRAKEN_ENABLED" \
    "$KRAKEN_TOP_GENERA" "${HPC_RUN_ID:-manual}" <<'PY'
import sys
import os
import glob
import shutil
import pandas as pd
import math
import csv
import json
import base64
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

outdir = sys.argv[1]
custom_dir = sys.argv[2]
samplesheet_path = sys.argv[3]
fastqc_enabled = sys.argv[4] == "yes"
mapping_enabled = sys.argv[5] == "yes"
insert_size_enabled = sys.argv[6] == "yes"
splice_junction_enabled = sys.argv[7] == "yes"
dropoff_enabled = sys.argv[8] == "yes"
kraken_enabled = sys.argv[9] == "yes"
kraken_top_genera = int(sys.argv[10])
run_id = sys.argv[11]
with open(samplesheet_path, newline="") as handle:
    sample_ids = [row["sample_id"] for row in csv.DictReader(handle, delimiter="\t") if row.get("sample_id")]
sample_order = {sample: i for i, sample in enumerate(sample_ids)}

os.makedirs(custom_dir, exist_ok=True)

def read_tsvs(pattern, current_only=False):
    files = glob.glob(pattern, recursive=True)
    dfs = []

    if current_only:
        files.sort(key=lambda f: sample_order.get(Path(f).parent.name, len(sample_order)))
    for f in files:
        try:
            if current_only:
                parent = Path(f).parent
                marker = parent / ".complete"
                if parent.name not in sample_order or not marker.is_file() or marker.read_text().strip() != run_id:
                    continue
            df = pd.read_csv(f, sep="\t")
            if not df.empty:
                dfs.append(df)
        except Exception as e:
            print(f"WARNING: could not read {f}: {e}")

    if dfs:
        return pd.concat(dfs, ignore_index=True)

    return pd.DataFrame()

def safe_num(x):
    if pd.isna(x):
        return None

    if isinstance(x, str):
        y = x.replace("%", "").strip()

        try:
            return float(y)
        except Exception:
            return x

    try:
        if math.isnan(float(x)):
            return None
    except Exception:
        pass

    return x

def table_column(col):
    """Presentation only: retain original TSVs and numeric precision."""
    labels = {
        "mean_q_last10bp": "Mean quality (last 10 bases)",
        "min_mean_q_anypos": "Lowest mean positional quality",
        "gc_peak_percent": "GC peak (%)", "max_adapter_percent": "Maximum adapter content (%)",
        "mapq_min": "Minimum MAPQ", "read": "Read mate",
        "fraction_spliced": "Spliced alignments (%)",
        "fraction_annotated": "Annotated junctions (%)",
        "fraction_novel": "Novel junctions (%)",
        "unmapped_percent_total": "Unmapped (% of total)",
        "classified_percent_unmapped": "Classified (% of unmapped)",
        "microbial_percent_total": "Microbial (% of total)",
        "microbial_percent_unmapped": "Microbial (% of unmapped)",
        "dropoff_exon_side_mean_norm": "Mean exon coverage (normalized)",
        "dropoff_intron_side_mean_norm": "Mean intron coverage (normalized)",
        "dropoff_near_exon_-25bp_norm": "Exon coverage at -25 bp (normalized)",
        "dropoff_near_intron_25bp_norm": "Intron coverage at +25 bp (normalized)",
        "dropoff_near_intron_to_exon_ratio": "Intron / exon coverage ratio",
    }
    millions = {"input_reads": "Input reads", "uniquely_mapped_reads": "Uniquely mapped reads",
                "total_unique_mapped_reads": "Qualifying mapped alignments",
                "spliced_reads": "Spliced alignments", "sum_unique_support": "Unique junction support",
                "sum_multi_support": "Multimapping junction support"}
    thousands = {"total_junctions", "annotated_junctions", "novel_junctions"}
    title = col.replace("_", " ").capitalize()
    factor, fmt = 1, "{:,.2f}"
    if col in millions:
        title, factor, fmt = millions[col] + " (millions)", 1e-6, "{:,.3f}"
    elif col in thousands:
        title, factor, fmt = title + " (thousands)", 1e-3, "{:,.3f}"
    elif "fraction" in col:
        title, factor = title.replace("Fraction", "Percentage").replace("fraction", "percentage") + " (%)", 100
    elif "percent" in col or col.endswith("_pct"):
        title = title.replace(" percent", "").replace(" pct", "") + " (%)"
    elif col.startswith("dropoff_"):
        fmt = "{:,.4f}"
    elif col == "mapq_min":
        fmt = "{:,.0f}"
    return labels.get(col, title), factor, fmt


def df_to_mqc_yaml(df, out_yaml, section_id, section_name, description):
    if df.empty:
        print(f"Skipping empty section: {section_name}")
        return

    df = df.copy()

    if "sample" not in df.columns:
        print(f"Skipping {section_name}: no sample column")
        return

    data = {}

    for _, row in df.iterrows():
        sample = str(row["sample"])
        # FASTQ mates are separate observations; do not overwrite R1 with R2.
        if "read" in df.columns:
            sample += " / " + str(row["read"])
        vals = {}

        for col in df.columns:
            if col == "sample":
                continue

            value = safe_num(row[col])
            vals[str(col)] = value * table_column(col)[1] if isinstance(value, (int, float)) else value

        data[sample] = vals

    with open(out_yaml, "w") as out:
        out.write(f'id: "{section_id}"\n')
        out.write(f'section_anchor: "{section_id}"\n')
        out.write(f'section_name: "{section_name}"\n')
        out.write(f'description: "{description}"\n')
        out.write('plot_type: "table"\n')
        if section_id == "custom_kraken_microbiome_summary":
            out.write('parent_id: "custom_kraken"\nparent_name: "Kraken2 Microbial Screening"\n')
        out.write("headers:\n")
        for col in df.columns:
            if col == "sample":
                continue
            title, _, fmt = table_column(col)
            out.write(f'  {json.dumps(str(col))}:\n')
            out.write(f'    title: {json.dumps(title)}\n')
            out.write(f'    format: {json.dumps(fmt)}\n')
        out.write("pconfig:\n")
        out.write(f'  id: "{section_id}_table"\n')
        out.write(f'  title: "{section_name}"\n')
        out.write("data:\n")

        for sample, vals in data.items():
            out.write(f'  "{sample}":\n')

            for k, v in vals.items():
                if v is None:
                    out.write(f'    "{k}": null\n')
                elif isinstance(v, (int, float)):
                    out.write(f'    "{k}": {v}\n')
                else:
                    vv = str(v).replace('"', '\\"')
                    out.write(f'    "{k}": "{vv}"\n')

    print(f"Wrote: {out_yaml}")

def write_splice_cohort(df):
    result_dir = Path(outdir) / "splice_junctions"
    result_dir.mkdir(parents=True, exist_ok=True)
    table = result_dir / "splice_read_fractions.tsv"
    summary_file = result_dir / "splice_read_fraction_cohort_summary.tsv"
    png_file = result_dir / "splice_read_fractions.png"
    pdf_file = result_dir / "splice_read_fractions.pdf"
    # Never display a plot or table from an earlier run after current failures.
    for path in (table, summary_file, png_file, pdf_file):
        path.unlink(missing_ok=True)
    if df.empty:
        print("WARNING: No current splice results; omitting cohort summaries.")
        return
    df.to_csv(table, sep="\t", index=False)
    valid = df["total_unique_mapped_reads"] > 0
    if not valid.all():
        excluded = ", ".join(df.loc[~valid, "sample"].astype(str))
        print(f"WARNING: Excluding samples with zero qualifying reads: {excluded}", file=sys.stderr)
        df = df.loc[valid].copy()

    if df.empty:
        print("WARNING: No samples with qualifying reads are available for plotting")
        return

    if df["condition"].isna().any() or df["condition"].astype(str).str.strip().eq("").any():
        print("WARNING: At least one sample has an empty condition; omitting splice plot")
        return

    df["condition"] = df["condition"].astype(str)
    conditions = list(dict.fromkeys(df["condition"]))
    groups = [
        100 * df.loc[df["condition"] == condition, "fraction_spliced"].to_numpy(dtype=float)
        for condition in conditions
    ]
    means = np.array([values.mean() for values in groups])

    summary = (
        df.groupby("condition", sort=False)["fraction_spliced"]
          .agg(n_samples="size", mean_fraction="mean", median_fraction="median",
               standard_deviation="std", minimum_fraction="min", maximum_fraction="max")
          .reset_index()
    )
    summary.to_csv(summary_file, sep="\t", index=False, float_format="%.6f")

    rng = np.random.default_rng(42)
    x_positions = np.arange(1, len(conditions) + 1)
    cmap = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(len(conditions))]

    fig_width = max(6.6, 1.35 * len(conditions) + 2.0)
    fig, ax = plt.subplots(figsize=(fig_width, 5.7))
    ax.set_facecolor("#EBEBEB")

    boxplot = ax.boxplot(
        groups,
        positions=x_positions,
        widths=0.55,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#000000", "linewidth": 1.5},
        whiskerprops={"color": "#4D4D4D", "linewidth": 1.1},
        capprops={"color": "#4D4D4D", "linewidth": 1.1},
        boxprops={"edgecolor": "#4D4D4D", "linewidth": 1.1},
    )

    for box, color in zip(boxplot["boxes"], colors):
        box.set_facecolor(color)
        box.set_alpha(0.55)

    for x, values, color in zip(x_positions, groups, colors):
        jitter = rng.uniform(-0.075, 0.075, size=len(values))
        ax.scatter(
            np.full(len(values), x) + jitter,
            values,
            s=34,
            color=color,
            alpha=0.95,
            edgecolors="#333333",
            linewidths=0.4,
            zorder=3,
        )

    label_heights = [min(107.5, values.max() + 5.5) for values in groups]
    for x, mean, label_y in zip(x_positions, means, label_heights):
        ax.text(x, label_y, f"{mean:.2f}%", ha="center", va="bottom", fontsize=15, fontweight="bold")

    ax.set_xlim(0.4, len(conditions) + 0.6)
    ax.set_ylim(0.0, 112)
    ax.set_xticks(x_positions, conditions)
    ax.set_yticks(np.arange(0, 101, 20))
    ax.set_yticks(np.arange(0, 101, 5), minor=True)
    ax.set_ylabel("Qualifying mapped alignments crossing splice junctions (%)")
    ax.set_xlabel("Condition")
    ax.spines[["top", "right"]].set_visible(False)
    for side in ["left", "bottom"]:
        ax.spines[side].set_visible(True)
        ax.spines[side].set_color("#000000")
        ax.spines[side].set_linewidth(1.1)

    ax.tick_params(axis="both", which="major", color="#000000", width=1.0)
    ax.tick_params(axis="y", which="minor", length=0)
    ax.grid(axis="y", which="major", color="#A8A8A8", linewidth=0.8)
    ax.grid(axis="y", which="minor", color="#CACACA", linewidth=0.45)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(png_file, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf_file, bbox_inches="tight", facecolor="white")
    plt.close(fig)

# -----------------------------
# Custom tables to include
# -----------------------------

mapping = read_tsvs(os.path.join(outdir, "mapping", "**", "*.mapping_summary.tsv"), current_only=True) if mapping_enabled else pd.DataFrame()
df_to_mqc_yaml(
    mapping,
    os.path.join(custom_dir, "custom_mapping_summary_mqc.yaml"),
    "custom_mapping_summary",
    "Custom Mapping Summary",
    "Mapping metrics parsed from STAR Log.final.out by the QC pipeline."
)

splice = read_tsvs(os.path.join(outdir, "splice_junctions", "**", "*.splice_junction_summary.tsv"), current_only=True) if splice_junction_enabled else pd.DataFrame()
df_to_mqc_yaml(
    splice,
    os.path.join(custom_dir, "custom_splice_junction_summary_mqc.yaml"),
    "custom_splice_junction_summary",
    "Splice Junction Summary",
    "STAR SJ.out.tab-derived splice junction summary, including annotated and novel junction fractions."
)

splice_read_fractions = read_tsvs(
    os.path.join(outdir, "splice_junctions", "**", "*.splice_read_fraction.tsv"), current_only=True
) if splice_junction_enabled else pd.DataFrame()
write_splice_cohort(splice_read_fractions)
df_to_mqc_yaml(
    splice_read_fractions,
    os.path.join(custom_dir, "custom_splice_read_fractions_mqc.yaml"),
    "custom_splice_read_fractions",
    "Spliced Read Percentages",
    "Percentage of primary, mapped, non-duplicate, QC-passing alignments with MAPQ at least 30 that cross one or more splice junctions. Counts are in millions."
)

splice_plot = os.path.join(outdir, "splice_junctions", "splice_read_fractions.png")
if splice_junction_enabled and os.path.isfile(splice_plot):
    multiqc_splice_plot = os.path.join(custom_dir, "splice_read_fractions_mqc.png")
    shutil.copy2(splice_plot, multiqc_splice_plot)
    print(f"Wrote: {multiqc_splice_plot}")
else:
    print("No splice read-fraction plot found; skipping MultiQC image.")

fastqc = read_tsvs(os.path.join(outdir, "fastqc", "**", "*.fastqc_parsed_metrics.tsv")) if fastqc_enabled else pd.DataFrame()
df_to_mqc_yaml(
    fastqc,
    os.path.join(custom_dir, "custom_fastqc_parsed_metrics_mqc.yaml"),
    "custom_fastqc_parsed_metrics",
    "Custom FastQC Parsed Metrics",
    "Extra FastQC metrics parsed by the QC pipeline, including last-10-bp quality, GC peak, and maximum adapter percentage."
)

insert_size = read_tsvs(
    os.path.join(outdir, "insert_size_distribution", "**", "*.insert_size_distribution_summary.tsv")
) if insert_size_enabled else pd.DataFrame()
df_to_mqc_yaml(
    insert_size,
    os.path.join(custom_dir, "custom_insert_size_distribution_summary_mqc.yaml"),
    "custom_insert_size_distribution_summary",
    "Insert Size Distribution Summary",
    "Insert-size summary metrics generated by the QC pipeline. Paired-end only."
)

# -----------------------------
# Dropoff compact summary table
# -----------------------------

drop_files = glob.glob(
    os.path.join(outdir, "dropoff", "**", "*.dropoff_profile.tsv"),
    recursive=True
) if dropoff_enabled else []

drop_rows = []

for f in drop_files:
    try:
        sample = os.path.basename(f).replace(".dropoff_profile.tsv", "")
        df = pd.read_csv(f, sep="\t")

        exon_mean = df.loc[df["dist_mid"] < 0, "norm_count"].mean()
        intron_mean = df.loc[df["dist_mid"] > 0, "norm_count"].mean()

        near_exon = df.loc[df["dist_mid"] == -25.0, "norm_count"]
        near_intron = df.loc[df["dist_mid"] == 25.0, "norm_count"]

        near_exon_val = near_exon.iloc[0] if len(near_exon) else None
        near_intron_val = near_intron.iloc[0] if len(near_intron) else None

        ratio = None

        if near_exon_val is not None and near_exon_val != 0 and near_intron_val is not None:
            ratio = near_intron_val / near_exon_val

        drop_rows.append({
            "sample": sample,
            "dropoff_exon_side_mean_norm": exon_mean,
            "dropoff_intron_side_mean_norm": intron_mean,
            "dropoff_near_exon_-25bp_norm": near_exon_val,
            "dropoff_near_intron_25bp_norm": near_intron_val,
            "dropoff_near_intron_to_exon_ratio": ratio,
        })

    except Exception as e:
        print(f"WARNING: could not summarize dropoff file {f}: {e}")

drop_summary = pd.DataFrame(drop_rows)

if not drop_summary.empty:
    drop_summary.to_csv(
        os.path.join(custom_dir, "dropoff_compact_summary.tsv"),
        sep="\t",
        index=False
    )

df_to_mqc_yaml(
    drop_summary,
    os.path.join(custom_dir, "custom_exon_intron_dropoff_summary_mqc.yaml"),
    "custom_exon_intron_dropoff_summary",
    "Exon-Intron Drop-off Summary",
    "Compact summary of exon-intron boundary drop-off profiles. Full TSVs are staged beside the report."
)

# -----------------------------
# Combined dropoff plot
# -----------------------------

if drop_files:
    n = len(drop_files)
    fig_height = max(3, 2.8 * n)

    fig, axes = plt.subplots(
        nrows=n,
        ncols=1,
        figsize=(7, fig_height),
        sharex=True
    )

    if n == 1:
        axes = [axes]

    for ax, f in zip(axes, sorted(drop_files)):
        sample = os.path.basename(f).replace(".dropoff_profile.tsv", "")

        try:
            df = pd.read_csv(f, sep="\t")
            df = df.sort_values("dist_mid")

            ax.plot(
                df["dist_mid"],
                df["norm_count"],
                marker="o"
            )

            ax.axvline(
                0,
                linestyle="--"
            )

            # Force y-axis to start at 0 and end at sample-specific max + 5% padding
            ymax = df["norm_count"].max()

            if pd.isna(ymax) or ymax <= 0:
                ymax = 1.0

            ax.set_ylim(0, ymax * 1.05)

            ax.set_title(sample)
            ax.set_ylabel("Norm. count")

        except Exception as e:
            ax.set_title(f"{sample} - failed to plot")
            ax.set_ylim(0, 1)
            ax.text(
                0.5,
                0.5,
                str(e),
                transform=ax.transAxes,
                ha="center",
                va="center"
            )

    axes[-1].set_xlabel("Distance from exon-intron boundary (bp)")

    fig.suptitle("Exon-Intron Drop-off Profiles", y=0.995)
    fig.tight_layout()

    out_png = os.path.join(custom_dir, "exon_intron_dropoff_profiles_mqc.png")
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"Wrote combined dropoff image: {out_png}")
else:
    print("No dropoff profile TSVs found; skipping combined dropoff plot.")

# -----------------------------
# Combined insert size distribution plot
# -----------------------------

insert_size_plot_files = glob.glob(
    os.path.join(outdir, "insert_size_distribution", "**", "*.insert_size_distribution_hist.png"),
    recursive=True
)

if insert_size_plot_files:
    n = len(insert_size_plot_files)
    fig_height = max(3, 2.8 * n)

    fig, axes = plt.subplots(
        nrows=n,
        ncols=1,
        figsize=(7, fig_height)
    )

    if n == 1:
        axes = [axes]

    for ax, png in zip(axes, sorted(insert_size_plot_files)):
        sample = os.path.basename(png).replace(".insert_size_distribution_hist.png", "")

        try:
            img = plt.imread(png)
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(sample)

        except Exception as e:
            ax.axis("off")
            ax.set_title(f"{sample} - failed to load")
            ax.text(
                0.5,
                0.5,
                str(e),
                transform=ax.transAxes,
                ha="center",
                va="center"
            )

    fig.suptitle("Insert Size Distribution Histograms", y=0.995)
    fig.tight_layout()

    out_png = os.path.join(custom_dir, "insert_size_distribution_histograms_mqc.png")
    fig.savefig(out_png, dpi=200)
    plt.close(fig)

    print(f"Wrote combined insert size distribution image: {out_png}")
else:
    print("No insert size distribution histogram PNGs found; skipping combined insert size distribution plot.")

# -----------------------------
# Kraken2 microbial-screening visualizations
# -----------------------------
# Per-sample classification is performed by Kraken.sh. Cohort summaries and
# figures are generated here because this script is called immediately before
# MultiQC stages its custom content.

if kraken_enabled:
    kraken_results = os.path.join(outdir, "kraken", "results")
    kraken_visualizations = os.path.join(outdir, "kraken", "visualizations")
    os.makedirs(kraken_visualizations, exist_ok=True)

    try:
        sample_metadata = pd.read_csv(samplesheet_path, sep="\t", dtype=str)
        required_metadata = {"sample_id", "layout", "condition"}
        if not required_metadata.issubset(sample_metadata.columns):
            missing = sorted(required_metadata.difference(sample_metadata.columns))
            raise ValueError("Samplesheet is missing Kraken metadata columns: " + ", ".join(missing))
        if sample_metadata["sample_id"].duplicated().any():
            raise ValueError("Samplesheet contains duplicate sample IDs")

        summary_frames = []
        taxa_frames = []
        unavailable = []
        for _, metadata in sample_metadata.iterrows():
            sample = str(metadata["sample_id"])
            sample_dir = os.path.join(kraken_results, sample)
            summary_file = os.path.join(sample_dir, sample + ".microbial_summary.tsv")
            taxa_file = os.path.join(sample_dir, sample + ".all_taxa.tsv")
            if not (os.path.isfile(summary_file) and os.path.getsize(summary_file) > 0
                    and os.path.isfile(taxa_file) and os.path.getsize(taxa_file) > 0):
                unavailable.append(sample)
                continue
            summary = pd.read_csv(summary_file, sep="\t")
            taxa = pd.read_csv(taxa_file, sep="\t")
            if len(summary) != 1:
                raise ValueError("Expected one Kraken summary row for " + sample)
            summary["sample"] = sample
            summary["layout"] = str(metadata["layout"])
            summary["condition"] = str(metadata["condition"])
            taxa["sample"] = sample
            taxa["layout"] = str(metadata["layout"])
            taxa["condition"] = str(metadata["condition"])
            summary_frames.append(summary)
            taxa_frames.append(taxa)

        if unavailable:
            print("WARNING: Kraken outputs unavailable for: " + ", ".join(unavailable))
        if not summary_frames:
            print("No completed Kraken outputs found; skipping Kraken custom content.")
        else:
            summary = pd.concat(summary_frames, ignore_index=True)
            taxa = pd.concat(taxa_frames, ignore_index=True)
            sample_order = summary["sample"].tolist()
            numeric_summary = ["total_input_fragments", "star_unmapped_fragments",
                               "kraken_classified_fragments", "kraken_unclassified_fragments",
                               "residual_human_fragments", "bacterial_fragments",
                               "archaeal_fragments", "viral_fragments", "fungal_fragments",
                               "other_classified_fragments"]
            numeric_taxa = ["clade_fragments", "direct_fragments", "minimizers",
                            "distinct_minimizers", "taxid"]
            for column in numeric_summary:
                summary[column] = pd.to_numeric(summary[column], errors="raise")
            for column in numeric_taxa:
                taxa[column] = pd.to_numeric(taxa[column], errors="raise")

            summary["star_mapped_fragments"] = (
                summary["total_input_fragments"] - summary["star_unmapped_fragments"]
            )
            categories = ["STAR mapped", "Residual human", "Bacteria", "Archaea",
                          "Viruses", "Fungi", "Other classified", "Unclassified"]
            colors = {"STAR mapped": "#4E79A7", "Residual human": "#9C755F",
                      "Bacteria": "#59A14F", "Archaea": "#B07AA1",
                      "Viruses": "#E15759", "Fungi": "#F28E2B",
                      "Other classified": "#A0A0A0", "Unclassified": "#D9D9D9"}
            domain_colors = {key: colors[key] for key in ["Bacteria", "Archaea", "Viruses", "Fungi"]}
            full_counts = pd.DataFrame({
                "STAR mapped": summary["star_mapped_fragments"].to_numpy(),
                "Residual human": summary["residual_human_fragments"].to_numpy(),
                "Bacteria": summary["bacterial_fragments"].to_numpy(),
                "Archaea": summary["archaeal_fragments"].to_numpy(),
                "Viruses": summary["viral_fragments"].to_numpy(),
                "Fungi": summary["fungal_fragments"].to_numpy(),
                "Other classified": summary["other_classified_fragments"].to_numpy(),
                "Unclassified": summary["kraken_unclassified_fragments"].to_numpy(),
            }, index=sample_order)
            full_percentages = full_counts.div(summary.set_index("sample")["total_input_fragments"], axis=0) * 100.0
            unmapped_percentages = full_counts.drop(columns=["STAR mapped"]).div(
                summary.set_index("sample")["star_unmapped_fragments"], axis=0
            ).replace([np.inf, -np.inf], 0).fillna(0) * 100.0

            taxa["fpm_total"] = taxa["clade_fragments"] / taxa["sample"].map(
                summary.set_index("sample")["total_input_fragments"]
            ) * 1_000_000.0
            genera = taxa.loc[(taxa["rank"] == "G") & taxa["microbial_domain"].isin(domain_colors)
                              & (taxa["clade_fragments"] > 0)].copy()
            top_count = kraken_top_genera
            selected = pd.DataFrame()
            genus_matrix = pd.DataFrame()
            if not genera.empty:
                metadata = (genera.groupby(["taxid", "name", "microbial_domain"], as_index=False)
                            .agg(total_fpm=("fpm_total", "sum"))
                            .sort_values(["total_fpm", "name"], ascending=[False, True]).head(top_count))
                metadata["taxon_label"] = metadata["name"]
                duplicate_names = metadata["name"].duplicated(keep=False)
                metadata.loc[duplicate_names, "taxon_label"] = (
                    metadata.loc[duplicate_names, "name"] + " [" +
                    metadata.loc[duplicate_names, "taxid"].astype(str) + "]")
                selected = genera.merge(
                    metadata[["taxid", "name", "microbial_domain", "total_fpm", "taxon_label"]],
                    on=["taxid", "name", "microbial_domain"], how="inner")
                selected = selected.groupby(["sample", "taxid", "name", "taxon_label", "microbial_domain", "total_fpm"], as_index=False).agg(
                    fpm_total=("fpm_total", "sum"), clade_fragments=("clade_fragments", "sum"),
                    distinct_minimizers=("distinct_minimizers", "sum"))
                ordered = metadata.sort_values(["total_fpm", "name"], ascending=[True, True])["taxon_label"].tolist()
                genus_matrix = selected.pivot_table(index="taxon_label", columns="sample", values="fpm_total", aggfunc="sum", fill_value=0.0).reindex(index=ordered, columns=sample_order, fill_value=0.0)

            summary.to_csv(os.path.join(kraken_visualizations, "kraken_summary.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
            taxa.to_csv(os.path.join(kraken_visualizations, "kraken_taxa_long.tsv"), sep="\t", index=False, quoting=csv.QUOTE_MINIMAL)
            genus_matrix.to_csv(os.path.join(kraken_visualizations, "kraken_top_genera_fpm_matrix.tsv"), sep="\t", index=True, index_label="genus")

            def draw_stacked(ax, frame, title):
                bottom = np.zeros(len(frame))
                for category in frame.columns:
                    values = frame[category].to_numpy(dtype=float)
                    ax.bar(np.arange(len(frame)), values, bottom=bottom, label=category,
                           color=colors[category], edgecolor="white", linewidth=0.4)
                    bottom += values
                ax.set_title(title, fontweight="bold")
                ax.set_ylim(0, 100)
                ax.set_ylabel("Percentage (%)")
                ax.set_xticks(np.arange(len(frame)))
                ax.set_xticklabels(sample_order, rotation=25, ha="right")
                ax.grid(axis="y", alpha=0.3)
                ax.set_axisbelow(True)
                ax.legend(frameon=False, fontsize=7, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1.0))

            markers = ["o", "s", "^", "D", "P", "X", "v", "<", ">", "*"]

            def draw_top_genera(ax, title):
                if selected.empty or genus_matrix.empty:
                    ax.text(0.5, 0.5, "No microbial genera were reported", ha="center", va="center")
                    ax.set_title(title, fontweight="bold")
                    ax.set_axis_off()
                    return
                taxa_order = list(genus_matrix.index)
                y_lookup = {taxon: index for index, taxon in enumerate(taxa_order)}
                offsets = [0.0] if len(sample_order) == 1 else np.linspace(-0.25, 0.25, len(sample_order))
                max_distinct = max(float(selected["distinct_minimizers"].max()), 1.0)
                for sample_index, sample in enumerate(sample_order):
                    sample_rows = selected.loc[selected["sample"] == sample].set_index("taxon_label")
                    for taxon in taxa_order:
                        if taxon in sample_rows.index:
                            row = sample_rows.loc[taxon]
                            if isinstance(row, pd.DataFrame):
                                row = row.iloc[0]
                            fpm, distinct, domain = float(row["fpm_total"]), float(row["distinct_minimizers"]), str(row["microbial_domain"])
                        else:
                            fpm, distinct = 0.0, 0.0
                            domain = str(selected.loc[selected["taxon_label"] == taxon, "microbial_domain"].iloc[0])
                        size = 28.0 + 92.0 * (math.log10(distinct + 1.0) / math.log10(max_distinct + 1.0))
                        ax.scatter(math.log10(fpm + 1.0), y_lookup[taxon] + offsets[sample_index],
                                   s=size, marker=markers[sample_index % len(markers)], color=domain_colors[domain],
                                   edgecolor="#222222", linewidth=0.5, alpha=0.85 if fpm > 0 else 0.18)
                ax.set_title(title, fontweight="bold")
                ax.set_xlabel("log10(genus fragments per million total + 1)")
                ax.set_yticks(range(len(taxa_order)))
                ax.set_yticklabels(taxa_order, fontsize=8)
                ax.grid(axis="x", alpha=0.3)
                ax.set_axisbelow(True)
                handles = [Patch(facecolor=color, edgecolor="#222222", label=domain) for domain, color in domain_colors.items()]
                handles += [Line2D([], [], linestyle="None", marker=markers[index % len(markers)], markerfacecolor="#B5B5B5", markeredgecolor="#222222", label=sample) for index, sample in enumerate(sample_order)]
                ax.legend(handles=handles, frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1.0))

            evidence = taxa.loc[taxa["rank"].isin(["G", "S"]) & taxa["microbial_domain"].isin(domain_colors) & (taxa["clade_fragments"] > 0)]

            def draw_evidence(ax, title):
                for domain, color in domain_colors.items():
                    for rank, marker in [("G", "o"), ("S", "^")]:
                        part = evidence.loc[(evidence["microbial_domain"] == domain) & (evidence["rank"] == rank)]
                        ax.scatter(np.log10(part["clade_fragments"] + 1), np.log10(part["distinct_minimizers"] + 1), color=color, marker=marker, s=32, alpha=0.7, edgecolor="#222222", linewidth=0.3)
                ax.set_title(title, fontweight="bold")
                ax.set_xlabel("log10(clade fragments + 1)")
                ax.set_ylabel("log10(distinct minimizers + 1)")
                ax.grid(alpha=0.3)
                ax.legend(handles=[Patch(facecolor=color, label=domain) for domain, color in domain_colors.items()] + [Line2D([], [], marker="o", linestyle="None", color="#555555", label="Genus"), Line2D([], [], marker="^", linestyle="None", color="#555555", label="Species")], frameon=False, fontsize=7)

            def save_kraken_figure(fig, stem):
                output = os.path.join(kraken_visualizations, stem + ".png")
                fig.savefig(output, dpi=300, bbox_inches="tight")
                plt.close(fig)
                # Only the combined overview belongs in the report. Individual
                # panels remain available as standalone source figures.
                if stem == "kraken_visualization_overview":
                    # HTML custom content supports shared parents in MultiQC
                    # 1.28; raw image sections are otherwise separate modules.
                    encoded = base64.b64encode(Path(output).read_bytes()).decode("ascii")
                    with open(os.path.join(custom_dir, stem + "_mqc.html"), "w") as handle:
                        handle.write('<!--\nid: kraken_visualization_overview\n'
                                     'section_anchor: kraken_visualization_overview\n'
                                     'parent_id: custom_kraken\n'
                                     'parent_name: "Kraken2 Microbial Screening"\n'
                                     'section_name: "Kraken2 overview"\n-->\n')
                        handle.write('<img alt="Kraken2 four-panel overview" '
                                     'style="max-width:100%;height:auto" '
                                     f'src="data:image/png;base64,{encoded}">\n')

            fig, axes = plt.subplots(2, 2, figsize=(19, 14), constrained_layout=True)
            draw_stacked(axes[0, 0], full_percentages, "A. Complete-library composition")
            draw_stacked(axes[0, 1], unmapped_percentages, "B. STAR-unmapped composition")
            draw_top_genera(axes[1, 0], "C. Top microbial genera")
            draw_evidence(axes[1, 1], "D. Fragment and minimizer evidence")
            fig.suptitle("Kraken2 screening of STAR-unmapped reads", fontsize=16, fontweight="bold")
            save_kraken_figure(fig, "kraken_visualization_overview")

            fig, ax = plt.subplots(figsize=(max(8.0, 1.35 * len(sample_order) + 4.0), 6.0), constrained_layout=True)
            draw_stacked(ax, full_percentages, "Complete-library composition")
            save_kraken_figure(fig, "01_total_library_composition")
            fig, ax = plt.subplots(figsize=(max(8.0, 1.35 * len(sample_order) + 4.0), 6.0), constrained_layout=True)
            draw_stacked(ax, unmapped_percentages, "Composition of STAR-unmapped reads")
            save_kraken_figure(fig, "02_star_unmapped_composition")
            fig, ax = plt.subplots(figsize=(11.5, max(6.0, 0.42 * max(len(genus_matrix.index), 1) + 2.0)), constrained_layout=True)
            draw_top_genera(ax, "Top microbial genera")
            save_kraken_figure(fig, "03_top_genera_dotplot")
            fig, ax = plt.subplots(figsize=(10.5, 7.0), constrained_layout=True)
            draw_evidence(ax, "Taxonomic support: fragments versus distinct minimizers")
            save_kraken_figure(fig, "04_taxon_evidence")

            with open(os.path.join(kraken_visualizations, "visualization_manifest.tsv"), "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["output", "description"])
                writer.writerows([
                    ("kraken_summary.tsv", "Combined sample-level Kraken QC metrics"),
                    ("kraken_taxa_long.tsv", "Combined long-format taxonomic table"),
                    ("kraken_top_genera_fpm_matrix.tsv", "Top-genus fragments-per-million matrix"),
                    ("01_total_library_composition.png", "Complete-library composition"),
                    ("02_star_unmapped_composition.png", "Composition within STAR-unmapped reads"),
                    ("03_top_genera_dotplot.png", "Top genera normalized per million total input fragments"),
                    ("04_taxon_evidence.png", "Clade-fragment versus distinct-minimizer evidence"),
                    ("kraken_visualization_overview.png", "Four-panel Kraken2 overview"),
                ])

            table_columns = ["sample", "condition", "layout", "unmapped_percent_total", "classified_percent_unmapped", "microbial_percent_total", "microbial_percent_unmapped", "top_genus", "top_species"]
            df_to_mqc_yaml(summary[table_columns], os.path.join(custom_dir, "custom_kraken_microbiome_summary_mqc.yaml"), "custom_kraken_microbiome_summary", "Kraken2 Microbial Screen Summary", "Kraken2 classification of STAR-unmapped reads. The overview plot and source tables are staged with this MultiQC report.")
            print("Wrote Kraken2 cohort visualizations: " + kraken_visualizations)
    except Exception as error:
        print("WARNING: Could not create Kraken2 MultiQC custom content: " + str(error))

print("Done creating MultiQC custom content.")
PY

echo "Custom MultiQC content created in:"
echo "$CUSTOM_MQC_DIR"
