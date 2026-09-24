#!/usr/bin/env python3
"""Local integration tests: real wrappers/modules with mocked HPC executables.

Run with a Python environment containing pandas and numpy to include bundle tests.
No jobs are submitted to a real scheduler and no repository artifacts are written.
"""

import csv
import importlib.util
import io
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

TOOL = Path(__file__).resolve().parents[1]
HAS_PANDAS = all(importlib.util.find_spec(name) for name in ("pandas", "numpy"))
HAS_PLOTS = HAS_PANDAS and importlib.util.find_spec("matplotlib") is not None
QC_SWITCHES = [
    "FASTQC", "MAPPING", "DUPLICATION", "INSERT_SIZE", "GENEBODY",
    "READ_DISTRIBUTION", "SPLICE_JUNCTION", "STRANDEDNESS", "DROPOFF", "KRAKEN",
]
MOCK_TOOL = r'''
import json, os, pathlib, sys
name = pathlib.Path(sys.argv[0]).name
args = sys.argv[1:]
with open(os.environ["MOCK_LOG"], "a") as out:
    out.write(json.dumps([name, args]) + "\n")
if name == "sbatch":
    counter = pathlib.Path(os.environ["MOCK_COUNTER"])
    job = int(counter.read_text()) if counter.exists() else 1000
    counter.write_text(str(job + 1))
    print(f"Submitted batch job {job}")
elif name == "STAR":
    if os.environ.get("FAIL_STAR"):
        sys.exit(7)
    prefix = args[args.index("--outFileNamePrefix") + 1]
    for suffix in ("Aligned.sortedByCoord.out.bam", "Log.final.out", "SJ.out.tab", "Unmapped.out.mate1", "Unmapped.out.mate2"):
        pathlib.Path(prefix + suffix).write_text("STAR fixture\n")
elif name == "samtools":
    if args[0] == "quickcheck" and os.environ.get("BAD_BAM"):
        sys.exit(1)
    if args[0] == "view" and "-H" in args:
        order = "queryname" if os.environ.get("UNSORTED_BAM") else "coordinate"
        print(f"@HD\tVN:1.6\tSO:{order}")
        if not os.environ.get("NO_SQ"):
            print("@SQ\tSN:chr1\tLN:1000")
    elif args[0] == "view" and "-c" in args:
        print(1000000 if ".downsampled.bam" in args[-1] else (2000000 if "S2" in args[-1] else 100))
    elif args[0] == "view" and "-b" in args:
        pathlib.Path(args[args.index("-o") + 1]).write_text("sampled fixture\n")
    elif args[0] == "view" and not os.environ.get("NO_READS"):
        for i in range(4):
            spliced = i < (2 if "S2" in args[-1] else 1)
            cigar = "25M100N25M" if spliced else "50M"
            print(f"read{i}\t0\tchr1\t1\t60\t{cigar}")
    if args[0] == "sort":
        pathlib.Path(args[args.index("-o") + 1]).write_text("sorted fixture\n")
elif name == "htseq-count":
    if os.environ.get("FAIL_HTSEQ"):
        print("partial\t1")
        sys.exit(8)
    if os.environ.get("EMPTY_HTSEQ"):
        sys.exit(0)
    if "S2" in args[-2]:
        print("geneB\t2\ngeneA\t3\n__no_feature\t4")
    else:
        print("geneA\t7\ngeneC\t8\n__no_feature\t9")
elif name == "java":
    for arg in args:
        if arg.startswith("O="):
            pathlib.Path(arg[2:]).write_text("Picard fixture\n")
        if arg.startswith("M="):
            pathlib.Path(arg[2:]).write_text("LIBRARY\tA\tB\tC\tD\tE\tF\tG\tPERCENT_DUPLICATION\nlib\t0\t0\t0\t0\t0\t0\t0\t0.2\n")
elif name == "geneBody_coverage.py":
    pathlib.Path(args[args.index("-o") + 1] + ".geneBodyCoverage.txt").write_text("coverage fixture\n")
elif name == "multiqc":
    if "--help" in args:
        print("--no-clean-up")
    elif "--version" in args:
        print("mock MultiQC")
    else:
        out = pathlib.Path(args[args.index("--outdir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / args[args.index("--filename") + 1]).write_text("<html>Mock MultiQC report</html>")
'''


class IntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cfrna_test_")
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name).resolve()
        self.tool = self.root / "HPC_tool"
        shutil.copytree(TOOL, self.tool, ignore=shutil.ignore_patterns("wrappers", "logs", "__pycache__"))
        self.bin = self.root / "bin"
        self.bin.mkdir()
        for name in ("sbatch", "STAR", "htseq-count", "samtools", "java", "geneBody_coverage.py", "multiqc"):
            script = self.bin / name
            script.write_text(f"#!{sys.executable}\n" + MOCK_TOOL)
            script.chmod(0o755)
        module = self.bin / "module"
        module.write_text("#!/bin/bash\nexit 0\n")
        module.chmod(0o755)
        # An extra symlink outside a venv can bypass its pyvenv.cfg discovery.
        for name in ("python", "python3"):
            shim = self.bin / name
            shim.write_text(f'#!/bin/bash\nexec {shlex.quote(sys.executable)} "$@"\n')
            shim.chmod(0o755)
        self.calls = self.root / "calls.jsonl"
        self.env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}",
                        MOCK_LOG=str(self.calls), MOCK_COUNTER=str(self.root / "counter"),
                        EBROOTPICARD=str(self.root), MPLCONFIGDIR=str(self.root / "mpl"))
        self.out = self.root / "out"
        self.sheet = self.root / "samples.tsv"
        self.config = self.root / "config.sh"
        self.index = self.root / "index"
        self.index.mkdir()
        self.gtf = self.root / "reference.gtf"
        self.gtf.touch()
        (self.root / "exons.bed").touch()
        self.settings = {
            "SAMPLESHEET": self.sheet, "OUTDIR": self.out, "GTF": self.gtf,
            "EXON_BED": self.root / "exons.bed", "STAR_INDEX": self.index,
            "CLUSTER_MODULE": "", "CLUSTER_ENV_MODULE": "",
            "DOWNSAMPLE_ENABLED": "no",
            **{f"{name}_ENABLED": "no" for name in QC_SWITCHES},
        }
        self.write_sheet()
        self.write_config()

    def write_sheet(self, minimal=False, ninth=False):
        columns = ["sample_id", "fastq_r1", "fastq_r2", "layout", "condition"] if minimal else [
            "sample_id", "fastq_r1", "fastq_r2", "bam", "star_log", "sj_tab", "layout", "condition"]
        if ninth:
            columns.append("transcriptome_bam")
        with self.sheet.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore")
            writer.writeheader()
            for sample, layout in (("S2", "PE"), ("S1", "SE")):
                paths = {"fastq_r1": f"{sample} read1.fastq.gz", "fastq_r2": f"{sample} read2.fastq.gz",
                         "bam": f"{sample}.bam", "star_log": f"{sample}.Log.final.out", "sj_tab": f"{sample}.SJ.out.tab"}
                for path in paths.values():
                    (self.root / path).write_text("input fixture\n")
                if layout == "SE":
                    paths["fastq_r2"] = ""  # Empty middle fields must not shift columns.
                for mate in (1, 2):
                    (self.root / f"{sample}.Unmapped.out.mate{mate}").touch()
                writer.writerow(dict(paths, sample_id=sample, layout=layout, condition="CONTROL" if sample == "S2" else "CASE",
                                     transcriptome_bam="NA"))

    def write_config(self, **updates):
        self.settings.update(updates)
        self.config.write_text("".join(f"{key}={shlex.quote(str(value))}\n" for key, value in self.settings.items()))

    def run_bash(self, script, *args, success=True, **env):
        result = subprocess.run(["bash", str(script), *map(str, args)], env=dict(self.env, **env),
                                cwd=self.root, text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def submit(self, success=True):
        result = self.run_bash(self.tool / "submit_HPC_QC.sh", self.config, success=success)
        if not success:
            return result
        self.runtime = max((self.out / "logs").glob("run.*/config.sh"), key=lambda p: p.stat().st_mtime_ns)
        with (self.out / "logs/submitted_jobs.tsv").open() as handle:
            self.jobs = {row["job_name"]: row for row in csv.DictReader(handle, delimiter="\t")}
        return result

    def commands(self, name):
        if not self.calls.exists():
            return []
        return [args for tool, args in map(json.loads, self.calls.read_text().splitlines()) if tool == name]

    def run_job(self, name, **env):
        # Slurm sets this from the wrapper directive; reproduce it locally.
        wrapper = Path(self.jobs[name]["wrapper_script"]).read_text()
        cpus = next(line.split("=", 1)[1] for line in wrapper.splitlines()
                    if line.startswith("#SBATCH --cpus-per-task="))
        env.setdefault("SLURM_CPUS_PER_TASK", cpus)
        return self.run_bash(self.jobs[name]["wrapper_script"], **env)

    def dependencies(self, name):
        value = self.jobs[name]["dependency"]
        return set(value.split(":")[1:]) if value != "none" else set()

    def count_path(self, sample):
        return self.out / "htseq" / sample / f"{sample}_htseq_counts.txt"

    def bundle(self):
        return zipfile.ZipFile(self.out / "summary/hpc_qc_transfer_bundle.zip")

    def test_switch_combinations_and_separate_jobs(self):
        for star in ("yes", "no"):
            for htseq in ("yes", "no"):
                with self.subTest(star=star, htseq=htseq):
                    self.write_config(STAR_ENABLED=star, HTSEQ_ENABLED=htseq)
                    self.submit()
                    for sample in ("S2", "S1"):
                        self.assertEqual(f"star_{sample}" in self.jobs, star == "yes")
                        self.assertEqual(f"htseq_{sample}" in self.jobs, htseq == "yes")
                        if htseq == "yes":
                            job = self.jobs[f"htseq_{sample}"]
                            self.assertEqual(job["extra_args"], sample)
                            expected = {self.jobs[f"star_{sample}"]["job_id"]} if star == "yes" else set()
                            self.assertEqual(self.dependencies(f"htseq_{sample}"), expected)
                            self.assertIn(job["job_id"], self.dependencies("multiqc"))
                    self.assertTrue(self.jobs["multiqc"]["dependency"].startswith("--dependency=afterany:"))
                    self.assertEqual(self.dependencies("aggregate"), {self.jobs["multiqc"]["job_id"]})

    def test_default_enabled_minimal_sheet_and_commands(self):
        self.write_sheet(minimal=True)
        original = self.sheet.read_bytes()
        # Legacy resource entries must not override allocations in the submitter.
        self.write_config(STAR_THREADS=3, STAR_MEM="40G", STAR_TIME="12:00:00",
                          HTSEQ_THREADS=7, HTSEQ_MEM="20G", HTSEQ_TIME="01:00:00",
                          HTSEQ_STRANDED="reverse")
        self.submit()
        self.assertEqual(self.sheet.read_bytes(), original)
        with (self.runtime.parent / "samplesheet.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows[0]), 9)
        self.assertEqual(rows[1]["fastq_r2"], "NA")
        self.assertEqual(rows[0]["bam"], str(self.out / "star/S2/S2.Aligned.sortedByCoord.out.bam"))
        for sample in ("S2", "S1"):
            self.run_job(f"star_{sample}")
            self.run_job(f"htseq_{sample}")
        star = self.commands("STAR")
        self.assertEqual(len(star), 2)
        for args, sample, mate_count in zip(star, ("S2", "S1"), (2, 1)):
            start = args.index("--readFilesIn") + 1
            self.assertEqual(len(args[start:args.index("--readFilesCommand")]), mate_count)
            for flag, value in {"--runThreadN": "8", "--twopassMode": "Basic", "--outMultimapperOrder": "Random",
                                "--outSAMmultNmax": "-1", "--outFilterMultimapNmax": "10",
                                "--outSAMprimaryFlag": "AllBestScore", "--outFilterScoreMinOverLread": "0.66",
                                "--outFilterMatchNminOverLread": "0.66", "--outFilterMatchNmin": "20",
                                "--outReadsUnmapped": "Fastx"}.items():
                self.assertEqual(args[args.index(flag) + 1], value)
            self.assertNotIn("--quantMode", args)
            self.assertIn(f"/{sample}/{sample}.", args[args.index("--outFileNamePrefix") + 1])
        htseq = self.commands("htseq-count")
        self.assertEqual(len(htseq), 2)
        for args, sample in zip(htseq, ("S2", "S1")):
            self.assertEqual(args[:8], ["--format", "bam", "--order", "pos", "--nonunique", "none", "--stranded", "reverse"])
            self.assertEqual(args[-2], str(self.out / f"star/{sample}/{sample}.Aligned.sortedByCoord.out.bam"))
            self.assertEqual(Path(f"{self.count_path(sample)}.complete").read_text().strip(), self.runtime.parent.name)
        wrapper = Path(self.jobs["star_S2"]["wrapper_script"]).read_text()
        self.assertIn("#SBATCH --cpus-per-task=8", wrapper)
        self.assertIn("#SBATCH --mem=60G", wrapper)
        self.assertIn("#SBATCH -t 36:00:00", wrapper)
        wrapper = Path(self.jobs["htseq_S2"]["wrapper_script"]).read_text()
        self.assertIn("#SBATCH --cpus-per-task=1", wrapper)
        self.assertIn("#SBATCH --mem=60G", wrapper)
        self.assertIn("#SBATCH -t 36:00:00", wrapper)

    def test_modules_use_slurm_cpu_allocation(self):
        self.write_config(STAR_THREADS=99, HTSEQ_THREADS=99, FASTQC_THREADS=99, DOWNSAMPLE_THREADS=99,
                          FASTQC_ENABLED="yes", DOWNSAMPLE_ENABLED="yes")
        submitter = self.tool / "submit_HPC_QC.sh"
        submitter.write_text(submitter.read_text().replace("STAR_THREADS=8", "STAR_THREADS=3")
                             .replace("HTSEQ_THREADS=1", "HTSEQ_THREADS=2"))
        self.submit()
        self.run_job("star_S2")
        self.run_job("htseq_S2", UNSORTED_BAM="yes")
        args = self.commands("STAR")[0]
        self.assertEqual(args[args.index("--runThreadN") + 1], "3")
        sort = next(args for args in self.commands("samtools") if args[0] == "sort")
        self.assertEqual(sort[sort.index("-@") + 1], "2")
        self.assertIn("#SBATCH --cpus-per-task=2", Path(self.jobs["fastqc_S2"]["wrapper_script"]).read_text())
        self.assertIn("#SBATCH --cpus-per-task=4", Path(self.jobs["downsample_S2"]["wrapper_script"]).read_text())

    def test_alignment_dependencies_and_transcriptome_selection(self):
        self.write_sheet(ninth=True)
        transcriptome = self.root / "transcriptome.bam"
        transcriptome.touch()
        self.sheet.write_text(self.sheet.read_text().replace("CONTROL\tNA", f"CONTROL\t{transcriptome}"))
        for star_enabled in ("yes", "no"):
            for downsample_enabled in ("yes", "no"):
                with self.subTest(star=star_enabled, downsample=downsample_enabled):
                    self.write_config(STAR_ENABLED=star_enabled, DOWNSAMPLE_ENABLED=downsample_enabled,
                                      **{f"{name}_ENABLED": "yes" for name in QC_SWITCHES})
                    self.submit()
                    gtf = {self.jobs["gtf_to_bed12"]["job_id"]}
                    bins = {self.jobs["make_dropoff_bins"]["job_id"]}
                    self.assertEqual(self.dependencies("gtf_to_bed12"), set())
                    self.assertEqual(self.dependencies("make_dropoff_bins"), set())
                    for sample in ("S2", "S1"):
                        star = {self.jobs[f"star_{sample}"]["job_id"]} if star_enabled == "yes" else set()
                        selected = {self.jobs[f"downsample_{sample}"]["job_id"]} if downsample_enabled == "yes" else star
                        expected = {
                            "fastqc": set(), "htseq": star, "mapping": star, "splice_junction": star,
                            "strandedness": star, "kraken": star, "duplication": selected,
                            "genebody": gtf | selected, "read_distribution": gtf | star,
                            "dropoff": bins | star,
                            "insert_size_distribution": set() if sample == "S2" else star,
                        }
                        if star_enabled == "yes":
                            expected["star"] = set()
                        if downsample_enabled == "yes":
                            expected["downsample"] = star
                        for module, deps in expected.items():
                            name = f"{module}_{sample}"
                            self.assertEqual(self.dependencies(name), deps, name)
                            self.assertEqual(self.jobs[name]["extra_args"], sample)
                    analysis_jobs = {job["job_id"] for name, job in self.jobs.items()
                                     if name not in ("multiqc", "aggregate")}
                    self.assertEqual(self.dependencies("multiqc"), analysis_jobs)
                    self.assertTrue(self.jobs["multiqc"]["dependency"].startswith("--dependency=afterany:"))
                    self.assertEqual(self.dependencies("aggregate"), {self.jobs["multiqc"]["job_id"]})
                    self.assertTrue(self.jobs["aggregate"]["dependency"].startswith("--dependency=afterany:"))
                    self.assertIn("Insert_Size_Distribution_Transcriptome.sh", Path(self.jobs["insert_size_distribution_S2"]["wrapper_script"]).read_text())
                    self.assertIn("Insert_Size_Distribution_Genomic.sh", Path(self.jobs["insert_size_distribution_S1"]["wrapper_script"]).read_text())
        for args in self.commands("sbatch"):
            if any(arg.startswith("--dependency=") for arg in args):
                self.assertIn("--kill-on-invalid-dep=yes", args)

    @unittest.skipUnless(HAS_PANDAS, "pandas/numpy required for real aggregation")
    def test_parallel_downsampling_manifests_and_consumers(self):
        self.write_config(STAR_ENABLED="no", HTSEQ_ENABLED="no", DOWNSAMPLE_ENABLED="yes",
                          DUPLICATION_ENABLED="yes", GENEBODY_ENABLED="yes")
        self.submit()
        # Run both sample jobs concurrently; no shared manifest exists yet.
        processes = [subprocess.Popen(["bash", self.jobs[f"downsample_{s}"]["wrapper_script"]],
                                      env=dict(self.env, SLURM_CPUS_PER_TASK="4"), cwd=self.root,
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                     for s in ("S2", "S1")]
        for process in processes:
            stdout, stderr = process.communicate(timeout=30)
            self.assertEqual(process.returncode, 0, stdout + stderr)
        manifests = self.out / "downsampled_bams/manifests" / self.runtime.parent.name
        records = {}
        for sample in ("S2", "S1"):
            with (manifests / f"{sample}.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["sample_id"], sample)
            records[sample] = rows[0]
        self.assertEqual(records["S2"]["status"], "downsampled")
        self.assertEqual(records["S2"]["retained_alignments"], "1000000")
        self.assertEqual(records["S1"]["status"], "not_downsampled")
        self.assertEqual(records["S1"]["selected_bam"], str(self.root / "S1.bam"))
        cohort = self.out / "downsampled_bams/downsampling_manifest.tsv"
        self.assertFalse(cohort.exists())
        annotation = self.out / "annotation"
        annotation.mkdir()
        bed = annotation / "genes.bed"
        bed.write_text("chr1\t0\t1000\n")
        (annotation / "BED12.path.txt").write_text(str(bed))
        self.run_job("duplication_S2")
        self.run_job("genebody_S2")
        for args in self.commands("java"):
            self.assertIn(f'I={records["S2"]["selected_bam"]}', args)
        self.run_job("aggregate")
        with cohort.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["sample_id"] for row in rows], ["S2", "S1"])
        # A failed retry invalidates only its own manifest; S1 remains usable.
        self.run_job("downsample_S2", success=False, BAD_BAM="yes")
        self.assertFalse((manifests / "S2.tsv").exists())
        self.assertTrue((manifests / "S1.tsv").exists())
        self.run_job("aggregate")
        with cohort.open() as handle:
            self.assertEqual([row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")], ["S1"])
        # A new submission must never merge the previous run's records.
        self.submit()
        self.run_job("aggregate")
        with cohort.open() as handle:
            self.assertEqual(list(csv.DictReader(handle, delimiter="\t")), [])

    def test_sample_selection_and_manual_all_modes(self):
        self.write_config(STAR_ENABLED="no", HTSEQ_ENABLED="no", DOWNSAMPLE_ENABLED="yes",
                          MAPPING_ENABLED="yes", SPLICE_JUNCTION_ENABLED="yes")
        self.submit()
        for job in ("downsample", "mapping", "splice_junction"):
            self.run_job(f"{job}_S2")
        for directory in ("mapping", "splice_junctions"):
            self.assertTrue((self.out / directory / "S2/.complete").exists())
            self.assertFalse((self.out / directory / "S1").exists())
        self.assertFalse((self.out / "splice_junctions/splice_read_fractions.tsv").exists())
        for module in ("Downsample.sh", "Map.sh", "Splice_Junction.sh"):
            self.run_bash(self.tool / "modules" / module, self.runtime, "missing", success=False)
            self.run_bash(self.tool / "modules" / module, self.runtime)
        self.assertTrue((self.out / "mapping/S1/.complete").exists())
        self.assertTrue((self.out / "splice_junctions/S1/.complete").exists())

    @unittest.skipUnless(HAS_PLOTS, "pandas/numpy/matplotlib required for actual reporting")
    def test_splice_cohort_reporting_statistics_and_stale_exclusion(self):
        self.write_config(STAR_ENABLED="no", HTSEQ_ENABLED="no", MAPPING_ENABLED="yes",
                          SPLICE_JUNCTION_ENABLED="yes")
        self.sheet.write_text(self.sheet.read_text().replace("CASE", "CONTROL"))
        self.submit()
        for sample in ("S2", "S1"):
            self.run_job(f"mapping_{sample}")
            self.run_job(f"splice_junction_{sample}")
        self.run_job("multiqc")
        result = self.out / "splice_junctions"
        with (result / "splice_read_fractions.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["sample"] for row in rows], ["S2", "S1"])
        self.assertEqual([float(row["fraction_spliced"]) for row in rows], [0.5, 0.25])
        with (result / "splice_read_fraction_cohort_summary.tsv").open() as handle:
            summary = list(csv.DictReader(handle, delimiter="\t"))[0]
        self.assertEqual(summary["n_samples"], "2")
        self.assertEqual(float(summary["mean_fraction"]), 0.375)
        self.assertEqual(float(summary["median_fraction"]), 0.375)
        self.assertAlmostEqual(float(summary["standard_deviation"]), 0.176777, places=6)
        self.assertTrue((result / "splice_read_fractions.png").is_file())
        self.assertTrue((result / "splice_read_fractions.pdf").is_file())
        self.assertTrue((self.out / "multiqc/custom_content/splice_read_fractions_mqc.png").is_file())
        self.run_job("aggregate")
        with self.bundle() as bundle:
            self.assertIn("multiqc/hpc_qc_multiqc_report.html", bundle.namelist())
        # Old sample outputs still exist but must not enter any current report.
        (result / "S1/.complete").write_text("old-run")
        (self.out / "mapping/S1/.complete").write_text("old-run")
        self.run_job("multiqc")
        self.run_job("aggregate")
        with (result / "splice_read_fractions.tsv").open() as handle:
            self.assertEqual([row["sample"] for row in csv.DictReader(handle, delimiter="\t")], ["S2"])
        with (self.out / "summary/hpc_qc_summary.tsv").open() as handle:
            summary = {row["sample"]: row for row in csv.DictReader(handle, delimiter="\t")}
        self.assertEqual(summary["S1"]["mapped_pct"], "NA")
        self.assertEqual(summary["S1"]["total_junctions"], "NA")
        staged = self.out / "multiqc/multiqc_input"
        self.assertFalse((staged / "S1.Log.final.out").exists())
        self.assertFalse((staged / "custom_tables/S1.splice_read_fraction.tsv").exists())
        # Zero qualifying reads keep the sample table but omit the cohort plot.
        self.run_job("splice_junction_S2", NO_READS="yes")
        self.run_job("multiqc")
        self.assertFalse((result / "splice_read_fractions.png").exists())
        self.assertFalse((result / "splice_read_fraction_cohort_summary.tsv").exists())
        # A failed sample must not abort custom content or MultiQC reporting.
        self.run_job("splice_junction_S2", success=False, BAD_BAM="yes")
        self.run_job("multiqc")
        self.assertFalse((result / "splice_read_fractions.tsv").exists())
        self.assertTrue((self.out / "multiqc/hpc_qc_multiqc_report.html").is_file())

    def test_invalid_inputs_fail_before_submission(self):
        original = dict(self.settings)
        for updates, expected in [
            ({"STAR_ENABLED": "maybe"}, "STAR_ENABLED"), ({"HTSEQ_ENABLED": "maybe"}, "HTSEQ_ENABLED"),
            ({"HTSEQ_STRANDED": "invalid"}, "HTSEQ_STRANDED"),
            ({"STAR_INDEX": self.root / "absent"}, "STAR_INDEX"), ({"GTF": self.root / "absent"}, "GTF"),
        ]:
            with self.subTest(updates=updates):
                self.settings = dict(original)
                self.write_config(**updates)
                self.assertIn(expected, self.submit(success=False).stderr)
                self.assertEqual(self.commands("sbatch"), [])
        self.settings = original
        self.write_config()
        sheet = self.sheet.read_text()
        for bad, expected in [(sheet.replace("S1\t", "S2\t"), "Duplicate"),
                              (sheet.replace("\tPE\t", "\tXX\t"), "layout"),
                              (sheet.replace("S1\t", "../S1\t"), "Invalid sample_id"),
                              (sheet.replace("S2 read2.fastq.gz", "missing.gz"), "fastq_r2")]:
            self.sheet.write_text(bad)
            self.assertIn(expected, self.submit(success=False).stderr)
            self.assertEqual(self.commands("sbatch"), [])

    def test_external_bams_without_fastqs_and_without_star_references(self):
        self.write_config(STAR_ENABLED="no")
        self.index.rmdir()
        self.sheet.write_text("sample_id\tbam\tlayout\tcondition\nS2\tS2.bam\tPE\tCONTROL\nS1\tS1.bam\tSE\tCASE\n")
        self.submit()
        for sample in ("S2", "S1"):
            self.run_job(f"htseq_{sample}")
        self.assertEqual(self.commands("htseq-count")[0][-2], str(self.root / "S2.bam"))

    def test_plain_fastq_and_standalone_selection(self):
        self.write_sheet(minimal=True)
        self.sheet.write_text(self.sheet.read_text().replace("S1 read1.fastq.gz", "S1.fastq"))
        (self.root / "S1.fastq").touch()
        self.run_bash(self.tool / "modules/Star_Alignment.sh", self.config, "S1")
        self.assertEqual(len(self.commands("STAR")), 1)
        self.assertNotIn("--readFilesCommand", self.commands("STAR")[0])
        self.run_bash(self.tool / "modules/Star_Alignment.sh", self.config, "missing", success=False)
        self.assertEqual(len(self.commands("STAR")), 1)

    def test_sorting_and_failed_count_publication(self):
        self.write_config(STAR_ENABLED="no", HPC_RUN_ID="fixture")
        bam = self.root / "S2.bam"
        original = bam.read_bytes()
        module = self.tool / "modules/HTSeq_Counts.sh"
        self.run_bash(module, self.config, "S2", UNSORTED_BAM="yes")
        self.assertEqual(bam.read_bytes(), original)
        self.assertEqual(len([args for args in self.commands("samtools") if args[0] == "sort"]), 1)
        self.assertIn("sorted.bam", self.commands("htseq-count")[0][-2])
        old_counts = self.count_path("S2").read_bytes()
        self.run_bash(module, self.config, "S2", success=False, FAIL_HTSEQ="yes")
        self.assertEqual(self.count_path("S2").read_bytes(), old_counts)
        self.assertFalse(Path(f"{self.count_path('S2')}.complete").exists())
        self.assertEqual(list(self.count_path("S2").parent.glob(".htseq.*")), [])
        for failure in ("EMPTY_HTSEQ", "BAD_BAM", "NO_SQ"):
            self.run_bash(module, self.config, "S2", success=False, **{failure: "yes"})
            self.assertFalse(Path(f"{self.count_path('S2')}.complete").exists())

    def test_star_failure_propagates(self):
        self.submit()
        self.run_job("star_S2", success=False, FAIL_STAR="yes")
        self.assertEqual(self.dependencies("htseq_S2"), {self.jobs["star_S2"]["job_id"]})

    @unittest.skipUnless(HAS_PANDAS, "pandas/numpy required for real aggregation")
    def test_complete_matrix_and_bundle(self):
        self.submit()
        for sample in ("S2", "S1"):
            self.run_job(f"star_{sample}")
            self.run_job(f"htseq_{sample}")
        # An unrelated/stale sample must never add another matrix column.
        stale = self.out / "htseq/OLD"
        stale.mkdir()
        (stale / "OLD_htseq_counts.txt").write_text("wrong\t99\n")
        self.run_job("aggregate")
        with self.bundle() as bundle:
            self.assertIn("hpc_qc_summary.tsv", bundle.namelist())
            all_rows = list(csv.reader(io.StringIO(bundle.read("counts/htseq_counts_combined_all.tsv").decode()), delimiter="\t"))
            self.assertEqual(all_rows, [["gene_id", "S2", "S1"], ["geneB", "2", "0"], ["geneA", "3", "7"],
                                        ["__no_feature", "4", "9"], ["geneC", "0", "8"]])
            genes = bundle.read("counts/htseq_counts_combined_genes_only.tsv").decode()
            self.assertNotIn("__no_feature", genes)
            self.assertIn("geneC\t0\t8", genes)

    @unittest.skipUnless(HAS_PANDAS, "pandas/numpy required for real aggregation")
    def test_missing_malformed_and_stale_counts_keep_qc_bundle(self):
        self.write_config(STAR_ENABLED="no")
        self.submit()
        self.run_job("htseq_S2")
        count = self.count_path("S1")
        count.parent.mkdir(parents=True)
        marker = Path(f"{count}.complete")
        for bad in (None, "", "geneA\tbroken\n", "geneA\t-1\n", "geneA\t1\ngeneA\t2\n", "__no_feature\t1\n", "stale"):
            with self.subTest(counts=bad):
                if bad is None:
                    count.unlink(missing_ok=True)
                else:
                    count.write_text("geneA\t1\n" if bad == "stale" else bad)
                marker.write_text("old_run" if bad == "stale" else self.runtime.parent.name)
                for name in ("all", "genes_only"):
                    (self.out / f"htseq/htseq_counts_combined_{name}.tsv").write_text("stale matrix")
                self.run_job("aggregate")
                with self.bundle() as bundle:
                    self.assertIn("hpc_qc_summary.tsv", bundle.namelist())
                    self.assertFalse(any("combined" in name for name in bundle.namelist()))
                    self.assertIn("S1\tmissing_or_invalid", bundle.read("counts/htseq_counts_status.tsv").decode())
                self.assertEqual(list((self.out / "htseq").glob("*combined*.tsv")), [])

    @unittest.skipUnless(HAS_PANDAS, "pandas/numpy required for real aggregation")
    def test_htseq_disabled_needs_no_counts_and_excludes_old_results(self):
        self.write_config(STAR_ENABLED="no", HTSEQ_ENABLED="no")
        self.submit()
        self.assertFalse((self.out / "htseq").exists())
        self.run_job("aggregate")
        with self.bundle() as bundle:
            self.assertIn("hpc_qc_summary.tsv", bundle.namelist())
            self.assertFalse(any(name.startswith("counts/") for name in bundle.namelist()))
        (self.out / "htseq").mkdir()
        (self.out / "htseq/htseq_counts_combined_all.tsv").write_text("stale matrix")
        (self.out / "summary/htseq_counts_status.tsv").write_text("stale failure")
        self.run_job("aggregate")
        with self.bundle() as bundle:
            self.assertFalse(any(name.startswith("counts/") for name in bundle.namelist()))
        self.assertFalse((self.out / "summary/htseq_counts_status.tsv").exists())


if __name__ == "__main__":
    unittest.main(verbosity=2)
