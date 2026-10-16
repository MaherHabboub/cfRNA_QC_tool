"""Regression checks for BAM reuse; Slurm and alignment tools are mocked."""
import csv
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


SOURCE = Path(__file__).resolve().parents[1]
COLUMNS = ['sample_id', 'fastq_r1', 'fastq_r2', 'bam', 'star_log', 'sj_tab',
           'layout', 'condition', 'transcriptome_bam']


class StarBamReuse(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.tool = self.root / 'HPC_tool'
        shutil.copytree(SOURCE, self.tool, ignore=shutil.ignore_patterns('logs', 'wrappers', '__pycache__'))
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        self.env = dict(os.environ, PATH=f'{self.bin}:{os.environ["PATH"]}',
                        MOCK_CALLS=str(self.root / 'calls'))
        self.command('sbatch', 'echo "$*" >> "$MOCK_CALLS"\necho "Submitted batch job $(wc -l < "$MOCK_CALLS" | tr -d " ")"')
        self.command('module', ':')
        self.command('STAR', 'echo "STAR $*" >> "$MOCK_CALLS"')
        self.command('samtools', '''case "$1" in
view) printf '@HD\\tSO:coordinate\\n@SQ\\tSN:chr1\\tLN:100\\n' ;;
quickcheck) exit 0 ;;
*) exit 1 ;;
esac''')
        self.command('htseq-count', 'echo "HTSEQ $*" >> "$MOCK_CALLS"\nprintf "gene1\\t5\\n"')
        for name in ('annotation.gtf', 'exons.bed', 'provided.bam', 'provided.log',
                     'provided.sj', 'existing.Unmapped.out.mate1', 'reads.fastq'):
            (self.root / name).touch()
        (self.root / 'index').mkdir()
        self.sheet = self.root / 'samples.tsv'
        self.config = self.root / 'config.sh'

    def command(self, name, body):
        path = self.bin / name
        path.write_text('#!/bin/bash\nset -eu\n' + body + '\n')
        path.chmod(0o755)

    def write_inputs(self, rows, **settings):
        with self.sheet.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t')
            writer.writeheader()
            for row in rows:
                writer.writerow(dict(dict.fromkeys(COLUMNS, 'NA'), layout='SE', condition='CASE', **row))
        config = dict(OUTDIR=self.root / 'output', SAMPLESHEET=self.sheet,
                      GTF=self.root / 'annotation.gtf', EXON_BED=self.root / 'exons.bed',
                      STAR_INDEX=self.root / 'index', STAR_ENABLED='yes', FASTQC_ENABLED='no')
        config.update(settings)
        self.config.write_text(''.join(f'{key}="{value}"\n' for key, value in config.items()))

    def existing(self, **changes):
        row = dict(sample_id='existing', bam='provided.bam', star_log='provided.log', sj_tab='provided.sj')
        row.update(changes)
        return row

    def run_script(self, script, config=None, target=None, success=True):
        args = ['bash', str(self.tool / script), str(config or self.config)]
        if target:
            args.append(target)
        result = subprocess.run(args, cwd=self.root, env=self.env, text=True, capture_output=True)
        if success:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def submit(self):
        self.run_script('submit_HPC_QC.sh')
        run = next((self.root / 'output/logs').glob('run.*'))
        with (run / 'samplesheet.tsv').open() as handle:
            samples = {r['sample_id']: r for r in csv.DictReader(handle, delimiter='\t')}
        with (self.root / 'output/logs/submitted_jobs.tsv').open() as handle:
            jobs = {r['job_name']: r for r in csv.DictReader(handle, delimiter='\t')}
        return run, samples, jobs

    def test_existing_bam_skips_star_without_index_or_fastqs(self):
        self.write_inputs([self.existing()], STAR_INDEX='')
        run, samples, jobs = self.submit()
        self.assertFalse(any(name.startswith('star_') for name in jobs))
        self.assertEqual(samples['existing']['bam'], str(self.root / 'provided.bam'))
        self.assertEqual(samples['existing']['star_log'], str(self.root / 'provided.log'))
        self.assertEqual(samples['existing']['sj_tab'], str(self.root / 'provided.sj'))
        self.assertEqual(jobs['htseq_existing']['dependency'], 'none')
        self.run_script('modules/Star_Alignment.sh', run / 'config.sh', 'existing')
        self.run_script('modules/HTSeq_Counts.sh', run / 'config.sh', 'existing')
        calls = (self.root / 'calls').read_text()
        self.assertNotIn('STAR ', calls)
        self.assertIn(str(self.root / 'provided.bam'), calls)

    def test_mixed_sheet_aligns_only_missing_bam_and_keeps_dependencies(self):
        self.write_inputs([self.existing(), dict(sample_id='new', fastq_r1='reads.fastq')])
        run, samples, jobs = self.submit()
        self.assertIn('star_new', jobs)
        self.assertNotIn('star_existing', jobs)
        star_id = jobs['star_new']['job_id']
        self.assertEqual(jobs['htseq_new']['dependency'], f'--dependency=afterok:{star_id}')
        self.assertEqual(jobs['htseq_existing']['dependency'], 'none')
        self.assertEqual(samples['new']['bam'], str(self.root / 'output/star/new/new.Aligned.sortedByCoord.out.bam'))
        self.run_script('modules/Star_Alignment.sh', run / 'config.sh')
        calls = (self.root / 'calls').read_text().splitlines()
        stars = [line for line in calls if line.startswith('STAR ')]
        self.assertEqual(len(stars), 1)
        self.assertIn('/new/new.', stars[0])

    def test_fastq_only_pe_sheet_still_aligns_and_counts_generated_bam(self):
        self.write_inputs([])
        self.sheet.write_text('sample_id\tfastq_r1\tfastq_r2\tlayout\tcondition\n'
                              'paired\treads.fastq\treads.fastq\tPE\tCASE\n')
        run, samples, jobs = self.submit()
        self.assertIn('star_paired', jobs)
        self.run_script('modules/Star_Alignment.sh', run / 'config.sh', 'paired')
        self.assertIn(f'--readFilesIn {self.root}/reads.fastq {self.root}/reads.fastq',
                      (self.root / 'calls').read_text())
        generated = Path(samples['paired']['bam'])
        generated.parent.mkdir(parents=True, exist_ok=True)
        generated.touch()
        # Original FASTQ-only config still resolves the generated BAM for HTSeq.
        self.run_script('modules/HTSeq_Counts.sh', target='paired')
        self.assertIn(str(generated), (self.root / 'calls').read_text())

    def test_standalone_star_and_htseq_reuse_bam(self):
        self.write_inputs([self.existing()], STAR_INDEX='')
        self.run_script('modules/Star_Alignment.sh')
        self.assertFalse((self.root / 'calls').exists())
        self.run_script('modules/HTSeq_Counts.sh')
        self.assertIn(str(self.root / 'provided.bam'), (self.root / 'calls').read_text())

    def test_invalid_supplied_bam_fails_before_submission(self):
        self.write_inputs([self.existing(bam='missing.bam')])
        result = self.run_script('submit_HPC_QC.sh', success=False)
        self.assertIn('bam for existing not found or unreadable', result.stderr)
        self.assertFalse((self.root / 'calls').exists())
        self.run_script('modules/Star_Alignment.sh', success=False)
        self.assertFalse((self.root / 'calls').exists())

    def test_reused_alignment_still_requires_mapping_inputs(self):
        self.write_inputs([self.existing(star_log='NA')])
        result = self.run_script('submit_HPC_QC.sh', success=False)
        self.assertIn('star_log for existing', result.stderr)
        self.assertFalse((self.root / 'calls').exists())

    def test_fastqc_still_requires_reads_for_reused_bam(self):
        self.write_inputs([self.existing()], FASTQC_ENABLED='yes')
        result = self.run_script('submit_HPC_QC.sh', success=False)
        self.assertIn('fastq_r1 for existing', result.stderr)
        self.assertFalse((self.root / 'calls').exists())

    def test_disabled_star_keeps_existing_bam(self):
        self.write_inputs([self.existing()], STAR_ENABLED='no', STAR_INDEX='')
        _, samples, jobs = self.submit()
        self.assertFalse(any(name.startswith('star_') for name in jobs))
        self.assertEqual(samples['existing']['bam'], str(self.root / 'provided.bam'))


if __name__ == '__main__':
    unittest.main()
