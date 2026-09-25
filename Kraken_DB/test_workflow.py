#!/usr/bin/env python3
"""Offline workflow tests. Every stage runs in an isolated copy with tiny mock inputs."""
import csv
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest


MOCK = r'''
import gzip, hashlib, json, os, pathlib, subprocess, sys
name = pathlib.Path(sys.argv[0]).name
a = sys.argv[1:]
root = pathlib.Path(os.environ['MOCK_ROOT'])
with (root / 'commands.jsonl').open('a') as f:
    f.write(json.dumps([name, *a]) + '\n')
def option(flag): return a[a.index(flag) + 1]
if name == 'sbatch':
    count = root / 'job_count'
    n = int(count.read_text()) + 1 if count.exists() else 101
    count.write_text(str(n))
    if os.environ.get('FAIL_SUBMIT') == 'yes': sys.exit(1)
    print(str(n) + ';test_cluster')
elif name == 'kraken2': print('Kraken version 2.1.3')
elif name == 'kraken2-build':
    if '--help' in a: print('--skip-maps'); sys.exit()
    if '--version' in a: print('Kraken version 2.1.3'); sys.exit()
    db = pathlib.Path(option('--db'))
    if '--download-taxonomy' in a:
        (db / 'taxonomy').mkdir(exist_ok=True)
        for n in ['names.dmp', 'nodes.dmp']: (db / 'taxonomy' / n).write_text('taxonomy\n')
    elif '--download-library' in a:
        d = db / 'library' / option('--download-library'); d.mkdir(parents=True)
        (d / 'library.fna').write_text('>ref\nACGT\n')
    elif '--build' in a:
        if os.environ.get('FAIL_BUILD') == 'yes': sys.exit(1)
        for n in ['hash.k2d', 'opts.k2d', 'taxo.k2d']: (db / n).write_text(n + '\n')
elif name == 'wget':
    if '-O' in a:
        row = ['na'] * 20
        row[5], row[11], row[19] = '2', 'Complete Genome', 'ftp://example.invalid/ASM1'
        pathlib.Path(option('-O')).write_text('\t'.join(row) + '\n')
    else:
        d = pathlib.Path(next(x.split('=', 1)[1] for x in a if x.startswith('--directory-prefix=')))
        with gzip.open(d / 'ASM1_genomic.fna.gz', 'wt') as f: f.write('>genome\nACGTACGT\n')
elif name == 'scan_fasta_file.pl': print('TAXID\tgenome\t2')
elif name == 'mask_low_complexity.sh': (pathlib.Path(a[0]) / 'library.fna.masked').touch()
elif name == 'kraken2-inspect':
    for taxid in [2, 2157, 10239, 4751, 9606]:
        if str(taxid) != os.environ.get('MISSING_TAXON'): print(f'10\t1\t1\tD\t{taxid}\tTaxon')
elif name == 'sha256sum':
    if os.environ.get('FAIL_CHECKSUM') == 'yes': sys.exit(1)
    for n in a: print(hashlib.sha256(pathlib.Path(n).read_bytes()).hexdigest() + '  ' + n)
elif name == 'readlink': print(pathlib.Path(a[-1]).resolve())
elif name == 'du': print('1000\t' + a[-1])
elif name == 'mock_time': sys.exit(subprocess.call(a[1:]))
elif name == 'xargs':
    # Run the supplied workers against the miniature manifest, without network.
    command = a[a.index('bash'):]
    for line in sys.stdin:
        subprocess.run(command + line.split(), check=True)
elif name == 'module':
    if os.environ.get('FAIL_MODULE') == 'yes': sys.exit(1)
    if a == ['list']: print('env/software/doduo Kraken2/2.1.3-gompi-2023a')
'''


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='kraken bundle tests ')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.bundle = self.root / 'bundle with spaces'
        shutil.copytree(Path(__file__).parent, self.bundle)
        # macOS time lacks GNU -v; substitute only this utility in the test copy.
        for p in self.bundle.glob('*.sbatch'):
            p.write_text(p.read_text().replace('/usr/bin/time -v', 'mock_time -v'))
        self.bin = self.root / 'bin'
        self.bin.mkdir()
        for name in ['module', 'sbatch', 'kraken2', 'kraken2-build', 'kraken2-inspect',
                     'wget', 'rsync', 'dustmasker', 'sha256sum', 'readlink', 'du',
                     'free', 'mock_time', 'scan_fasta_file.pl', 'mask_low_complexity.sh', 'xargs']:
            p = self.bin / name
            p.write_text('#!' + sys.executable + '\n' + MOCK)
            p.chmod(0o755)
        self.config = self.root / 'custom settings.sh'
        self.storage = self.root / 'database storage'
        self.config.write_text((self.bundle / 'config.sh').read_text().replace(
            '"/path/to/kraken_storage"', shlex.quote(str(self.storage))))
        self.db = self.storage / 'database/cfrna_k2_bacteria_archaea_viral_human_fungi_full'
        self.env = dict(os.environ, PATH=str(self.bin) + os.pathsep + os.environ['PATH'],
                        MOCK_ROOT=str(self.root))

    def run_script(self, name, success=True, extra_env=None, path=None):
        p = subprocess.run(['bash', str(path or self.bundle / name), str(self.config), str(self.bundle)],
                           env=dict(self.env, **(extra_env or {})), text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if success:
            self.assertEqual(p.returncode, 0, p.stdout)
        else:
            self.assertNotEqual(p.returncode, 0, p.stdout)
        return p

    def commands(self, name):
        p = self.root / 'commands.jsonl'
        return [r for r in map(json.loads, p.read_text().splitlines()) if r[0] == name] if p.exists() else []

    def download(self):
        self.run_script('00_setup_directories.sh')
        self.run_script('02_download_database.sbatch')

    def build(self):
        self.download()
        self.run_script('03_build_database.sbatch')

    def test_submission_paths_dependencies_and_spooled_script(self):
        self.run_script('submit_all.sh')
        jobs = self.commands('sbatch')
        self.assertEqual(len(jobs), 3)
        for i, args in enumerate(jobs):
            self.assertEqual(args[-2:], [str(self.config), str(self.bundle)])
            self.assertIn('--output=' + str(self.storage / 'logs' / self.db.name) +
                          '/' + ['02_download_database', '03_build_database', '04_validate_database'][i] + '_%j.out', args)
            deps = [x for x in args if x.startswith('--dependency=')]
            self.assertEqual(deps, [] if i == 0 else ['--dependency=afterok:' + str(100 + i)])
            self.assertEqual('--kill-on-invalid-dep=yes' in args, i > 0)
        spooled = self.root / 'slurm_script'
        shutil.copyfile(self.bundle / '02_download_database.sbatch', spooled)
        self.run_script('', path=spooled)
        self.assertTrue((self.db / '.taxonomy_download_complete').exists())

    def test_invalid_config_and_environment_prevent_submission(self):
        original = self.config.read_text()
        for old, new in [('KRAKEN_KMER_LENGTH=35', 'KRAKEN_KMER_LENGTH=0'),
                         ('KRAKEN_SHARED_GROUP=""', 'KRAKEN_SHARED_GROUP="nonexistent-kraken-test-group"'),
                         ('(archaea viral fungi human bacteria)', '(unsupported)'),
                         (shlex.quote(str(self.storage)), '"/path/to/example"')]:
            self.config.write_text(original.replace(old, new))
            self.run_script('submit_all.sh', success=False)
            self.assertEqual(self.commands('sbatch'), [])
        self.config.write_text(original)
        self.run_script('submit_all.sh', success=False, extra_env={'FAIL_MODULE': 'yes'})
        self.assertEqual(self.commands('sbatch'), [])

    def test_failed_submission_stops_chain(self):
        self.run_script('submit_all.sh', success=False, extra_env={'FAIL_SUBMIT': 'yes'})
        self.assertEqual(len(self.commands('sbatch')), 1)

    def test_default_config_and_resource_directives(self):
        (self.bundle / 'config.sh').write_text(self.config.read_text())
        result = subprocess.run(['bash', str(self.bundle / 'submit_all.sh')],
                                env=self.env, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(self.commands('sbatch')[0][-2], str(self.bundle / 'config.sh'))
        for stage, cpus, memory, hours in [('02_download_database', 16, '32G', '3-00:00:00'),
                                         ('03_build_database', 32, '200G', '3-00:00:00'),
                                         ('04_validate_database', 4, '200G', '1-00:00:00')]:
            content = (self.bundle / (stage + '.sbatch')).read_text()
            for directive in [f'--cpus-per-task={cpus}', '--mem=' + memory, '--time=' + hours]:
                self.assertIn('#SBATCH ' + directive, content)

    def test_library_subset(self):
        self.config.write_text(self.config.read_text().replace(
            '(archaea viral fungi human bacteria)', '(bacteria)'))
        self.build()
        self.run_script('04_validate_database.sbatch', extra_env={'MISSING_TAXON': '9606'})
        self.assertTrue((self.db / '.database_validation_complete').exists())

    def test_download_build_resume_and_provenance(self):
        self.build()
        self.assertIn('>kraken:taxid|2|genome', (self.db / 'library/bacteria/library.fna').read_text())
        args = next(a for a in self.commands('kraken2-build') if '--build' in a)
        for flag, value in [('--threads', '32'), ('--kmer-len', '35'),
                            ('--minimizer-len', '31'), ('--minimizer-spaces', '7')]:
            self.assertEqual(args[args.index(flag) + 1], value)
        self.assertNotIn('--max-db-size', args)
        before = {p.name: p.read_bytes() for p in (self.db / 'provenance').iterdir()}
        downloads = len(self.commands('wget'))
        builds = len([a for a in self.commands('kraken2-build') if '--build' in a])
        self.run_script('02_download_database.sbatch')
        self.run_script('03_build_database.sbatch')
        self.assertEqual(downloads, len(self.commands('wget')))
        self.assertEqual(builds, len([a for a in self.commands('kraken2-build') if '--build' in a]))
        self.assertEqual(before, {p.name: p.read_bytes() for p in (self.db / 'provenance').iterdir()})
        self.run_script('04_validate_database.sbatch')
        self.assertTrue((self.db / '.database_validation_complete').exists())
        rows = dict(list(csv.reader((self.db / 'database_manifest.tsv').read_text().splitlines(), delimiter='\t'))[1:])
        for stage in ['download', 'build', 'validation']:
            self.assertRegex(rows[stage + '_completed_utc'], r'^\d{4}-.*Z$')
        self.assertIn('2.1.3', rows['kraken2_version'])
        self.assertIn('library/bacteria/https_manifest.tsv', (self.db / 'reference_metadata_checksums.sha256').read_text())
        self.assertNotIn(str(self.storage), (self.db / 'checksums.sha256').read_text())
        self.assertEqual((self.db / 'hash.k2d').stat().st_mode & 0o777, 0o644)
        stamp = (self.db / 'provenance/validation_completed_utc').read_bytes()
        self.run_script('04_validate_database.sbatch')
        self.assertEqual(stamp, (self.db / 'provenance/validation_completed_utc').read_bytes())

    def test_failed_build_and_validation_remove_success_markers(self):
        self.download()
        self.run_script('03_build_database.sbatch', success=False, extra_env={'FAIL_BUILD': 'yes'})
        self.assertFalse((self.db / '.database_build_complete').exists())
        self.assertFalse((self.db / 'provenance/build_completed_utc').exists())
        self.run_script('03_build_database.sbatch')
        for env in [{'MISSING_TAXON': '9606'}, {'FAIL_CHECKSUM': 'yes'}, {'FAIL_MODULE': 'yes'}]:
            (self.db / '.database_validation_complete').touch()
            self.run_script('04_validate_database.sbatch', success=False, extra_env=env)
            self.assertFalse((self.db / '.database_validation_complete').exists())
        self.run_script('04_validate_database.sbatch')


if __name__ == '__main__':
    unittest.main(verbosity=2)
