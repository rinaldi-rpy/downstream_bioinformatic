"""Single-entry launcher tests; inference tools are replaced by test stubs."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class LauncherTests(unittest.TestCase):
    def test_help_without_dependencies(self):
        result = subprocess.run(['bash', str(ROOT/'run_all.sh'), '--help'], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('bootstrap', result.stdout)
        self.assertIn('treemix', result.stdout)

    def test_launcher_selects_every_stage_and_preserves_arguments(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder/'shell').mkdir()
            (folder/'run_all.sh').write_text((ROOT/'run_all.sh').read_text())
            # Isolate only launcher behavior, not the scientific calculations.
            (folder/'shell/common.sh').write_text('need(){ :; }\nrequire_file(){ :; }\nrequire_plink19(){ :; }\n')
            recorder = folder/'run_pipeline.sh'
            recorder.write_text('printf "%s\\n" "$STEPS" "$1" "$2" > "$RECORD"\n')
            rscript = folder/'Rscript'
            rscript.write_text('#!/usr/bin/env bash\nexit 0\n')
            rscript.chmod(0o755)
            env = dict(os.environ, PYTHON=sys.executable, BCFTOOLS='bcftools', VCFTOOLS='vcftools',
                       ADMIXTURE='admixture', TREEMIX='treemix', RECORD=str(folder/'record'),
                       PATH=str(folder)+os.pathsep+os.environ['PATH'], STEPS='pca')
            result = subprocess.run(['bash', str(folder/'run_all.sh'), 'input with spaces.vcf.gz', 'metadata.csv'],
                                    env=env, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual((folder/'record').read_text().splitlines(),
                             ['qc,ld,pca,admixture,fst,bootstrap,treemix', 'input with spaces.vcf.gz', 'metadata.csv'])

    def test_rejects_missing_arguments(self):
        result = subprocess.run(['bash', str(ROOT/'run_all.sh')], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Usage:', result.stderr)

if __name__ == '__main__':
    unittest.main()
