"""Orchestration tests; external PCA execution is replaced with a command recorder."""
import csv
import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from test_core import ROOT

class WorkflowTests(unittest.TestCase):
    def test_partial_pca_run_records_effective_pc_and_preserves_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); prior=p/'prior/02_ld';prior.mkdir(parents=True)
            prefix=prior/'merged.pruned'
            Path(str(prefix)+'.bed').write_bytes(b'fake-command-test-only')
            Path(str(prefix)+'.bim').write_text('1 a 0 1 A G\n1 b 0 2 A G\n')
            Path(str(prefix)+'.fam').write_text('s1 s1 0 0 0 -9\ns2 s2 0 0 0 -9\ns3 s3 0 0 0 -9\n')
            meta=p/'meta.csv';meta.write_text('ID,Group Label\ns1,A\ns2,B\ns3,C\n')
            raw=p/'raw.vcf.gz';raw.write_bytes(b'data')
            tool=p/'plink';tool.write_text('#!/usr/bin/env bash\nif [[ ${1:-} == --version ]]; then echo \"PLINK v1.90 synthetic command recorder\"; exit 0; fi\nprintf "%s\\n" "$@" >> "${COMMAND_LOG}"\n');tool.chmod(0o755)
            run=p/'new'
            env=dict(os.environ,STEPS='pca',INPUT_RUN=str(p/'prior'),RUN_DIR=str(run),PLINK=str(tool),PYTHON=sys.executable,COMMAND_LOG=str(p/'commands'))
            result=subprocess.run(['bash',str(ROOT/'run_pipeline.sh'),str(raw),str(meta)],env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr+result.stdout)
            self.assertIn('253\t2\t3\t2',(run/'03_pca/pca_dimensions.tsv').read_text())
            self.assertEqual(Path(str(prefix)+'.bed').read_bytes(),b'fake-command-test-only')
            self.assertTrue((run/'provenance.json').is_file())

    def test_resume_summarizes_existing_runs_without_reanalysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'runs').mkdir()
            (p/'settings.txt').write_text('MAX_M=1\nREPS=2\n')
            (p/'population_order.txt').write_text('A\nB\nC\n')
            for m in range(2):
                for rep in range(1,3):
                    prefix=p/'runs'/f'treemix_m{m}_rep{rep}'
                    Path(str(prefix)+'.llik').write_text(f'Starting ln(likelihood): -100\nExiting ln(likelihood): {-30+m*10+rep}\niterations: 1000\n')
                    with gzip.open(str(prefix)+'.treeout.gz','wt') as handle:handle.write('(A,(B,C));\n')
            before={q.name:q.read_bytes() for q in (p/'runs').iterdir()}
            result=subprocess.run(['bash',str(ROOT/'summarize_treemix.sh'),str(p)],env=dict(os.environ,PYTHON=sys.executable),capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr+result.stdout)
            outputs=list(p.glob('postprocess_*'));self.assertEqual(len(outputs),1)
            with (outputs[0]/'best_runs.csv').open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([r['Replicate'] for r in rows],['2','2'])
            self.assertEqual([float(r['LogLikelihood']) for r in rows],[-28,-18])
            self.assertEqual(before,{q.name:q.read_bytes() for q in (p/'runs').iterdir()})

    def test_metadata_population_cli_subsets_real_samples(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'meta.csv').write_text('ID,Group Label\n001,A\n002,B\n003,C\n')
            (p/'samples').write_text('001\n002\n')
            result=subprocess.run([sys.executable,str(ROOT/'python/metadata.py'),'populations','--meta',str(p/'meta.csv'),'--samples',str(p/'samples'),'--out',str(p/'pop')],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(sorted(q.name for q in (p/'pop').glob('*.txt')),['A.txt','B.txt'])
            self.assertEqual((p/'pop/A.txt').read_text(),'001\n')

if __name__=='__main__':unittest.main()
