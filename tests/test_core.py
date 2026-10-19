"""Regression tests for sample identity, SNP order, and statistical estimators."""
import csv
import gzip
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'python'))

class CoreTests(unittest.TestCase):
    def test_metadata_keeps_leading_zero_and_rejects_collision(self):
        from metadata import load_metadata
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)/'meta.csv'
            p.write_text('ID,Group Label\n001,A B\n002,A_B\n')
            with self.assertRaises(ValueError):
                load_metadata(p, 'Group Label')
            p.write_text('ID,Group Label\n001,A B\n002,C\n')
            self.assertEqual(load_metadata(p, 'Group Label')[0]['001'], 'A_B')

    def test_treemix_preserves_order_and_uses_exact_counts(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'population_sizes.csv').write_text('Population,N_individuals\nA,2\nB,2\n')
            (p/'autosomal.frq.strat').write_text(
                'CHR SNP CLST A1 A2 MAF MAC NCHROBS\n'
                '1 z A A G .5 2 4\n1 z B A G .25 1 4\n'
                '1 a A A G .25 1 4\n1 a B A G .75 3 4\n')
            env=dict(os.environ,OUT_DIR=tmp,BLOCK_SIZE='1')
            subprocess.run([sys.executable,str(ROOT/'python/treemix_convert.py')],env=env,check=True,capture_output=True)
            self.assertEqual((p/'retained_snps.tsv').read_text().splitlines()[1:],['1\tz','1\ta'])
            with gzip.open(p/'treemix_input.gz','rt') as f:
                self.assertEqual(f.read(),'A B\n2,2 1,3\n1,3 3,1\n')

    def test_wc_complete_separation_and_no_structure(self):
        import numpy as np
        from fst_analysis import wc_components
        n=np.array([[2.,2.],[2.,2.]])
        a,d,valid=wc_components(n,np.array([[0.,4.],[2.,2.]]),np.array([[0.,0.],[2.,2.]]),np.array([0]),np.array([1]),2)
        self.assertTrue(valid.all())
        self.assertAlmostEqual((a/d)[0,0],1.)
        self.assertAlmostEqual((a/d)[1,0],0.)

    def test_bootstrap_reproducible(self):
        import numpy as np
        from fst_analysis import joint_bootstrap
        a=np.array([[1.,2.],[2.,3.],[3.,4.],[4.,5.]])
        d=np.full_like(a,10.)
        blocks=[('1',0),('1',1),('2',0),('2',1)]
        first=joint_bootstrap(a,d,blocks,20,42)
        second=joint_bootstrap(a,d,blocks,20,42)
        np.testing.assert_array_equal(first[0],second[0])
        self.assertEqual(first[2],0)

    def test_help_does_not_need_data(self):
        result=subprocess.run(['bash',str(ROOT/'run_pipeline.sh'),'--help'],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)
        self.assertIn('STEPS',result.stdout)

if __name__=='__main__':
    unittest.main()
