import csv
import gzip
import itertools
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import numpy as np
from test_core import ROOT
sys.path.insert(0,str(ROOT/'python'))
from metadata import load_metadata,match_samples
from likelihood import read_likelihood
from fst_analysis import wc_components,read_blocks,aggregate_genotypes,joint_bootstrap
from qc_summary import summarize

class ValidationTests(unittest.TestCase):
    def test_bad_metadata(self):
        for text in ['ID,Group Label\nx,A\nx,B\n','ID,Group Label\nx,\n','Other,Group Label\nx,A\n','ID,Group Label\nx y,A\n','ID,Group Label\n']:
            with self.subTest(text=text),tempfile.TemporaryDirectory() as tmp:
                p=Path(tmp)/'meta.csv';p.write_text(text)
                with self.assertRaises(ValueError):load_metadata(p)

    def test_fam_uses_actual_fid_and_fails_unmatched(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.fam';p.write_text('family 001 0 0 0 -9\n')
            self.assertEqual(match_samples(p,{'001':'A'}),[('family','001','A')])
            with self.assertRaises(ValueError):match_samples(p,{'002':'A'})

    def test_likelihood_final_line_and_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.llik'
            p.write_text('Starting ln(likelihood): -999\nExiting ln(likelihood): -123.5\niterations: 1000\n')
            self.assertEqual(read_likelihood(p),-123.5)
            for value in ['Exiting ln(likelihood): nan\n','m=3\niterations: 10\n','Exiting ln(likelihood): inf\n']:
                p.write_text(value)
                with self.assertRaises(ValueError):read_likelihood(p)

    def test_qc_info_depth_is_not_format_depth(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'x.vcf'
            p.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\n'
                         '1\t1\t.\tA\tG\t20\tPASS\tXDP=999;DP=12\tGT:DP\t0/0:5\t0/1:7\n'
                         '1\t2\t.\tA\tAT\t.\tPASS\t.\tGT\t0/0\t0/1\n')
            records,snps,stats=summarize(p)
            self.assertEqual((records,snps),(2,1))
            self.assertEqual(stats['INFO_DP'].mean,12)
            self.assertEqual(stats['FORMAT_DP'].mean,6)
            self.assertAlmostEqual(stats['FORMAT_DP'].values()[-1],2**.5)

    def test_negative_fst_and_missing_are_retained_or_excluded_correctly(self):
        a,d,valid=wc_components(np.array([[2.,2.],[1.,2.]]),np.array([[2.,2.],[1.,2.]]),np.zeros((2,2)),np.array([0]),np.array([1]),2)
        self.assertAlmostEqual(a[0,0]/d[0,0],-1.)
        self.assertFalse(valid[1,0]);self.assertEqual(a[1,0],0)

    def test_packed_bed_chunk_boundaries_match_direct_calculation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);prefix=str(p/'x')
            # 8 individuals, 4 populations; complete diploid BED codes.
            codes=np.array([[0,0,3,3,0,2,2,3],[0,2,2,3,0,0,3,3],
                            [3,3,0,0,2,2,0,0],[0,2,2,3,2,3,0,2],
                            [0,0,2,3,3,3,0,2],[2,2,0,0,0,2,3,3]],dtype=np.uint8)
            packed=np.zeros((len(codes),2),dtype=np.uint8)
            for i in range(8):packed[:,i//4] |= codes[:,i] << (2*(i%4))
            Path(prefix+'.bed').write_bytes(b'\x6c\x1b\x01'+packed.tobytes())
            subpops=[np.array([i,i+1]) for i in range(0,8,2)]
            ids=np.array([0,0,1,1,2,2]);pairs=list(itertools.combinations(range(4),2))
            first=aggregate_genotypes(prefix,8,subpops,ids,pairs,3,1,2)
            second=aggregate_genotypes(prefix,8,subpops,ids,pairs,3,4,2)
            for one,two in zip(first,second):np.testing.assert_allclose(one,two)
            # Independent scalar variance-component calculation from diploid allele dosages.
            a1_count=np.where(codes==0,2,np.where(codes==2,1,0))
            expected_a=np.zeros((3,6));expected_d=expected_a.copy()
            for snp,row in enumerate(codes):
                for k,(i,j) in enumerate(pairs):
                    p1=a1_count[snp,subpops[i]].sum()/4
                    p2=a1_count[snp,subpops[j]].sum()/4
                    mean=(p1+p2)/2;s2=(p1-mean)**2+(p2-mean)**2
                    h=((row[subpops[i]]==2).sum()+(row[subpops[j]]==2).sum())/4
                    a=s2-(mean*(1-mean)-s2/2-h/4)
                    b=2*(mean*(1-mean)-s2/2-3*h/8)
                    denominator=a+b+h/2
                    if denominator>0:expected_a[ids[snp],k]+=a;expected_d[ids[snp],k]+=denominator
            np.testing.assert_allclose(first[0],expected_a)
            np.testing.assert_allclose(first[1],expected_d)

    def test_fst_python_complete_fixture(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp); prefix=str(p/'x');out=p/'out';out.mkdir()
            (p/'meta.csv').write_text('ID,Group Label\n'+''.join(f's{i},P{i//2}\n' for i in range(8)))
            Path(prefix+'.fam').write_text(''.join(f'fid{i} s{i} 0 0 0 -9\n' for i in range(8)))
            Path(prefix+'.bim').write_text(''.join(f'NC_059157.1 v{i} 0 {i*1000000+1} A G\n' for i in range(4)))
            Path(prefix+'.bed').write_bytes(b'\x6c\x1b\x01'+bytes([0xf0,0xaa,0xaa,0xf0,0x0f,0xaa,0xaa,0x0f]))
            env=dict(os.environ,OUT_DIR=str(out),ANALYSIS_BFILE=prefix,META_PATH=str(p/'meta.csv'),FST_DIR=tmp,
                     AUTOSOME_MAP=str(ROOT/'config/autosomes.tsv'),BLOCK_BP='1000000',N_BOOT='20',SEED='42',CHUNK_SNPS='3',MIN_CALLED='2')
            result=subprocess.run([sys.executable,str(ROOT/'python/fst_analysis.py')],env=env,capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            with (out/'pairwise_fst_summary.csv').open() as handle: rows=list(csv.DictReader(handle))
            self.assertEqual(len(rows),6)
            self.assertTrue(all(int(r['N_bootstrap'])==20 for r in rows))
            self.assertTrue((out/'block_components.npz').is_file())

    def test_treemix_rejects_inconsistent_orientation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/'population_sizes.csv').write_text('Population,N_individuals\nA,2\nB,2\n')
            (p/'autosomal.frq.strat').write_text('CHR SNP CLST A1 A2 MAF MAC NCHROBS\n1 x A A G .5 2 4\n1 x B G A .5 2 4\n')
            result=subprocess.run([sys.executable,str(ROOT/'python/treemix_convert.py')],env=dict(os.environ,OUT_DIR=tmp,BLOCK_SIZE='1'),capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('orientation',result.stderr)

    def test_admixture_summary_preserves_replicate_cv(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            for rep,cv in [(1,.6),(2,.4)]:
                d=p/f'K2_rep{rep}';d.mkdir()
                (d/'admixture.log').write_text(f'CV error (K=2): {cv}\n')
                (d/'input.2.Q').write_text(str(rep));(d/'input.2.P').write_text(str(rep))
            subprocess.run([sys.executable,str(ROOT/'python/admixture_summary.py'),tmp,'2','2','2'],check=True)
            self.assertEqual((p/'best/K2.Q').read_text(),'2')
            with (p/'cv_all_replicates.csv').open() as f:self.assertEqual(len(list(csv.DictReader(f))),2)

if __name__=='__main__':unittest.main()
