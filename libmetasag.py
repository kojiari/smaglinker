###############################################################################
# class list
#   Checkm
#   MaFasta
#   MaClassifier
#   RefFastaSet
###############################################################################
from collections import defaultdict
from distutils.dir_util import copy_tree
import csv
import glob
import gzip
import os
import pysam
import re
import shutil
import subprocess
import sys

class Checkm:
    def __init__(self, odir):
        self.checkm_dir = os.path.join(odir, 'checkm')
        self.checkm_res = os.path.join(odir, 'checkm.result')
        self.nrsag_comp = 0.0
        self.sgbin_comp = 0.0

    def remove_intermediate(self):
        shutil.rmtree(self.checkm_dir)

    def run(self, idir, cpu=1):
        cmds = ['checkm', 'lineage_wf', '-r', '--tab_table',
                '-t', str(cpu),
                '-x', 'fasta',
                '-f', self.checkm_res,
                idir,
                self.checkm_dir]
        subprocess.run(cmds)
        with open(self.checkm_res) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                self.__dict__[row['Bin Id'] + '_comp'] = float(row['Completeness'])

class GenomeRna:
    def __init__(self, odir):
        self.fasta = os.path.join(odir, 'final.contigs.fasta')

    def _rrna_and_trna(self, gff):
        with open(gff) as fh:
            for line in fh:
                if line[0] == '>':
                    break
                if line[0] == '#':
                    continue
                cols = line.rstrip().split('\t')
                if cols[2] != 'rRNA' and cols[2] != 'tRNA':
                    continue
                attrs = dict(list(map(lambda v: v.split('='), cols[8].split(';'))))
                yield(cols[0], attrs['product'])

    def salvage(self, mfna, mgff, sfna, sgff):
        main_prod = {}
        sub_contig = {}
        for contig_id, prod in self._rrna_and_trna(mgff):
            main_prod[prod] = True
        for contig_id, prod in self._rrna_and_trna(sgff):
            if prod not in main_prod:
                sub_contig[contig_id] = True
        with open(self.fasta, 'w') as fo:
            with open(mfna) as fh:
                print(fh.read(), file=fo)
            for fa_id, seqs in MaFasta.each_fasta(sfna):
                if fa_id in sub_contig:
                    print('>' + fa_id, file=fo)
                    print('\n'.join(seqs), file=fo)

class Hm2:
    def __init__(self, odir, hm2_root='HaploMerger2_20180603'):
        self.hm2_lastz = os.path.join(hm2_root, 'lastz_1.02.00_centOS6')
        self.hm2_chainnet = os.path.join(hm2_root, 'chainNet_jksrc20100603_centOS6')
        self.hm2_bin = os.path.join(hm2_root, 'bin')
        self.hm2_template = os.path.join(hm2_root, 'project_template')
        self.merge_dir = os.path.join(odir, 'hm2')
        self.merge_bin = os.path.join(self.merge_dir, 'bin')
        self.merge_work = os.path.join(self.merge_dir, 'work')
        self.asm_name = 'assembly'
        self.asm = os.path.join(self.merge_work, '{}.fa.gz'.format(self.asm_name))
        self.merge_fa = os.path.join(self.merge_work, '{}_A_ref_D.fa.gz'.format(self.asm_name))
        self.merge_fasta = os.path.join(self.merge_dir, 'merge.fasta')

    def remove_intermediate(self):
        shutil.rmtree(self.merge_dir)

    def run(self, gzfasta, cpu=1):
        copy_tree(self.hm2_bin, self.merge_bin)
        copy_tree(self.hm2_template, self.merge_work)
        shutil.copy(gzfasta, self.asm)
        cdir = os.getcwd()
        os.chdir(self.merge_work)
        #-- confiture --#
        cmds = ['sed', '-i', '-e', 's/threads=1/threads={}/g'.format(cpu)]
        cmds.extend(glob.glob('hm.batch*'))
        subprocess.run(cmds)
        subprocess.run(['sed', '-i', '-e', 's#^PATH=.*#PATH={}:{}:\$PATH#g'.format(self.hm2_lastz, self.hm2_chainnet), 'run_all.batch'])
        subprocess.run(['sed', '-i', '-e', 's/your_assembly_name/{}/g'.format(self.asm_name), 'run_all.batch'])
        subprocess.run(['sed', '-i', '-e', 's#./hm.batchC.*##g', 'run_all.batch'])
        subprocess.run(['sed', '-i', '-e', 's#./hm.batchE.*##g', 'run_all.batch'])
        subprocess.run(['sed', '-i', '-e', 's/{}_A_ref_C/{}_A_ref/g'.format(self.asm_name, self.asm_name), 'run_all.batch'])
        subprocess.run('./run_all.batch')
        os.chdir(cdir)
        if os.path.exists(self.merge_fa):
            with gzip.open(self.merge_fa, 'rb') as fh:
               with open(self.merge_fasta, 'wb') as fo:
                   shutil.copyfileobj(fh, fo)

class MaFasta:
    def __init__(self, odir):
        self.fasta = os.path.join(odir, 'ma.fasta')
        self.contigs = {}

    def load(self, path, mlen):
        for faid, seqs in MaFasta.each_fasta(path, mlen):
            self.contigs[faid] = seqs
        with open(self.fasta, 'w') as fo:
            for faid, seqs in self.contigs.items():
                print('>{}'.format(faid), file=fo)
                print('\n'.join(seqs), file=fo)

    def remove_intermediate(self):
        os.remove(self.fasta)

    @classmethod
    def each_fasta(self, path, mlen=0, xlen=0):
        with open(path) as fh:
            faid = None
            seqs = []
            for line in fh:
                if line[0] == '>':
                    seqlen = len(''.join(seqs))
                    if faid is not None and seqlen > mlen and (xlen == 0 or seqlen < xlen):
                        yield(faid, seqs)
                    faid = line.strip().split(' ')[0].replace('>', '')
                    seqs = []
                else:
                    seqs.append(line.strip())
            if faid is not None and seqlen > mlen and (xlen == 0 or seqlen < xlen):
                yield(faid, seqs)

class MaClassifier:
    def __init__(self, odir, mdir='intermediate'):
        self.pref = os.path.join(mdir, 'ref')
        self.bam = os.path.join(mdir, 'ma2ref.bam')
        self.flt = os.path.join(mdir, 'ma2ref.flt.bam')
        self.classification = os.path.join(odir, 'classification')
        self.bindir = os.path.join(odir)
        self.bins = defaultdict(list)

    def classify(self, ma, ref, cpu=1, ident=0.99, min_match=200):
        self._index(ref.fasta)
        self._mapping(ma.fasta, cpu=cpu)
        self._filter(ident=ident)
        self._classify(ma.contigs, min_ident=ident, min_match=min_match)
        self._binning(ma.contigs)

    def remove_intermediate(self):
        os.remove('{}.amb'.format(self.pref))
        os.remove('{}.ann'.format(self.pref))
        os.remove('{}.bwt'.format(self.pref))
        os.remove('{}.pac'.format(self.pref))
        os.remove('{}.sa'.format(self.pref))
        os.remove(self.bam)
        os.remove(self.flt)

    def _filter(self, ident=0.99):
        with pysam.AlignmentFile(self.bam, 'rb') as fh:
            fo = pysam.AlignmentFile(self.flt, "wb", template=fh)
            for read in fh:
                match = sum(list(map(lambda v: v[1], filter(lambda v: v[0] == 0, read.cigartuples))))
                align = len(list(filter(lambda v: v[2].isupper(), read.get_aligned_pairs(matches_only=True, with_seq=True))))
                if float(align) / float(match) >= ident:
                    fo.write(read)
            fo.close()

    def _index(self, fasta):
        cmds = ['bwa', 'index', '-p', self.pref, fasta]
        subprocess.run(cmds)

    def _mapping(self, fasta, cpu=1):
        tmp = '{}.tmp'.format(self.bam)
        cmd1 = ['bwa', 'mem', '-t', str(cpu), self.pref, fasta]
        cmd2 = ['samtools', 'view', '-F', '2052', '-b']
        cmd3 = ['samtools', 'sort', '-@', str(cpu), '-o', self.bam, '-T', tmp]
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stdin=p1.stdout)
        subprocess.run(cmd3, stdin=p2.stdout)

    def _classify(self, contigs, min_ident=0.99, min_match=200):
        with open(self.classification, 'w') as fo:
            classified = {}
            with pysam.AlignmentFile(self.bam, 'rb') as fh:
                for read in fh:
                    match = sum(list(map(lambda v: v[1], filter(lambda v: v[0] == 0, read.cigartuples))))
                    align = len(list(filter(lambda v: v[2].isupper(), read.get_aligned_pairs(matches_only=True, with_seq=True))))
                    ident = float(align) / float(match)
                    rec = [read.query_name, len(''.join(contigs[read.query_name])), match, ident]
                    if ident < min_ident or match < min_match:
                        rec.append('|'.join(['unclassified', read.reference_name]))
                        self.bins['unclassified'].append(read.query_name)
                    else:
                        rec.append(read.reference_name)
                        self.bins[read.reference_name.split('|')[0]].append(read.query_name)
                    classified[read.query_name] = True
                    print('\t'.join(list(map(lambda v: str(v), rec))), file=fo)
            for maid, seqs in contigs.items():
                if maid in classified:
                    continue
                rec = [maid, len(''.join(seqs)), 0, 0, 'unclassified']
                self.bins['unclassified'].append(maid)
                print('\t'.join(list(map(lambda v: str(v), rec))), file=fo)

    def _binning(self, contigs):
        os.makedirs(self.bindir, exist_ok=True)
        for refid, maid_lst in self.bins.items():
            fname = os.path.join(self.bindir, '{}.fasta'.format(refid))
            with open(fname, 'w') as fo:
                for maid in maid_lst:
                    print('>{}'.format(maid), file=fo)
                    print('\n'.join(contigs[maid]), file=fo)

class Prokka:
    def __init__(self, odir):
        self.prokka_dir = os.path.join(odir, 'prokka')
        self.pref_main = 'main'
        self.pref_sub = 'sub'
        self.main_gff = os.path.join(self.prokka_dir, self.pref_main + '.gff')
        self.main_fna = os.path.join(self.prokka_dir, self.pref_main + '.fna')
        self.sub_gff = os.path.join(self.prokka_dir, self.pref_sub + '.gff')
        self.sub_fna = os.path.join(self.prokka_dir, self.pref_sub + '.fna')

    def remove_intermediate(self):
        shutil.rmtree(self.prokka_dir)

    def run(self, main, sub, cpu=1):
        for v in ((self.pref_main, main), (self.pref_sub, sub)):
            cmds = ['prokka', '--force', '--rawproduct', '--cpus', str(cpu),
                    '--prefix', v[0], '--locustag', v[0], '--centre', 'bB',
                    '--mincontiglen', '200', '--outdir', self.prokka_dir, v[1]]
            subprocess.run(cmds)

class RefFastaSet:
    def __init__(self, odir):
        self.fasta = os.path.join(odir, 'ref.fasta')

    def collect(self, refs):
        lines = []
        for path in refs:
            fname = os.path.basename(path)
            name = os.path.splitext(fname)[0]
            with open(path) as fh:
                for line in fh:
                    line = re.sub('^>', '>{}|'.format(name), line.strip())
                    lines.append(line)
        with open(self.fasta, 'w') as fo:
            for line in lines:
                print(line, file=fo)

    def remove_intermediate(self):
        os.remove(self.fasta)

class TwoFasta:
    def __init__(self, odir, nrsag, sgbin):
        self.fasta_dir = os.path.join(odir, 'fasta')
        self.nrsag = os.path.join(self.fasta_dir, 'nrsag.fasta')
        self.sgbin = os.path.join(self.fasta_dir, 'sgbin.fasta')
        self.main = os.path.join(self.fasta_dir, 'main.fasta')
        self.sub = os.path.join(self.fasta_dir, 'sub.fasta')
        os.makedirs(self.fasta_dir, exist_ok=True)
        shutil.copy(nrsag, self.nrsag)
        shutil.copy(sgbin, self.sgbin)
        self.nrsag_len = 0
        for sid, seqs in MaFasta.each_fasta(self.nrsag):
            self.nrsag_len += len(''.join(seqs))
        self.sgbin_len = 0
        for sid, seqs in MaFasta.each_fasta(self.sgbin):
            self.sgbin_len += len(''.join(seqs))
        self.master = self.nrsag
        self.slave = self.sgbin
        self.master_name = 'nrSAG'

    def create_main_sub(self, master_min, slave_min):
        with open(self.main, 'w') as fo:
            for fa_id, seqs in MaFasta.each_fasta(self.master, master_min):
                print('>' + fa_id, file=fo)
                print('\n'.join(seqs), file=fo)
            for fa_id, seqs in MaFasta.each_fasta(self.slave, slave_min):
                print('>' + fa_id, file=fo)
                print('\n'.join(seqs), file=fo)
        with open(self.sub, 'w') as fo:
            for fa_id, seqs in MaFasta.each_fasta(self.master, xlen=master_min + 1):
                print('>' + fa_id, file=fo)
                print('\n'.join(seqs), file=fo)
            for fa_id, seqs in MaFasta.each_fasta(self.slave, xlen=slave_min + 1):
                print('>' + fa_id, file=fo)
                print('\n'.join(seqs), file=fo)

    def remove_intermediate(self):
        shutil.rmtree(self.fasta_dir)

    def set_master_slave(self, nrsag_comp, sgbin_comp):
        self.master = self.nrsag
        self.slave = self.sgbin
        self.master_name = 'nrSAG'
        if (nrsag_comp == sgbin_comp and self.nrsag_len < self.sgbin_len) or (nrsag_comp <= sgbin_comp):
            self.master = self.sgbin
            self.slave = self.nrsag
            self.master_name = 'sgBin'

