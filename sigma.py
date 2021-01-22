#!/usr/bin/env python

#-----------------------------------------------------------------------------#
# requirements
#   bwa
#   samtools
#   prokka
#   haplomerger2
#   python lib
#     pysam
#-----------------------------------------------------------------------------#
import argparse
import os
import shutil
import sys

from version import __version__
from libsigma import Checkm, GenomeRna, Hm2, MaClassifier, MaFasta, Prokka, RefFastaSet, TwoFasta

DEFAULT_BDIR = 'sgbin'
DEFAULT_GDIR = 'sigma'
DEFAULT_MDIR = 'intermediate'
DEFAULT_HM2 = None
DEFAULT_CPU = 1
DEFAULT_MLEN = 1000
DEFAULT_SLEN = 10000
DEFAULT_IDENT = 0.99
DEFAULT_HM2_ENV = 'HAPLOMERGER2_PATH'

#-- binning process --#
def proc_binning(args):
    mdir = os.path.join(args.odir, DEFAULT_MDIR)
    os.makedirs(mdir, exist_ok=True)
    ref = RefFastaSet(mdir)
    ref.collect(args.refs)
    ma = MaFasta(mdir)
    ma.load(args.ma, args.mlen)
    classifier = MaClassifier(args.odir, mdir)
    classifier.classify(ma, ref, cpu=args.cpu, ident=args.ident)
    if not args.keep:
        ref.remove_intermediate()
        ma.remove_intermediate()
        classifier.remove_intermediate()

#-- merge process --#
def proc_merge(args):
    os.makedirs(args.odir, exist_ok=True)
    #-- create objects --#
    tf = TwoFasta(args.odir, args.nrsag, args.sgbin)
    checkm = Checkm(args.odir)
    hm2 = Hm2(args.odir, os.environ.get(DEFAULT_HM2_ENV, args.hm2))
    prokka = Prokka(args.odir)
    grna = GenomeRna(args.odir)
    #-- run --#
    checkm.run(tf.fasta_dir, args.cpu)
    tf.set_master_slave(checkm.nrsag_comp, checkm.sgbin_comp)
    tf.create_main_sub(args.mlen, args.slen)
    hm2.run(tf.main, args.cpu)
    main_fa = hm2.merge_fasta if os.path.exists(hm2.merge_fasta) else tf.master
    prokka.run(main_fa, tf.sub)
    grna.salvage(prokka.main_fna, prokka.main_gff, prokka.sub_fna, prokka.sub_gff)
    fasta = os.path.join(args.odir, 'mgsag.fasta' if tf.master_name == 'nrSAG' else 'sgmag.fasta')
    shutil.move(grna.fasta, fasta)
    if not args.keep:
        tf.remove_intermediate()
        checkm.remove_intermediate()
        hm2.remove_intermediate()
        prokka.remove_intermediate()

#-- positional arguments --#
def opt_ma_fasta(op):
    op.add_argument('ma', metavar='ma_fasta',
                    type=str, help='metagenome assembly (fasta)')
def opt_nrsag_fasta(op):
    op.add_argument('nrsag', metavar='nrsag', type=str, help='nrsag (fasta)')
def opt_refs(op):
    op.add_argument('refs', metavar='sag_fasta', nargs='+',
                    type=str, help='reference SAGs (fasta)')
def opt_sgbin_fasta(op):
    op.add_argument('sgbin', metavar='sgbin', type=str, help='sgbin (fasta)')

#-- optional arguments --#
def opt_haplomerger2(op):
    op.add_argument('--haplomerger2', action='store', dest='hm2',
                    default=DEFAULT_HM2, type=str,
                    metavar='DIR', help='haplomerger2 directory, default: {}'.format(DEFAULT_HM2))
def opt_keep(op):
    op.add_argument('--keep-intermediate', action='store_true', dest='keep',
                    help='do not remove intermediate files')
def opt_lower_ident(op):
    op.add_argument('-l', '--lower-identity', action='store', dest='ident',
                    default=DEFAULT_IDENT, type=float,
                    metavar='FLOAT', help='lower identity of alignment, default: {}'.format(DEFAULT_IDENT))
def opt_min_length(op):
    op.add_argument('-m', '--min-length', action='store', dest='mlen',
                    default=DEFAULT_MLEN, type=int,
                    metavar='INT', help='minimum length of metagenome assembly, default: {}'.format(DEFAULT_MLEN))
def opt_min_master_length(op):
    op.add_argument('--min-master', action='store', dest='mlen',
                    default=DEFAULT_MLEN, type=int,
                    metavar='INT', help='minimum length of master assembly, default: {}'.format(DEFAULT_MLEN))
def opt_min_slave_length(op):
    op.add_argument('--min-slave', action='store', dest='slen',
                    default=DEFAULT_SLEN, type=int,
                    metavar='INT', help='minimum length of slave assembly, default: {}'.format(DEFAULT_SLEN))
def opt_output(op, default_dir):
    op.add_argument('-o', '--output-dir', action='store', dest='odir',
                    default=default_dir, type=str,
                    metavar='DIR', help='output directory, default: {}'.format(default_dir))
def opt_output_binning(op):
    opt_output(op, DEFAULT_BDIR)
def opt_output_merge(op):
    opt_output(op, DEFAULT_GDIR)
def opt_threads(op):
    op.add_argument('-t', '--threads', action='store', dest='cpu',
                    default=DEFAULT_CPU, type=int,
                    metavar='INT', help='number of threads, default: {}'.format(DEFAULT_CPU))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--version', action='version', version='%(prog)s ' + __version__)
    sp = ap.add_subparsers()

    #-- binning --#
    sp_binning = sp.add_parser('binning', help='classify metagenome assembly into sgbin')
    for opt_func in (opt_ma_fasta, opt_refs, opt_keep, opt_lower_ident, opt_min_length, opt_output_binning, opt_threads):
        opt_func(sp_binning)
    sp_binning.set_defaults(handler=proc_binning)

    #-- merge --#
    sp_merge = sp.add_parser('merge', help='merge nrsag and sgbin')
    for opt_func in (opt_nrsag_fasta, opt_sgbin_fasta, opt_haplomerger2, opt_keep, opt_min_master_length, opt_min_slave_length, opt_output_merge, opt_threads):
        opt_func(sp_merge)
    sp_merge.set_defaults(handler=proc_merge)

    args = ap.parse_args()
    if hasattr(args, 'handler'):
        args.handler(args)
    else:
        ap.print_help()

if __name__ == '__main__':
    main()
