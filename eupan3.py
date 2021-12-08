#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/11/1 8:30 PM
    @Usage: python3 eupan3.py command [options]
"""
import argparse
import os
import sys

DESCRIPTION = """
EUPAN3: EUkaryote Pan-genome ANalysis toolkit for 3rd generation sequencing
                                                                                                         
"""
VERSION = "Version 0.1.0"

if __name__ == "__main__":
    dir_name, file_name = os.path.split(os.path.abspath(sys.argv[0]))
    parser = argparse.ArgumentParser(description=DESCRIPTION + VERSION)
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    sub_parser = parser.add_subparsers(dest='subcmd', title='sub commands')
    # AssemSta
    parser_as = sub_parser.add_parser("assemsta", help='''Use Quast to get assemblies statistics''')
    parser_as.add_argument('-a', '--assembly_path', metavar='<assembly.fa>',
                           help='Path of contigs or scaffolds', type=str, required=True)
    parser_as.add_argument('-r', '--reference_path', metavar='<reference.fa>',
                           help='Path of reference genome', type=str, required=True)
    parser_as.add_argument('-o', '--output_dir', metavar='<output_directory>',
                           help='Path of output of Quast', type=str, required=True)
    parser_as.add_argument('-q', '--quast_path', metavar='<quast.py>',
                           help='Path of quast.py (default: quast.py in $PATH)', type=str, default='quast.py')
    parser_asq = parser_as.add_argument_group('quast parameters')
    parser_asq.add_argument('-t', '--thread', metavar='<int>',
                            help='number of threads using minimap2 in Quast (default: 1)', type=int, default=1)
    parser_asq.add_argument('-l', '--length_threshold', metavar='<int>',
                            help='Min length of contigs/scaffolds to consider (default: 500, if --large_genome, '
                                 'default: 3000)',
                            type=int, default=-1)
    parser_asq.add_argument('-g', '--gff_path', metavar='<reference.gff>',
                            help='Anotations of reference genome (default: None)', type=str, default=None)
    parser_asq.add_argument('-i', '--identity', metavar='<int>',
                            help='Min alignment identity of sequences (choices: 80/90/95. default: 90)', type=int,
                            choices=[80, 90, 95],
                            default=90)
    parser_asq.add_argument('--large_genome', default=False, action='store_true',
                            help='Use optimal parameters for evaluation of large genomes (typically > 100 Mbp)')

    # UnalnBlockSeq
    parser_ubs = sub_parser.add_parser("unalnbseq", help='''Get partial/full unaligned block sequences''')
    parser_ubs.add_argument('-a', '--assembly_path', metavar='<assembly.fa>',
                            help='Path of contigs or scaffolds', type=str, required=True)
    parser_ubs.add_argument('-u', '--unaln_path', metavar='<unaligned.info>',
                            help='Path of unaligned table (at quast_output/contigs_reports/contigs_report_xxx.unaligned'
                                 '.info)',
                            type=str, required=True)
    parser_ubs.add_argument('-o', '--output_path', metavar='<output.fa>',
                            help='Path of output unaln sequences', type=str, required=True)
    parser_ubs.add_argument('-l', '--length_filter', metavar='<int>',
                            help='Min length of block sequences to consider (default: 500)', type=int, default=500)
    parser_ubs.add_argument('-k', '--kind_unaln', metavar='<str>',
                            help='Use full/partial/all unaligned sequences (choices: full/partial/all. default: all)',
                            choices=['full', 'partial', 'all'], type=str, default='all')
    parser_ubs.add_argument('-s', '--sample_tag', metavar='<str>',
                            help='Add sample tag before each contig/scaffold \
                            (default: None, e.g. -s Sample1 , "_" will be used, Chr1 -> Sample1_Chr1)', type=str,
                            default='')
    parser_ubs.add_argument('-n', '--nbase_ignore', metavar='<int>',
                            help='Max percentage of N bases in block sequences to ignore (default: 100)',
                            type=int, default=100)
    parserr_ubsr = parser_ubs.add_argument_group('realign parameters')
    parserr_ubsr.add_argument('-rr', '--realign_reference', metavar='<reference.fa>',
                              help='Using minimap2 to realign  to reference genome/mitochondrion/plastid and drop high '
                                   'similar unaligned block sequences',
                              type=str, default='')
    parserr_ubsr.add_argument('-ri', '--realign_identity', metavar='<int>',
                              help='[Only use when -rr on] Min alignment identity of sequences in realign step ('
                                   'choices: '
                                   '80/90/95. default: 90)',
                              type=int, choices=[80, 90, 95],
                              default=90)
    parserr_ubsr.add_argument('-rc', '--realign_coverage', metavar='<int>',
                              help='[Only use when -rr on] Min alignment coverage of sequences (default: 80)',
                              type=int,
                              default=80)
    parserr_ubsr.add_argument('-rm', '--realign_minimap2', metavar='<minimap2_path>',
                              help='[Only use when -rr on] Path of minimap2 (default: minimap2 in $PATH)', type=str,
                              default='minimap2')
    parserr_ubsr.add_argument('-rt', '--realign_thread', metavar='<int>',
                              help='[Only use when -rr on] Number of threads using minimap2 (default: 1)', type=int,
                              default=1)
    parserr_ubsr.add_argument('-rd', '--realign_dir', metavar='<str>',
                              help='[Only use when -rr on] Temp directory to realign (default: temp_dir)', type=str,
                              default='temp_dir_[Time]')
    # MergebSeq
    parser_ms = sub_parser.add_parser("mergebseq", help='''Merge the block sequences from many files''')
    parser_ms.add_argument('-o', '--output', metavar='<output.fa>',
                           help='Path of output fasta', type=str, default='merged_unaln.fa')
    parser_ms.add_argument('--allow_samectg', action='store_true',
                           help='Allow the sample contig/scaffold name', default=False)
    parser_ms.add_argument('input_fa', nargs='+',
                           help='Path of input fasta', type=str)

    # RmRedunant
    parser_rr = sub_parser.add_parser("rmredundant", help='''Cluster the sequences and remove redundant sequences''')
    parser_rr.add_argument('-i', '--input_fa', metavar='<input.fa>',
                           help='Path of input fasta', type=str, required=True)
    parser_rr.add_argument('-o', '--output_dir', metavar='<output_dir>',
                           help='Path of output directory', type=str, default='Rmredundant_output')
    parser_rr.add_argument('-c', '--sequence_identity', metavar='<int>',
                           help='Sequence identity threshold (default: 90)', type=int, default=90)
    parser_rr.add_argument('-m', '--method_path', metavar='<cluter_method_path>',
                           help='Path of cluster method, support gclust/cd-hit-est/blastn (default: gclust in $PATH)',
                           type=str, default='gclust')
    parser_rr.add_argument('-t', '--thread', metavar='<int>',
                           help=' Number of threads when clustering (default: 1)', type=int,
                           default=1)
    # rmCtm
    parser_rc = sub_parser.add_parser("rmctm", help='''Detect and discard the contaminated sequences''')
    parser_rc.add_argument('-i', '--input_fa', metavar='<input.fa>',
                           help='Path of input fasta', type=str, required=True)
    parser_rc.add_argument('-o', '--output_dir', metavar='<output_dir>',
                           help='Path of output directory', type=str, default='Rmctm_output')
    parser_rc.add_argument('--rewrite', default=False, action='store_true',
                           help='Rewrite the blastout in output directory')
    parsert_rc = parser_rc.add_argument_group('database parameters')
    parsert_rc.add_argument('-nt', '--nt', metavar='<nt_path>',
                            help='Path of nt (download and decompress from '
                                 'https://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz)',
                            required=True)
    parsert_rc.add_argument('-at', '--at', metavar='<accession2taxid_path>',
                            help='Path of accession2taxid (download and decompress from '
                                 'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)',
                            required=True)
    parsert_rc.add_argument('-rl', '--rl', metavar='<rankedlineage_path>',
                            help='Path of rankedlineage (download and decompress from '
                                 'https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)', required=True)
    parserf_rc = parser_rc.add_argument_group('filter parameters')
    parserf_rc.add_argument('-wl', '--white_list', metavar='<str>',
                            help='legal tax, others will be considered as contamination (default: Viridiplantae)',
                            type=str, default='Viridiplantae')
    parserf_rc.add_argument('-ai', '--alignment_identity', metavar='<int>',
                            help='Alignment identity threshold (default: 90)', type=int, default=90)
    parserf_rc.add_argument('-al', '--alignment_length', metavar='<int>',
                            help='Alignment length threshold (default: 100)', type=int, default=100)
    parserf_rc.add_argument('--force_remove', default=False, action='store_true',
                            help='Force remove all contamination even if same regions of sequences map to white list')
    parserb_rc = parser_rc.add_argument_group('blastn parameters')
    parserb_rc.add_argument('-m', '--blastn_path', metavar='<blastn_path>',
                            help='Path of blastn (default: blastn in $PATH)',
                            type=str, default='blastn')
    parserb_rc.add_argument('-e', '--evalue', metavar='<float>',
                            help='E-value threshold of blastn (default: 1e-5)',
                            type=float, default=1e-5)
    parserb_rc.add_argument('-t', '--thread', metavar='<int>',
                            help=' Number of threads when blastn, recommend less than 4 (default: 1)', type=int,
                            default=1)

    # FastaSta
    parser_fs = sub_parser.add_parser("fastasta", help='''Calculate statistics of fasta file (DNA)''')
    parser_fs.add_argument('-i', '--input_path', metavar='<input.fa>',
                           help='Path of input fasta file', type=str, required=True)
    parser_fs.add_argument('-o', '--output_path', metavar='<output.fa>',
                           help='Path of output fasta file', type=str, required=True)
    parser_fs.add_argument('-l', '--length_filter', metavar='<int>',
                           help='Min length of sequences to consider (default: 500)', type=int, default=500)

    # pTpG
    parser_ptpg = sub_parser.add_parser("ptpg", help='''Extract the longest transcript elements from genes''')
    parser_ptpg.add_argument('-i', '--input', metavar='<input.gff/gtf>',
                             help='Path of input gff/gtf', type=str, required=True)
    parser_ptpg.add_argument('-r', '--region', metavar='<str>',
                             help='CDS or exon (default: CDS)', type=str,
                             choices=['CDS', 'exon'], default='CDS')
    parser_ptpg.add_argument('-o', '--output', metavar='<output.gff/gtf>',
                             help='Path of output gff/gtf', type=str, required=True)

    # GeneCov
    parser_gc = sub_parser.add_parser("genecov", help='''Compute gene coverage and depth''')
    parser_gc.add_argument('-a', '--annotation', metavar='<input.gff/gtf>', help='Path of input gff/gtf', type=str,
                           required=True)
    parser_gc.add_argument('-b', '--bed', metavar='<input.bed>', help='bed of of coverage from bedtools', type=str,
                           required=True)
    parser_gc.add_argument('-r', '--region', metavar='<str>', help='CDS or exon (default: CDS)', type=str,
                           choices=['CDS', 'exon'], default='CDS')
    parser_gc.add_argument('-o', '--output', metavar='<output.cov>', help='Path of output cov', type=str, required=True)
    parser_gc.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
    parser_gc.add_argument('-m', '--min_depth', metavar='<int>', help='Min depth (default: 1)', type=int, default=1)

    # EleCov
    parser_ec = sub_parser.add_parser("elecov", help='''Compute genomic element coverage and depth''')
    parser_ec.add_argument('-a', '--annotation', metavar='<anotation.bed>', help='Path of anotation.bed', type=str,
                           required=True)
    parser_ec.add_argument('-b', '--bed', metavar='<input.bed>', help='bed of of coverage from bedtools', type=str,
                           required=True)
    parser_ec.add_argument('-o', '--output', metavar='<output.cov>', help='Path of output cov', type=str, required=True)
    parser_ec.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
    parser_ec.add_argument('-m', '--min_depth', metavar='<int>', help='Min depth (default: 1)', type=int, default=1)

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        exit()

    args = vars(parser.parse_args())
    command = dir_name + '/' + "lib/{subcmd}.py ".format(subcmd=args["subcmd"]) + ' '.join(sys.argv[2:])
    os.system(command)
