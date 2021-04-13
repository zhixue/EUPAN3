#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 2:30 PM
    @Usage: python3 eupan3.py command [options]
"""
import argparse
import os
import sys

DESCRIPTION = """
Third Generation Sequencing EUkaryote Pan-genome ANalysis toolkit (TGSEUPAN)
                                                                                                         
"""
VERSION = "Version 0.1.0"

if __name__ == "__main__":
    dir_name, file_name = os.path.split(os.path.abspath(sys.argv[0]))
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument('-v', '--version', action='version', version=VERSION)

    sub_parser = parser.add_subparsers(dest='subcmd', title='sub commands')
    # AssemSta
    parser_as = sub_parser.add_parser("assemsta", help='Use Quast to get assembled contigs/scaffolds to reference to '
                                                       'check the statistics  (The script will call QUAST program, '
                                                       'so the directory where quast.py locates is needed.)')
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
                            help='Min length of sequences to consider (default: None)', type=str, default=None)
    parser_asq.add_argument('-i', '--identity', metavar='<int>',
                            help='Min alignment identity of sequences (choices: 80/90/95. default: 90)', type=int,
                            choices=[80, 90, 95],
                            default=90)
    parser_asq.add_argument('--large_genome', default=False, action='store_true',
                            help='Use optimal parameters for evaluation of large genomes (typically > 100 Mbp)')

    # UnalnBlockSeq
    parser_ubs = sub_parser.add_parser("unalnsseq",
                                       help="Get partial/full/all unaligned block sequences of contigs/scaffolds from "
                                            "Quast")
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
    parser_ms = sub_parser.add_parser("rmredundant", help='''Cluster the sequences and remove redundant sequences''')
    parser_ms.add_argument('-i', '--input_fa', metavar='<input.fa>',
                           help='Path of input fasta', type=str, required=True)
    parser_ms.add_argument('-o', '--output_dir', metavar='<output_dir>',
                           help='Path of output directory', type=str, default='Rmredundant_output')
    parser_ms.add_argument('-c', '--sequence_identity', metavar='<int>',
                           help='Sequence identity threshold (default: 90)', type=int, default=90)
    parser_ms.add_argument('-m', '--method_path', metavar='<cluter_method_path>',
                           help='Path of cluster method, support gclust/cd-hit-est/blastn (default: gclust in $PATH)',
                           type=str, default='gclust')
    parser_ms.add_argument('-t', '--thread', metavar='<int>',
                           help=' Number of threads when clustering (default: 1)', type=int,
                           default=1)
    if len(sys.argv[1:]) == 0:
        parser.print_help()
        exit()

    args = vars(parser.parse_args())
    command = dir_name + '/' + "lib/{subcmd}.py ".format(subcmd=args["subcmd"]) + ' '.join(sys.argv[2:])
    os.system(command)
