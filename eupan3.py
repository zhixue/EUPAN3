#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2022/1/17 14:30 PM
    @Usage: python3 eupan3.py subcommand [options]
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
    parser_as = sub_parser.add_parser("assemsta", help='''Use Quast to get assemblies statistics''', prefix_chars=':')
    parser_as.add_argument('integers', metavar='options',nargs='*') 

    # UnalnBlockSeq
    parser_ubs = sub_parser.add_parser("unalnbseq", help='''Get partial/full unaligned block sequences''', prefix_chars=':')
    parser_ubs.add_argument('integers', metavar='options',nargs='*') 

    # MergebSeq
    parser_ms = sub_parser.add_parser("mergebseq", help='''Merge the block sequences from many files''', prefix_chars=':')
    parser_ms.add_argument('integers', metavar='options',nargs='*')

    # RmRedunant
    parser_rr = sub_parser.add_parser("rmredundant", help='''Cluster the sequences and remove redundant sequences''', prefix_chars=':')
    parser_rr.add_argument('integers', metavar='options',nargs='*')

    # rmCtm
    parser_rc = sub_parser.add_parser("rmctm", help='''Detect and discard the contaminated sequences''', prefix_chars=':')
    parser_rc.add_argument('integers', metavar='options',nargs='*')

    # FastaSta
    parser_fs = sub_parser.add_parser("fastasta", help='''Calculate statistics of fasta file (DNA)''', prefix_chars=':')
    parser_fs.add_argument('integers', metavar='options',nargs='*')

    # pTpG
    parser_ptpg = sub_parser.add_parser("ptpg", help='''Extract the longest transcript elements from genes''', prefix_chars=':')
    parser_ptpg.add_argument('integers', metavar='options',nargs='*')

    # GeneCov
    parser_gc = sub_parser.add_parser("genecov", help='''Compute gene coverage and depth''', prefix_chars=':')
    parser_gc.add_argument('integers', metavar='options',nargs='*')

    # EleCov
    parser_ec = sub_parser.add_parser("elecov", help='''Compute genomic element coverage and depth''', prefix_chars=':')
    parser_ec.add_argument('integers', metavar='options',nargs='*')

    if len(sys.argv[1:]) == 0:
        parser.print_help()
        exit()

    args = vars(parser.parse_args())
    command = sys.executable + ' ' + dir_name + '/' + "lib/{subcmd}.py ".format(subcmd=args["subcmd"]) + ' '.join(sys.argv[2:])
    os.system(command)
