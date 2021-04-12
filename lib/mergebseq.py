#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/12 4:17 PM
    @Usage: python3 mergebseq.py -o merged_unaln.fa sample1.fa sample2.fa sample3.fa ...
"""
import argparse
from tlog import *
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Merge the sequences''')
    parser.add_argument('-o', '--output', metavar='<output.fa>',
                        help='Path of output fasta', type=str, default='merged_unaln.fa')
    parser.add_argument('--allow_samectg', action='store_true',
                        help='Allow the sample contig/scaffold name', default=False)
    parser.add_argument('input_fa', nargs='+',
                        help='Path of input fasta', type=str)
    args = vars(parser.parse_args())
    try:
        if os.path.isdir(args['output']):
            logging.warning("# {output} exists! It will be rewrited!".format(output=args['output']))
        visited_ctg = dict()
        open_in_file_n = 0
        with open(args['output'], 'w') as fout:
            for in_file in args['input_fa']:
                with open(in_file) as fin:
                    open_in_file_n += 1
                    for line in fin:
                        if line.startswith('>'):
                            ctg_tag = line[1:].rstrip()
                            ctg_name = ctg_tag.split()[0]
                            if ctg_name in visited_ctg:
                                if not args['allow_samectg']:
                                    raise ValueError("Same contig/scaffold name exists: %s !" % ctg_name)
                            visited_ctg[ctg_name] = 1
                        fout.write(line)
        logging.info("# Read {file_n} input files, with {ctg_n} block sequences.".format(file_n=open_in_file_n,
                                                                                         ctg_n=len(visited_ctg)))
    except Exception as e:
        logging.error(e)
