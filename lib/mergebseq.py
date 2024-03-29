#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/12 4:17 PM
    @Usage: python3 mergebseq.py -o merged_unaln.fa sample1.fa sample2.fa sample3.fa ...
"""
import argparse
from tlog import *
import os
import glob

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Merge the block sequences from many files''')
    parser.add_argument('-o', '--output_fa', metavar='<output.fa>',
                        help='Path of output fasta', type=str, default='merged_unaln.fa')
    parser.add_argument('--ignore_samename', action='store_true',
                        help='Ignore the sample sequence name', default=False)
    parser.add_argument('input_fa', nargs='+',
                        help='Path of inputs with fasta', type=str)
    args = vars(parser.parse_args())
    try:
        if os.path.isfile(args['output_fa']):
            logging.warning("# {output} exists! It will be rewrited!".format(output=args['output_fa']))
        visited_ctg = dict()
        in_files = args['input_fa']
        if len(in_files) == 1:
            in_files = glob.glob(in_files[0])
        open_in_file_n = 0
        with open(args['output_fa'], 'w') as fout:
            for in_file in in_files:
                with open(in_file) as fin:
                    ctg_n_in_this_file = 0
                    write_flag = 1
                    open_in_file_n += 1
                    for line in fin:
                        if line.startswith('>'):
                            write_flag = 1
                            ctg_tag = line[1:].rstrip()
                            ctg_name = ctg_tag.split()[0]
                            if ctg_name in visited_ctg:
                                if not args['ignore_samename']:
                                    raise ValueError("# Same contig/scaffold name exists: %s !" % ctg_name)
                                write_flag = 0
                            visited_ctg[ctg_name] = 1
                            ctg_n_in_this_file += 1
                        if write_flag:
                            fout.write(line)
                    logging.info("# Load {file} with {n} block sequences.".format(file=in_file, n=ctg_n_in_this_file))
        logging.info("# Total load {file_n} input files, with {ctg_n} block sequences.".format(file_n=open_in_file_n,
                                                                                               ctg_n=len(visited_ctg)))
    except Exception as e:
        logging.error(e)
