#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/6/22 9:40 AM
    @Usage: python3 assemsta.py [options]
"""
import argparse
from tlog import *
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Use Quast to get assembled contigs/scaffolds to reference and to 
    check the statistics (The script will call Quast program, so the directory where quast.py locates is needed.)''')
    parser.add_argument('-a', '--assembly_path', metavar='<assembly.fa>',
                        help='Path of contigs or scaffolds', type=str, required=True)
    parser.add_argument('-r', '--reference_path', metavar='<reference.fa>',
                        help='Path of reference genome', type=str, required=True)
    parser.add_argument('-o', '--output_dir', metavar='<output_directory>',
                        help='Path of output of Quast', type=str, required=True)
    parser.add_argument('-q', '--quast_path', metavar='<quast.py>',
                        help='Path of quast.py (default: quast.py in $PATH)', type=str, default='quast.py')

    parserq = parser.add_argument_group('quast parameters')
    parserq.add_argument('-t', '--thread', metavar='<int>',
                         help='Number of threads using minimap2 in Quast (default: 1)', type=int, default=1)
    parserq.add_argument('-l', '--length_threshold', metavar='<int>',
                         help='Min length of contigs/scaffolds to consider (default: 500, if --large_genome, default: '
                              '3000)',
                         type=int, default=-1)
    parserq.add_argument('-g', '--gff_path', metavar='<reference.gff>',
                         help='Anotations of reference genome (default: None)', type=str, default=None)
    parserq.add_argument('-i', '--identity', metavar='<int>',
                         help='Min alignment identity of sequences (choices: 80/90/95. default: 90)', type=int,
                         choices=[80, 90, 95],
                         default=90)
    parserq.add_argument('--large_genome', default=False, action='store_true',
                         help='Use optimal parameters for evaluation of large genomes (typically > 100 Mbp)')

    args = vars(parser.parse_args())

    # start
    if args["large_genome"]:
        large_par = "--large"
        if args["length_threshold"] == -1:
            length_threshold_par = '-m 3000'
        else:
            length_threshold_par = '-m ' + str(args["length_threshold"])
    else:
        large_par = ''
        if args["length_threshold"] == -1:
            length_threshold_par = '-m 500'
        else:
            length_threshold_par = '-m ' + str(args["length_threshold"])

    if args["gff_path"]:
        gff_par = "-g " + args["gff_path"]
    else:
        gff_par = ''

    command = "{quast} -t {thread} {length_threshold_par} {large_par} -r {ref} {gff_par} -o {output_dir} " \
              "--no-icarus --no-snps " \
              "--min-identity {min_indentity} {assembly}".format(quast=args["quast_path"],
                                                                 thread=args["thread"],
                                                                 length_threshold_par=length_threshold_par,
                                                                 large_par=large_par,
                                                                 ref=args["reference_path"],
                                                                 gff_par=gff_par,
                                                                 output_dir=args["output_dir"],
                                                                 min_indentity=args["identity"],
                                                                 assembly=args["assembly_path"]
                                                                 )
    logging.info("# Start running Quast using this command:\n{cmd}".format(cmd=command))
    os.system(command)
    logging.info("# Finish running Quast.")
