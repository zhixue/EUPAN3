#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/24 10:54 PM
    @Usage: python3 fastasta.py -i input.fa -o out fasta_summary.tsv
"""
import argparse
import sys
from tlog import *
import os
from collections import OrderedDict


def seq_statistics(seq_m):
    all_base = len(seq_m)
    a_base = seq_m.count('A') + seq_m.count('a')
    t_base = seq_m.count('T') + seq_m.count('t')
    c_base = seq_m.count('C') + seq_m.count('c')
    g_base = seq_m.count('G') + seq_m.count('g')
    n_base = seq_m.count('N') + seq_m.count('n')
    nonn_base = a_base + t_base + c_base + g_base + n_base
    other_base = all_base - nonn_base
    if nonn_base == 0:
        gc_content = 0
    else:
        gc_content = (c_base + g_base) / nonn_base
    return all_base, nonn_base, a_base, t_base, c_base, g_base, other_base, round(gc_content, 4)


def process_fasta(fasta, min_len):
    seq_summary = OrderedDict()
    ctg = ''
    seq = ''
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                if seq != '':
                    if len(seq) >= min_len:
                        seq_summary[ctg] = seq_statistics(seq)
                ctg = line.rstrip()[1:].split()[0]
                seq = ''
            else:
                seq += line.rstrip()
    return seq_summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Calculate statistics of fasta file (DNA)''')
    parser.add_argument('-i', '--input_path', metavar='<input.fa>',
                        help='Path of input fasta file', type=str, required=True)
    parser.add_argument('-o', '--output_path', metavar='<output.fa>',
                        help='Path of output fasta file', type=str, required=True)
    parser.add_argument('-l', '--length_filter', metavar='<int>',
                        help='Min length of sequences to consider (default: 500)', type=int, default=500)
    args = vars(parser.parse_args())

    fasta_ctg_summary = process_fasta(args['input_path'], args['length_filter'])
    logging.info("# Summarize the statistics of fasta.")

    # check output dir
    if os.path.isdir(args["output_path"]):
        logging.warning("# {output_path} exists! It will be rewrited!".format(output_path=args["output_path"]))
    with open(args["output_path"], 'w') as fout:
        # command
        fout.write('# ' + ' '.join(sys.argv) + '\n')
        # total summary
        total_n = 0
        total_base = 0
        total_nonnbase = 0
        total_gc_content = 0
        for key in fasta_ctg_summary:
            total_n += 1
            total_base += fasta_ctg_summary[key][0]
            total_nonnbase += fasta_ctg_summary[key][1]
            total_gc_content += fasta_ctg_summary[key][1] * fasta_ctg_summary[key][-1]
        total_gc_content = round(total_gc_content / total_nonnbase, 4)
        summary_log = "# Total {n} sequences, with {base} bases ({nonnbase} bases are not N). " \
                      "Mean length = {lth}, GC content = {gc}".format(n=total_n,
                                                                      base=total_base,
                                                                      nonnbase=total_nonnbase,
                                                                      lth=round(total_base / total_n, 1),
                                                                      gc=total_gc_content)
        logging.info(summary_log)
        summary_header = "#" + '\t'.join(["Seq", "Bases", "NonNBases", "A", "T", "C", "G", "Others", "GCcontent"])
        temp_strings = [summary_log, summary_header]
        for key in fasta_ctg_summary:
            temp_strings += ['\t'.join([key] + [str(x) for x in fasta_ctg_summary[key]])]
        fout.write('\n'.join(temp_strings))
