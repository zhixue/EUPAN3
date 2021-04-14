#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/14 10:36 AM
    @Usage: python3 rmCtm.py [options]
"""
import argparse
from tlog import *
import os
from getTax import gettax, detectcontamination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Detect and discard the potentail contaminated sequences''')
    parser.add_argument('-i', '--input_fa', metavar='<input.fa>',
                        help='Path of input fasta', type=str, required=True)
    parser.add_argument('-o', '--output_dir', metavar='<output_dir>',
                        help='Path of output directory', type=str, default='Rmctm_output')
    parser.add_argument('--rewrite', default=False, action='store_true',
                        help='Rewrite the blastout in output directory')

    parsert = parser.add_argument_group('database parameters')
    parsert.add_argument('-nt', '--nt', metavar='<nt_path>',
                         help='Path of nt (download and decompress from https://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz)',
                         required=True)
    parsert.add_argument('-at', '--at', metavar='<accession2taxid_path>',
                         help='Path of accession2taxid (download and decompress from '
                              'https://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)',
                         required=True)
    parsert.add_argument('-rl', '--rl', metavar='<rankedlineage_path>',
                         help='Path of rankedlineage (download and decompress from '
                              'https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz)', required=True)

    parserf = parser.add_argument_group('filter parameters')
    parserf.add_argument('-wl', '--white_list', metavar='<str>',
                         help='legal tax, others will be considered as contamination (default: Viridiplantae)',
                         type=str, default='Viridiplantae')
    parserf.add_argument('-ai', '--alignment_identity', metavar='<int>',
                         help='Alignment identity threshold (default: 90)', type=int, default=90)
    parserf.add_argument('-al', '--alignment_length', metavar='<int>',
                         help='Alignment length threshold (default: 100)', type=int, default=100)

    parserb = parser.add_argument_group('blastn parameters')
    parserb.add_argument('-m', '--blastn_path', metavar='<blastn_path>',
                         help='Path of blastn (default: blastn in $PATH)',
                         type=str, default='blastn')
    parserb.add_argument('-e', '--evalue', metavar='<float>',
                         help='E-value threshold of blastn (default: 1e-5)',
                         type=float, default=1e-5)
    parserb.add_argument('-t', '--thread', metavar='<int>',
                         help=' Number of threads when blastn, recommend less than 4 (default: 1)', type=int,
                         default=1)

    args = vars(parser.parse_args())

    # check output dir
    if os.path.isdir(args["output_dir"]):
        logging.warning("# {output_path} exists! It will be rewrited!".format(output_path=args["output_dir"]))
    else:
        os.mkdir(args["output_dir"])
        logging.info("# Output path is {output_path}.".format(output_path=args["output_dir"]))

    # blast
    blast_out = args["output_dir"] + '/' + "uniquebseq2nt_blastout.txt"
    command = "{blastn} " \
              "-query {in_fa} " \
              "-db {dbidx} " \
              "-out {blastout} " \
              "-evalue 1e-5 " \
              "-outfmt '6 qacc sacc qlen slen length qstart qend sstart send pident evalue' " \
              "-num_threads {thread}".format(blastn=args["blastn_path"],
                                             in_fa=args["input_fa"],
                                             dbidx=args["nt"],
                                             blastout=blast_out,
                                             thread=args["thread"])

    if os.path.isfile(blast_out) and (not args["rewrite"]):
        logging.info("# Skip blast block sequences to nt.")
    else:
        logging.info("# Blast block sequences to nt.")
        os.system(command)

    # get tax from blast
    acc2tax_out = args["output_dir"] + '/' + "uniquebseq2nt_taxout.txt"
    load_acc_n, load_tax_n, successacc2taxid_n = gettax(blast_out,
                                                        args["at"],
                                                        args["rl"],
                                                        acc2tax_out)
    logging.info("# Read {accn} accessions, "
                 "{taxn} taxes, "
                 "{acc2taxn} accessions successfully transformed to taxes".format(accn=load_acc_n,
                                                                                  taxn=load_tax_n,
                                                                                  acc2taxn=successacc2taxid_n))
    # filter tax
    filtered_seq_path = args["output_dir"] + '/' + "nocontamination.fa"
    droped_seq_path = args["output_dir"] + '/' + "contamination.fa"
    droped_seq_inf_path = args["output_dir"] + '/' + "contamination_information.txt"
    detectcontamination(blast_out,
                        acc2tax_out,
                        args["alignment_identity"],
                        args["alignment_length"],
                        args["white_list"],
                        args["input_fa"],
                        filtered_seq_path,
                        droped_seq_path,
                        droped_seq_inf_path)
    logging.info("# Finish removing contamination.")
