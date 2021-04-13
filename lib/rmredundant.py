#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/12 3:28 PM
    @Usage: python3 rmredundant.py [options]
"""
import argparse
from tlog import *
import os
import sys
import Genome_Interval as gi

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Cluster the sequences and remove redundant sequences''')
    parser.add_argument('-i', '--input_fa', metavar='<input.fa>',
                        help='Path of input fasta', type=str, required=True)
    parser.add_argument('-o', '--output_dir', metavar='<output_dir>',
                        help='Path of output directory', type=str, default='Rmredundant_output')
    parser.add_argument('-c', '--sequence_identity', metavar='<int>',
                        help='Sequence identity threshold (default: 90)', type=int, default=90)
    parser.add_argument('-m', '--method_path', metavar='<cluter_method_path>',
                        help='Path of cluster method, support gclust/cd-hit-est/blastn (default: gclust in $PATH)',
                        type=str, default='gclust')
    parser.add_argument('-t', '--thread', metavar='<int>',
                        help=' Number of threads when clustering (default: 1)', type=int,
                        default=1)

    args = vars(parser.parse_args())

    if os.path.isdir(args["output_dir"]):
        logging.warning("# {dir_path} exists! It will be rewrited!".format(dir_path=args["output_dir"]))
    else:
        os.mkdir(args["output_dir"])
        logging.info("# Output directory path is {dir_path}".format(dir_path=args["output_dir"]))

    clust_fa = args["output_dir"] + '/' + "unique_bseq.fa"
    clust_out = args["output_dir"] + '/' + "unique_bseq.fa.clstr"

    if args['method_path'].split('/')[-1] == "gclust":
        gclust_dir = args['method_path'].rstrip('gclust').rstrip('/')
        # sort
        gclust_sort = gclust_dir + '/' + "script/sortgenome.pl"
        sorted_file = args["output_dir"] + '/' + 'sorted_temp.fa'
        command = "{gsort} " \
                  "--genomes-file {rawfile} " \
                  "--sortedgenomes-file {sortfile}".format(gsort=gclust_sort,
                                                           rawfile=args["input_fa"],
                                                           sortfile=sorted_file)
        logging.info("# Sort block sequences with command: {command}".format(command=command))
        os.system(command)
        # cluster
        gclust_log = args["output_dir"] + '/' + 'gclust.log'
        command = "{gclust} " \
                  "-loadall " \
                  "-minlen 20 " \
                  "-both -nuc " \
                  "-threads {thread} " \
                  "-ext 1 " \
                  "-sparse 2 " \
                  "-memiden {idt} " \
                  "{sortfile} > {gclust_out} 2>{gclust_log}".format(gclust=args["method_path"],
                                                                    thread=args["thread"],
                                                                    idt=args["sequence_identity"],
                                                                    sortfile=sorted_file,
                                                                    gclust_out=clust_out,
                                                                    gclust_log=gclust_log)
        logging.info("# Cluster block sequences with command: {command}".format(command=command))
        os.system(command)
        # gclust2fa
        clust_n = gi.gclust2fa(sorted_file, clust_out, clust_fa)
        logging.info("# Total {n} cluster block sequences are got".format(n=clust_n))
    elif args['method_path'].split('/')[-1] == "cd-hit-est":
        command = "{cdhitest} " \
                  "-M 0 " \
                  "-T {thread} " \
                  "-i {in_fa} " \
                  "-o {clust_fa} " \
                  "-c {idt}".format(cdhitest=args["method_path"],
                                    thread=args["thread"],
                                    in_fa=args["input_fa"],
                                    clust_fa=clust_fa,
                                    idt=float(args["sequence_identity"])/100)
        logging.info("# Cluster block sequences with command: {command}".format(command=command))
        os.system(command)
        logging.info("# Finish clustering block.")
    elif args['method_path'].split('/')[-1] == "blastn":
        blastn_dir = args['method_path'].rstrip('blastn').rstrip('/')
        makedb = blastn_dir + '/' + "makeblastdb"
        dbidx = args["output_dir"] + '/' + "temp_blast_idx"
        # make blast index
        command = "{makeblastdb} " \
                  "-dbtype nucl " \
                  "-in {in_fa} " \
                  "-out {dbidx}".format(makeblastdb=makedb,
                                        in_fa=args["input_fa"],
                                        dbidx=dbidx)
        os.system(command)
        logging.info("# Finish building index.")
        # blastn
        blast_out = args["output_dir"] + '/' + "temp_blastout.txt"
        command = "{blastn} " \
                  "-query {in_fa} " \
                  "-db {dbidx} " \
                  "-out {blastout} " \
                  "-evalue 1e-5 " \
                  "-outfmt '6 qseqid sseqid qlen slen length qstart qend sstart send pident evalue' " \
                  "-max_target_seqs 1000 -num_threads {thread}".format(blastn=args["method_path"],
                                                                       in_fa=args["input_fa"],
                                                                       dbidx=dbidx,
                                                                       blastout=blast_out,
                                                                       thread=args["thread"])
        logging.info("# Blast block sequences with command: {command}".format(command=command))
        os.system(command)
        logging.info("# Finish blast.")
        # run blastcluster
        dir_name, file_name = os.path.split(os.path.abspath(sys.argv[0]))
        command = "perl " + dir_name + '/' + "lib/blastCluster.pl " + ' '.join([args["input_fa"],
                                                                                str(args["sequence_identity"]),
                                                                                blast_out,
                                                                                clust_fa])
        logging.info("# Finish clustering.")
