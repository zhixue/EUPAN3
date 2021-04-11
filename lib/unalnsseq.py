#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 2:30 PM
    @Usage:
"""
import argparse
import re
import os
from tlog import *
import Genome_Interval as gi


def trans_contig_name(string_a, char='_'):
    return re.sub('[^0-9a-zA-Z]+', char, string_a)


def read_unaln_table(unalign_table, kind_unalign, min_len):
    unalign_dict = {}
    with open(unalign_table) as f:
        # contig, total_length, unaligned_length, type, unaligned_parts
        i = 0
        for line in f:
            i += 1
            if i == 1:
                continue
            temp = line.rstrip().split('\t')
            unaln_type = temp[3]
            interval_list = gi.GIntervalList(temp[4])
            filtered_interval_list = gi.filter_length(interval_list, min_len=min_len)
            if filtered_interval_list.isempty():
                continue
            if unaln_type in kind_unalign:
                unalign_dict[trans_contig_name(temp[0].rstrip())] = tuple(
                    (temp[1], temp[2], temp[3], filtered_interval_list))
    return unalign_dict


def write_interval_seq(fasta_file, chrn_intervals, output_fasta, sample_tag, nbase_ignore, nbase='N'):
    # init
    input_record_num = 0
    output_record_num = 0
    output_subseq_num = 0
    seq_name = ''
    seq_seq = ''

    with open(output_fasta, 'w') as fout:
        with open(fasta_file) as fin:
            for line in fin:
                if line.startswith('>'):
                    input_record_num += 1
                    # check record
                    if seq_name != '':
                        if seq_name in chrn_intervals:
                            output_record_num += 1
                            for interval in chrn_intervals[seq_name][3].intervals:
                                sub_seq = interval.getsubseq(seq_seq)
                                if sub_seq.count(nbase)/len(sub_seq) * 100 > nbase_ignore:
                                    continue
                                fout.write('>' + sample_tag + seq_name + ':' + str(
                                    interval.lower_bound) + '-' + str(interval.upper_bound) + ' ' +
                                           chrn_intervals[seq_name][2] + '\n')
                                fout.write(sub_seq + '\n')
                                output_subseq_num += 1
                    seq_name = trans_contig_name(line.rstrip()[1:])
                    seq_seq = ''
                else:
                    seq_seq += line.rstrip()
            # final one
            if seq_name != '':
                if seq_name in chrn_intervals:
                    output_record_num += 1
                    for interval in chrn_intervals[seq_name][3].intervals:
                        sub_seq = interval.getsubseq(seq_seq)
                        if sub_seq.count('N') / len(sub_seq) * 100 > nbase_ignore:
                            continue
                        fout.write('>' + sample_tag + seq_name + ':' + str(
                            interval.lower_bound) + '-' + str(interval.upper_bound) + ' ' + chrn_intervals[seq_name][
                                       2] + '\n')
                        fout.write(sub_seq + '\n')
                        output_subseq_num += 1
    return input_record_num, output_record_num, output_subseq_num


def drop_seq_from_paf(paf_path, least_coverage, out_seq_list_path):
    # paf format
    # 1 string, query name
    # 2 int, query length
    # 3 int, query start (0 base,closed)
    # 4 int, query end (0 base,open)
    # 5 char, strand (+/-)
    # 6 string, target name
    # 7 int, target length
    # 8 int, target start
    # 9 int, target end
    # 10 int, number of residue matched
    # 11 int, alignment block length
    # 12 int, mapping quality (0~255,255:missing)
    with open(paf_path) as fin:
        with open(out_seq_list_path, 'w') as fout:
            output_set = set()
            for line in fin:
                temp = line.rstrip().split('\t')
                if float(temp[9]) / float(temp[1]) * 100 >= least_coverage:
                    output_set.add(temp[0])
            fout.write('\n'.join(output_set))
    return len(output_set)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='''Get partial/full/all unaligned sub sequences of contigs/scaffolds from Quast''')
    parser.add_argument('-a', '--assembly_path', metavar='<assembly.fa>',
                        help='Path of contigs or scaffolds', type=str, required=True)
    parser.add_argument('-u', '--unaln_path', metavar='<unaligned.info>',
                        help='Path of unaligned table (at quast_output/contigs_reports/contigs_report_xxx.unaligned'
                             '.info)',
                        type=str, required=True)
    parser.add_argument('-o', '--output_path', metavar='<output.fa>',
                        help='Path of output unaln sequences', type=str, required=True)
    parser.add_argument('-l', '--length_filter', metavar='<int>',
                        help='Min length of sub sequences to consider (default: 500)', type=int, default=500)
    parser.add_argument('-k', '--kind_unaln', metavar='<str>',
                        help='Use full/partial/all unaligned sequences (choices: full/partial/all. default: all)',
                        choices=['full', 'partial', 'all'], type=str, default='all')
    parser.add_argument('-s', '--sample_tag', metavar='<str>',
                        help='Add sample tag before each contig/scaffold \
                        (default: None, e.g. -s Sample1 , "_" will be used, Chr1 -> Sample1_Chr1)', type=str,
                        default='')
    parser.add_argument('-n', '--nbase_ignore', metavar='<int>',
                        help='Max percentage of N bases in sub sequences to ignore (default: 100)',
                        type=int, default=100)

    parserr = parser.add_argument_group('realign parameters')
    parserr.add_argument('-rr', '--realign_reference', metavar='<reference.fa>',
                         help='Using minimap2 to realign  to reference genome/mitochondrion/plastid and drop high '
                              'similar unaligned sub sequences',
                         type=str, default='')
    parserr.add_argument('-ri', '--realign_identity', metavar='<int>',
                         help='[Only use when -rr on] Min alignment identity of sequences in realign step (choices: '
                              '80/90/95. default: 90)',
                         type=int, choices=[80, 90, 95],
                         default=90)
    parserr.add_argument('-rc', '--realign_coverage', metavar='<int>',
                         help='[Only use when -rr on] Min alignment coverage of sequences (default: 80)',
                         type=int,
                         default=80)
    parserr.add_argument('-rm', '--realign_minimap2', metavar='<minimap2_path>',
                         help='[Only use when -rr on] Path of minimap2 (default: minimap2 in $PATH)', type=str,
                         default='minimap2')
    parserr.add_argument('-rt', '--realign_thread', metavar='<int>',
                         help='[Only use when -rr on] Number of threads using minimap2 (default: 1)', type=int,
                         default=1)
    parserr.add_argument('-rd', '--realign_dir', metavar='<str>',
                         help='[Only use when -rr on] Temp directory to realign (default: temp_dir)', type=str,
                         default='temp_dir')
    args = vars(parser.parse_args())
    assembly_file = args['assembly_path']
    unalign_table = args['unaln_path']
    output_fasta = args['output_path']
    length_cutoff = args['length_filter']
    kind_use = args['kind_unaln']
    sample_tag = args['sample_tag']
    nbase_ignore = args['nbase_ignore']
    if kind_use == 'all':
        kind_use = ['full', 'partial']
    else:
        kind_use = [kind_use]
    if sample_tag:
        sample_tag = sample_tag + '_'

    if args["realign_reference"]:
        realign_step = 1
        logging.info("# Realign step is ON.")
    else:
        realign_step = 0
        logging.info("# Realign step is OFF.")

    # start
    unalign_dict = read_unaln_table(unalign_table, kind_use, length_cutoff)
    logging.info(
        "# Load {record_n} chromosomes/contigs/scaffolds from unalign.info".format(record_n=len(unalign_dict.keys())))
    input_record_num, output_record_num, output_subseq_num = write_interval_seq(assembly_file, unalign_dict,
                                                                                output_fasta, sample_tag, nbase_ignore)
    logging.info(
        "# Load {inseq_n} chromosomes/contigs/scaffolds from fasta, write {outseq_n} sequences with {suboutseq_n} sub "
        "sequences.".format(
            inseq_n=input_record_num,
            outseq_n=output_record_num,
            suboutseq_n=output_subseq_num
        ))
    if realign_step:
        # check temp dir
        if os.path.isdir(args["realign_dir"]):
            logging.warning("# {temp_path} exists! It will be rewrited!".format(temp_path=args["realign_dir"]))
        else:
            os.mkdir(args["realign_dir"])
            logging.info("# Realign path is {temp_path}.".format(temp_path=args["realign_dir"]))

        # map with minimap2
        logging.info("# Start realigning sub sequences to references.")
        if args["realign_identity"] == 90:
            minimap2_idt_par = "-x asm10"
        elif args["realign_identity"] == 95:
            minimap2_idt_par = "-x asm5"
        elif args["realign_identity"] == 80:
            minimap2_idt_par = "-x asm20"
        else:
            minimap2_idt_par = ""
        temp_paf = "sseqmap2ref.paf"
        command = '{minimap2} -t {thread} {idt_par} {ref} {query} > {temp_dir}/{temp_paf}'.format(
            minimap2=args["realign_minimap2"],
            thread=args["realign_thread"],
            idt_par=minimap2_idt_par,
            ref=args["realign_reference"],
            query=output_fasta,
            temp_dir=args["realign_dir"],
            temp_paf=temp_paf
        )
        os.system(command)
        logging.info("# Finish realigning sub sequences to references.")
        # filter
        temp_out_seq_list_path = "drop_sseq.txt"
        temp_out_seq_filtered_path = "remain_sseq.fa"
        drop_n = drop_seq_from_paf(args["realign_dir"] + '/' + temp_paf, args["realign_coverage"],
                                   args["realign_dir"] + '/' + temp_out_seq_list_path)
        fsr_return_status = gi.fa_some_record(output_fasta, args["realign_dir"] + '/' + temp_out_seq_list_path,
                                              args["realign_dir"] + '/' + temp_out_seq_filtered_path, exclude=True)
        logging.info("# Remove {n} sub sequences similiar to references.".format(n=drop_n))
        # mv and rm temp dir
        logging.info("# Update sub sequences and remove temp files.")
        os.system('mv {temp_dir}/{temp_fa} {output_fa}'.format(temp_dir=args["realign_dir"],
                                                               temp_fa=temp_out_seq_filtered_path,
                                                               output_fa=output_fasta))
        os.system('rm -rf {temp_dir}'.format(temp_dir=args["realign_dir"]))
        logging.info("# Realign step is finished.")
