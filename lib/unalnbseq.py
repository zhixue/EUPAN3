#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 2:30 PM
    @Usage: python3 unalnbseq.py [options]
"""
import argparse
import re
import os
from tlog import *
import Genome_Interval as gi


def trans_contig_name(string_a, char='_'):
    return re.sub('[^0-9a-zA-Z]', char, string_a)


def read_unaln_table(unalign_table_path, kind_unalign, min_len):
    unalign_inf_dict = {}
    with open(unalign_table_path) as f:
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
                unalign_inf_dict[temp[0].rstrip()] = tuple(
                    (temp[1], temp[2], temp[3], filtered_interval_list))
    return unalign_inf_dict


def write_interval_seq(fasta_file, chrn_intervals, output_fa, sample_tagname, nbase_ignore_precent=100, nbase='N'):
    # init
    in_record_num = 0
    out_record_num = 0
    out_blockseq_num = 0
    seq_name = ''
    seq_name_trans = ''
    seq_seq = ''

    with open(output_fa, 'w') as fout:
        with open(fasta_file) as fin:
            for line in fin:
                if line.startswith('>'):
                    in_record_num += 1
                    # check record
                    if seq_name_trans != '':
                        if seq_name_trans in chrn_intervals:
                            out_record_num += 1
                            for interval in chrn_intervals[seq_name_trans][3].intervals:
                                block_seq = interval.getsubseq(seq_seq)
                                if block_seq.count(nbase) / len(block_seq) * 100 > nbase_ignore_precent:
                                    continue
                                out_blockseq_num += 1
                                bseq_anno = '>' + sample_tagname + '_' + 'bseq' + str(out_blockseq_num) + \
                                    ' ' + seq_anno + \
                                    ', ' + 'contiglength_' + chrn_intervals[seq_name_trans][0] + \
                                    ':' + str(interval.lower_bound) + '-' + str(interval.upper_bound) + \
                                    ', ' + 'bseqlength_' + str(interval.length) + \
                                    ', ' + chrn_intervals[seq_name_trans][2]
                                fout.write(bseq_anno + '\n')
                                fout.write(block_seq + '\n')

                    seq_anno = line.rstrip()[1:]
                    seq_name = seq_anno.split()[0]
                    seq_name_trans = trans_contig_name(seq_name)
                    seq_seq = ''
                else:
                    seq_seq += line.rstrip()
            # final one
            if seq_name_trans != '':
                if seq_name_trans in chrn_intervals:
                    out_record_num += 1
                    for interval in chrn_intervals[seq_name_trans][3].intervals:
                        block_seq = interval.getsubseq(seq_seq)
                        if block_seq.count('N') / len(block_seq) * 100 > nbase_ignore_precent:
                            continue
                        out_blockseq_num += 1
                        bseq_anno = '>' + sample_tagname + '_' + 'bseq' + str(out_blockseq_num) + \
                                    ' ' + seq_anno + \
                                    ', ' + 'contiglength_' + chrn_intervals[seq_name_trans][0] + \
                                    ':' + str(interval.lower_bound) + '-' + str(interval.upper_bound) + \
                                    ', ' + 'bseqlength_' + str(interval.length) + \
                                    ', ' + chrn_intervals[seq_name_trans][2]
                        fout.write(bseq_anno + '\n')
                        fout.write(block_seq + '\n')
    return in_record_num, out_record_num, out_blockseq_num


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
        description='''Get partial/full/all unaligned block sequences of contigs/scaffolds from Quast''')
    parser.add_argument('-a', '--assembly_path', metavar='<assembly.fa>',
                        help='Path of contigs or scaffolds', type=str, required=True)
    parser.add_argument('-u', '--unaln_path', metavar='<unaligned.info>',
                        help='Path of unaligned table (at quast_output/contigs_reports/contigs_report_xxx.unaligned'
                             '.info)',
                        type=str, required=True)
    parser.add_argument('-o', '--output_path', metavar='<output.fa>',
                        help='Path of output unaln sequences', type=str, required=True)
    parser.add_argument('-l', '--length_filter', metavar='<int>',
                        help='Min length of block sequences to consider (default: 500)', type=int, default=500)
    parser.add_argument('-k', '--kind_unaln', metavar='<str>',
                        help='Use full/partial/all unaligned sequences (choices: full/partial/all. default: all)',
                        choices=['full', 'partial', 'all'], type=str, default='all')
    parser.add_argument('-s', '--sample_tag', metavar='<str>',
                        help='Add sample tag before each contig/scaffold \
                        (default: None, e.g. -s Sample1 , "_" will be used, Chr1 -> Sample1_Chr1)', type=str,
                        default='')
    parser.add_argument('-n', '--nbase_ignore', metavar='<int>',
                        help='Max percentage of N bases in block sequences to ignore (default: 100)',
                        type=int, default=100)

    parserr = parser.add_argument_group('realign parameters')
    parserr.add_argument('-rr', '--realign_reference', metavar='<reference.fa>',
                         help='Using minimap2 to realign  to reference genome/mitochondrion/plastid and drop high '
                              'similar unaligned block sequences',
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
                         help='[Only use when -rr on] Temp directory to realign (default: unalnbseq_temp_[Time])',
                         type=str, default='unalnbseq_temp')
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
    input_record_num, output_record_num, output_blockseq_num = write_interval_seq(assembly_file, unalign_dict,
                                                                                  output_fasta, sample_tag,
                                                                                  nbase_ignore)
    logging.info(
        "# Load {inseq_n} chromosomes/contigs/scaffolds from fasta, write {outseq_n} sequences with {blockoutseq_n} "
        "block "
        "sequences.".format(
            inseq_n=input_record_num,
            outseq_n=output_record_num,
            blockoutseq_n=output_blockseq_num
        ))
    if realign_step:
        # check temp dir
        temp_dir = args["realign_dir"]  # + '_' + time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
        if os.path.isdir(args["realign_dir"]):
            logging.warning("# {temp_path} exists! It will be rewrited!".format(temp_path=temp_dir))
        else:
            os.mkdir(args["realign_dir"])
            logging.info("# Realign path is {temp_path}.".format(temp_path=temp_dir))

        # map with minimap2
        logging.info("# Start realigning block sequences to references.")
        if args["realign_identity"] == 90:
            minimap2_idt_par = "-x asm10"
        elif args["realign_identity"] == 95:
            minimap2_idt_par = "-x asm5"
        elif args["realign_identity"] == 80:
            minimap2_idt_par = "-x asm20"
        else:
            minimap2_idt_par = ""
        temp_paf = "bseq2ref.paf"
        command = '{minimap2} -t {thread} {idt_par} {ref} {query} > {temp_dir}/{temp_paf}'.format(
            minimap2=args["realign_minimap2"],
            thread=args["realign_thread"],
            idt_par=minimap2_idt_par,
            ref=args["realign_reference"],
            query=output_fasta,
            temp_dir=temp_dir,
            temp_paf=temp_paf
        )
        os.system(command)
        logging.info("# Finish realigning block sequences to references.")
        # filter
        temp_out_seq_list_path = "mapped_bseq.txt"
        temp_out_seq_filtered_path = "remain_bseq.fa"
        drop_n = drop_seq_from_paf(temp_dir + '/' + temp_paf, args["realign_coverage"],
                                   temp_dir + '/' + temp_out_seq_list_path)
        fsr_return_status = gi.fa_some_record(output_fasta, temp_dir + '/' + temp_out_seq_list_path,
                                              temp_dir + '/' + temp_out_seq_filtered_path, exclude=True)
        logging.info("# Remove {n} block sequences similiar to references.".format(n=drop_n))
        # mv and rm temp dir
        logging.info("# Update block sequences and remove temp files.")
        os.system('mv {temp_dir}/{temp_fa} {output_fa}'.format(temp_dir=temp_dir,
                                                               temp_fa=temp_out_seq_filtered_path,
                                                               output_fa=output_fasta))
        # os.system('rm -rf {temp_dir}'.format(temp_dir=temp_dir))
        logging.info("# Realign step is finished.")
