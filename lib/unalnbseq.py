#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2023/2/17 10:38 AM
    @Usage: python3 unalnbseq.py [options]
"""
import argparse
import re
import os
from tlog import *
import Genome_Interval as gi
import time


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
                unalign_inf_dict[trans_contig_name(temp[0].rstrip())] = tuple(
                    (temp[1], temp[2], temp[3], filtered_interval_list))
    return unalign_inf_dict


def write_interval_seq(fasta_file, chrn_intervals, output_fa,
                       unalign_path, sample_tagname, nbase_ignore_precent=100, nbase='N'):
    # init
    in_record_num = 0
    out_record_num = 0
    out_blockseq_num = 0
    seq_name = ''
    previous_seq_name = ''
    seq_name_trans = ''
    seq_seq = ''
    bseq_prefix = sample_tagname + '_' + 'bs'
    gff_string = ''

    with open(output_fa, 'w') as fout:
        with open(fasta_file) as fin:
            for line in fin:
                if line.startswith('>'):
                    in_record_num += 1
                    # check record
                    if seq_name != '':
                        if seq_name_trans in chrn_intervals:
                            out_record_num += 1
                            for interval in chrn_intervals[seq_name_trans][3].intervals:
                                block_seq = interval.getsubseq(seq_seq)
                                if block_seq.count(nbase) / len(block_seq) * 100 >= nbase_ignore_precent:
                                    continue
                                out_blockseq_num += 1
                                bseq_anno = '>' + bseq_prefix + str(out_blockseq_num) + \
                                            ' ' + seq_anno + \
                                            ', ' + 'contiglength_' + chrn_intervals[seq_name_trans][0] + \
                                            ':' + str(interval.lower_bound) + '-' + str(interval.upper_bound) + \
                                            ', ' + 'bseqlength_' + str(interval.length) + \
                                            ', ' + chrn_intervals[seq_name_trans][2]
                                gff_cols = [seq_name, "eupan3", "unalignblock", str(interval.lower_bound),
                                            str(interval.upper_bound), '.', '.', '.',
                                            gi.dict2string({"ID": bseq_prefix + str(out_blockseq_num),
                                                            "length": str(interval.length)})]
                                if seq_name != previous_seq_name:
                                    gff_string += '##sequence-region {chrn} 1 {end}\n'.format(chrn=seq_name,
                                                                                              end=chrn_intervals[
                                                                                                  seq_name_trans][0])
                                    previous_seq_name = seq_name
                                gff_string += '\t'.join(gff_cols) + '\n'
                                fout.write(bseq_anno + '\n')
                                fout.write(block_seq + '\n')

                    seq_anno = line.rstrip()[1:]
                    seq_name = seq_anno.split()[0]
                    seq_name_trans = trans_contig_name(seq_name)
                    seq_seq = ''
                else:
                    seq_seq += line.rstrip()
            # final one
            if seq_name != '':
                if seq_name_trans in chrn_intervals:
                    out_record_num += 1
                    for interval in chrn_intervals[seq_name_trans][3].intervals:
                        block_seq = interval.getsubseq(seq_seq)
                        if block_seq.count('N') / len(block_seq) * 100 >= nbase_ignore_precent:
                            continue
                        out_blockseq_num += 1
                        bseq_anno = '>' + bseq_prefix + str(out_blockseq_num) + \
                                    ' ' + seq_anno + \
                                    ', ' + 'contiglength_' + chrn_intervals[seq_name_trans][0] + \
                                    ':' + str(interval.lower_bound) + '-' + str(interval.upper_bound) + \
                                    ', ' + 'bseqlength_' + str(interval.length) + \
                                    ', ' + chrn_intervals[seq_name_trans][2]
                        gff_cols = [seq_name, "eupan3", "unalignblock", str(interval.lower_bound),
                                    str(interval.upper_bound), '.', '.', '.',
                                    gi.dict2string({"ID": bseq_prefix + str(out_blockseq_num),
                                                    "length": str(interval.length)})]
                        gff_string += '\t'.join(gff_cols) + '\n'
                        fout.write(bseq_anno + '\n')
                        fout.write(block_seq + '\n')
    # write gff
    with open(output_fa + '.gff3', 'w') as fout:
        fout.write('##gff-version 3\n')
        fout.write("# Generate from {file1} and {file2}.\n".format(
            file1=fasta_file,
            file2=unalign_path
        ) + gff_string)
    return in_record_num, out_record_num, out_blockseq_num


def drop_seq_from_paf(paf_path, least_idt, least_coverage, filtered_length,
                      out_seq_list_path, out_seq_stat_path, segment_mode=1):
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

    # init
    seq_cov_dict = dict()
    output_set = set()

    if segment_mode != 0:
        query_length = dict()
        query_map_region = dict()
        query_cov_region = dict()
        query_uncov_region = dict()
        with open(paf_path) as fin:
            for line in fin:
                temp = line.rstrip().split('\t')
                # like paf
                if float(temp[9]) / float(temp[10]) * 100 < least_idt:
                    continue
                if temp[0] not in query_length:
                    query_length[temp[0]] = int(temp[1])
                    # query end open -> query end closed
                if temp[0] not in query_map_region:
                    query_map_region[temp[0]] = [(int(temp[2]), int(temp[3]) - 1)]
                else:
                    if (int(temp[2]), int(temp[3]) - 1) not in query_map_region[temp[0]]:
                        query_map_region[temp[0]] += [(int(temp[2]), int(temp[3]) - 1)]

        # compute coverage and write
        query_map_length = dict()
        fout1 = open(out_seq_list_path, 'w')
        fout2 = open(out_seq_stat_path, 'w')
        fout2.write('\t'.join(['#query', 'query_length',
                               'query_aligned_length',
                               'query_aligned_cov',
                               'query_aligned_region',
                               'query_cov_region',
                               'query_uncov_region']) + '\n')
        for key in query_map_region:
            query_map_region[key] = gi.GIntervalList(sorted(query_map_region[key]))
            query_cov_region[key] = gi.union_interval(query_map_region[key])
            query_map_length[key] = gi.union_length(query_cov_region[key])
            query_uncov_region[key] = gi.complementary_interval(query_cov_region[key],
                                                                gi.GInterval([0, query_length[key]]),
                                                                min_len=filtered_length)
            qcov = query_map_length[key] / query_length[key]
            seq_cov_dict[key] = qcov
            fout2.write('\t'.join([str(x) for x in [key, query_length[key],
                                                    query_map_length[key],
                                                    seq_cov_dict[key],
                                                    query_map_region[key],
                                                    query_cov_region[key],
                                                    query_uncov_region[key]]]) + '\n')
            if qcov * 100 > least_coverage:
                output_set.add(key)
        fout1.write('\n'.join(output_set))
        fout1.close()
        fout2.close()
        return len(output_set), seq_cov_dict

    else:
        # max aligned length of hsp
        with open(paf_path) as fin:
            with open(out_seq_list_path, 'w') as fout:
                for line in fin:
                    temp = line.rstrip().split('\t')
                    # qcov 0 ~ 100
                    qcov = float(temp[9]) / float(temp[1]) * 100
                    # write seq cov 0 ~ 1
                    if temp[0] in seq_cov_dict:
                        seq_cov_dict[temp[0]] = max(qcov / 100, seq_cov_dict[temp[0]])
                    else:
                        seq_cov_dict[temp[0]] = qcov / 100
                    if qcov > least_coverage:
                        output_set.add(temp[0])
                fout.write('\n'.join(output_set))
        return len(output_set), seq_cov_dict


def gff_add_cov(rawgff, newgff, cov_dict):
    with open(rawgff) as fin:
        with open(newgff, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    fout.write(line)
                    continue
                temp = line.rstrip().split('\t')
                attr = gi.string2dict(temp[8])
                seqid = attr['ID']
                if 'qcov' in attr:
                    continue
                if seqid in cov_dict:
                    fout.write(line.rstrip() + "qcov={i};\n".format(i=str(round(cov_dict[seqid], 6))))
                else:
                    fout.write(line.rstrip() + "qcov={i};\n".format(i=str(0)))
    return len(cov_dict)


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
                         help='Using minimap2 to realign to reference genome/mitochondrion/plastid and drop high '
                              'similar unaligned block sequences',
                         type=str, default='')
    parserr.add_argument('-ri', '--realign_identity', metavar='<int>',
                         help='[Only use when -rr on] Min alignment identity of sequences in realign step (choices: '
                              '80/90/95. default: 90)',
                         type=int, choices=[80, 90, 95],
                         default=90)
    parserr.add_argument('-rc', '--realign_coverage', metavar='<int>',
                         help='[Only use when -rr on] Min alignment coverage of dropped sequences (default: 80)',
                         type=int,
                         default=80)
    parserr.add_argument('-rs', '--realign_segmentmode', metavar='<int>',
                         help='[Only use when -rr on] Compute coverage with all segment parts (default: 0)',
                         type=int, choices=[0, 1],
                         default=0)
    parserr.add_argument('-rm', '--realign_minimap2', metavar='<minimap2_path>',
                         help='[Only use when -rr on] Path of minimap2 (default: minimap2 in $PATH)', type=str,
                         default='minimap2')
    parserr.add_argument('-rt', '--realign_thread', metavar='<int>',
                         help='[Only use when -rr on] Number of threads using minimap2 (default: 1)', type=int,
                         default=1)
    parserr.add_argument('-rd', '--realign_dir', metavar='<str>',
                         help='[Only use when -rr on] Temp directory to realign (default: unalnbseq_temp_[Time])',
                         type=str, default='unalnbseq_temp')
    parserr.add_argument('-rk', '--realign_keepdir', metavar='<int>', choices=[0, 1],
                         help='[Only use when -rr on] Keep temp directory (default: 0)',
                         type=int, default=0)
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
                                                                                  output_fasta, unalign_table,
                                                                                  sample_tag, nbase_ignore)
    logging.info(
        "# Load {inseq_n} chromosomes/contigs/scaffolds from fasta, write {blockoutseq_n} blocks "
        "from {outseq_n} chromosomes/contigs/scaffolds.".format(
            inseq_n=input_record_num,
            outseq_n=output_record_num,
            blockoutseq_n=output_blockseq_num
        ))

    # check if blocks are got
    if output_blockseq_num == 0:
        logging.error("# No blocks are get, please check names in {af} and {ut} are same.".format(
            af=assembly_file,
            ut=unalign_table
        ))
        exit(1)

    if realign_step:
        # check temp dir
        if args["realign_dir"] == "unalnbseq_temp":
            temp_dir = args["realign_dir"] + '_' + time.strftime("%Y%m%d%H%M%S", time.localtime(time.time()))
        else:
            # self defined dir
            temp_dir = args["realign_dir"]
        if os.path.isdir(temp_dir):
            logging.warning("# {temp_path} exists! It will be rewrited!".format(temp_path=temp_dir))
        else:
            os.mkdir(temp_dir)
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

        if args["realign_segmentmode"]:
            add_parameter = '--no-long-join'
        else:
            add_parameter = ''

        command = '{minimap2} {addp} -c -t {thread} {idt_par} {ref} {query} > {temp_dir}/{temp_paf}'.format(
            addp=add_parameter,
            minimap2=args["realign_minimap2"],
            thread=args["realign_thread"],
            idt_par=minimap2_idt_par,
            ref=args["realign_reference"],
            query=output_fasta,
            temp_dir=temp_dir,
            temp_paf=temp_paf
        )
        logging.info("# Realigning block sequences to references: {cmd}".format(cmd=command))
        os.system(command)
        logging.info("# Finish realigning block sequences to references.")
        # scan paf
        temp_out_seq_list_path = "mapped_bseq.txt"
        temp_out_seq_stat_path = "mapped_bseq.stat"
        logging.info("# Compute block sequence coverage. \
        # Segment mode = {smode}.".format(smode=args["realign_segmentmode"]))
        drop_n, seq_cov = drop_seq_from_paf(temp_dir + '/' + temp_paf,
                                            args["realign_identity"],
                                            args["realign_coverage"],
                                            length_cutoff,
                                            temp_dir + '/' + temp_out_seq_list_path,
                                            temp_dir + '/' + temp_out_seq_stat_path,
                                            args["realign_segmentmode"])
        # copy raw bsq.gff
        temp_out_gff_raw_path = "unalign_bseq.gff3"
        os.system('cp {output_gff} {temp_dir}/{temp_gff}'.format(temp_dir=temp_dir,
                                                                 temp_gff=temp_out_gff_raw_path,
                                                                 output_gff=output_fasta + ".gff3"))
        # add query cov in gff
        aqc_return_status = gff_add_cov(temp_dir + '/' + temp_out_gff_raw_path,
                                        temp_dir + '/' + 'temp.gff3',
                                        seq_cov)
        os.system('mv {newgff} {oldgff}'.format(newgff=temp_dir + '/' + 'temp.gff3',
                                                oldgff=temp_dir + '/' + temp_out_gff_raw_path))
        # filter bseq.fa and bseq.gff
        temp_out_seq_filtered_path = "remain_bseq.fa"
        temp_out_gff_filtered_path = "remain_bseq.gff3"
        fsr_return_status = gi.fa_some_record(output_fasta, temp_dir + '/' + temp_out_seq_list_path,
                                              temp_dir + '/' + temp_out_seq_filtered_path, exclude=True)
        gsr_return_status = gi.gff_some_record(temp_dir + '/' + temp_out_gff_raw_path,
                                               temp_dir + '/' + temp_out_seq_list_path,
                                               temp_dir + '/' + temp_out_gff_filtered_path, key="ID", exclude=True)
        logging.info("# Remove {n} block sequences similiar to references.".format(n=drop_n))
        # mv / rm temp files
        logging.info("# Update block sequences.")
        os.system('cp {temp_dir}/{temp_fa} {output_fa}'.format(temp_dir=temp_dir,
                                                               temp_fa=temp_out_seq_filtered_path,
                                                               output_fa=output_fasta))
        os.system('cp {temp_dir}/{temp_gff} {output_gff}'.format(temp_dir=temp_dir,
                                                                 temp_gff=temp_out_gff_filtered_path,
                                                                 output_gff=output_fasta + ".gff3"))
        if not args["realign_keepdir"]:
            logging.info("# Remove temp files.")
            os.system('rm -rf {temp_dir}'.format(temp_dir=temp_dir))
        logging.info("# Realign step is finished.")
