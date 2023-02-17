#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2023/1/`6 9:44 AM
    @Usage: python3 minimap2clust.py -i raw_sorted.fa -p xx.paf -g idt -o seqfa
"""
import argparse


def read_fa(sorted_fa_file, read_len=False):
    # read fasta
    ctg_ids = dict()
    ctg_len = dict()
    i = 1
    if read_len:
        with open(sorted_fa_file) as f:
            for line in f:
                if line.startswith('>'):
                    ctg = line.rstrip().split()[0][1:]
                    ctg_ids[ctg] = i
                    ctg_len[ctg] = 0
                    i += 1
                else:
                    ctg_len[ctg] += len(line.rstrip())
    else:
        with open(sorted_fa_file) as f:
            for line in f:
                if line.startswith('>'):
                    ctg = line.rstrip().split()[0][1:]
                    ctg_ids[ctg] = i
                    i += 1
    return ctg_ids, ctg_len


def line2temp(string, dataformat='paf'):
    temp = string.rstrip().split('\t')
    if dataformat == 'paf':
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
        query_name = temp[0]
        target_name = temp[5]
        query_length = float(temp[1])
        target_length = float(temp[6])
        matchedlength = float(temp[9])
    else:
        query_name = temp[0]
        target_name = temp[1]
        block_len = int(temp[4])
        block_idt = float(temp[9])
        query_length = int(temp[2])
        target_length = int(temp[3])
        matchedlength = block_len * block_idt / 100
    return query_name, target_name, query_length, target_length, matchedlength


def minimap2clust(sorted_fasta_path,
                  paf_path,
                  global_identity,
                  output_fa,
                  output_clust):
    if global_identity > 1:
        global_identity = global_identity / 100

    if paf_path.endswith('paf'):
        map_format = 'paf'
    else:
        map_format = 'txt'

    # load fasta record order
    ctgids, ctglen = read_fa(sorted_fasta_path)

    # read map results
    tree_dict = dict()  # {branch: [leaf1, leaf2] } , len(branch) >= len(leaf)
    with open(paf_path) as fin:
        for fline in fin:
            qn, tn, ql, tl, ml = line2temp(fline, map_format)
            if qn == tn:
                continue
            if ml / min(ql, tl) >= global_identity:
                qn_order = ctgids[qn]
                tn_order = ctgids[tn]
                if tn_order < qn_order:
                    branch, leaf = tn, qn
                else:
                    branch, leaf = qn, tn
                if branch not in tree_dict:
                    tree_dict[branch] = set()
                tree_dict[branch].add(leaf)

    skip_set = set()
    group_id = 1
    # work in python3.6+ , output clust out
    fout_c = open(output_clust, 'w')
    for name in ctgids:
        if name in skip_set:
            continue
        members_string = ''
        if name in tree_dict:
            members = tree_dict[name]
            filtered_members = set()
            for member in members:
                if member not in skip_set:
                    filtered_members.add(member)
                    skip_set.add(member)
            if filtered_members:
                members_string = ' ' + ' '.join(filtered_members)
        output_line = 'Group {i}: '.format(i=group_id) + name + members_string
        fout_c.write(output_line + '\n')
        group_id += 1
    fout_c.close()
    # output fasta out
    write_flag = 1
    with open(output_fa, 'w') as fout_f:
        with open(sorted_fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    name = line.rstrip().split()[0][1:]
                    if name in skip_set:
                        write_flag = 0
                    else:
                        write_flag = 1
                if write_flag:
                    fout_f.write(line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Cluster the sequences using minimap2 output''')
    parser.add_argument('-i', '--input_fa', metavar='<input.fa>',
                        help='Path of input fasta', type=str, required=True)
    parser.add_argument('-p', '--paf', metavar='<minimap2out.paf>',
                        help='Path of minimap2 paf', type=str, required=True)
    parser.add_argument('-g', '--global_identity', metavar='<int>',
                        help='global identity threshold (default: 90)', type=int, default=90)
    parser.add_argument('-o', '--output_fa', metavar='<output.fa>',
                        help='Path of output fasta', type=str, required=True)
    parser.add_argument('-oc', '--output_clust', metavar='<output.clust>',
                        help='Path of output clust (default: <output.fa>.clstr)', type=str, default='')
    args = vars(parser.parse_args())

    sortedfastapath = args["input_fa"]
    mappath = args["paf"]
    globalidentity = args["global_identity"]
    outfastapath = args["output_fa"]
    if args["output_clust"] == '':
        outclustpath = args["output_fa"] + '.clstr'
    else:
        outclustpath = args["output_clust"]

    minimap2clust(sortedfastapath,
                  mappath,
                  globalidentity,
                  outfastapath,
                  outclustpath)
