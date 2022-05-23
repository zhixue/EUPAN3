#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2022/5/23 1:26 PM
    @Usage: python3 elongatedbseq.py -n novel_bseq.gff3 -a annotation.gff3 -p elongated_bseq.bed
"""

import argparse
import os
from tlog import *
from Genome_Interval import GFFElement, overlap


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Elongate blocks with annotations to keep complete genes''')
    parser.add_argument('-n', '--novel_bseq_gff_path', metavar='<novel_bseq.gff3>',
                        help='Path of novel bseq gff', type=str, required=True)
    parser.add_argument('-a', '--annotation_gff_path', metavar='<annotation.gff3>',
                        help='Path of annotation gff', type=str, required=True)
    parser.add_argument('-p', '--prefix_output', metavar='<elongated_bseq>',
                        help='Prefix of output (.table/.fa/.gff)', type=str, required=True)
    parser.add_argument('-f', '--fa', metavar='<input.fa>', help='Path of input fa', type=str, required=True)
    parser.add_argument('-t', '--tool_bedtools', metavar='<str>', help='Path of bedtools (Default: bedtools)', type=str,
                        default='bedtools')

    args = vars(parser.parse_args())
    novel_bseq_gff_path = os.path.abspath(args['novel_bseq_gff_path'])
    annotation_gff_path = os.path.abspath(args['annotation_gff_path'])
    input_fa = os.path.abspath(args['fa'])
    prefix_output = args['prefix_output']
    output_table_path = prefix_output + '.table'
    output_fa_path = prefix_output + '.fa'
    output_gff_path = prefix_output + '.gff'
    bedtools = args['tool_bedtools']

    # load gene gff
    gene_dict = dict()  # {'chr1':{'geneA':''}}
    ele_dict = dict()  # {'geneA':(start,end)}
    region_dict = dict()  # {'chr1':[(start,end),]}
    with open(annotation_gff_path) as fin:
        for line in fin:
            if line.startswith('#') or line.rstrip() == '':
                continue
            temp = line.rstrip().split('\t')
            if temp[2] == "gene":
                ele = GFFElement(line)
                ele_id = ele.get_id()
                if ele.chrn not in gene_dict:
                    gene_dict[ele.chrn] = dict()
                    region_dict[ele.chrn] = []
                gene_dict[ele.chrn][ele_id] = ''
                ele_region = (ele.start, ele.end)
                ele_dict[ele_id] = ele_region
                if ele_region not in region_dict[ele.chrn]:
                    region_dict[ele.chrn] += [ele_region]
    for key in region_dict:
        region_dict[key] = sorted(region_dict[key])
    logging.info('# Finish loading {n1} chrs, {n2} genes from gene annotations.'.format(n1=len(gene_dict),
                                                                                        n2=len(ele_dict)))


    # load bseq gff, find overlap and write table
    current_anno_idx = 0
    current_record = []  # ["Chr", "EStart", "EEnd", "Eid", "BseqIDs", "GeneIDs"]
    last_record = []
    table_fout = open(output_table_path, 'w')

    with open(novel_bseq_gff_path) as fin:
        for line in fin:
            if line.startswith('#') or line.rstrip() == '':
                continue
            cblock = GFFElement(line)

            last_record = current_record
            current_record = [cblock.chrn, cblock.start, cblock.end, cblock.get_id() + 'E', cblock.get_id(), '.']
            # skip chromosome without annotations
            if cblock.chrn not in gene_dict:
                # no elongation and write
                if last_record:
                    last_record[1] -= 1
                    table_fout.write('\t'.join([str(x) for x in last_record]) + '\n')
            else:
                # update current_anno_idx
                while current_anno_idx < len(region_dict[cblock.chrn]) - 1:
                    if region_dict[cblock.chrn][current_anno_idx][1] < cblock.start:
                        current_anno_idx += 1
                    else:
                        break

                # overlap with last record, update current record left
                if not last_record:
                    continue
                if last_record[0] == current_record[0] \
                    and overlap((last_record[1], last_record[2]),
                                (current_record[1], current_record[2])):
                    current_record[1] = min(current_record[1], last_record[1])
                    current_record[2] = max(current_record[2], last_record[2])
                    current_record[3] = last_record[3]
                    current_record[4] = last_record[4] + ',' + current_record[4]
                    current_record[5] = last_record[5]
                else:
                    last_record[1] -= 1
                    table_fout.write('\t'.join([str(x) for x in last_record]) + '\n')
                    if last_record[0] != current_record[0]:
                        current_anno_idx = 0

                # overlap with gene, update current record left/right
                while overlap((current_record[1], current_record[2]), region_dict[cblock.chrn][current_anno_idx]):
                    current_record[1] = min(current_record[1], region_dict[cblock.chrn][current_anno_idx][0])
                    current_record[2] = max(current_record[2], region_dict[cblock.chrn][current_anno_idx][1])
                    if current_record[5] == '.':
                        current_record[5] = tuple(gene_dict[cblock.chrn].keys())[current_anno_idx]
                    else:
                        current_record[5] += ',' + tuple(gene_dict[cblock.chrn].keys())[current_anno_idx]
                    if current_anno_idx == len(gene_dict[cblock.chrn]) - 1:
                        break
                    current_anno_idx += 1
        # final one
        last_record = current_record
        last_record[1] -= 1
        table_fout.write('\t'.join([str(x) for x in last_record]) + '\n')
    table_fout.close()
    logging.info('# Finish elongating and writing table.')

    # write elongated fasta (use bedtools getfasta)
    cmd2 = bedtools + ' getfasta -fi ' + input_fa + ' -bed ' + output_table_path + ' -fo ' + output_fa_path + ' -nameOnly'
    os.system(cmd2)
    # os.system('rm -rf ' + temp_bed_path)
    logging.info('# Finish writing elongated fasta.')

    # transform gene annotation to elongated fasta (chrn, location change)
    # load used genes
    used_genes_dict = dict()
    with open(output_table_path) as fin:
        for line in fin:
            temp = line.rstrip().split()
            if len(temp) < 6:
                continue
            if temp[5] == '' or temp[5] == '.':
                continue
            genes = temp[5].split(',')
            local_0_position = int(temp[1])
            local_chrn_id = temp[3]
            for gene in genes:
                used_genes_dict[gene] = (local_chrn_id, local_0_position)

    gff_fout = open(output_gff_path, 'w')
    # read gff and write
    with open(annotation_gff_path) as fin:
        for line in fin:
            if line.startswith('#') or line.rstrip() == '':
                continue
            try:
                ele = GFFElement(line)
            except:
                print(line.split())
                exit()
            if ele.type == 'gene':
                if ele.get_id() in used_genes_dict:
                    # transform columns 1, 5, 6
                    ele.chrn = used_genes_dict[ele.get_id()][0]
                    ele.start -= used_genes_dict[ele.get_id()][1]
                    ele.end -= used_genes_dict[ele.get_id()][1]
                    gff_fout.write(ele.tostring() + '\n')
            else:
                parentid = ele.find_parentid()
                if parentid in used_genes_dict:
                    if ele.get_id() not in used_genes_dict:
                        used_genes_dict[ele.get_id()] = used_genes_dict[parentid]
                    # transform columns 1, 5, 6
                    ele.chrn = used_genes_dict[ele.get_id()][0]
                    ele.start -= used_genes_dict[ele.get_id()][1]
                    ele.end -= used_genes_dict[ele.get_id()][1]
                    gff_fout.write(ele.tostring() + '\n')
    gff_fout.close()
    logging.info('# Finish writing elongated gff.')




