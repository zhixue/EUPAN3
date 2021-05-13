#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/28 3:32 PM
    @Usage: python3 ptpg.py -i x.gtf/x.gff3 -r cds -o x_pTpG.gtf/x_pTpG.gff
"""
import logging
import os
import argparse
from Genome_Interval import string2dict
from tlog import *


def ptpg_gtf(ingtf, region, outgtf):
    current_gene = ''
    block_string = dict()
    block_length = dict()

    with open(ingtf) as f:
        with open(outgtf, 'w') as fout:
            for line in f:
                temp = line.rstrip().split('\t')
                # like transcript_id "LOC_Os01g01010.1"; gene_id "LOC_Os01g01010";
                attr = string2dict(temp[8], eq=' ', rm_quote=True)
                gene_id = attr['gene_id']
                if trans_id in attr:
                    trans_id = attr['transcript_id']
                else:
                    continue
                if current_gene != gene_id and current_gene != '':
                    # check
                    max_len = 0
                    max_len_trans = ''
                    if block_length:
                        for trans in block_length:
                            if block_length[trans] > max_len:
                                max_len = block_length[trans]
                                max_len_trans = trans
                        fout.write(block_string[max_len_trans])
                    # init
                    block_string = dict()
                    block_length = dict()

                current_gene = gene_id
                if trans_id not in block_string:
                    block_string[trans_id] = line
                else:
                    block_string[trans_id] += line
                if temp[2] == region:
                    start = int(temp[3])
                    end = int(temp[4])
                    if trans_id not in block_length:
                        block_length[trans_id] = 0
                    block_length[trans_id] += end - start + 1

            # last one
            if current_gene != '':
                # check
                max_len = 0
                max_len_trans = ''
                for trans in block_length:
                    if block_length[trans] > max_len:
                        max_len = block_length[trans]
                        max_len_trans = trans
                fout.write(block_string[max_len_trans])


def ptpg_gff(ingff, region, outgff):
    gene_length_dict = dict()  # {'gene1':{'mRNA1':89,'mRNA2':81}}
    transcripts_dict = dict()  # {'mRNA1':[(120,200),(202,209)],'mRNA2':[(120,200)]}

    with open(ingff) as f:
        for line in f:
            # pass comment and blank
            if line.startswith('#'):
                continue
            if line.rstrip() == '':
                continue

            temp = line.rstrip().split('\t')
            line_type = temp[2]

            lineid = string2dict(temp[-1])['ID']
            if line_type == 'gene':
                gene_length_dict[lineid] = dict()
            elif line_type in ('transcript', 'mRNA'):
                geneid = string2dict(temp[-1])['Parent']
                transcriptid = lineid
                gene_length_dict[geneid][transcriptid] = 0
                transcripts_dict[transcriptid] = []
            elif line_type == region:
                transcriptid = string2dict(temp[-1])['Parent']
                # avoid same id of exon/CDS
                # e_id = lineid
                e_start = int(temp[3])
                e_end = int(temp[4])
                transcripts_dict[transcriptid] += [(e_start, e_end)]
                gene_length_dict[geneid][transcriptid] += e_end - e_start + 1
    logging.info("# Load {n1} genes, {n2} transcripts.".format(n1=len(gene_length_dict), n2=len(transcripts_dict)))

    # get max
    for gene in gene_length_dict:
        maxtranscript = ''
        maxlen = 0
        if len(gene_length_dict[gene]) >= 2:
            for transcript in gene_length_dict[gene]:
                if gene_length_dict[gene][transcript] > maxlen:
                    maxtranscript = transcript
                    maxlen = gene_length_dict[gene][transcript]
            gene_length_dict[gene] = {maxtranscript: maxlen}

    # write
    write_gene_n = 0
    with open(ingff) as f:
        with open(outgff, 'w') as fout:
            for line in f:
                # pass comment and blank
                if line.startswith('#'):
                    continue
                if len(line.rstrip()) == 0:
                    continue

                temp = line.rstrip().split('\t')
                line_type = temp[2]
                attr = string2dict(temp[-1])
                lineid = attr['ID']
                # may not write genes, transcripts with no exon/CDS
                if line_type == 'gene':
                    geneline = line
                elif line_type in ('transcript', 'mRNA'):
                    geneid = attr['Parent']
                    transcriptid = lineid
                    if transcriptid in gene_length_dict[geneid]:
                        transcriptline = line
                else:
                    transcriptid = attr['Parent']
                    if transcriptid in gene_length_dict[geneid]:
                        if gene_length_dict[geneid][transcriptid] > 0:
                            if geneline != '' and transcriptline != '':
                                fout.write(geneline + transcriptline)
                                geneline = ''
                                transcriptline = ''
                                write_gene_n += 1
                            fout.write(line)
    logging.info("# Write {n1} genes, {n1} transcripts.".format(n1=write_gene_n))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Extract the longest transcript elements' records for each gene, 
    according to CDS or exon length''')
    parser.add_argument('-i', '--input', metavar='<input.gff/gtf>',
                        help='Path of input gff/gtf', type=str, required=True)
    parser.add_argument('-r', '--region', metavar='<str>',
                        help='CDS or exon (default: CDS)', type=str,
                        choices=['CDS', 'exon'], default='CDS')
    parser.add_argument('-o', '--output', metavar='<output.gff/gtf>',
                        help='Path of output gff/gtf', type=str, required=True)

    args = vars(parser.parse_args())
    input_path = os.path.abspath(args['input'])
    output_path = os.path.abspath(args['output'])
    selected_region = args['region']

    # error check
    if input_path == output_path:
        logging.error('# Error: input is same as output gff!')
        exit(1)

    if input_path.endswith('gtf'):
        ptpg_gtf(input_path, selected_region, output_path)
    elif input_path.endswith('gff') or input_path.endswith('gff3'):
        ptpg_gff(input_path, selected_region, output_path)
