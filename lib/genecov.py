#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/28 5:41 PM
    @Usage:
"""
import argparse
import logging

from Genome_Interval import *
from tlog import *
import os

def readgff(gff, ele_select="CDS"):
    gene_dict = dict()  # {'chr1':{'geneA':{'transcript1':{'CDS1':''}}}}
    ele_dict = dict()  # {'geneA':(start,end)}
    region_dict = dict()  # {'chr1':[(start,end)]}
    current_gene_id = ''
    current_transcript_id = ''
    with open(gff) as f:
        for line in f:
            if line.startswith('#') or line.rstrip() == '':
                continue
            temp = line.rstrip().split('\t')
            if temp[2] not in ["gene", "mRNA", "transcript", ele_select]:
                continue
            ele = GFFElement(line)
            ele_chrn = ele.chrn
            if ele_chrn not in gene_dict:
                gene_dict[ele_chrn] = dict()
                region_dict[ele_chrn] = []
            ele_type = ele.type
            if ele_type == "gene":
                # ele = GFFGene(line)
                current_gene_id = ele.get_id()
                if current_gene_id not in gene_dict[ele_chrn] and current_gene_id != '':
                    gene_dict[ele_chrn][current_gene_id] = dict()
            elif ele_type == "mRNA" or ele_type == "transcript":
                # ele = GFFTranscript(line)
                current_transcript_id = ele.get_id()
                if current_transcript_id not in gene_dict[ele_chrn][current_gene_id]:
                    gene_dict[ele_chrn][current_gene_id][current_transcript_id] = dict()
            elif ele_type == ele_select:
                # ele = GFFElement(line)
                gene_dict[ele.chrn][current_gene_id][current_transcript_id][ele.get_id()] = ''
            ele_region = (ele.start, ele.end)
            ele_dict[ele.get_id()] = ele_region
            if ele_region not in region_dict[ele_chrn]:
                region_dict[ele_chrn] += [ele_region]
    return gene_dict, ele_dict, region_dict


def readgtf(gtf, ele_select="CDS"):
    gene_dict = dict()  # {'chr1':{'geneA':{'transcript1':{'CDS1':''}}}}
    ele_dict = dict()  # {'geneA':(1,2)}
    region_dict = dict()  # {'chr1':{(1,2):['xx']}]}
    current_gene_id = ''
    current_transcript_id = ''
    current_ele_id = ''
    current_eles = []
    current_transcripts = []
    with open(gtf) as f:
        for line in f:
            if line.startswith('#') or line.rstrip() == '':
                continue
            temp = line.rstrip().split('\t')
            if temp[2] != ele_select:
                continue
            ele_chrn = temp[0]
            if ele_chrn not in gene_dict:
                gene_dict[ele_chrn] = dict()
                region_dict[ele_chrn] = dict()
            ele_type = temp[2]
            if ele_type == ele_select:
                ele = GTFElement(line)
                # new transcript/gene
                if current_transcript_id != '':
                    if ele.details["transcript_id"] != current_transcript_id:
                        transcript_ele = GTFVirtual(current_eles)
                        transcript_region = (transcript_ele.start, transcript_ele.end)
                        current_transcripts += [transcript_ele]
                        ele_dict[current_transcript_id] = transcript_ele
                        if transcript_region not in region_dict[ele_chrn]:
                            region_dict[ele_chrn][transcript_region] = []
                        region_dict[ele_chrn][transcript_region] += [current_transcript_id]
                        current_eles = []
                    if ele.details["gene_id"] != current_gene_id:
                        gene_ele = GTFVirtual(current_transcripts)
                        gene_region = (gene_ele.start, gene_ele.end)
                        ele_dict[current_gene_id] = gene_ele
                        if gene_region not in region_dict[ele_chrn]:
                            region_dict[ele_chrn][gene_region] = []
                        region_dict[ele_chrn][gene_region] += [current_gene_id]
                        current_transcripts = []
                current_gene_id = ele.details["gene_id"]
                current_transcript_id = ele.details["transcript_id"]
                if current_gene_id not in gene_dict[ele_chrn]:
                    gene_dict[ele_chrn][current_gene_id] = dict()
                if current_transcript_id not in gene_dict[ele_chrn][current_gene_id]:
                    gene_dict[ele_chrn][current_gene_id][current_transcript_id] = dict()
                current_eles += [ele]
                current_ele_id = ele.get_id()
                gene_dict[ele_chrn][current_gene_id][current_transcript_id][current_ele_id] = ''
                ele_region = (ele.start, ele.end)
                if ele_region not in region_dict[ele_chrn]:
                    region_dict[ele_chrn][ele_region] = current_ele_id
        # last one
        transcript_ele = GTFVirtual(current_eles)
        transcript_region = (transcript_ele.start, transcript_ele.end)
        current_transcripts += [transcript_ele]
        ele_dict[current_transcript_id] = transcript_ele
        if transcript_region not in region_dict[ele_chrn]:
            region_dict[ele_chrn][transcript_region] = []
        region_dict[ele_chrn][transcript_region] += [current_transcript_id]

        gene_ele = GTFVirtual(current_transcripts)
        gene_region = (gene_ele.start, gene_ele.end)
        ele_dict[current_gene_id] = gene_ele
        if gene_region not in region_dict[ele_chrn]:
            region_dict[ele_chrn][gene_region] = []
        region_dict[ele_chrn][gene_region] += [current_gene_id]
    return gene_dict, ele_dict, region_dict


def compute_cov(cov_region, annotation_dict, chrn, sample_tag):
    outresults = []
    cov_region_intervallist = GIntervalList(cov_region)
    if not chrn in annotation_dict[0]:
        return outresults
    for gene in annotation_dict[0][chrn]:
        # gene
        gene_region_intervallist = GIntervalList([annotation_dict[1][gene]])
        gene_cov_interval = intersect_interval(gene_region_intervallist, cov_region_intervallist)
        gene_cov = union_length(gene_cov_interval) / union_length(gene_region_intervallist)
        # transcript
        for transcript in annotation_dict[0][chrn][gene].keys():
            transcript_region = annotation_dict[1][transcript]
            transcript_region_intervallist = GIntervalList([transcript_region])
            if transcript_region_intervallist == gene_region_intervallist:
                transcript_cov = gene_cov
            else:
                transcript_cov_interval = intersect_interval(transcript_region_intervallist, cov_region_intervallist)
                transcript_length = union_length(transcript_region_intervallist)
                if transcript_length == 0:
                    transcript_cov = 0
                else:
                    transcript_cov = union_length(transcript_cov_interval) / transcript_length
            # elements
            elements = annotation_dict[0][chrn][gene][transcript].keys()
            elements_region = [annotation_dict[1][ele] for ele in elements]
            element_region_intervallist = GIntervalList(elements_region)
            if element_region_intervallist == gene_region_intervallist:
                element_cov = transcript_cov
            else:
                element_cov_interval = intersect_interval(element_region_intervallist, cov_region_intervallist)
                element_length = union_length(element_region_intervallist)
                if element_length == 0:
                    element_cov = 0
                else:
                    element_cov = union_length(element_cov_interval) / element_length
            outresults += ['\t'.join([str(x) for x in [gene, transcript, sample_tag, gene_cov, transcript_cov, element_cov]])]
    return outresults


def scan_bed(bedfile, annotation_dict, output, sample_tag='', at_least_depth=1, scan_depth=False):
    current_chrn = ''
    current_chrn_covregion = []
    current_chrn_depthregion = dict()
    fout = open(output, 'w')
    with open(bedfile) as f:
        for line in f:
            temp = line.rstrip().split('\t')
            if temp[0] != current_chrn:
                # compute
                temp_results = compute_cov(current_chrn_covregion, annotation_dict, current_chrn, sample_tag)
                if not temp_results:
                    fout.write('\n'.join([str(x) for x in temp_results]) + '\n')
                    fout.flush()
                current_chrn_covregion = []
            cov = float(temp[3])  # int error if 1.1e6
            current_chrn = temp[0]
            start_pos = int(temp[1]) + 1
            end_pos = int(temp[2])
            if cov < at_least_depth:
                continue
            # update chrn_covregion
            if not current_chrn_covregion:
                current_chrn_covregion += [[start_pos, end_pos]]
            elif start_pos == current_chrn_covregion[-1][1]+1:
                current_chrn_covregion[-1][1] = end_pos
            else:
                current_chrn_covregion += [[start_pos, end_pos]]
        # final one, compute
        temp_results = compute_cov(current_chrn_covregion, annotation_dict, current_chrn, sample_tag)
        if not temp_results:
            fout.write('\n'.join([str(x) for x in temp_results]))
            fout.flush()
    fout.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Compute gene coverage''')
    parser.add_argument('-a', '--annotation', metavar='<input.gff/gtf>', help='Path of input gff/gtf', type=str, required=True)
    parser.add_argument('-b', '--bed', metavar='<input.bed>', help='bed of of coverage from bedtools', type=str, required=True)
    parser.add_argument('-r', '--region', metavar='<str>', help='CDS or exon (default: CDS)', type=str, choices=['CDS', 'exon'], default='CDS')
    parser.add_argument('-o', '--output', metavar='<output.cov>', help='Path of output cov', type=str, required=True)
    parser.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
    parser.add_argument('-m', '--min_depth', metavar='<int>', help='Min depth', type=int, default=1)

    args = vars(parser.parse_args())
    annotation_path = os.path.abspath(args['annotation'])
    output_path = os.path.abspath(args['output'])
    bed_path = os.path.abspath(args['bed'])
    sample_tag = args['sample_name']
    region = args['region']
    min_depth = args['min_depth']

    if annotation_path.endswith('gtf'):
        annotation = readgtf(annotation_path, region)
    elif annotation_path.endswith('gff') or annotation_path.endswith('gff3'):
        annotation = readgff(annotation_path, region)

    chrn_n = 0
    gene_n = 0
    trans_n = 0
    ele_n = 0
    for chrn in annotation[0]:
        chrn_n += 1
        for gene in annotation[0][chrn]:
            gene_n += 1
            for transcript in annotation[0][chrn][gene]:
                trans_n += 1
                ele_n += len(annotation[0][chrn][gene][transcript])

    logging.info('# Load {chrn_n} chromosomes, {gene_n} genes, {tran_n} transcripts, {ele_n} {region}.'.format(
        chrn_n=chrn_n,
        gene_n=gene_n,
        tran_n=trans_n,
        ele_n=ele_n,
        region=region
    ))
    scan_bed(bed_path, annotation, output_path, sample_tag, min_depth, scan_depth=False)





