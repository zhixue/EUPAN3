#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/28 5:41 PM
    @Usage: bedtools genomecov -bga -split -ibam xx.bam > xx.bed;
            python3 genecov.py -a ref_pTpG.gff -o xx.cov -n xx -b xx.bed
"""
import argparse
import os
import sys
from Genome_Interval import *
from tlog import *


def readgff(gff, ele_select="CDS"):
    gene_dict = dict()  # {'chr1':{'geneA':{'transcriptA_1':{'CDSA_1_1':''}}}}
    ele_dict = dict()  # {'geneA':(start,end),'transcriptA_1':(start,end),'CDSA_1_1':(start,end)}
    region_dict = dict()  # {'chr1':[(start,end)]}
    current_gene_id = ''
    current_transcript_id = ''
    exist_same_element = False
    with open(gff) as f:
        for line in f:
            if line.startswith('#') or line.rstrip() == '':
                continue
            temp = line.rstrip().split('\t')
            if temp[2] not in ["gene", "mRNA", "transcript", ele_select]:
                continue
            ele = GFFElement(line)
            ele_id = ele.get_id()
            if ele.chrn not in gene_dict:
                gene_dict[ele.chrn] = dict()
                region_dict[ele.chrn] = []
            if ele.type == "gene":
                current_gene_id = ele_id
                if current_gene_id not in gene_dict[ele.chrn] and current_gene_id != '':
                    gene_dict[ele.chrn][current_gene_id] = dict()
            elif ele.type == "mRNA" or ele.type == "transcript":
                current_transcript_id = ele_id
                if current_transcript_id not in gene_dict[ele.chrn][current_gene_id] and current_gene_id != '':
                    gene_dict[ele.chrn][current_gene_id][current_transcript_id] = dict()
                    # init element time
                    uniq_element_exist_time = 1
            elif ele.type == ele_select:
                # avoid same ID of exon/CDS
                if ele_id in gene_dict[ele.chrn][current_gene_id][current_transcript_id]:
                    exist_same_element = True
                    ele_part = ele_id.split(':')
                    if ele_part[-1] == str(uniq_element_exist_time):
                        ele_id = ':'.join(ele_part[:-1]) + ':' + str(uniq_element_exist_time + 1)
                    else:
                        ele_id += ':' + str(uniq_element_exist_time + 1)
                    uniq_element_exist_time += 1
                gene_dict[ele.chrn][current_gene_id][current_transcript_id][ele_id] = 1
            ele_region = (ele.start, ele.end)
            ele_dict[ele_id] = ele_region
            if ele_region not in region_dict[ele.chrn]:
                region_dict[ele.chrn] += [ele_region]
    if exist_same_element:
        logging.warning("# Exist same ID in " + ele_select + '! Add ":[number]" (e.g. xx:2,...) after the same ID!')
    return gene_dict, ele_dict, region_dict


def readgtf(gtf, ele_select="CDS"):
    gene_dict = dict()  # {'chr1':{'geneA':{'transcript1':{'CDS:xx:1':''}}}}
    ele_dict = dict()  # {'geneA':(1,2)}
    region_dict = dict()  # {'chr1':[(1,2),]}
    current_gene_id = ''
    current_transcript_id = ''
    current_eles = []
    current_transcripts = []
    # store the gene/transcript region if the record exists
    gene_obj_dict = dict()
    transcript_obj_dict = dict()
    no_gene_record_flag = 1
    no_transcript_record_flag = 1
    with open(gtf) as f:
        for line in f:
            if line.startswith('#') or line.rstrip() == '':
                continue
            temp = line.rstrip().split('\t')
            if temp[2] not in ("gene", "transcript", "mRNA", ele_select):
                continue
            ele_chrn = temp[0]
            if ele_chrn not in gene_dict:
                gene_dict[ele_chrn] = dict()
                region_dict[ele_chrn] = []
            ele_type = temp[2]
            ele = GTFElement(line)
            if ele_type == "gene":
                gene_obj_dict[ele.details["gene_id"]] = ele
                no_gene_record_flag = 0
            elif ele_type in ("transcript", "mRNA"):
                transcript_obj_dict[ele.details["transcript_id"]] = ele
                no_transcript_record_flag = 0
            elif ele_type == ele_select:
                # new transcript/gene
                if current_transcript_id != '':
                    # new transcipt
                    if ele.details["transcript_id"] != current_transcript_id:
                        if current_transcript_id in transcript_obj_dict:
                            transcript_ele = transcript_obj_dict[current_transcript_id]
                        else:
                            transcript_ele = GTFVirtual(current_eles)
                        current_transcripts += [transcript_ele]
                        ele_dict[current_transcript_id] = transcript_ele.region
                        if transcript_ele.region not in region_dict[ele_chrn]:
                            region_dict[ele_chrn] += [transcript_ele.region]
                        current_eles = []
                    # new gene
                    if ele.details["gene_id"] != current_gene_id:
                        if current_gene_id in gene_obj_dict:
                            gene_ele = gene_obj_dict[current_gene_id]
                        else:
                            gene_ele = GTFVirtual(current_transcripts)
                        ele_dict[current_gene_id] = gene_ele.region
                        if gene_ele.region not in region_dict[ele_chrn]:
                            region_dict[ele_chrn] += [gene_ele.region]
                        current_transcripts = []
                # update current
                current_gene_id = ele.details["gene_id"]
                current_transcript_id = ele.details["transcript_id"]
                if current_gene_id not in gene_dict[ele_chrn]:
                    gene_dict[ele_chrn][current_gene_id] = dict()
                if current_transcript_id not in gene_dict[ele_chrn][current_gene_id]:
                    gene_dict[ele_chrn][current_gene_id][current_transcript_id] = dict()
                current_eles += [ele]
                current_ele_id = ele.get_id()
                gene_dict[ele_chrn][current_gene_id][current_transcript_id][current_ele_id] = 1
                ele_dict[current_ele_id] = ele.region
                if ele.region not in region_dict[ele_chrn]:
                    region_dict[ele_chrn] += [ele.region]
        # last one
        if current_transcript_id in transcript_obj_dict:
            transcript_ele = transcript_obj_dict[current_transcript_id]
        else:
            transcript_ele = GTFVirtual(current_eles)
        current_transcripts += [transcript_ele]
        ele_dict[current_transcript_id] = transcript_ele.region
        if transcript_ele.region not in region_dict[ele_chrn]:
            region_dict[ele_chrn] += [transcript_ele.region]

        if current_gene_id in gene_obj_dict:
            gene_ele = gene_obj_dict[current_gene_id]
        else:
            gene_ele = GTFVirtual(current_transcripts)
        ele_dict[current_gene_id] = gene_ele.region
        if gene_ele.region not in region_dict[ele_chrn]:
            region_dict[ele_chrn] += [gene_ele.region]

    if no_gene_record_flag:
        logging.warning("# No gene records in gtf, use transcripts of gene to infer the start and end!")
    if no_transcript_record_flag:
        logging.warning("# No transcript/mRNA records in gtf, use elements of transcript to infer the start and end!")
    return gene_dict, ele_dict, region_dict


def compute_cov(annotation_dicts, anno_list_object, used_region, chrn, sample_tag):
    outresults = []
    for gene in annotation_dicts[0][chrn]:
        # gene
        gene_region = annotation_dicts[1][gene]
        gene_obj = anno_list_object.get_GI(gene_region)
        if gene_obj:
            gene_cov = gene_obj.get_cov()
            gene_depth = gene_obj.get_depth()
        else:
            gene_cov = 0
            gene_depth = 0

        # transcript
        for transcript in annotation_dicts[0][chrn][gene]:
            if not annotation_dicts[0][chrn][gene][transcript]:
                continue
            transcript_region = annotation_dicts[1][transcript]
            if transcript_region == gene_region:
                transcript_cov = gene_cov
                transcript_depth = gene_depth
            else:
                transcript_obj = anno_list_object.get_GI(transcript_region)
                if transcript_obj:
                    transcript_cov = transcript_obj.get_cov()
                    transcript_depth = transcript_obj.get_depth()
                else:
                    transcript_cov = 0
                    transcript_depth = 0
            transcript_ele_sumcov = 0
            transcript_ele_sumdepth = 0
            transcript_ele_sumlength = 0
            transcript_ele_region = []
            transcript_ele = []
            # elements
            for ele in annotation_dicts[0][chrn][gene][transcript]:
                level = used_region
                transcript_ele += [ele]
                element_region = annotation_dicts[1][ele]
                element_obj = anno_list_object.get_GI(element_region)
                transcript_ele_region += [element_region]
                transcript_ele_sumlength += element_obj.length
                transcript_ele_sumcov += element_obj.sumcov
                transcript_ele_sumdepth += element_obj.sumdepth
                # if element_region == (99171, 99793):
                #    print(element_obj)
                if element_region == transcript_region:
                    element_cov = transcript_cov
                    element_depth = transcript_depth
                else:
                    if element_obj:
                        element_cov = element_obj.get_cov()
                        element_depth = element_obj.get_depth()
                    else:
                        element_cov = 0
                        element_depth = 0
                outresults += ['\t'.join([str(x) for x in [sample_tag,
                                                           chrn,
                                                           level,
                                                           gene,
                                                           gene_region,
                                                           round(gene_cov, 4),
                                                           round(gene_depth, 4),
                                                           transcript,
                                                           transcript_region,
                                                           round(transcript_cov, 4),
                                                           round(transcript_depth, 4),
                                                           ele,
                                                           element_region,
                                                           round(element_cov, 4),
                                                           round(element_depth, 4)
                                                           ]])]
            # sum at transcript level
            level = 'mRNA'
            if transcript_ele_sumlength != 0:
                transcript_ele_cov = transcript_ele_sumcov / transcript_ele_sumlength
                transcript_ele_depth = transcript_ele_sumdepth / transcript_ele_sumlength
            else:
                transcript_ele_cov = 0
                transcript_ele_depth = 0
            if transcript_ele:
                outresults += ['\t'.join([str(x) for x in [sample_tag,
                                                           chrn,
                                                           level,
                                                           gene,
                                                           gene_region,
                                                           round(gene_cov, 4),
                                                           round(gene_depth, 4),
                                                           transcript,
                                                           transcript_region,
                                                           round(transcript_cov, 4),
                                                           round(transcript_depth, 4),
                                                           transcript_ele,
                                                           transcript_ele_region,
                                                           round(transcript_ele_cov, 4),
                                                           round(transcript_ele_depth, 4)
                                                           ]])]
    return outresults


def scan_bed(bedfile, annotation_dicts, output, used_region, comd, sample_tag='', at_least_depth=1):
    current_chrn = ''
    current_anno_idx = 0
    visited_chrns = set()
    if os.path.isfile(output):
        logging.warning("# {output} exists!".format(output=output))
        with open(output) as f:
            for line in f:
                if line.startswith('#') or line.startswith('\n'):
                    continue
                temp = line.split('\t')
                if len(temp) >= 2:
                    visited_chrns.add(line.split('\t')[1])
        logging.warning("# Skip {n} chromosomes!".format(n=len(visited_chrns)))
        fout = open(output, 'a')
    else:
        fout = open(output, 'w')
        # write '#', command & colnames
        header_cols = ["# Sample", "Chr", "Level", "Gene", "Gene_Region", "Gene_Cov", "Gene_Depth",
                       "mRNA", "mRNA_Region", "mRNA_Cov", "mRNA_Depth",
                       "Element", "Element_Region", "Element_Cov", "Element_Depth"]
        fout.write('#' + comd + '\n' + '\t'.join(header_cols) + '\n')
    with open(bedfile) as f:
        for line in f:
            temp = line.rstrip().split('\t')
            # skip chromosome without annotations
            if temp[0] not in annotation_dicts[0]:
                continue
            # skip visited chromosome
            if temp[0] in visited_chrns:
                continue
            if temp[0] != current_chrn:
                # get depth, cov of last chromosome
                if current_chrn != '':
                    temp_results = compute_cov(annotation_dicts, anno_list_obj, used_region, current_chrn, sample_tag)
                    fout.write('\n'.join([str(x) for x in temp_results]) + '\n')
                # init for a new chromosome
                if temp[0] in annotation_dicts[2]:
                    anno_list_obj = GIntervalList(annotation_dicts[2][temp[0]], min_depth=at_least_depth)
                    anno_list_obj.sort()
                    current_anno_idx = 0
            # bed format, 0-based, right open
            current_chrn = temp[0]
            # [0-based, right open] to [1-based, right closed]
            # update scan pos
            end_pos = int(temp[2])
            if end_pos < anno_list_obj.intervals[current_anno_idx].lower_bound:
                continue

            depth = int(float(temp[3]))  # int error if 1.1e6
            # ignore 0 depth record
            if depth == 0:
                continue

            # no any scan in current chromosome
            start_pos = int(temp[1]) + 1
            if current_anno_idx == anno_list_obj.count - 1 and \
                    start_pos > anno_list_obj.intervals[current_anno_idx].upper_bound:
                continue

            temp_scan = GInterval((start_pos, end_pos),
                                  sdepth=depth,
                                  min_depth=at_least_depth)
            # update anno pos
            while start_pos > anno_list_obj.intervals[current_anno_idx].upper_bound:
                if current_anno_idx >= anno_list_obj.count - 1:
                    break
                current_anno_idx += 1

            # add sum
            for i in range(current_anno_idx, anno_list_obj.count):
                add_flag = anno_list_obj.intervals[i].update_sum(temp_scan)
                if add_flag == 0 and anno_list_obj.intervals[i].lower_bound > end_pos:
                    break

        # final one, compute
        if current_chrn != '':
            temp_results = compute_cov(annotation_dicts, anno_list_obj, used_region, current_chrn, sample_tag)
            fout.write('\n'.join([str(x) for x in temp_results]) + '\n')
    fout.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Compute gene coverage and depth''')
    parser.add_argument('-a', '--annotation', metavar='<input.gff/gtf>', help='Path of input gff/gtf', type=str,
                        required=True)
    parser.add_argument('-b', '--bed', metavar='<input.bed>', help='bed of of coverage from bedtools', type=str,
                        required=True)
    parser.add_argument('-r', '--region', metavar='<str>', help='CDS or exon (default: CDS)', type=str,
                        choices=['CDS', 'exon'], default='CDS')
    parser.add_argument('-o', '--output', metavar='<output.cov>', help='Path of output cov', type=str, required=True)
    parser.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
    parser.add_argument('-m', '--min_depth', metavar='<int>', help='Min depth (default: 1)', type=int, default=1)

    args = vars(parser.parse_args())
    annotation_path = os.path.abspath(args['annotation'])
    output_path = os.path.abspath(args['output'])
    bed_path = os.path.abspath(args['bed'])
    sample_tag = args['sample_name']
    used_region = args['region']
    min_depth = args['min_depth']
    # command
    logging.info('# ' + ' '.join(sys.argv))

    if annotation_path.endswith('gtf'):
        annotation = readgtf(annotation_path, used_region)
    elif annotation_path.endswith('gff') or annotation_path.endswith('gff3'):
        annotation = readgff(annotation_path, used_region)
    else:
        annotation = ''
        logging.error("# No gff/gff3/gtf!")
        exit()

    chrn_n = 0
    gene_n = 0
    trans_n = 0
    ele_n = 0
    chrn_used_n = 0
    gene_used_n = 0
    trans_used_n = 0
    chrn_used = 0
    gene_used = 0
    for chrn in annotation[0]:
        chrn_n += 1
        chrn_used = 0
        for gene in annotation[0][chrn]:
            gene_n += 1
            gene_used = 0
            for transcript in annotation[0][chrn][gene]:
                trans_n += 1
                temp_elen = len(annotation[0][chrn][gene][transcript])
                if temp_elen:
                    trans_used_n += 1
                    chrn_used = 1
                    gene_used = 1
                    ele_n += temp_elen
            if gene_used:
                gene_used_n += 1
        if chrn_used:
            chrn_used_n += 1

    logging.info('# Load {chrn_n} chromosome(s), {gene_n} gene(s), {tran_n} transcript(s), {ele_n} {region}(s).'.format(
        chrn_n=chrn_n,
        gene_n=gene_n,
        tran_n=trans_n,
        ele_n=ele_n,
        region=used_region
    ))
    if trans_n != trans_used_n:
        logging.warning('# Drop Chromosome(s)/gene(s)/transcript(s) without {region}(s)!'.format(
            region=used_region
        ))
    logging.info('# Use {chrn_n} chromosome(s), {gene_n} gene(s), {tran_n} transcript(s), {ele_n} {region}(s).'.format(
        chrn_n=chrn_used_n,
        gene_n=gene_used_n,
        tran_n=trans_used_n,
        ele_n=ele_n,
        region=used_region
    ))
    scan_bed(bed_path, annotation, output_path, used_region, ' '.join(sys.argv), sample_tag, min_depth)
    logging.info('# Finish computing coverage and depth.')
