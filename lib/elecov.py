#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/11/8 4:49 PM
    @Usage: bedtools genomecov -bga -split -ibam xx.bam > xx.bed;
            python3 elecov.py -a ele.bed -o xx.cov -n xx -b xx.bed
"""
import argparse
from Genome_Interval import *
from tlog import *
import os
import sys


def readbed(bed):
    region_dict = dict()  # {'chr1':[(start,end)]}
    rep_region_n = 0
    with open(bed) as f:
        for line in f:
            if line.startswith('#') or line.rstrip() == '':
                continue
            ele = BEDElement(line)
            if ele.chrn not in region_dict:
                region_dict[ele.chrn] = set()
            if (ele.start, ele.end) in region_dict[ele.chrn]:
                rep_region_n += 1
                continue
            region_dict[ele.chrn].add((ele.start, ele.end))
        for key in region_dict:
            region_dict[key] = tuple(region_dict[key])
    if rep_region_n:
        logging.warning('# Remove {rep_region_n} same region elements.'.format(rep_region_n=rep_region_n))
    return dict(), dict(), region_dict


def compute_cov(annotation_dicts, anno_list_object, chrn, sample_tag):
    outresults = []
    # elements
    for element_region in annotation_dicts[2][chrn]:
        element_obj = anno_list_object.get_GI(element_region)
        element_cov = element_obj.get_cov()
        element_depth = element_obj.get_depth()
        outresults += ['\t'.join([str(x) for x in [sample_tag,
                                                   chrn,
                                                   element_region,
                                                   round(element_cov, 4),
                                                   round(element_depth, 4)
                                                   ]])]
    return outresults


def scan_bed(bedfile, annotation_dicts, output, sample_tag='', at_least_depth=1):
    current_chrn = ''
    current_anno_idx = 0
    fout = open(output, 'w')
    header_cols = ["# Sample", "Chr", "Element_Region", "Element_Cov", "Element_Depth"]
    fout.write('\t'.join(header_cols) + '\n')
    with open(bedfile) as f:
        for line in f:
            temp = line.rstrip().split('\t')
            # skip chromosome without annotations
            if temp[0] not in annotation_dicts[2]:
                continue
            if temp[0] != current_chrn:
                # get depth, cov of last chromosome
                if current_chrn != '':
                    temp_results = compute_cov(annotation_dicts, anno_list_obj, current_chrn, sample_tag)
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
            temp_results = compute_cov(annotation_dicts, anno_list_obj, current_chrn, sample_tag)
            fout.write('\n'.join([str(x) for x in temp_results]) + '\n')
    fout.close()
    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Compute gene coverage''')
    parser.add_argument('-a', '--annotation', metavar='<anotation.bed>', help='Path of anotation.bed', type=str,
                        required=True)
    parser.add_argument('-b', '--bed', metavar='<input.bed>', help='bed of of coverage from bedtools', type=str,
                        required=True)
    parser.add_argument('-o', '--output', metavar='<output.cov>', help='Path of output cov', type=str, required=True)
    parser.add_argument('-n', '--sample_name', metavar='<str>', help='Name of sample', type=str, required=True)
    parser.add_argument('-m', '--min_depth', metavar='<int>', help='Min depth (default: 1)', type=int, default=1)

    args = vars(parser.parse_args())
    annotation_path = os.path.abspath(args['annotation'])
    output_path = os.path.abspath(args['output'])
    bed_path = os.path.abspath(args['bed'])
    sample_tag = args['sample_name']

    min_depth = args['min_depth']
    # command
    logging.info('# ' + ' '.join(sys.argv))

    if annotation_path.endswith('bed'):
        annotation = readbed(annotation_path)
    else:
        logging.error("# No bed!")
        exit()

    chrn_n = 0
    ele_n = 0
    for chrn in annotation[2]:
        chrn_n += 1
        for ele in annotation[2][chrn]:
            ele_n += 1

    logging.info('# Load {chrn_n} chromosomes, {ele_n} elements.'.format(
        chrn_n=chrn_n,
        ele_n=ele_n
    ))
    scan_bed(bed_path, annotation, output_path, sample_tag, min_depth)
    logging.info('# Finish computing coverage and depth.')
