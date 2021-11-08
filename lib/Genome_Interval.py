#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 10:17 AM
    @Usage:
"""


class GInterval(object):
    # region legal: list, tuple, string
    def __init__(self, region, chrn='', sdepth=0, min_depth=1, base_0=False, right_closed_interval=True):
        if isinstance(region, list) or isinstance(region, tuple):
            self.lower_bound = region[0]
            self.upper_bound = region[1]
        elif isinstance(region, str):
            self.lower_bound = int(float(region[1:-1].split(',')[0].strip()))
            self.upper_bound = int(float(region[1:-1].split(',')[1].strip()))
        else:
            raise Exception("Unacceptable data type", region)
        if self.lower_bound > self.upper_bound:
            self.upper_bound, self.lower_bound = self.lower_bound, self.upper_bound
        if not right_closed_interval:
            self.upper_bound -= 1
        if base_0:
            self.lower_bound += 1
            self.upper_bound += 1
        self.region = [self.lower_bound, self.upper_bound]
        self.length = self.upper_bound - self.lower_bound + 1
        self.chrn = chrn
        self.sumdepth = self.length * sdepth
        if sdepth >= min_depth:
            self.sumcov = self.length
        else:
            self.sumcov = 0

    def __str__(self):
        return str(self.chrn) + ':' + str(self.region) + \
               ', cov=' + str(self.get_cov()) + ', depth=' + str(self.get_depth())

    def getsubseq(self, seq, maxchar_in_one_line=0):
        subseq = seq[self.lower_bound - 1:self.upper_bound]
        if maxchar_in_one_line == 0:
            return subseq
        else:
            return add_newline(subseq, maxchar_in_one_line, endnewline=False)

    def overlap_length(self, gi_object):
        if self.lower_bound > gi_object.upper_bound or self.upper_bound < gi_object.lower_bound:
            return 0
        else:
            return min(self.upper_bound, gi_object.upper_bound) - max(self.lower_bound, gi_object.lower_bound) + 1

    def update_sum(self, gi_object):
        region_overlap_length = self.overlap_length(gi_object)
        if region_overlap_length == 0:
            return 0
        else:
            self.sumdepth += region_overlap_length * gi_object.get_depth()
            if gi_object.sumcov > 0:
                self.sumcov += region_overlap_length
            return 1

    def get_cov(self):
        return self.sumcov / self.length

    def get_depth(self):
        return self.sumdepth / self.length


class GIntervalList(object):
    def __init__(self, regions_string, chrn_string='', depth_tuple=[], min_depth=1,
                 sep_intervals=',', sep_bounds='-', base_0=False, right_closed_interval=True, min_len=1):
        # 1,2|3,5
        coordinates = []
        if isinstance(regions_string, str):
            if regions_string != '':
                coordinates = list([coord for coord in regions_string.split(sep_intervals)])
                coordinates = list([tuple([int(y) for y in x.split(sep_bounds)]) for x in coordinates])
        else:
            coordinates = regions_string
        self.count = len(coordinates)
        intervals = []
        if len(depth_tuple) == self.count:
            for i in range(self.count):
                temp_interval = GInterval(coordinates[i],
                                          chrn=chrn_string,
                                          sdepth=depth_tuple[i],
                                          min_depth=min_depth,
                                          base_0=base_0,
                                          right_closed_interval=right_closed_interval)
                if temp_interval.length >= min_len:
                    intervals += [temp_interval]
        else:
            # using global min depth information
            for i in range(self.count):
                temp_interval = GInterval(coordinates[i],
                                          chrn=chrn_string,
                                          min_depth=min_depth,
                                          base_0=base_0,
                                          right_closed_interval=right_closed_interval)
                if temp_interval.length >= min_len:
                    intervals += [temp_interval]

        self.intervals = intervals
        self.update_count()
        self.sorted = False
        self.sep_intervals = sep_intervals
        self.sep_bounds = sep_bounds

    def __str__(self):
        return self.sep_intervals.join([
            self.sep_bounds.join(
                [str(x.lower_bound), str(x.upper_bound)]
            ) for x in self.intervals])

    def update_count(self):
        self.count = len(self.intervals)

    def sort(self, reverse=False):
        temp = self.intervals
        sorted_temp = sorted(enumerate([x.region for x in temp]), key=lambda y: y[1], reverse=reverse)
        indices = [y[0] for y in sorted_temp]
        self.intervals = [temp[idx] for idx in indices]
        self.sorted = True

    def isempty(self):
        return self.count == 0

    def get_GI(self, region):
        for i in range(self.count):
            if self.intervals[i].lower_bound == region[0] and self.intervals[i].upper_bound == region[1]:
                return self.intervals[i]
        return None


def filter_length(interval_list, min_len=500):
    filtered_string = ''
    for i in range(interval_list.count):
        if interval_list.intervals[i].length >= min_len:
            filtered_string += str(interval_list.intervals[i].lower_bound) + \
                               interval_list.sep_bounds + \
                               str(interval_list.intervals[i].upper_bound) + \
                               interval_list.sep_intervals
    return GIntervalList(filtered_string.rstrip(interval_list.sep_intervals))


def union_interval(interval_list):
    if not interval_list.sorted:
        interval_list.sort()
    if interval_list.count == 0:
        return None
    elif interval_list.count == 1:
        return interval_list
    else:
        union_string = ''
        last_interval = interval_list.intervals[0]
        last_start = last_interval.lower_bound
        last_end = last_interval.upper_bound
        for i in range(1, interval_list.count):
            current_interval = interval_list.intervals[i]
            current_start = current_interval.lower_bound
            current_end = current_interval.upper_bound
            if last_end < current_start:
                # xxx---
                # ----xx
                union_string += str(last_start) + interval_list.sep_bounds + str(last_end) + interval_list.sep_intervals
                last_start = max(last_start, current_start)
            else:
                last_start = min(last_start, current_start)
            last_end = max(last_end, current_end)
        # final one
        union_string += str(last_start) + interval_list.sep_bounds + str(last_end)
        return GIntervalList(union_string)


def union_length(interval_list, overlap=True):
    if not interval_list.sorted:
        interval_list.sort()
    if not overlap:
        return sum([x.length for x in interval_list.intervals])
    count_length = 0
    if interval_list.count == 0:
        return count_length
    for i in range(1, interval_list.count):
        if interval_list.intervals[i - 1].upper_bound > interval_list.intervals[i].upper_bound:
            interval_list.intervals[i] = GInterval(
                [interval_list.intervals[i].lower_bound, interval_list.intervals[i - 1].upper_bound]
            )
        count_length += \
            interval_list.intervals[i - 1].length - \
            max(0, interval_list.intervals[i - 1].upper_bound - interval_list.intervals[i].lower_bound + 1)
    # final one
    count_length += interval_list.intervals[-1].length
    return count_length


def complementary_interval(interval_list, universal_region, min_len=500):
    rest_string = ''
    p_start = universal_region.lower_bound
    if not interval_list.sorted:
        interval_list.sort()
    for interval in interval_list.intervals:
        if universal_region.upper_bound < interval.lower_bound:
            continue
        p_end = interval.lower_bound - 1
        if p_end - min_len >= p_start:
            rest_string += str(p_start) + interval_list.sep_bounds + str(p_end) + interval_list.sep_intervals
        p_start = interval.upper_bound + 1
    # final one
    p_end = universal_region.upper_bound
    if p_end - min_len >= p_start:
        rest_string += str(p_start) + interval_list.sep_bounds + str(p_end)
    return GIntervalList(rest_string)


def intersect_interval(interval_list1, interval_list2):
    i, j = 0, 0
    intersect_string = ''
    while i < interval_list1.count and j < interval_list2.count:
        a1, a2 = interval_list1.intervals[i].lower_bound, interval_list1.intervals[i].upper_bound
        b1, b2 = interval_list2.intervals[j].lower_bound, interval_list2.intervals[j].upper_bound
        if b2 >= a1 and a2 >= b1:
            intersect_string += str(max(a1, b1)) + \
                                interval_list1.sep_bounds + str(min(a2, b2)) + \
                                interval_list1.sep_intervals
        if b2 < a2:
            j += 1
        else:
            i += 1
    return GIntervalList(intersect_string.rstrip(interval_list1.sep_intervals))


# sequence and fasta filter
def add_newline(seq, maxchar_in_one_line=60, endnewline=True):
    newseq = ''
    seqlen = len(seq)
    for i in range(0, seqlen, maxchar_in_one_line):
        newseq += seq[i:min(i + maxchar_in_one_line, seqlen)] + '\n'
    if endnewline:
        return newseq
    else:
        return newseq.rstrip('\n')


def fa_some_record(infa_file, ctg_list_file, outfa_file, exclude=False):
    infa_record_n, ctg_list_recond_n, outfa_record_n = 0, 0, 0
    if infa_file == outfa_file:
        return 0, 0, 0
    # read ctg_list_file
    ctg_list = []
    with open(ctg_list_file) as flist:
        for line in flist:
            if line.startswith('#'):
                continue
            ctg_list += [line.rstrip().split()[0]]
            ctg_list_recond_n += 1
    # read fasta
    write_flag = exclude
    with open(outfa_file, 'w') as fout:
        with open(infa_file) as fin:
            for line in fin:
                if line.startswith('>'):
                    infa_record_n += 1
                    if line[1:].rstrip().split()[0] in ctg_list:
                        if exclude:
                            write_flag = 0
                        else:
                            write_flag = 1
                            outfa_record_n += 1
                    else:
                        if exclude:
                            write_flag = 1
                            outfa_record_n += 1
                        else:
                            write_flag = 0
                if write_flag:
                    fout.write(line)
    return infa_record_n, ctg_list_recond_n, outfa_record_n


def gclust2fa(fasta_file, clust_file, out_fa):
    if fasta_file == out_fa:
        return 0
    representative_ctgs = dict()
    representative_n = 0
    with open(clust_file) as f:
        for line in f:
            if line.startswith('>'):
                representative_n += 1
            else:
                temp = line.rstrip().split()
                if temp[-1] == '*':
                    representative = 1
                else:
                    representative = 0
                if representative == 1:
                    ctgname = temp[2].rstrip('.')
                    representative_ctgs[ctgname] = ''
    out_flag = 0
    with open(out_fa, 'w') as fout:
        with open(fasta_file) as f:
            for line in f:
                if line.startswith('>'):
                    if line.rstrip().split()[0] in representative_ctgs:
                        out_flag = 1
                    else:
                        out_flag = 0
                if out_flag == 1:
                    fout.write(line)
    return representative_n


def string2dict(long_string, sep=';', eq='=', rm_quote=False):
    if rm_quote:
        long_string = long_string.replace('"', '').replace("'", '')
    long_string = long_string.replace('; ', ';')
    out_dict = dict()
    tmp = long_string.rstrip(sep).split(sep)
    for i in tmp:
        key, value = i.split(eq)
        out_dict[key] = value
    return out_dict


# gff
class GFFElement(object):
    def __init__(self, line_string):  # https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
        self.string = line_string  # including '\n'
        # type only support: gene, mRNA/transcript, exon, CDS, five_prime_UTR, three_prime_UTR, start_codon, end_codon
        self.chrn, self.source, self.type, self.start, self.end, \
        self.score, self.strand, self.phase, self.attributes = line_string.rstrip().split('\t')[:9]
        self.details = string2dict(self.attributes)
        # update start, end
        self.start = int(self.start)
        self.end = int(self.end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        # GInterval.__init__(self, [self.start, self.end])

    def find_parentid(self):
        parent_id = ''
        if 'Parent' in self.details:
            parent_id = self.details['Parent']
            if parent_id.find(','):  # many parents, return first one
                return parent_id.split(',')[0]
        return parent_id

    def get_id(self):
        return self.details["ID"]

    def get_length(self):
        return self.end - self.start + 1


class GFFTranscript(GFFElement):
    def __init__(self, transcript_line):
        GFFElement.__init__(self, transcript_line)
        self.children = []
        self.parent = []
        self.cdslength = 0
        self.exonlength = 0

    def add_parent(self, gene):
        self.parent += [gene]

    def add_child(self, element):
        self.children += [element]

    def get_key_length(self, key='CDS'):
        sum_len = 0
        for ele in self.children:
            if ele.type == key:
                sum_len += ele.get_length()
        return sum_len

    def update_length(self):
        for ele in self.children:
            if self.type == 'exon':
                self.exonlength += ele.get_length()
            elif self.type == 'CDS':
                self.CDSlength += ele.get_length()


class GFFGene(GFFElement):
    def __init__(self, gene_line):
        GFFElement.__init__(self, gene_line)
        self.children = []

    def add_child(self, element):
        self.children += [element]


# gtf
class GTFElement(object):
    def __init__(self, line_string):  # https://mblab.wustl.edu/GTF22.html
        self.string = line_string  # including '\n'
        # type only support: CDS, start_codon and stop_codon (required)
        # 5UTR, 3UTR, inter, inter_CNS, intron_CNS and exon (optional)
        self.chrn, self.source, self.type, self.start, self.end, \
        self.score, self.strand, self.frame, self.attributes = line_string.rstrip().split('\t')[:9]
        self.details = string2dict(self.attributes, eq=' ', rm_quote=True)
        # update start, end
        self.start = int(self.start)
        self.end = int(self.end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        self.region = (self.start, self.end)
        # GInterval.__init__(self, [self.start, self.end])

    def find_parentid(self):
        return self.details['transcript_id']

    def get_id(self):
        return self.details["transcript_id"] + ':' + self.type + '_' + str(self.start) + '-' + str(self.end)

    def get_length(self):
        return self.end - self.start + 1


class GTFVirtual(object):
    def __init__(self, gtf_elements):
        self.chrn = gtf_elements[0].chrn
        self.strand = gtf_elements[0].strand
        self.cdslength = 0
        self.exonlength = 0
        # start ,end
        self.start = gtf_elements[0].start
        self.end = gtf_elements[0].end
        for ele in gtf_elements:
            if self.end < ele.end:
                self.end = ele.end
            if self.start > ele.start:
                self.start = ele.start
        self.region = (self.start, self.end)
        self.details = gtf_elements[0].details
        if "gene_id" in gtf_elements[0].details:
            self.details = {"gene_id": gtf_elements[0].details["gene_id"]}
        if "transcript_id" in gtf_elements[0].details:
            self.details["transcript_id"] = gtf_elements[0].details["transcript_id"]
        # self.children = gtf_elements

    def get_length(self):
        return self.end - self.start + 1

# bed
class BEDElement(object):
    def __init__(self, line_string):  # https://m.ensembl.org/info/website/upload/bed.html
        self.string = line_string  # including '\n'
        self.chrn,  self.start, self.end = line_string.rstrip().split('\t')[:3]
        # update start, end
        # [0-based, right open] to [1-based, right closed]
        self.start = int(self.start) + 1
        self.end = int(self.end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        self.region = (self.start, self.end)

