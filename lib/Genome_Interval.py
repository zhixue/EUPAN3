#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 10:17 AM
    @Usage:
"""


class GInterval(object):
    # region legal: list, tuple, string
    def __init__(self, region, base_0=False, right_closed_interval=True):
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
            self.upper_bound = self.upper_bound - 1
        if base_0:
            self.lower_bound += 1
            self.upper_bound += 1
        self.region = [self.lower_bound, self.upper_bound]
        self.length = self.upper_bound - self.lower_bound + 1

    def __str__(self):
        return str(self.region)

    def getsubseq(self, seq, maxchar_in_one_line=0):
        subseq = seq[self.lower_bound-1:self.upper_bound]
        if maxchar_in_one_line == 0:
            return subseq
        else:
            return add_newline(subseq, maxchar_in_one_line, endnewline=False)


class GIntervalList(object):
    def __init__(self, regions_string, sep_intervals=',', sep_bounds='-', base_0=False, right_closed_interval=True, min_len=2):
        # 1,2|3,5
        coordinates = []
        if regions_string != '':
            coordinates = list([coord for coord in regions_string.split(sep_intervals)])
            coordinates = list([tuple([int(y) for y in x.split(sep_bounds)]) for x in coordinates])
        self.count = len(coordinates)
        intervals = []
        for i in range(self.count):
            temp_interval = GInterval(coordinates[i], base_0, right_closed_interval)
            if temp_interval.length >= min_len:
                intervals += [temp_interval]
        self.intervals = intervals
        self.update_count()
        self.sorted = False
        self.sep_intervals = sep_intervals
        self.sep_bounds = sep_bounds

    def __str__(self):
        return self.sep_intervals.join([self.sep_bounds.join([str(x.lower_bound), str(x.upper_bound)]) for x in self.intervals])

    def update_count(self):
        self.count = len(self.intervals)

    def sort(self, reverse=False):
        temp = self.intervals
        sorted_temp = sorted(enumerate([x.region for x in temp]), reverse=reverse)
        indices = [y[0] for y in sorted_temp]
        self.intervals = [temp[idx] for idx in indices]
        self.sorted = True

    def isempty(self):
        return self.count == 0


def filter_length(interval_list, min_len=500):
    filtered_string = ''
    for i in range(interval_list.count):
        if interval_list.intervals[i].length >= min_len:
            filtered_string += str(interval_list.intervals[i].lower_bound) + interval_list.sep_bounds + str(interval_list.intervals[i].upper_bound) + interval_list.sep_intervals
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
        if interval_list.intervals[i-1].upper_bound > interval_list.intervals[i].upper_bound:
            interval_list.intervals[i] = GInterval([interval_list.intervals[i].lower_bound, interval_list.intervals[i-1].upper_bound])
        count_length += interval_list.intervals[i-1].length - max(0, interval_list.intervals[i-1].upper_bound - interval_list.intervals[i].lower_bound + 1)
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
            intersect_string += str(max(a1, b1)) + interval_list1.sep_bounds + str(min(a2, b2)) + interval_list1.sep_intervals
        if b2 < a2:
            j += 1
        else:
            i += 1
    return GIntervalList(intersect_string.rstrip(interval_list1.sep_intervals))


## sequence and fasta filter

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
    Write_Flag = exclude
    with open(outfa_file, 'w') as fout:
        with open(infa_file) as fin:
            for line in fin:
                if line.startswith('>'):
                    infa_record_n += 1
                    if line[1:].rstrip().split()[0] in ctg_list:
                        if exclude:
                            Write_Flag = 0
                        else:
                            Write_Flag = 1
                            outfa_record_n += 1
                    else:
                        if exclude:
                            Write_Flag = 1
                            outfa_record_n += 1
                        else:
                            Write_Flag = 0
                if Write_Flag:
                    fout.write(line)
    return infa_record_n, ctg_list_recond_n, outfa_record_n