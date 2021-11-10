#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/11/8 4:49 PM
    @Usage: python3 dat.py xxx.trfoutput.dat > xxx.trf.bed
"""

import sys
import re


with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('Sequence:'):
            chrn = line.split()[1]
        if line == '\n':
            continue
        if re.match(r'^\d', line):
            temp = line.split()
            start = int(temp[0]) - 1
            end = int(temp[1])
            rp_str = temp[13]
            rp_cn = temp[3]
            print('\t'.join([str(x) for x in [chrn, start, end, rp_str, rp_cn]]))
