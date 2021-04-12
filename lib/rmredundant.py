#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/12 3:28 PM
    @Usage: python3 rmredundant.py [options]
"""
import argparse
from tlog import *
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Cluster the sequences and remove redundant sequences''')
    parser.add_argument('-i', '--input', metavar='<input.fa>',
                        help='Path of input fasta', type=str, required=True)
    parser.add_argument('-o', '--output', metavar='<output.fa>',
                        help='Path of output fasta', type=str, default='non_redundant.fa')
    parser.add_argument('-c', '--', metavar='<output.fa>',
                        help='Path of output fasta', type=str, default='non_redundant.fa')
