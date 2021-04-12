#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2021/4/8 10:07 PM
    @Usage:
"""
import logging

logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s %(filename)s %(levelname)s] %(message)s',
                    datefmt='%a %d %b %Y %H:%M:%S',
                    # filename='{subcmd}.log'.format(subcmd=sys.argv[0].split('/')[-1].split('.')[0]),
                    # filemode='w'
                    )
