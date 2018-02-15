#!/usr/bin/env python3
# -*- coding: utf-8 -*-

""" spike.py: Abstracts some of the logging functions called in cactus.py and spike.py """
from time import gmtime, strftime

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m' 

def write(msg, msgType = bcolors.HEADER):
    print(msgType + str('[') + strftime("%Y-%m-%d %H:%M:%S") + str('] ') + str(msg) + bcolors.ENDC)
