#!/usr/bin/env python3

import argparse

from src.main import *


parser = argparse.ArgumentParser(
    description='Smaps A tool to extends contigs to reduce gaps with unmapped reads')

parser.add_argument('-read1', action='store', dest='read1',
                    required=True, help='first fastq file (required)')
parser.add_argument('-read2', action='store', dest='read2',
                    required=False, default=None, help='second fastq file (optional)')
parser.add_argument('-output', action='store', dest='output', required=True,
                    help='output for results files (please, give the full path for the output)')
parser.add_argument('-gff', action='store', dest='gff',
                    required=False, default=None, help='gff file (optional)')
parser.add_argument('-reference', action='store', dest='reference',
                    required=False, default=None, help='reference file (optional)')
parser.add_argument('-sspace', action='store', dest='sspace', required=False,
                    default=5, help='number of runs of sspace software (optional, default=5)')
parser.add_argument('-minreads', action='store', dest='minreads', required=False, default=5,
                    help='minimum number of unmapped reads to extend a base (optional, default=5)')
parser.add_argument('--version', action='version', version='Smaps 1.1')

results = parser.parse_args()

read1 = results.read1
read2 = results.read2
output = results.output
gff = results.gff
reference = results.reference
sspace = results.sspace
minreads = results.minreads

Main(read1, output, read2, gff, reference, sspace, minreads)