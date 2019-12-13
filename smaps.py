#!/usr/bin/env python3

import argparse
from main import *

parser = argparse.ArgumentParser(
    description='Smaps - A tool to extends contigs to reduce gaps with unmapped reads')

parser.add_argument('-read1', action='store', dest='read1',
                    required=True, help='first fastq file (required)')
parser.add_argument('-read2', action='store', dest='read2',
                    required=False, default=None, help='second fastq file (optional)')
parser.add_argument('-reference', action='store', dest='reference',
                    required=False, default=None, help='reference genome (optional)')
parser.add_argument('-gff', action='store', dest='gff',
                    required=False, default=None, help='gff (optional)')
parser.add_argument('-project', action='store', dest='project',
                    required=True, help='name of project (required)')
parser.add_argument('-sspace', action='store', dest='sspace', type=int, default=5, required=False,
                    help='minimum number of reads needed to extend a base (optional)')
parser.add_argument('--version', action='version', version='Smaps 1.0')
results = parser.parse_args()

read1 = results.read1
read2 = results.read2
project = results.project
sspace = results.sspace
gff = results.gff
reference = results.reference

smaps(read1, project, sspace, read2, reference, gff)
