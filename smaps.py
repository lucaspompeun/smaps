import argparse
import main

# args structure
parser = argparse.ArgumentParser(description = 'Smaps - A tool to extends contigs to reduce gaps with unmapped reads')

parser.add_argument('-read1', action='store', dest='read1', required = True, help='first fastq file (required)')
parser.add_argument('-read2', action='store', dest='read2', required = True, help='second fastq file (required)')
parser.add_argument('-project', action='store', dest='project', required = True, help='name of project (required)')
parser.add_argument('-o', action='store', dest='number', type = int, default = 5, required = False, help='minimum number of reads needed to call a base during an extension (optional)')
parser.add_argument('--version', action='version', version='Smaps 1.0')
results = parser.parse_args()