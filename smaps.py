import argparse
import main

parser = argparse.ArgumentParser(description = 'Smaps - A tool to extends contigs to reduce gaps with unmapped reads')

parser.add_argument('-read1', action='store', dest='read1', required = True, help='first fastq file (required)')
parser.add_argument('-read2', action='store', dest='read2', required = True, help='second fastq file (required)')
parser.add_argument('-project', action='store', dest='project', required = True, help='name of project (required)')
parser.add_argument('-o', action='store', dest='o', type = int, default = 5, required = False, help='minimum number of reads needed to extend a base (optional)')
parser.add_argument('--version', action='version', version='Smaps 1.0')
results = parser.parse_args()

read1 = results.read1
read2 = results.read2
project = results.project
o = results.o

main(read1, read2, project, o)
