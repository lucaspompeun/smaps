from Bio import Entrez, SeqIO
import os

# read fastQ
def readFastq(filename):
    with open(filename) as fastq:
        