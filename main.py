from Bio import Entrez, SeqIO
import os
import requests
import time


# insert all gcbias here (output needed: scaffold.fasta, read1.fastq and read2.fastq)


# global variable to directory
global path
path = os.path.dirname(os.path.realpath(__file__))
print(path)

# write out
def write_file(filename, data, mode="w"):
    with open(filename, mode) as out:
        out.write(data)


# G C  B I A S  H E R E (RETURN CONTIGS.FASTA)
"""
# uploading to gcbias
def upload_gcbias(url, read1, read2 = None, reference = None):
    files = {'read1': open(read1, 'rb')}
    if read2:
        files.update({'read2': open(read2, 'rb')})
    if reference:
        files.update({'reference': open(reference, 'rb')})

    r = requests.post(url, files=files)

    return r # 200 é que deu ok
"""


# prokka void function (ffn - nucleotideo(quantidade de genes) e faa - proteina)
def prokka(filename, project):
    out = project + 'prokka/'
    if not os.path.exists(out):
        os.mkdir(out)

    prokka = 'prokka --outdir ' + out + ' --prefix ' + project + ' ' + filename
    os.system(prokka)

def get_prokka(project, file): # file requires the extension from files (.gbk, .ffn, faa)
    folder = project + 'prokka/'

    return folder + '*.' + file


# mapping sequences with bowtie2 (reference should receive 'contigs.fasta' from spades function)
def bowtie2(read1, reference, project, read2 = None, N = '1', L = '22', threads = '16'):
    out = project + 'bowtie/'
    database = out + 'database'
    if not os.path.exists(out):
        os.mkdir(out)

    os.system('bowtie2-build ' + reference + ' ' + database + ' > ' + out + 'database.log')

    cmd = 'bowtie2 -p ' + threads + ' -x ' + database

    if read2:
        cmd += ' - 1 ' + read1 + ' -2 ' + read2
    else:
        cmd += ' -U ' + read1

    cmd += ' -N ' + N + ' -L ' + L + ' -S ' + out + 'output.sam > ' + out + 'execution.log'

    write_file(out + 'commandline.txt', cmd)

    os.system(cmd)

    return out + 'output.sam'


# from sam file to sorted bam file
def samtools(samfile, project):
    out = project + 'samtools/'
    if not os.path.exists(out):
        os.mkdir(out)

    outbam = out + 'output.bam'
    sortedbam = out + 'output_sorted.bam'

    os.system(path + '/samtools view -Sb ' + samfile + ' > ' + outbam)
    os.system(path + '/samtools sort ' + outbam + ' -o ' + sortedbam)
    os.system(path + '/samtools index ' + sortedbam)

    return sortedbam 


# unmapped reads
def unmappedreads(bamfile, project):
    out = project + 'unmappedreads/'
    if not os.path.exists(out):
        os.mkdir(out)

    unmapped_sam = 'unmapped.sam'
    unmapped_bam = 'unmapped.bam'
    unmapped_header = 'unmapped.header'
    unmapped_header_sam  = 'unmapped_header.sam'

    try:
        unmapped = 'samtools view -f4 ' + bamfile + ' > ' + unmapped_sam
        views = 'samtools view -Sb ' + unmapped_sam + ' > ' + unmapped_bam

        os.system(path + '/unmappedreads ' + unmapped)
        os.system(path + '/unmappedreads ' + views)

    except:
        unmapped = 'samtools view -f4 ' + bamfile + ' > ' + unmapped_sam
        get_header = 'samtools view -H ' + unmapped_bam + ' > ' + unmapped_header
        add_header = 'cat ' + unmapped_header + unmapped_sam + ' > ' + unmapped_header_sam
        views = 'samtools view -Sb ' + unmapped_sam + ' > ' + unmapped_bam

        os.system(path + '/unmappedreads ' + unmapped)
        os.system(path + '/unmappedreads ' + get_header)
        os.system(path + '/unmappedreads ' + add_header)
        os.system(path + '/unmappedreads ' + views)

    sam_to_fastq = 'java -jar SamToFastq.jar I=out_with_header.sam F=out_with_header.fastqSamToFastq.jar I=' + unmapped_bam + ' F=unmapped_read_1.fastq F2=unmapped_read_2.fastq FU=unmapped_unpaired.fastq'
    os.system(sam_to_fastq)

def get_unmapped_fastq(project, value):
    folder = project + 'unmappedreads/'

    return folder + 'unmapped_read_' + value + '.fastq'


# SSPACE function
def sspace(project, contigs, fastq):


    return
