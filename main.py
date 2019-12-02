#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Smaps - A web tool to extends contigs to reduce gaps with unmapped reads
Author: Lucas Pompeu
"""

import os
import sys
import time


global path
path = os.path.dirname(os.path.realpath(__file__))


def write_file(filename, data, mode="w"):
    with open(filename, mode) as out:
        out.write(data)


def spades_default(read1, project, read2=None, trusted_contigs=None, S=None, threads='16', untrusted_contigs=None):

    print("\n\n\n<><> SPADES DEFAULT <><>\n\n\n")

    out = 'projects/' + project + '/spades_default/'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd = path + '/spades/bin/spades.py -o ' + out + ' -t ' + threads + ' '
    if read2:
        cmd += '-1 ' + read1 + ' -2 ' + read2 + ' '
        if S:
            cmd += '-s ' + S + ' '
    else:
        cmd += '-s ' + read1

    if trusted_contigs:
        cmd += ' --trusted-contigs ' + trusted_contigs + ' '

    if untrusted_contigs:
        cmd += ' --untrusted-contigs ' + untrusted_contigs + ' '

    os.system(cmd)

    return out + 'contigs.fasta'


def spades_unmapped(read1, project, read2=None, trusted_contigs=None, S=None, threads='16', untrusted_contigs=None):

    print("\n\n\n<><> SPADES UNMMAPED <><>\n\n\n")

    out = 'projects/' + project + '/spades_unmmaped/'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd = path + '/spades/bin/spades.py -o ' + out + ' -t ' + threads + ' '
    if read2:
        cmd += '-1 ' + read1 + ' -2 ' + read2 + ' '
        if S:
            cmd += '-s ' + S + ' '
    else:
        cmd += '-s ' + read1

    if trusted_contigs:
        cmd += ' --trusted-contigs ' + trusted_contigs + ' '

    if untrusted_contigs:
        cmd += ' --untrusted-contigs ' + untrusted_contigs + ' '

    os.system(cmd)

    return out + 'contigs.fasta'


def prokka(filename, project):

    print('\n\n\n<><> PROKKA <><>\n\n\n')

    out = 'projects/' + project + '/prokka/'

    prokka = 'prokka --outdir ' + out + ' --prefix ' + project + ' ' + filename
    os.system(prokka)


def bowtie2(read1, reference, project, read2=None, N='1', L='22', threads='16'):

    print('\n\n\n>>> BOWTIE2 <<<\n\n\n')

    out = 'projects/' + project + '/' + 'bowtie/'
    database = out + 'database'
    if not os.path.exists(out):
        os.mkdir(out)

    os.system('bowtie2-build ' + reference + ' ' + database)

    cmd = 'bowtie2 -p ' + threads + ' -x ' + database

    if read2:
        cmd += ' -1 ' + read1 + ' -2 ' + read2
    else:
        cmd += ' -U ' + read1

    cmd += ' -N ' + N + ' -L ' + L + ' -S ' + out + 'output.sam'

    os.system(cmd)

    return out + 'output.sam'


def samtools(samfile, project):

    print('\n\n\n<><> SAMTOOLS <><>\n\n\n')

    out = 'projects/' + project + '/samtools/'
    if not os.path.exists(out):
        os.mkdir(out)

    outbam = out + 'output.bam'
    sortedbam = out + 'output_sorted.bam'

    os.system(path + '/samtools view -Sb ' + samfile + ' > ' + outbam)
    os.system(path + '/samtools sort ' + outbam + ' -o ' + sortedbam)
    os.system(path + '/samtools index ' + sortedbam)

    return sortedbam


def unmappedreads(bamfile, project):

    print('\n\n\n<><> SAMTOOLS and SAMTOFASTQ <><>\n\n\n')

    out = 'projects/' + project + '/unmappedreads/'
    if not os.path.exists(out):
        os.mkdir(out)

    unmapped_sam = out + 'unmapped.sam'
    unmapped_bam = out + 'unmapped.bam'

    unmapped = 'samtools view -f4 ' + bamfile + ' > ' + unmapped_sam
    views = 'samtools view -Sb ' + unmapped_sam + ' > ' + unmapped_bam
    sam_to_fastq = 'java -jar SamToFastq.jar I=' + unmapped_bam + ' F=' \
        + out + 'unmapped_read_1.fastq F2=' + out + 'unmapped_read_2.fastq FU=' \
        + out + 'unmapped_unpaired.fastq 2>&1 | tee ' + out + 'log_1.txt'

    os.system(unmapped)
    os.system(views)
    os.system(sam_to_fastq)

    x = open(out + 'unmapped_read_1.fastq', 'r')
    x = x.read()
    if len(x) == 0:
        unmapped_sam = out + 'unmapped.sam'
        unmapped_header = out + 'unmapped.header'
        unmapped_header_sam = out + 'unmapped_header.sam'
        unmapped_header_bam = out + 'unmapped_header.bam'

        get_header = 'samtools view -H ' + bamfile + ' > ' + unmapped_header
        add_header = 'cat ' + unmapped_header + ' ' + \
            unmapped_sam + ' > ' + unmapped_header_sam
        view = views = 'samtools view -Sb ' + \
            unmapped_header_sam + ' > ' + unmapped_header_bam
        sam_to_fastq = 'java -jar SamToFastq.jar I=' + unmapped_header_bam + ' F=' + out \
            + 'unmapped_read_1.fastq F2=' + out + 'unmapped_read_2.fastq FU=' + out \
            + 'unmapped_unpaired.fastq'

        os.system(get_header)
        os.system(add_header)
        os.system(view)
        os.system(sam_to_fastq)


def get_unmapped_fastq(project, value):

    print('\n\n\n<><> GETTING UMMAPED FASTQ <><>\n\n\n')

    folder = 'projects/' + project + '/unmappedreads/'

    return folder + 'unmapped_read_' + str(value) + '.fastq'


def sspace(project, contigs, fastq1, o, fastq2=None):
    out = 'projects/' + project + '/sspace'
    out1 = 'projects/' + project

    if fastq2:
        data = 'Lib1 bowtie ' + path + '/' + fastq1 + \
            ' ' + path + '/' + fastq2 + ' 400 0.25 FR'
        write_file('projects/' + project + '/library.txt', data)
    else:
        data = 'unpaired bowtie ' + path + '/' + fastq1
        write_file('projects/' + project + '/library.txt', data)

    sspace = path + '/SSPACE/SSPACE.pl -l ' + path + '/projects/' + project + '/library.txt -s ' \
        + path + '/' + contigs + ' -x 1 -o ' + str(o) + ' -T 8 -p 1 -b sspace'
    os.system('cd ' + out1 + ' && ' + sspace)

    return out + '/' + 'sspace.final.scaffolds.fasta'


def quast(contig_list, project, reference=None, gff=None):
    
    print('\n\n\n<><> QUAST <><>\n\n\n')

    out = 'projects/' + project + '/quast/'
    if not os.path.exists(out):
        os.mkdir(out)
    
    cmd = "quast.py" + " -o " + out + " "
    for contig in contig_list:
        if contig:
            cmd += contig + ' '

    if reference:
        cmd += '-r ' + reference

    if gff:
        cmd += " -G " + gff

    os.system(cmd)

    return out


def gaa(query, target, project):
    """
    Query: unmapped contig
    Target: extended contig
    """
    print('\n\n\n<><> GAA <><>\n\n\n')

    out = 'projects/' + project + '/gaa/'
    if not os.path.exists(out):
        os.mkdir(out)
    
    cmd = '/gaa/gaa.pl -t ' + target + ' -q ' + query + ' -o ' + out
    os.system(cmd)


def smaps(read1, project, o, read2=None, reference=None, gff=None):
    print('\n\n\n<><> SMAPS <><>\n\n\n')
    print('Thanks for using Smaps, please cite us.\n\n')
    time.sleep(10)

    out = 'projects/' + project
    if not os.path.exists(out):
        os.mkdir(out)

    if read2:
        contigs_spades = spades_default(read1, project, read2)
    else:
        contigs_spades = spades_default(read1, project)

    sam_file = bowtie2(read1, contigs_spades, project, read2)
    sorted_bam = samtools(sam_file, project)

    unmappedreads(sorted_bam, project)
    unmapped_fastq1 = get_unmapped_fastq(project, 1)
    unmapped_fastq2 = get_unmapped_fastq(project, 2)

    print('\n\n\n<><> SSPACE <><>\n\n\n')
    if read2:
        extended_contigs = sspace(project, contigs_spades,
                                 read1, o, read2)
    else:
        extended_contigs = sspace(project, contigs_spades,
                                 read1, o)

    
    unmmaped_contigs = spades_unmapped(unmapped_fastq1, project, unmapped_fastq2)

    
    gaa(unmmaped_contigs, extended_contigs, project)

    """prokka(gapblaster_contigs, project)

    contig_list = [contigs_spades, gapblaster_contigs]
    if reference and gff:
        quast(contig_list, project, reference, gff)
    else:
        if reference:
            quast(contig_list, project, reference)
        else:
            quast(contig_list, project, gff)"""
