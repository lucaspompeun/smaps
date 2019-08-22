# -*- coding: utf-8 -*-

import os
import time
import sys


global path
path = os.path.dirname(os.path.realpath(__file__))


def write_file(filename, data, mode="w"):
    with open(filename, mode) as out:
        out.write(data)


def gcbias(read1, read2, project):
    out = project + 'gcbias/'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd = 'gcbias read1 ' + read1 + ' read2 ' + read2 + ' project ' + project + 'gcbias' + ' outcontig ' + path + '/' + out + '/contigs_gcbias.fasta'
    os.system(cmd)

    return out + '/contigs_gcbias.fasta'


def prokka(filename, project):
    out = project + 'prokka/'

    prokka = 'prokka --outdir ' + out + ' --prefix ' + project + ' ' + filename
    os.system(prokka)



def bowtie2(read1, reference, project, read2 = None, N = '1', L = '22', threads = '16'):
    out = project + 'bowtie/'
    database = out + 'database'
    if not os.path.exists(out):
        os.mkdir(out)

    print('bowtie2-build ' + reference + ' ' + database)

    cmd = 'bowtie2 -p ' + threads + ' -x ' + database

    if read2:
        cmd += ' -1 ' + read1 + ' -2 ' + read2
    else:
        cmd += ' -U ' + read1

    cmd += ' -N ' + N + ' -L ' + L + ' -S ' + out + 'output.sam'

    write_file(out + 'commandline.txt', cmd)

    print(cmd)

    return out + 'output.sam'


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
    

def unmappedreads(bamfile, project):
    out = project + 'unmappedreads/'
    if not os.path.exists(out):
        os.mkdir(out)

    unmapped_sam = out + 'unmapped.sam'
    unmapped_bam = out + 'unmapped.bam'

    unmapped ='samtools view -f4 ' + bamfile + ' > ' + unmapped_sam
    views = 'samtools view -Sb ' + unmapped_sam + ' > ' + unmapped_bam
    sam_to_fastq = 'java -jar SamToFastq.jar I=' + unmapped_bam + ' F=' + out + 'unmapped_read_1.fastq F2=' + out + 'unmapped_read_2.fastq FU=' + out + 'unmapped_unpaired.fastq 2>&1 | tee ' + out + 'log_1.txt'
    
    os.system(unmapped)
    os.system(views)
    os.system(sam_to_fastq)

    x = open(out + 'unmapped_unpaired.fastq', 'r')
    x = x.read()
    if len(x) == 0:
        unmapped_sam = out + 'unmapped.sam'
        unmapped_bam = out + 'unmapped.bam'
        unmapped_header = out + 'unmapped.header'
        unmapped_header_sam = out + 'unmapped_header.sam'
        unmapped_header_bam = out + 'unmapped_header.bam'

        get_header = 'samtools view -H ' + bamfile+ ' > ' + unmapped_header
        add_header = 'cat ' + unmapped_header + ' ' + unmapped_sam + ' > ' + unmapped_header_sam
        view = views = 'samtools view -Sb ' + unmapped_header_sam + ' > ' + unmapped_header_bam
        sam_to_fastq = 'java -jar SamToFastq.jar I=' + unmapped_header_bam + ' F=' + out +'unmapped_read_1.fastq F2=' + out +'unmapped_read_2.fastq FU=' + out + 'unmapped_unpaired.fastq 2>&1 | tee ' + out + 'log_2.txt'

        os.system(get_header)
        os.system(add_header)
        os.system(view)
        os.system(sam_to_fastq)

def get_unmapped_fastq(project, value):
    folder = project + 'unmappedreads/'

    return folder + 'unmapped_read_' + str(value) + '.fastq'


def sspace(project, contigs, fastq1, fastq2, o = 5):
    out = project + 'sspace/'

    data = 'Lib1 bowtie ' + fastq1 + ' ' + fastq2 + ' 400 0.25 FR'
    write_file('library.txt', data)

    sspace = path + 'perl /SSPACE/SSPACE.pl -l library.txt -s ' + contigs + ' -x 1 -o ' + str(o) + ' -T 8 -p 1 -b ' + out
    os.system(sspace)

    return out + 'scaffold.fasta'


def quast(contig_list, ref_fasta, project):
	out = 'projects/' + project + 'quast/'
	if not os.path.exists(out):
		os.mkdir(out)
	cmd = "quast.py" + " -o " + out + " "
	for contig in contig_list:
		if contig:
			cmd += contig + ' '

	cmd += " > " + out + "quast.log"

	os.system(cmd)

	return out


def main(read1, read2, project):

    contigs_gcbias = gcbias(read1, read2, project)

    sam_file = bowtie2(read1, contigs_gcbias, project, read2)
    sorted_bam = samtools(sam_file, project)

    unmappedreads(sorted_bam, project)
    unmapped_fastq1 = get_unmapped_fastq(project, 1)
    unmapped_fastq2 = get_unmapped_fastq(project, 2)

    scaffolds_fasta = sspace(project, contigs_gcbias, unmapped_fastq1, unmapped_fastq2)

    prokka(scaffolds_fasta, project)

main(sys.argv[1], sys.argv[2], sys.argv[3])


