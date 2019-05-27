from Bio import Entrez, SeqIO
import os


# global variable to directory
global path
path = os.path.dirname(os.path.realpath(__file__))


# write out
def write_file(filename, data, mode="w"):
    with open(filename, mode) as out:
        out.write(data)


# assembly with spades
def spades(read1, project, prefix, read2 = None, trusted_contigs = None, threads = '16', untrusted_contigs = None):
    out = project + prefix + 'spades/'
    if not os.path.exists(out):
        os.mkdir(out)

    cdm = path + '/spades/bin/spades.py -o ' + out + ' -t ' + threads + ' '

    if read2:
        cmd += '-1 ' + read1 + ' -2 ' + read2 + ' '
        if S:
            cmd += '-s ' + read1
    else:prokka --kindgom Bacteria --outdir prokka_SRR1424625 --genus Escherichia --locustag SRR1424625 /home/lucas/pipeline_teste1/SRR1424625/montagem_default/contigs.fasta
        cmd += '-s ' + read1

    if trusted_contigs:
        cmd += ' --trusted-contigs ' + trusted_contigs + ' '

    if untrusted_contigs:
        cmd += ' --untrusted-contigs ' + untrusted_contigs + ' '

    cmd += ' > ' + out + 'execution.log'

    write_file(out + 'commandline.txt', cmd) # write execution log
    os.system(cmd)

    return out + 'contigs.fasta'


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

    unmapped = 'samtools view -f4 ' + bamfile + ' > ' + unmapped_sam
    ################################### adicionar o cabeçalho ##########################################
    views = 'samtools view -Sb ' + unmapped_sam + ' > ' + unmapped_bam
    #output = 'samtools fastq ' + unmapped_bam + ' > output.fastq'
    # unmapped bam vai pro sam to fastq

    os.system(path + '/unmappedreads ' + unmapped)
    os.system(path + '/unmappedreads ' + views)
    os.system(path + '/unmappedreads ' + output)
    
    return out + 'unmapped.fastq'


# SSPACES function


# prokka function (prokka --outdir pastasaida --genus genero --locustag valordecabecalho contigs.fasta)
def prokka(filename, project):
    out = project + 'prokka/'

