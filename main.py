from Bio import Entrez, SeqIO
import os


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
    else:
        cmd += '-s ' + read1

    if trusted_contigs:
        cmd += ' --trusted-contigs ' + trusted_contigs + ' '

    if untrusted_contigs:
        cmd += ' --untrusted-contigs ' + untrusted_contigs + ' '

    cmd += ' > ' + out + 'execution.log'

    write_file(out + 'commandline.txt', cmd) # write execution log
    os.system(cmd)

    return out + 'contigs.fasta'


# mapping sequences with bowtie2
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

