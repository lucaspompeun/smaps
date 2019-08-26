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

unmappedreads('data/output_sorted.bam', 'teste_unmapped')







"""import argparse
import main

# args structure
parser = argparse.ArgumentParser(description = 'Smaps - A tool to extends contigs to reduce gaps with unmapped reads')

parser.add_argument('-read1', action='store', dest='read1', required = True, help='first fastq file (required)')
parser.add_argument('-read2', action='store', dest='read2', required = True, help='second fastq file (required)')
parser.add_argument('-project', action='store', dest='project', required = True, help='name of project (required)')
parser.add_argument('-o', action='store', dest='number', type = int, default = 5, required = False, help='minimum number of reads needed to call a base during an extension (optional)')
parser.add_argument('--version', action='version', version='Smaps 1.0')
results = parser.parse_args()

"""