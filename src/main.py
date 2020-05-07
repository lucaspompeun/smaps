import time
import os


global localPath
localPath = os.path.dirname(os.path.realpath(__file__))


def WriteTimeLog(state, startTime, output):
    elapsedTime = time.time() - startTime
    finalTime = time.strftime("%H:%M:%S", time.gmtime(elapsedTime))

    with open(output + '/log_time.txt', 'a+') as log:
        log.write(state + ': ' + finalTime + '\n')


def WriteLog(output, data):
    with open(output + '/log.txt', 'a+') as log:
        log.write(data + '\n')


def Spades(read1, output, read2=None):
    timeSpades = time.time()
    WriteTimeLog('Spades - Start: ', timeSpades, output)

    out = output + '/spades_default'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd = localPath + '/softwares/spades/bin/spades.py -o ' + out

    if read2:
        cmd += ' -1 ' + read1 + ' -2 ' + read2
    else:
        cmd += ' -s ' + read1

    WriteLog(output, 'CLI Spades: ' + cmd)
    os.system(cmd)

    WriteTimeLog('Spades - End: ', timeSpades, output)

    return out + '/contigs.fasta'


def UnmappedAssembly(read1, read2, output):
    timeUnmappedAssembly = time.time()
    WriteTimeLog('Unmapped Assembly - Start: ', timeUnmappedAssembly, output)

    out = output + '/spades_unmapped'
    if not os.path.exists(out):
        os.mkdir(out)

    cmd = localPath + '/softwares/spades/bin/spades.py -o ' + \
        out + ' -1 ' + read1 + ' -2 ' + read2
    WriteLog(output, 'CLI Unmapped Reads Assembly: ' + cmd)
    os.system(cmd)

    WriteTimeLog('Unmapped Assembly - End: ', timeUnmappedAssembly, output)

    return out + '/contigs.fasta'


def bowtie2(read1, reference, output, read2=None):
    timeBowtie = time.time()
    WriteTimeLog('Bowtie2 - Start: ', timeBowtie, output)

    out = output + '/bowtie'
    database = out + '/database'
    if not os.path.exists(out):
        os.mkdir(out)

    databaseCommandLine = 'bowtie2-build ' + reference + ' ' + database
    WriteLog(output, 'CLI Database Bowtie2: ' + databaseCommandLine)
    os.system(databaseCommandLine)

    cmd = 'bowtie2 -p 16 -x ' + database
    if read2:
        cmd += ' -1 ' + read1 + ' -2 ' + read2
    else:
        cmd += ' -U ' + read1

    cmd += ' -S ' + out + '/output.sam'

    WriteLog(output, 'CLI Bowtie2: ' + cmd)

    os.system(cmd)

    WriteTimeLog('Bowtie2 - End: ', timeBowtie, output)

    return out + '/output.sam'


def Samtools(samFile, output):
    timeSamtools = time.time()
    WriteTimeLog('Samtools - Start: ', timeSamtools, output)

    out = output + '/samtools'
    if not os.path.exists(out):
        os.mkdir(out)

    samtools = localPath + '/bin/samtools '

    cmdView = samtools + 'view -Sb ' + samFile + ' > ' + out + '/output.bam'
    WriteLog(output, 'CLI Samtools View: ' + cmdView)
    os.system(cmdView)

    cmdSort = samtools + 'sort ' + out + \
        '/output.bam -o ' + out + '/output_sorted.bam'
    WriteLog(output, "CLI Samtools Sort: " + cmdSort)
    os.system(cmdSort)

    cmdIndex = samtools + 'index ' + out + '/output_sorted.bam'
    WriteLog(output, 'CLI Samtools Index: ' + cmdIndex)
    os.system(cmdIndex)

    WriteTimeLog('Samtools - End: ', timeSamtools, output)

    return out + '/output_sorted.bam'


def UnmappedReads(bamFile, output):
    """
    The return of this function returns two files of unmapped reads, you have to save both variables in distincts variables.
    """

    timeUnmappedReads = time.time()
    WriteTimeLog('Unmapped Reads - Start: ', timeUnmappedReads, output)

    out = output + '/unmapped_reads'
    if not os.path.exists(out):
        os.mkdir(out)

    unmappedSam = out + '/unmapped.sam'
    unmappedBam = out + '/unmapped.bam'

    samtools = localPath + '/bin/samtools '

    unmappedView = samtools + ' view -f4 ' + bamFile + ' > ' + unmappedSam
    WriteLog(output, "CLI Unmapped Reads View: " + unmappedView)
    os.system(unmappedView)

    unmappedViews = samtools + 'view -Sb ' + unmappedSam + ' > ' + unmappedBam
    WriteLog(output, "CLI Unmapped Reads Views: " + unmappedViews)
    os.system(unmappedViews)

    samToFastq = 'java -jar ' + localPath + '/bin/SamToFastq.jar I=' + unmappedBam + ' F=' + out + \
        '/unmapped_read_1.fastq F2=' + out + \
        '/unmapped_read_2.fastq FU=' + out + '/unmapped_unpaired.fastq'
    WriteLog(output, 'CLI SamToFastq: ' + samToFastq)
    os.system(samToFastq)

    unmappedFile = open(out + '/unmapped_read_1.fastq', 'r')
    afile = unmappedFile.read()
    if len(afile) == 0:
        unmappedSam = out + '/unmapped.sam'
        unmappedHeader = out + '/unmapped.header'
        unmappedHeaderSam = out + '/unmapped_header.sam'
        unmappedHeaderBam = out + '/unmapped_header.bam'
        samtools = localPath + '/bin/samtools '

        cmdGetHeader = samtools + 'view -H ' + bamFile + ' > ' + unmappedHeader
        os.system(cmdGetHeader)

        cmdAddHeader = 'cat ' + unmappedHeader + ' ' + \
            unmappedSam + ' > ' + unmappedHeaderSam
        os.system(cmdAddHeader)

        cmdViews = samtools + 'view -Sb ' + unmappedHeaderSam + ' > ' + unmappedHeaderBam
        os.system(cmdViews)

        cmdSamToFastq = 'java -jar ' + localPath + '/bin/SamToFastq.jar I=' + unmappedHeaderBam + ' F=' + \
            out + '/unmapped_read_1.fastq F2=' + out + \
            '/unmapped_read_2.fastq FU=' + out + '/unmapped_unpaired.fastq'
        os.system(cmdSamToFastq)

    WriteTimeLog('Unmapped Reads - End: ', timeUnmappedReads, output)

    return out + '/unmapped_read_1.fastq', out + '/unmapped_read_2.fastq'


def WriteLibrary(filename, data):
    with open(filename, 'w') as lib:
        lib.write(data)


def Sspace(output, contig, read1, minReadsExtension, read2=None):
    timeSspace = time.time()
    WriteTimeLog('Sspace - Start: ', timeSspace, output)

    out = output + '/sspace'
    # the sspace create the output folder

    if read2:
        libData = 'Lib1 bowtie ' + read1 + ' ' + read2 + ' 400 0.25 FR'
        WriteLibrary(output + '/library.txt', libData)
    else:
        libData = 'unpaired bowtie ' + read1
        WriteLibrary(output + '/library.txt', libData)

    # depois verificar se adiciono o -b que é o nome final do arquivo
    cmdSspace = localPath + '/softwares/sspace/SSPACE.pl -l ' + output + \
        '/library.txt -s ' + contig + ' -x 1 -o ' + \
        str(minReadsExtension) + ' -T 8 -b sspace'
    WriteLog(output, 'CLI Sspace: ' + cmdSspace)
    os.system('cd ' + output + ' && ' + cmdSspace)

    WriteTimeLog('Sspace - End: ', timeSspace, output)

    return out + '/sspace.final.scaffolds.fasta'


def Awk(query, target):
    contigs = [query, target]
    for contig in contigs:
        os.system("awk '/^>/{print " + '">Contig0."' + " ++i; next}{print}' < " +
                  contig + " > " + contig + ".mod | mv " + contig + ".mod " + contig)


def Gaa(query, target, output):
    """
    Query: unmapped contig
    Target: extended contig
    """

    timeGaa = time.time()
    WriteTimeLog('GAA - Start: ', timeGaa, output)

    out = output + '/gaa'
    if not os.path.exists(out):
        os.mkdir(out)

    cmdGaa = 'perl ' + localPath + '/softwares/gaa/gaa.pl -t ' + \
        target + ' -q ' + query + ' -o ' + out
    WriteLog(output, 'CLI GAA: ' + cmdGaa)
    os.system(cmdGaa)

    WriteTimeLog('GAA - End: ', timeGaa, output)

    return out + '/*.fa'


def Prokka(contig, output):
    timeProkka = time.time()
    WriteTimeLog('Prokka - Start: ', timeProkka, output)

    out = output + '/prokka'
    if not os.path.exists(out):
        os.mkdir(out)

    cmdProkka = 'prokka --outdir ' + out + ' --prefix smaps_annotation ' + contig
    WriteLog(output, 'CLI Prokka: ' + cmdProkka)
    os.system(cmdProkka)


def Quast(output, contigList, reference=None, gff=None):
    timeQuast = time.time()
    WriteTimeLog('Quast - Start: ', timeQuast, output)

    out = output + '/quast'
    if not os.path.exists(out):
        os.mkdir(out)

    cmdQuast = 'quast -o ' + out + ' '
    for contig in contigList:
        cmdQuast += contig + ' '

    if reference:
        cmdQuast + '-r ' + reference

    if gff:
        cmdQuast += ' -G ' + gff

    WriteLog(output, 'CLI Quast: ' + cmdQuast)
    os.system(cmdQuast)

    WriteTimeLog('Quast - End: ', timeQuast, output)


def Main(read1, output, read2=None, gff=None, reference=None, sspace=5, minReadsExtension=5):
    timeMain = time.time()

    if not os.path.exists(output):
        os.mkdir(output)

    WriteTimeLog('Smaps - Start: ', timeMain, output)

    print('\n\n\nThanks for using Smaps, please cite us.\n\n\n')
    time.sleep(5)

    if not os.path.exists(output):
        os.mkdir(output)

    # Default assembly and alignment
    if read2:
        contigSpades = Spades(read1, output, read2)
        samFile = bowtie2(read1, contigSpades, output, read2)
    else:
        contigSpades = Spades(read1, output)
        samFile = bowtie2(read1, contigSpades, output)

    sortedBamFile = Samtools(samFile, output)

    unmappedRead1, unmappedRead2 = UnmappedReads(sortedBamFile, output)

    # Scaffolding contigs
    if read2:
        preContigSspace = Sspace(output, contigSpades,
                                 read1, minReadsExtension, read2)

        for _ in range(int(sspace)):
            preContigSspace = Sspace(
                output, preContigSspace, read1, minReadsExtension, read2)

        extendedContig = Sspace(output, preContigSspace,
                                read1, minReadsExtension, read2)

    else:
        preContigSspace = Sspace(output, contigSpades,
                                 read1, minReadsExtension)

        for _ in range(int(sspace)):
            preContigSspace = Sspace(
                output, preContigSspace, read1, minReadsExtension)

        extendedContig = Sspace(output, preContigSspace,
                                read1, minReadsExtension)

    # Assemblying unmapped reads
    unmappedContig = UnmappedAssembly(unmappedRead1, unmappedRead2, output)

    # Graph accordance assembly
    Awk(unmappedContig, extendedContig)

    gaaContig = Gaa(unmappedContig, extendedContig, output)

    # Characterization and quality measure of contigs
    Prokka(gaaContig, output)

    contigList = [contigSpades, gaaContig]

    if reference and gff:
        Quast(output, contigList, reference, gff)
    else:
        if reference:
            Quast(output, contigList, reference)
        else:
            Quast(output, contigList, gff=gff)

    WriteTimeLog('Smaps - End: ', timeMain, output)
