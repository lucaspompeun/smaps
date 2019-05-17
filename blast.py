from Bio import SeqIO
from Bio.Blast import NCBIWWW
import sys, os

# blast de nucleotideos com o NCBI para um contig
def blastn(filename):
    query = SeqIO.read(filename, format="fasta") 
    blast = NCBIWWW.qblast("blastn", "nt", query.seq)
    blastFile = open("blastn.xml", "w") 
    blastFile.write(blast.read())
    blastFile.close()
    blast.close()

    return


# blast de nucleotideos com o NCBI para varios contigs
def multiblastn(filename):
    for record in SeqIO.parse(filename, "fasta"):
        blast = NCBIWWW.qblast('blastn', 'nt', record.seq)
        blastFile = open(record.id + '-blastn.xml', 'w')
        blastFile.write(blast.read())
        blastFile.close()
        blast.close()

    return


# blast local de 2 sequencias que tem o tabular 6 no output
def blastnlocal6(query, subject):
    out = 'blastn.outfmt6'
    cmd = 'blastn -query ' + str(query) + ' -subject ' + str(subject) + ' -outfmt 6 > ' + out
    os.system(cmd)

    return


#### funcao para buscar as sequencias pelo entrez no NCBI
