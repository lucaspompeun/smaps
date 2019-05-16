from Bio import SeqIO
from Bio.Blast import NCBIWWW
import sys, os

# blast de nucleotideos com o NCBI
def blastn(filename):
    query = SeqIO.read(filename, format="fasta") 
    blast = NCBIWWW.qblast("blastn", "nt", query.seq)
    blastFile = open("my_blast.xml", "w") 
    blastFile.write(blast.read())
    blastFile.close()
    blast.close()

def blastx(query):
    cmd = 'blastx -query ' + str(query) + '-db nt -out result.xml -evalue 0.001 -outfmt 6'
    os.system(cmd)




#blastn(sys.argv[1])
blastx(sys.argv[1])

