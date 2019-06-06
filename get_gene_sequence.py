from Bio import Entrez, SeqIO
Entrez.email = 'Your.Name.Here@example.org'  # alterar o email depois

# lista de genes, com o nome, locus tag e protein id
x = ['Cp258_1020', 'Cp258_1649', 'fbaA', 'rpsB', 'Cp258_1028', 'rplP', 'ccdA', 'corA', 'Cp258_0318', 'Cp258_0559', 'Cp258_1264', 'Cp258_1429', 'entD', 'galE', 'gyrB', 'oppCD1', 'purC', 'secF', 'sigH', 'yceG', 'Cp258_0847', 'Cp258_1314', 'Cp258_2095', 'metG', 'pimT', 'cobM', 'engA', 'ppgK', 'Cp258_1926', 'sufR', 'Cp258_1483', 'glgE', 'terC', 'dctA', 'tcsS3', 'AY851612', 'AY851612', 'AY851612']

def findSeq(genelist):
    alist = ''

    for gene in genelist:
        sterm = gene
        handle = Entrez.esearch(db="gene", retmode="xml", term=sterm)
        record = Entrez.read(handle)
        IDArray = record["IdList"]

        if len(IDArray) == 0:
            sterm = gene
            handle = Entrez.esearch(db="nuccore", term=sterm)
            record = Entrez.read(handle)
            IDArray = record["IdList"]
            toString = str(IDArray[0])
            gene_id = toString
            id = '>' + gene + '\n'
            alist += id
        else:
            toString = str(IDArray[0])
            gene_id = toString
            id = '>' + gene + '\n'
            alist += id

        try:
            handle = Entrez.efetch(db="nuccore", id=gene_id, rettype="gb", retmode="fasta")
            whole_sequence = SeqIO.read(handle, "genbank")
            seq = whole_sequence.seq
            alist += seq

        except:
            print('gene ' + gene + ' not found')


    with open('genelist.fasta','w') as o:
        o.write(str(alist))

findSeq(x)