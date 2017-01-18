# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:31:20 2016

@author: grupo5
"""

from Bio import SeqIO

record = SeqIO.read('sequenceGenbank.gb', 'genbank')
accession = 'NC_002942.5'
f = open('tabelaNCBI.csv', 'w')
f.write('Accession;Start;Stop;Strand;GeneID;Locus;Locus tag;Protein Id;Length;Protein Name;Function\n')
for i in range(len(record.features)):
    if record.features[i].type == "CDS":
        geneID = record.features[i].qualifiers['db_xref'][0].replace('GeneID:','')
        if ('gene' in record.features[i].qualifiers):
            geneName = record.features[i].qualifiers['gene'][0]
        else:
            geneName = '-'
        locustag = record.features[i].qualifiers['locus_tag'][0]
        start = str(record.features[i].location.start+998651) 
        stop = str(record.features[i].location.end+998651)
        strand = record.features[i].location.strand
        proteinID = record.features[i].qualifiers['protein_id'][0]
        proteinName = record.features[i].qualifiers['product'][0]
        aminoNumber = len(record.features[i].qualifiers['translation'][0])
        if ('function' in record.features[i].qualifiers):
            function = record.features[i].qualifiers['function'][0]
        else:
            function = '-'
        geneID = record.features[i].qualifiers['db_xref'][0].replace('GeneID:','')
        f.write(accession+';'+str(start)+';'+str(stop)+';'+str(strand)+';'+geneID+';'+geneName+';'+locustag+';'+proteinID+';'+str(aminoNumber)+';'+proteinName+';'+function+'\n')
f.close()

