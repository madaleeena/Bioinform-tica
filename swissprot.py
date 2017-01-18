# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 00:58:56 2017

@author: danie
"""

from Bio import SwissProt
from Bio import SeqIO
from Bio.Blast import NCBIWWW


def open_swissprot():
    handle = open("uniprot-proteome.txt")
    swiss_record = [record for record in SwissProt.parse(handle)]
    return swiss_record
    
    
#Procura as features das proteínas associadas aos nossos locus tag
def select_features():
    seq_record = SeqIO.read('seq_record.gb', 'genbank')
    f = seq_record.features
    features = list()
    for i in range(925,1155):
        for ft in f:
            if 'locus_tag' in ft.qualifiers.keys():
                if '0'*(4-int(len(str(i))))+str(i) in ft.qualifiers['locus_tag'][0]:
                    features.append(ft)
    return features
             
#Guarda o registo das diferentes entries obtidas da swissprot 
def select_swissprot_entries():
    features = select_features()
    entries = list()
    for ft in features:
        if ft.type == 'CDS':
            for entry in open_swissprot():
                for ref in entry.cross_references:
                    if ref[0] == 'RefSeq':
                        if 'protein_id' in ft.qualifiers.keys():
                            if ref[1] == ft.qualifiers['protein_id'][0] and ref[1] not in entries:
                                entries.append(entry)
    return entries
    
#Seleciona as "reviewed entries" dentro da lista de entries obtida na função anterior
def swissprot_reviewed():
    sse = select_swissprot_entries()
    r_entries = list()
    for entry in sse:
        if entry.data_class == 'Reviewed':
            r_entries.append(entry)
    return r_entries

"""    
Rev = swissprot_reviewed()
fich_rev = open("Prot_Reviewed.txt", "w")
for record in Rev:
    for ref in record.references:
        fich_rev.write(str("Nome:%s" %record.entry_name))
        fich_rev.write(str("\nAccession:%s" %record.accessions[0]))
        fich_rev.write(str("\nNome do gene:%s\n\n" %record.gene_name))
fich_rev.close()
"""

# Com a lista obtida, é possível obter uma lista dos ID's das proteínas que será útil para o BLAST
def interest_proteins():
    IDs = []
    Rev = swissprot_reviewed()
    for ss in Rev:
        for ref in ss.cross_references:
            if ref[0] == 'RefSeq':
                IDs.append(ref[1])
    return IDs
    
# print(interest_proteins())
    
trad = interest_proteins() 
ficheiro = open("my_blast.xml","w")

for prot in trad:
    handle = NCBIWWW.qblast("blastp","nr", prot)
    ficheiro.write(handle.read())
ficheiro.close()
handle.close()   
