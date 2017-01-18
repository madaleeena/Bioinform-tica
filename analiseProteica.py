# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 13:31:20 2016

@author: grupo5
"""

from Bio import SeqIO
from Bio import SwissProt 

record = SeqIO.read('sequenceGenbank.gb', 'genbank')
handle = open("uniprot-proteome.txt")
swiss_record = [record for record in SwissProt.parse(handle)]
entries=[]
accession = 'NC_002942.5'
fich_rev = open('AnaliseProteinas.txt', 'w')
fich_rev.write('Análise das Propriedades das Proteínas')
for i in range(len(record.features)):
    if record.features[i].type == "CDS":
        if ('gene' in record.features[i].qualifiers):
            geneName = record.features[i].qualifiers['gene'][0]
        else:
            geneName = '-'
        locustag = record.features[i].qualifiers['locus_tag'][0]
        start = str(record.features[i].location.start+998651) 
        stop = str(record.features[i].location.end+998651)
        proteinName = record.features[i].qualifiers['product'][0]
        if ('function' in record.features[i].qualifiers):
            function = record.features[i].qualifiers['function'][0]
        else:
            function = '-'
        geneID = record.features[i].qualifiers['db_xref'][0].replace('GeneID:','')
        for entry in swiss_record:
            for ref in entry.cross_references:
                if ref[0] == 'RefSeq':
                    if 'protein_id' in record.features[i].qualifiers.keys():
                        if ref[1] == record.features[i].qualifiers['protein_id'][0] and ref[1] not in entries:
                            nomegene = locustag+'('+geneName+')'
                            fich_rev.write(str("\n\nNome do gene: %s" %nomegene))
                            location = str(start) + '-' + str(stop)
                            fich_rev.write(str("\nZona do genoma: %s" %location))
                            fich_rev.write(str("\nStatus: %s" %entry.data_class))
                            uniprotID = entry.entry_name.replace('_LEGPH','')
                            fich_rev.write(str("\nUniProt ID: %s" %uniprotID))
                            fich_rev.write(str("\nNome da proteína: %s" %proteinName))
                            fich_rev.write(str("\nFunção molecular: %s" %function))
                            localizacao = None
                            for cr in entry.cross_references:
                                if cr[0] == "GO":
                                    (type, ids, cl, pis) = cr
                                    if type == "GO":
                                        cls = str(cl).split(":")
                                        if cls[0] == 'C':
                                            localizacao = cls[1]
                            if localizacao is None:
                                fich_rev.write(str("\nLocalização celular: -"))
                            else:
                                fich_rev.write(str("\nLocalização celular: %s" %localizacao))
fich_rev.close()
