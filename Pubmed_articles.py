# -*- coding: utf-8 -*-
"""
Created on Fri Dec 16 14:14:22 2016

@author: grupo5
"""

#Ref: Tutorial Biopython (9.14.1)
from Bio import Entrez
from Bio import Medline

Entrez.email = "pg32939@alunos.uminho.pt"   
handle = Entrez.egquery(term="legionella pneumophila subsp. pneumophila str. Philadelphia 1")
record = Entrez.read(handle)

#numero de artigos
for row in record["eGQueryResult"]:
    if row["DbName"]=="pubmed":
        nr = row["Count"]

#download dos identificadores do Pubmed dos artigos
handle = Entrez.esearch(db="pubmed", term="Legionella pneumophila subsp. pneumophila str. Philadelphia 1", retmax=nr)
record = Entrez.read(handle)
idlist = record["IdList"]


handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
records = list(Medline.parse(handle))

#criacao de fichero txt com os artigos(titulo,autor(es), source)
result = open("qqLegionellaP_literature.txt","w")
for record in records:
    title = ("Title: " + str(record.get("TI", "?")))
    authors = ("Authors: " + str(record.get("AU", "?")))
    source = ("Source: " + str(record.get("SO", "?")))
    space = ("\n")
    result.write(title)
    result.write(space)
    result.write(authors)
    result.write(space)
    result.write(source)
    result.write(space)
    result.write(space)
result.close()


