# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 22:02:25 2017

@author: danie
"""

from Bio.Blast import NCBIXML

def ler_blast():
    result = open('my_blast.xml') 
    blast_record = NCBIXML.parse(result)
    for record in blast_record: 
        for alignment in record.alignments:
            print('\n\nAcession:',alignment.accession)
            print('ID:',alignment.hit_id)
            print('Description:',alignment.hit_def)
            print('Length:',alignment.length)
            for hsp in alignment.hsps:
                print('\ne-value:',hsp.expect)
                print('Score:',hsp.score)
                print('Match:',hsp.match)
            
ler_blast()
