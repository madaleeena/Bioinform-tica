from Bio import SeqIO

 
#recolherfuncoes
#traduzir os genes na sequencia


record = SeqIO.read("sequenceGenbank.gb", "gb")
funcao=[]
for f in record.features:
    if f.type == "CDS":
        if 'function' in f.qualifiers:
            funcao.append((f.qualifiers["function"][0]))

#auxilio grupo 14

#print de INF importante

print("#######" + "  Numero de genes com funcao definida e desconhecida  " + "#######" + "\n")
print ("Numero de genes com funcao definida: " + str(len(func)))
print ("Numero de genes com funcao desconhecida: " + str(230-len(func))) #230 genes no total
