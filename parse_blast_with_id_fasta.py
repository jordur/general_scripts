#! /usr/bin/python3

"""
Script per mapejar un identificador fasta en un resultat blast sepatrat per tabuladors
Dividim el resultat blast pels diferents camps i busquem el primer d'ells en l'identificador d'un arxiu fasta
En aquest cas el vaig usar per buscar les toxines de Giulia en uns resultats blast que vaig fer
"""

from Bio import SeqIO
from compiler.pyassem import DONE
toxines = "/home/bec2-jcalvete/Feina_Jordi/Libia/Dispholidus/boomslang_toxins_translated.fas.fasta"
blast = open("/home/bec2-jcalvete/Feina_Jordi/Libia/Dispholidus/all_toxins_Dispholidus.txt", "r")
file_out=open('/home/bec2-jcalvete/Feina_Jordi/Libia/Dispholidus/output.fasta','w')

for line in blast:
    name = line.split()[0].lower()
    desc1 = line.split()[1]
    desc2 = line.split()[2] 
#    print desc1
    for seq_record in SeqIO.parse( toxines, "fasta"):
        if seq_record.id in name:
           file_out.write(">%s\t%s\t%s\n%s\n" % (seq_record.id,desc1,desc2,seq_record.seq)) 
           #file_out.write(out)
           #print(repr(seq_record.seq))
           #print(len(seq_record))
print "Done"     
