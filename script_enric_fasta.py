#!/usr/bin/python
#De Sergio para Enric with love. 21/6/2017 --> Me debes una comida.
import sys
archivo_fasta = sys.argv[1]
archivo_nombres = sys.argv[2]

not_take = []
with open(archivo_nombres,"r") as file2:
	for lines in file2:
		not_take.append(lines.rstrip())

		
dic = {}
n = 0
m = 0
list_headers = []
with open(archivo_fasta,"r") as file1:
	for line in file1:
		line = line.rstrip()
		if line.startswith(">"):
			header = line
			if header in dic:
				print header
			dic[header] = ""
			list_headers.append(header)
		else:
			dic[header] += line
		n += 1

#print m,n
#print len(dic)

print len(not_take)

output = open("fasta_for_enric.fa","w")
for x,y in dic.items():
	if x[1:] not in not_take:
		output.write(x+"\n"+y+"\n")
	