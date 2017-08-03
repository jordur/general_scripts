#!/bin/bash

#---------------------------------
#sg_unificar_intervalos_captura.sh
#Sistemas Genomicos, 23 Julio 2011
#Sheila
#---------------------------------


i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

for j in ${i[@]}
do
	awk '{if ($1=="chr'$j'") print $0}' $1 > tmp_chr$j
	sort -k2n tmp_chr$j > tmp_chr$j\_ordenado
	sg_unificar_intervalos_de_captura_solapantes.pl tmp_chr$j\_ordenado > exoma50_chr$j	
done

cat exoma50_chr* > exoma50_intervalos_captura.bed
rm tmp_*
