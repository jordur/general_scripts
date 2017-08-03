#/bin/bash

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

for j in ${i[@]}
do

	awk '{if ($1=="'$j'") print $0}' $1 | sort -u > xxx_$j
	sg_convertir_intervalo_N_columna_posiciones.pl xxx_$j | sort -u -k2n > xxxX_$j
	wc -l xxxX_$j >> o
done
awk '{print $1}' o > oo
sg_sumatorio_columna.pl oo
rm o oo xxx*
