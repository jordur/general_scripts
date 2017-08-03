#!/bin/bash
#El input de este script es la ultima columna del output de mirDeep2, columna que describe la poscion de la horquilla y su orientacion
#Ejemplo: chr15:83424731..83424839:+
#El input es un unico fichero generado a partir de la concatenacion de esta ultima columna para todas las muestras incluidas en el estudio.
#El objetivo del script es obtener un unico fichero que contenga todas las horquillas identificadas en todas las muestras para poder asi utilizar estas coordenadas para extraer la cuantificacion por miRNA nuevo utilizando BedTools 

chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

#Formateo el fichero de input para poder reutilizar el script sg_unificar_intervalos_miRNAs.pl. Separo los miRNAs tambien por cadena, positiva o negativa
sed 's/:/\t/g' $1 | sed 's/\.\./\t/g' | awk '{if ($4=="+") print $0}' > $1\_pos
sed 's/:/\t/g' $1 | sed 's/\.\./\t/g' | awk '{if ($4=="-") print $0}' > $1\_neg
for cromo in ${chrom[@]}
do
	awk '{if ($1=="chr'$cromo'") print $0}' $1\_pos | sort -k2n > $cromo\_pos
	awk '{if ($1=="chr'$cromo'") print $0}' $1\_neg | sort -k2n > $cromo\_neg
done

for nome in ${chrom[@]}
do
	sg_unificar_intervalors_miRNAs.pl $nome\_pos | awk '{print $0"\t+"}' > $nome\_pos_unif
	sg_unificar_intervalors_miRNAs.pl $nome\_neg | awk '{print $0"\t-"}' > $nome\_neg_unif		
	cat $nome\_pos_unif $nome\_neg_unif | awk '{print $1"\t"$2"\t"$3"\t"$1"_"$2"_"$3"_"$4"\t1\t"$4}'  > $nome\_OK_bed
done
cat *_OK_bed > $1\_final.bed
rm *_pos *_pos_unif *_OK_bed *_neg_unif *_neg

