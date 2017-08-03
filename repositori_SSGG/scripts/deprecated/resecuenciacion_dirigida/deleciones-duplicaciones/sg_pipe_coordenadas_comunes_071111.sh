#!/bin/bash

#El input del script son ficheros pileup. Es importante poner los ficheros pileup en un mismo directorio y dar la entrada como *pileup
echo $(date)
#b=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
b=(1)
mkdir tmp
cd tmp
for sample in $@
do
	for i in ${b[@]}
	do
		#Separa por cromosomas
		awk '{if ($1=="chr'$i'") print $0}' ../$sample  > $sample\_bases_$i
		#Saca posiciones en rango. 
		#-------------> Cambiar intervalos de captura
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/cardio/chr$i\_intervalos_sondas $sample\_bases_$i > $sample\_bases_N_rango_$i
		#Coge posiciones cubiertas y coverage correspondiente
		awk '{print $2"\t"$8}'  $sample\_bases_N_rango_$i > $sample\_$i
		#Saca posiciones que han de cubrirse a partir de los intervalos de captura. Convierte los intervalos a posiciones unicas en una sola columna. Trabaja por cromosoma
		#------------->  Cambiar intervalos de captura
		/share/apps/scripts/resecuenciacion_dirigida/deleciones-duplicaciones/sg_convertir_intervalo_N_columna_posiciones.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/cardio/chr$i\_intervalos_sondas | sort -u -n >  posiciones_$i
		
	done

done

rm *_bases_*

echo "#Chr Positions" $@ > o && sed 's/ /\t/gi' o > oo
 

arr=()
for x in ${b[@]}
do
	for  samples in $@
	do
		#El script sg_coordenadas_comunes_011111.pl coge como input un array de nombres y un array de posiciones de referencia. Cada muestra del array se compara con las posiciones de referencia. 
		#Meto en la variable arr el nombre de la muestra y el cromosoma. 
		arr=("${arr[@]}" "$samples"_"$x")
	done
	#Chequea si el fichero posiciones_$x esta vacio o no
	if [[ -s posiciones_$x ]] ; then
		#El script coge las posiciones de referencia que han de capturarse y un array o listado de ficheros de muestras que contienen dos columnas, posicion y coverage. Si la posicion de referencia y de la muestra coinciden, es decir, la base esta cubierta en la muestra, el script anyade' el valor del coverage, si la base no se cubre en la muestra el script anyade' un 0.
		/share/apps/scripts/resecuenciacion_dirigida/deleciones-duplicaciones/sg_coordenadas_comunes_011111.pl ${arr[@]} posiciones_$x
		
	fi
	#Una vez que termina el script vacio la variable para utilizarla para cada uno de los cromosomas
	arr=()
done

for fin in ${b[@]}
do
	#La salida del script sg_coordenadas_comunes_011111.pl genera ficheros con extension comp que en este paso se van a utilizar para generar un fichero por cromosoma en el que consten las posiciones de la referencia y el valor del coverage para cada una de las muestras
	paste posiciones_$fin *_$fin\comp > cromosoma_$fin
	
	if [[ -s cromosoma_$fin ]] ; then
		awk '{print "chr'$fin'\t"$0}' cromosoma_$fin > matriz_comparacion_entre_muestras_chr$fin.tab  
		/share/apps/scripts/resecuenciacion_dirigida/deleciones-duplicaciones/sg_intervalos_comunes_071111.pl matriz_comparacion_entre_muestras_chr$fin.tab covered > intervalos_$fin
	fi
done

cat oo intervalos_* > ../intervalos_todas_muestras


cd ..

cat tmp/matriz_comparacion_entre_muestras_chr*.tab > matriz_comparacion_entre_muestras


num_col=$(head -n 1 matriz_comparacion_entre_muestras | awk '{ print NF}')
#Normaliza cada columna de la matriz por la mediana
cp matriz_comparacion_entre_muestras moc
for ((a=3; a <= $num_col ; a++))
do
	valor="0"
	awk -v temp=$a '{print $temp}' matriz_comparacion_entre_muestras > col
	out=$(/share/apps/scripts/resecuenciacion_dirigida/deleciones-duplicaciones/sg_calcular_mediana_columna.pl col)
	awk -v temp=$a '{if ($temp==0) {print $valor} else {$temp=$temp/'$out'; print}}' moc > t
	mv t moc
done
mv moc  tmp/matriz_comparacion_entre_muestras_todos_cromosomas_normalizado_por_cols_tmp && sed 's/ /\t/gi' tmp/matriz_comparacion_entre_muestras_todos_cromosomas_normalizado_por_cols_tmp > tmp/matriz_comparacion_entre_muestras_todos_cromosomas_normalizado_por_cols
/share/apps/scripts/resecuenciacion_dirigida/deleciones-duplicaciones/sg_calcular_mediana_N_fila.pl tmp/matriz_comparacion_entre_muestras_todos_cromosomas_normalizado_por_cols > tmp/o
cat tmp/oo tmp/o > matriz_comparacion_entre_muestras_normalizado
rm -rf tmp col
echo $(date)
