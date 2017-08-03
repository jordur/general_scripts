#!/bin/bash
#Con este script se pretende identificar las bases que nunca se capturan con el kit de enriquecimiento utilizado de manera que podamos utilizarlas con dos fines: 1) ver que mutaciones se localizan en estas bases, 2) performance del sistema de captura
#Input: ficheros pileup y los intervalos de sondas /data/results/Solid0065/info_paneles/cardio.... Para que el programa coja los invervalos de las sondas se ha de generar un enlace simbolico al directorio donde estemos trabajando utilizando el comando ln -s /data/results/Solid0065/info_paneles/cardio/* . (cambiamos cardio por el panel de genes que vayamos a usar).

b=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
mkdir tmp
cd tmp

#Separa cada uno de los ficheros pileup por cromosoma
for sample in $@
do
	for i in ${b[@]}
	do
		awk '{if ($1=="chr'$i'") print $2}' ../$sample  | sort -n > $sample\_$i
	done
done
#echo "Positions" $@ > oo 

#Calcula las posiciones unicas que se cubren en todas las muestras
for j in  ${b[@]}
do
#	cat *_$j | sort -u -n > posiciones_$j
	sg_bed2pileup_modificado.pl ../chr$j\_intervalos_sondas > posiciones_$j
done

#Ejecuta el programa principal para cada muestra y por cromosoma
arr=()
for x in ${b[@]}
do
	for  samples in $@
	do
		arr=("${arr[@]}" "$samples"_"$x")
	done
	#Chequea si posiciones_$x esta vacio o no
	if [[ -s posiciones_$x ]] ; then
		#echo "../sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango.pl ${arr[@]} posiciones_$x"
		sg_calcular_bases_q_no_se_cubren_desde_pileup_en_rango.pl ${arr[@]} posiciones_$x
	fi
	arr=()
done

for ok in ${b[@]}
do
	#Como se pone
	#if [-f *_$ok."comp"] ; then
		#echo "paste posiciones_$ok *_$ok\comp > cromosoma_$ok"
		paste posiciones_$ok *_$ok\comp > cromosoma_$ok
	#fi
done 

#Anade el nombre del cromosoma al que pertenece cada fichero con el objetivo de juntar todos los resultados en un unico fichero final
for fin in ${b[@]}
do
	if [[ -s cromosoma_$fin ]] ; then
		awk '{print "chr'$fin'\t"$0}' cromosoma_$fin > ../matriz_comparacion_entre_muestras_chr$fin.tab   
	fi
done

cd ..
sg_calcular_rangos_d_posiciones_no_cubiertas_a_partir_d_matriz.pl matriz_comparacion_entre_muestras_chr*.tab 

#cd ..

#Concatena todos los resultados obtenidos que se encuentran separados por cromosoma
#cat matriz_comparacion_entre_muestras_chr*.tab > matriz_comparacion_entre_muestras_todos_cromosomas

#Limpia todo el directorio dejando solo el resultado final matriz_comparacion_entre_muestras_todos_cromosomas
#rm -rf tmp  matriz_comparacion_entre_muestras_chr*.tab 
#rm -rf tmp
for crom in ${b[@]}
do
	sg_miRNA_extract_seq.pl matriz_comparacion_entre_muestras_chr$crom.tabrangos /data/results/Solid0065/referencia_genoma_humano/GRch37.58/maskedVersion/Homo_sapiens.GRCh37.62.dna_rm.chromosome.$crom.fa > tmp/seq_interval_noCubiertos_chr$crom.fa
	/usr/local/bin/scripts/resecuenciacion_dirigida/pipelineSNV_Biobase/sg_noRepeat_intervals.pl   tmp/seq_interval_noCubiertos_chr$crom.fa > tmp/seq_noCubiertos_Nopeat_NolowComplexy_chr$crom.fa
	grep ">" tmp/seq_noCubiertos_Nopeat_NolowComplexy_chr$crom.fa | sed 's/>/chr'$crom'\t/' > tmp/interval_noCubiertos_Nopeat_NolowComplexy.$crom.bed
        /usr/local/bin/scripts/resecuenciacion_dirigida/pipelineSNV_Biobase/sg_extraer_SNPs_fichero_Biobase.pl tmp/interval_noCubiertos_Nopeat_NolowComplexy.$crom.bed /data/results/Solid0065/info_paneles/HGMD-Biobase/bioBaseSortUniques_chr$crom.txt > ALL_SNV_Biobase_noCubiertos_$crom

done

cat tmp/seq_noCubiertos_Nopeat_NolowComplexy_chr*.fa > zonas_genoma_no_cubiertas_no_repetitivas.fasta

rm -rf tmp matriz_comparacion_entre_muestras_chr*

