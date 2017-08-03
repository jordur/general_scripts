#/bin/bash

############################################
#sg_unificar_intervalos_de_captura_desde_BED.sh v1
#Creado por Sheila. 20/07/12
############################################


EXPECTED_ARGS=1
E_BADARGS=666

if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'DESCRIPCION: Genera un fichero con los intervalos de captura no solapantes a partir de un fichero BED con el siguiente formato\n'
        echo ""
	echo browser position chr1:762097-762217
        echo track db="hg19" name="Human All Exon V4 plus UTRs" description="" visibility=2 color=0,128,0 useScore=1
        echo chr1    762097      762217   A_37_B00586175  1000    +
        echo chr1    762200      762320   A_37_B00705005  1000    +
	echo ""
        echo 'INPUT: Fichero BED con la ruta completa donde se encuentra'
        echo 'OUTPUT: Fichero tabulado con tres columnas Ej: chr1<tab>1200<tab>1320'
        exit $E_BADARGS
fi


i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
mkdir tmp-unificar-intervalos
cd tmp-unificar-intervalos
grep -v position $1 | grep -v track | awk '{print $1"\t"$2"\t"$3}' | sed 's/chr//' | sed 's/\t0*/\t/gi' > fichero-inicio-formateado
for j in ${i[@]}
do
	awk '{if ($1=="'$j'") print $0}' fichero-inicio-formateado > spl_tmp1_$j
	sg_convertir_intervalo_en_pileup.pl spl_tmp1_$j | sort -u -k2n  > spl_tmp2_$j
	awk '{print "chr"$1"\t"$2"\t"$8}' spl_tmp2_$j > spl_tmp3_$j
	if [ -s "spl_tmp3_$j" ]; 
	then
		sg_obtener_intervalos_de_pileup.pl spl_tmp3_$j "covered"   > ../nuevos-intervalos_$j
        fi
	 
done
cd ..
rm -rf tmp-unificar-intervalos

