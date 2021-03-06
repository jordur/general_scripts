#!/bin/bash

echo "Inicio $(date)"

preprocess ()
{
	CORRECTED=-1
	COUNT=0
	INPUT="$1"
	INPUT2=(${INPUT//,/ })
	for i in "${INPUT2[@]}"
	do
		echo "Iniciando  $i"
		mkdir "primary_stats_$COUNT"

		SCALAR=$(grep -c ">" $i)
		SCALAR2=$(grep  -v "#" $i | grep  -v ">" | wc -l )
		if [ "$SCALAR" != "$SCALAR2" ]; then
			echo "Control del fichero $i No OK!, Secuencias de calidad en diferentes lineas"
			echo "corrigiendo ***"
			/data/results/solid0065/jmrosa/muestras_junior/calidad_datos_primarios.pl $i 1 null null > primary_stats_$COUNT/ghost.txt ##CAMBIAR ESTO CUANDO ESTE EN EL PATH EL SCRIPT
			mv primary_stats_$COUNT/ghost.txt primary_stats_$COUNT/$COUNT\_STEP1
			let CORRECTED=1
			echo "*** corregido"
		else
			echo "Control del fichero $i OK!"
			let CORRECTED=0
		fi

		if [ $SCALAR -ge 40000000 ]; then
			RAND_LINE=10000000
		else
			let RAND_LINE=$SCALAR/3
		fi
		echo "Calculando el n�mero de lecturas aleatorias para $SCALAR"
		echo "Calculado $RAND_LINE lecturas"

		RAND_LINE2=$[$RAND_LINE/$2]
		echo "Calculado un n�mero $RAND_LINE2 lecturas en $2 threads"

		cd primary_stats_$COUNT/

		if [ "$CORRECTED" -eq "0" ]; then
			FILES=$(split -l $RAND_LINE $i)
			templates=$(ls x*)
			for tpl in $templates
			do
				SCALAR_SUB=$(grep -c ">" $tpl)
				if [ "$SCALAR_SUB" -gt "500" ];##AQUI
				then
					if [ "$RAND_LINE2" -ge "$SCALAR_SUB" ];
					then
						let SCALAR_SUB=$RAND_LINE2
					fi
					/data/results/solid0065/jmrosa/muestras_junior/calidad_datos_primarios.pl $tpl 2 $RAND_LINE2 $SCALAR_SUB $3&
				fi
			done
			wait
		else
			FILES=$(split -l $RAND_LINE $COUNT\_STEP1)
			templates=$(ls x*)

			for tpl in $templates
			do
				SCALAR_SUB=$(grep -c ">" $tpl)
				if [ "$SCALAR_SUB" -gt "50" ];##AQUI
				then
					if [ "$RAND_LINE2" -ge "$SCALAR_SUB" ];
					then
						let SCALAR_SUB=$RAND_LINE2
					fi
					/data/results/solid0065/jmrosa/muestras_junior/calidad_datos_primarios.pl $tpl 2 $RAND_LINE2 $SCALAR_SUB $3 &
				fi
			done
			wait
		fi
		
		/data/results/solid0065/jmrosa/muestras_junior/calidad_datos_primarios.pl null 3 null

		Rscript /data/results/solid0065/jmrosa/muestras_junior/primary_stats.R sequence_positions_final.stats  qual_per_positions_final.stats  qual_mean_final.stats  indt_per_positions_final.stats indt_number_final.stats $COUNT

		mv *.pdf ../

		cd ..
		rm -rf primary_stats_$COUNT/
		let COUNT++
	done
}
# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

EXPECTED_ARGS=3 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


##PIPELINE CANCER $1 BAM OUTPUT FROM BIOSCOPE PAIRED END $2 BAM OUTPUT FROM BIOSCOPE FRAGMENT  $3 Output Bioscope Small indel (*gff) in Paired End $4 Output Bioscope Small Indels in Fragment (*.gff) $5 Type Analysis (Exoma38 Exoma50 Cardio Cancer)

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Name: primary_Stats.sh'
        echo 'Description: Process and Visualization of Quality Stats of primary Data'
        echo 'Mandatory parameters:'
        echo '       $1 : One o more input *.qual. Separated by comma. Complete Path'
        echo '       $2 : Number of threads'
	echo '       $3 : Junior yes/no'
        exit $E_BADARGS
fi

preprocess $1 $2 $3
