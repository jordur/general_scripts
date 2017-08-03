#!/bin/bash

echo "INFO: Starting at $(date)"

preprocess ()
{
	prefix=$4
	CORRECTED=-1
	COUNT=0
	INPUT="$1"
	INPUT2=(${INPUT//,/ })
	for i in "${INPUT2[@]}"; do
		mkdir "primary_stats_$COUNT"

		PATTERN=${i##*/}
		
		echo "Input file: $PATTERN"
		
		if [ ! -z "$prefix" ]; then
			PATTERN=$prefix
		fi
		
		SCALAR=$(grep -c ">" $i)
		SCALAR2=$(grep  -v "#" $i | grep  -v ">" | wc -l )
		if [ "$SCALAR" != "$SCALAR2" ]; then
			sg_primary_quality.pl  $i 1 null null > primary_stats_$COUNT/ghost.txt
			mv primary_stats_$COUNT/ghost.txt primary_stats_$COUNT/$COUNT\_STEP1
			let CORRECTED=1
		else
			let CORRECTED=0
		fi

		if [ $SCALAR -ge 40000000 ]; then
			RAND_LINE=10000000
		else
			let RAND_LINE=$SCALAR/3
		fi

		RAND_LINE2=$[$RAND_LINE/$2]

		#cd primary_stats_$COUNT/

		if [ "$CORRECTED" -eq "0" ]; then
			FILES=$(split -l $RAND_LINE $i)
			templates=$(ls x*)
			mv x* primary_stats_$COUNT/
			cd primary_stats_$COUNT/

			for tpl in $templates
			do
				SCALAR_SUB=$(grep -c ">" $tpl)
				if [ "$SCALAR_SUB" -gt "500" ];##AQUI
				then
					if [ "$RAND_LINE2" -ge "$SCALAR_SUB" ];
					then
						let SCALAR_SUB=$RAND_LINE2
					fi
					sg_primary_quality.pl  $tpl 2 $RAND_LINE2 $SCALAR_SUB $3&
				fi
			done
			wait
		else
			cd primary_stats_$COUNT/

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
					sg_primary_quality.pl $tpl 2 $RAND_LINE2 $SCALAR_SUB $3 &
				fi
			done
			wait
		fi
		
		sg_primary_quality.pl  null 3 null
		sg_primary_quality_stats.R sequence_positions_final.stats qual_per_positions_final.stats  qual_mean_final.stats  indt_per_positions_final.stats indt_number_final.stats $PATTERN

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

# It will be checked for the proper number of arguments
if [ ${#} -lt $EXPECTED_ARGS ]; then
	echo 'Name: primary_Stats.sh'
	echo 'Description: Process and Visualization of Quality Stats of primary Data'
	echo 'Mandatory parameters:'
	echo '       $1 : One o more input *.qual. Separated by comma'
	echo '       $2 : Number of threads'
	echo '       $3 : Junior yes/no'
	echo '       $4 : Output prefix (optional, default is qual file name)'
	exit $E_BADARGS
fi

preprocess $1 $2 $3 $4