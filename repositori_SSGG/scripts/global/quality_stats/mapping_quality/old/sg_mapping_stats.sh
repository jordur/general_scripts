#!/bin/bash

main()
{
	# Definitions
	chrs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y )
	bam=`readlink -f $1`
	coverage_threshold=$2
	offset=$3
	slot_limit=$4
	output=`readlink -f $5`
	priority=$6
	target=`readlink -f $7`
	cwd=`pwd`

	echo "Inicio $(date)"
	RAND1=$(date +%N)
	
	# Generate temporal and logs folders
	mkdir "temp.$RAND1"
	if [ ! -d stats/mapping ]; then
		mkdir -p stats/mapping
	fi
	cd stats/mapping

	echo "cd $cwd/temp.$RAND1
bamtools split -in $bam -reference -stub split
for chr in ${chrs[@]}; do awk -v chr=\$chr '{if(\$1==\"chr\"chr) print \$0}' $target > chr\${chr}_target_intervals; done" > split_$RAND1.job
	job_1=`submit_job -p $6 split_$RAND1.job | awk '{print $3}'`

	#CREATE THE JOBS ARRAY FILE
	echo "cd $cwd/temp.$RAND1
echo \"Inicio $(date)\" 
chromosomes=( ${chrs[@]} );

SGE_TASK_ID=\$((SGE_TASK_ID-1))

FILE=chr\${chromosomes[\$SGE_TASK_ID]}_target_intervals

samtools view split.REF_chr\${chromosomes[\$SGE_TASK_ID]}.bam > split.REF_chr\${chromosomes[\$SGE_TASK_ID]}.sam
samtools view -bT /share/references/genomes/human/hg19/reference/chr\${chromosomes[\$SGE_TASK_ID]} split.REF_chr\${chromosomes[\$SGE_TASK_ID]}.sam > split.REF_chr\${chromosomes[\$SGE_TASK_ID]}_OK.bam
samtools index split.REF_chr\${chromosomes[\$SGE_TASK_ID]}_OK.bam

if [ -s \"\$FILE\" ]; then
	TEQC.R split.REF_chr\${chromosomes[\$SGE_TASK_ID]}_OK.bam \$FILE chr\${chromosomes[\$SGE_TASK_ID]} $coverage_threshold $offset
	echo \${chromosomes[\$SGE_TASK_ID]}  >> chromosome
	grep -v V1 coverage_threshold.txt > coverage_threshold_OK.txt

	if [ -s \"reads_target_offSet\" ]; then
		paste reads_target reads_target_offSet | grep -v x | awk '{print \$2\"\t\"\$3\"\t\"\$6}' > specificity.txt
	else
		 grep -v x reads_target | awk '{print \$2\"\t\"\$3}' > specificity.txt 
	fi
fi" > array_mapping_stats_$RAND1.job

	#LAUNCH THE JOBS ARRAY:
	submit_job -hold_jid $job_1 -t 1-24 -tc $slot_limit -p $priority -N sort_$RAND1 array_mapping_stats_$RAND1.job

    echo "cd $cwd/temp.$RAND1
	sed 's/\"//g' chromosome | sort -k1n  > chromosome_Sort.txt
	sed 's/chr//' coverage_targetCoverage.txt | sed 's/\"//g' | sort -k2n -k3n -k4n > coverage_targetCoverage_Sort.txt
	sed 's/chr//' coverage_threshold_OK.txt | sed 's/\"//g' | sort -k2n -k3n > coverage_threshold_OK_Sort.txt
	sed 's/chr//' specificity.txt | sed 's/\"//g' | sort -k1n > specificity_Sort.txt

	TEQC_Graphics.R coverage_targetCoverage_Sort.txt coverage_threshold_OK_Sort.txt chromosome_Sort.txt specificity_Sort.txt $offset $coverage_threshold
	if [ ! -d $output ]; then
		mkdir $output
	fi
	mv *.pdf $output
	mv coverage_targetCoverage_Sort.txt $output
	mv coverage_threshold_OK_Sort.txt $output
	mv specificity_Sort.txt $output
	cd $cwd
	#rm -rf temp.$RAND1" > sort_matrix_$RAND1.job

	submit_job -hold_jid sort_$RAND1 -p $6  sort_matrix_$RAND1.job 
	echo " Final $(date)"
}

# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

EXPECTED_ARGS=7 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments

##COVERAGE STATS PIPELINE
# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Name: TEQC.sh'
        echo 'Description: Coverage Stats Pipeline'
	echo 'Mandatory parameters:'
	echo '       $1 : Bam file:'
	echo '       $2 : Coverage Thresholds: (Ejm 1,10,20), please no more than 10 thresholds'
	echo '       $3 : Offset: The pipeline compute the target reads and Offset taregts reads'
	echo '       $4 : Slot Limit:'
	echo '       $5 : Name output directory'
	echo '       $6 : Priority submit job'
	echo '       $7 : Path to target intervals files'
        exit $E_BADARGS
fi

main $1 $2 $3 $4 $5 $6 $7
