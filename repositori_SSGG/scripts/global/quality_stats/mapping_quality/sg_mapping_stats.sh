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
	queue=$6
	target=`readlink -f $7`
	priority=$8
	job_to_wait=$9
	cwd=`pwd`
	RAND1=$(date +%N)

	# Generate temporal and logs folders
	mkdir "temp.$RAND1"
	
	echo "cd $cwd/temp.$RAND1
bamtools split -in $bam -reference -stub split.$RAND1
for i in ${chrs[@]}; do
	chr=\"chr\"
	end=\"\W\"
	chr_complet=\$chr\$i\$end
	grep \"\$chr_complet\"  $target > chr\$i\_target_intervals; 
done" > split_$RAND1.job
	if [ "$job_to_wait" == "" ]; then
		job_1=`submit_job -q $queue -p $priority split_$RAND1.job | awk '{print $3}'`
	else 
		job_1=`submit_job -q $queue -p $priority -hold_jid $job_to_wait split_$RAND1.job | awk '{print $3}'`
	fi

	#CREATE THE JOBS ARRAY FILE
	echo "cd $cwd/temp.$RAND1
echo \"Inicio $(date)\"
chromosomes=( ${chrs[@]} )

SGE_TASK=\$((\$SGE_TASK_ID - 1))
FILE=chr\${chromosomes[\$SGE_TASK]}_target_intervals

samtools view split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}.bam > split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}.sam
samtools view -bT /share/references/genomes/human/hg19/reference/chr\${chromosomes[\$SGE_TASK]} split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}.sam > split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}_OK.bam
samtools index split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}_OK.bam

if [ -s "\$FILE" ]; then
	echo \"Starting TEQC for split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}_OK.bam\"
	TEQC.R split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}_OK.bam \$FILE chr\${chromosomes[\$SGE_TASK]} $coverage_threshold $offset
	echo \"Finishing TEQC for split.$RAND1.REF_chr\${chromosomes[\$SGE_TASK]}_OK.bam\"
	echo \${chromosomes[\$SGE_TASK]}  >> chromosome
	grep -v V1 chr\${chromosomes[\$SGE_TASK]}_coverage_threshold.txt > chr\${chromosomes[\$SGE_TASK]}_coverage_threshold_OK.txt

	if [ -s \"chr\${chromosomes[\$SGE_TASK]}_reads_target_offSet\" ]; then
		paste chr\${chromosomes[\$SGE_TASK]}_reads_target chr\${chromosomes[\$SGE_TASK]}_reads_target_offSet | grep -v x | awk '{print \$2\"\t\"\$3\"\t\"\$6}' > chr\${chromosomes[\$SGE_TASK]}_specificity.txt
	fi
else
	grep -v x chr\${chromosomes[\$SGE_TASK]}_reads_target | awk '{print \$2\"\t\"\$3}' > chr\${chromosomes[\$SGE_TASK]}_specificity.txt
fi" > array_mapping_stats_$RAND1.job

	job_2=`submit_job -q $queue -p $priority -hold_jid split_$RAND1.job -t 1-24 -tc $slot_limit -N sort_$RAND1 array_mapping_stats_$RAND1.job | awk '{print $3}' | sed 's/\./ /' | awk '{print $1}'`

	echo "cd $cwd/temp.$RAND1
echo \"Inicio $(date)\"
chrs=( ${chrs[@]} );

	cat chr*_coverage_targetCoverage.txt >  global_targetCoverage.txt
	sed 's/chrX/chr23/' global_targetCoverage.txt | sed 's/\\\"//g' > global_targetCoverage_ghost.txt
	mv global_targetCoverage_ghost.txt global_targetCoverage.txt
	sed 's/chrY/chr24/' global_targetCoverage.txt > global_targetCoverage_ghost.txt
	mv global_targetCoverage_ghost.txt global_targetCoverage.txt
	sed 's/chr//g' global_targetCoverage.txt > global_targetCoverage_ghost.txt
	mv global_targetCoverage_ghost.txt global_targetCoverage.txt
	sort -k2n -k3n -k4n global_targetCoverage.txt  > global_targetCoverage_ghost.txt 
	mv global_targetCoverage_ghost.txt global_targetCoverage.txt

	cat chr*_coverage_threshold_OK.txt > global_coverage_threshold_OK.txt
	sed 's/chrX/chr23/'  global_coverage_threshold_OK.txt | sed 's/\\\"//g' > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	sed 's/chrY/chr24/' global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt 
	sed 's/chr//g' global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	sort -k2n -k3n -k4n global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	
	awk '{print \"chr\"\$0}' global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	sed 's/chr23/chrX/' global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	sed 's/chr24/chrY/' global_coverage_threshold_OK.txt > global_coverage_threshold_OK_ghost.txt
	mv global_coverage_threshold_OK_ghost.txt global_coverage_threshold_OK.txt
	

	for i in ${chrs[@]}; do
		cat chr\$i\_*_coverage_bases.txt | sed 's/\"//g' |  sort -k3n  > chr\$i.pileup
	done

	cat chr*.pileup > All_chr.pileup
	sed 's/chr//' All_chr.pileup > All_chr_ghost.pileup
	sort -k2n -k3n All_chr_ghost.pileup > All_chr.pileup
	rm All_chr_ghost.pileup
	

	cat chr*_specificity.txt > specificity.txt
	sed 's/chrX/chr23/' specificity.txt | sed 's/\"//g' > specificity_ghost.txt
	mv specificity_ghost.txt specificity.txt
	sed 's/chrY/chr24/' specificity.txt > specificity_ghost.txt 
 	mv specificity_ghost.txt specificity.txt
	sed 's/chr//g'  specificity.txt > specificity_ghost.txt
	mv specificity_ghost.txt specificity.txt
	sort -k1n specificity.txt > specificity_ghost.txt
	mv specificity_ghost.txt specificity.txt
	
	sed 's/X/23/' chromosome | sed 's/\"//g' > chromosome_ghost
	mv chromosome_ghost chromosome

	sed 's/Y/24/' chromosome > chromosome_ghost
	mv chromosome_ghost chromosome
	sort -k1n chromosome > chromosome_ghost
	mv chromosome_ghost chromosome 
	sed 's/23/X/' chromosome > chromosome_ghost
	mv chromosome_ghost chromosome
	sed 's/24/Y/' chromosome > chromosome_ghost
	mv chromosome_ghost chromosome" > global_file_parsing_$RAND1.job
	 
	job_3=`submit_job -q $queue -p $priority -hold_jid sort_$RAND1 -N global_$RAND1 global_file_parsing_$RAND1.job | awk '{print $3}'`

	echo "cd $cwd/temp.$RAND1
echo \"Inicio $(date)\"

	TEQC_Graphics.R global_targetCoverage.txt global_coverage_threshold_OK.txt chromosome specificity.txt $offset $coverage_threshold All_chr.pileup

	mkdir $output
	mkdir $output/images/
	mv *.pdf $output/images/
	mkdir $output/pileup/
	mv *.pileup $output/pileup/
	mkdir $output/tables/
	mv global_targetCoverage.txt $output/tables/
	mv global_coverage_threshold_OK.txt $output/tables/
	mv specificity.txt $output/tables/
	mkdir $output/chr/
	mv chromosome $output/chr/
	cd $cwd
	rm -rf $output/pileup/
	rm -rf $cwd/temp.$RAND1" > graphics_final_$RAND1.job

	submit_job -q $queue -p $priority -hold_jid global_$RAND1  graphics_final_$RAND1.job
}

# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

EXPECTED_ARGS=7 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments

##COVERAGE STATS PIPELINE
# It will be checked for the proper number of arguments
if [ ${#} -lt $EXPECTED_ARGS ]; then
	echo 'Name: TEQC.sh'
	echo 'Description: Coverage Stats Pipeline'
	echo 'Mandatory parameters:'
	echo '       $1 : Bam file'
	echo '       $2 : Coverage Thresholds: (Ejm 1,10,20), please no more than 10 thresholds'
	echo '       $3 : Offset: The pipeline computes the target reads and offset target reads'
	echo '       $4 : Slot Limit'
	echo '       $5 : Output folder name'
	echo '       $6 : Submit queue (low,med,high,mem)'
	echo '       $7 : Path to target intervals file'
	echo '       $8 : Job priority'
	echo '       $9 : Job ID to wait before start'
	exit $E_BADARGS
fi

main $1 $2 $3 $4 $5 $6 $7 $8 $9
