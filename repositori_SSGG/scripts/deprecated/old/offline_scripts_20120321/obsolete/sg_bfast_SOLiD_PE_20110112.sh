#!/bin/bash

EXPECTED_ARGS=4 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo 'Name: sg_bfast_SOLiD_PE.sh'
	echo 'Description: Pipeline for carrying out the alignment on PE reads from SOLiD by means of bfast.'
	echo '	This scripts splits the files in different files so that paralellize the process, and therefore make it faster.'
	echo '	The reference parameter is hardcoded (that is, it is included in the script code, in environment parameters section.'
	echo 'Mandatory parameters:'
	echo '	$1 : fasta_file_1'
	echo '	$2 : qual_file_1'
	echo '	$3 : fasta_file_2'
	echo '	$4 : qual_file_2'
	exit $E_BADARGS
fi


# ---------------------
# ---- definitions ----
# ---------------------

fasta_file_1=$1
fasta_file_2=$3
qual_file_1=$2
qual_file_2=$4

# Environment parameters:
reference="/data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa"


# ----------------------------------------
# ------- Main body of the script --------
# ----------------------------------------

export exec_time=`date`
echo Starting on $exec_time

# Check the amount of reads in the primary files:
echo Proccessing file $fasta_file_1
export number_reads_1=`grep -c ">" $fasta_file_1` #104132987
echo Number of reads: $number_reads_1
echo Proccessing file $fasta_file_2
export number_reads_2=`grep -c ">" $fasta_file_2` #104132987
echo Number of reads: $number_reads_2

# If the primary files have more than 10M reads, they get splitted:
export exec_time=`date`
echo "Splitting files in 10M reads files, if needed (for acceleranting computing process). Started on $exec_time"

if [ $number_reads_1 -gt 10000000 ]; then
	let "number_files_1=$number_reads_1 / 10000000"
	echo "Number of files for fasta1: $number_files_1"
	sg_split_fasta_file_en_varios_ficheros.pl $fasta_file_1 10000000
	for ((k=1;k<${number_files_1}+2;k+=1)); do
		mv resultado$k F3_$k.csfasta
	done
	sg_split_fasta_file_en_varios_ficheros.pl $qual_file_1 10000000
	for ((k=1;k<${number_files_1}+2;k+=1)); do
		mv resultado$k F3_$k.qual
	done
else
	cp $fasta_file_1 F3_1.csfasta
	cp $qual_file_1 F3_1.qual
fi

if [ $number_reads_2 -gt 10000000 ]; then
        let "number_files_2=$number_reads_2 / 10000000"
	echo "Number of files for fasta2: $number_files_2"
        sg_split_fasta_file_en_varios_ficheros.pl $fasta_file_2 10000000
        for ((k=1;k<${number_files_2}+2;k+=1)); do
                mv resultado$k F5_$k.csfasta
        done
        sg_split_fasta_file_en_varios_ficheros.pl $qual_file_2 10000000
        for ((k=1;k<${number_files_2}+2;k+=1)); do
                mv resultado$k F5_$k.qual
        done
else
	cp $fasta_file_2 F5_1.csfasta
	cp $qual_file_2 F5_1.qual
fi


# The read files get transformed into fastq files:
export exec_time=`date`
echo Transformation of read files into fastq files started on $exec_time

for ((k=1;k<${number_files_1}+2;k+=1)); do
	solid2fastq F3_$k.csfasta F3_$k.qual -o F3_${k}_bfast
done

for ((k=1;k<${number_files_2}+2;k+=1)); do
        solid2fastq F5_$k.csfasta F5_$k.qual -o F5_${k}_bfast
done


# Mapping and search for indexes in read files
if [ -e mapping.txt ]; then
	rm mapping.txt
fi

export exec_time=`date`
echo Starting mapping of F3_bfast.fastq with bfast match on $exec_time
for ((k=1;k<${number_files_1}+2;k+=1)); do
	echo "/share/apps/bfast+bwa-0.6.5a/bin/bfast match -n 1 -t -f $reference -A 1 -r F3_${k}_bfast.fastq > F3_${k}_matches_bfast.bmf &" >> mapping.txt
done

export exec_time=`date`
echo Starting mapping of F5_bfast.fastq with bfast bwaaln on $exec_time
for ((k=1;k<${number_files_2}+2;k+=1)); do
        echo "/share/apps/bfast+bwa-0.6.5a/bin/bfast bwaaln -c $reference F5_${k}_bfast.fastq -f F5_${k}_matches_bfastbwaaln.bmf &" >> mapping.txt
done

sg_commands_nodes_sequencer.sh 1 mapping.txt

# Matches files get concatenated
if [ -e F3_matches_bfast.bmf ]; then
	rm F3_matches_bfast.bmf
fi
for ((k=1;k<${number_files_1}+2;k+=1)); do
	cat F3_${k}_matches_bfast.bmf >> F3_matches_bfast.bmf
done

if [ -e F5_matches_bfastbwaaln.bmf ]; then
        rm F5_matches_bfastbwaaln.bmf
fi
for ((k=1;k<${number_files_2}+2;k+=1)); do
        cat F5_${k}_matches_bfastbwaaln.bmf >> F5_matches_bfastbwaaln.bmf
done

# Alineamientos locales y juntar las alineaciones F3 y F5
# bfast localalign se inicia a las 14:40h. Finaliza a las 22:30h. Duración: 8 horas
export exec_time=`date`
echo Starting local alignment of data with bfast localalign on $exec_time
echo "/share/apps/bfast+bwa-0.6.5a/bin/bfast localalign -f $reference -1 F3_matches_bfast.bmf -2 F5_matches_bfastbwaaln.bmf -A 1 -t -U -n 4 >PE_F3-bfast_F5-bfastbwaaln.baf 2>>nohup.out" > localalign.txt
sg_commands_nodes_sequencer.sh 1 localalign.txt

# Filtrado de alineaciones y creación de .sam
# bfast postprocess se inicia a las 9:43h del 13/07/2011. Finaliza a las 9:50h. Duración: 10 minutos
export exec_time=`date`
echo Starting alignment filtering and sam creation with bfast postprocess on $exec_time
echo "/share/apps/bfast+bwa-0.6.5a/bin/bfast postprocess -A 1 -R -n 1 -f /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -i PE_F3-bfast_F5-bfastbwaaln.baf >PE_F3-bfast_F5-bfastbwaaln.sam 2>>nohup.out" > postprocess.txt
sg_commands_nodes_sequencer.sh 1 postprocess.txt

# Creación de .bam
# samtools view se inicia a las 9:51h. Finaliza a las 10:05h. Duración: 15 minutos
export exec_time=`date`
echo Starting .bam creation with samtools view on $exec_time
/share/apps/samtools-0.1.16/samtools view -b -S PE_F3-bfast_F5-bfastbwaaln.sam > PE_F3-bfast_F5-bfastbwaaln.bam

export exec_time=`date`
echo Finished on $exec_time

