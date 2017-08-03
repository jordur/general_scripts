#!/bin/bash

# -------------------
# - GATK_indels_SNPs -
# -------------------

GATK_indels_SNPs ()
# Description: sequence of actions to carry out on the files for obtaining the mapping metrics
# Parameters:
#       $1 : root directory (ending with /)
#       $2 : sample name
{
        # ---------------
        # - Definitions -
        # ---------------


        # ------------------------
        # - body of the pipeline -
        # ------------------------

	# -N SAET_AZ1_F3_18_11_10
	#PBS -o SAET_AZ1_F3_18_11_10.out
	#PBS -e SAET_AZ1_F3_18_11_10.err
	#PBS -q pipeline
	
	# The read group has to be added to each single read, since it's not present in the original .sam.
	# Please be aware that the RG:Z: value added to the reads (20110221001232913) has to be the same of the value in the header (cab)
	awk '{print $0,"\tRG:Z:20110221001232913"}' $1$2/454Reads.$2.sam > $1$2/454Reads.$2_RG.sam
	
	#The header will be added to the .sam
	cat $1cab $1$2/454Reads.$2_RG.sam > $1$2/454Reads.$2_headed.sam
	
	#Eliminate unmaped reads and convert to bam, maintaining the header
	samtools view -F 4 -h -S $1$2/454Reads.$2_headed.sam -b > $1$2/454Reads.$2_headed.bam
	
	#Es necesario ordenar el bam antes de indexarlo
	samtools sort $1$2/454Reads.$2_headed.bam $1$2/454Reads.$2_headed_sorted
	#Este comando sirve para indexar el bam
	samtools index $1$2/454Reads.$2_headed_sorted.bam
	
	#Calibrar GATK
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I $1$2/454Reads.$2_headed_sorted.bam -o $1$2/intervals.intervals -B:dbsnp,dbsnp $1datos_partida/selected131_sort.rod
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T IndelRealigner -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I $1$2/454Reads.$2_headed_sorted.bam -targetIntervals $1$2/intervals.intervals -o $1$2/realigner.bam -B:dbsnp,dbsnp $1datos_partida/selected131_sort.rod
	#Es necesario ordenar el nuevo bam 
	samtools sort $1$2/realigner.bam $1$2/realigner_sorted
	#Es necesario indexar el nuevo bam
	samtools index $1$2/realigner_sorted.bam
	
	#Calibrar GATK. Se definen los parámetros por defecto --default_read_group y --default_platform, por no venir incluidos en el fichero
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T CountCovariates -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa --DBSNP $1datos_partida/selected131_sort.rod -I $1$2/realigner_sorted.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $1$2/recal.csv --default_read_group $2 --default_platform 454
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T TableRecalibration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I $1$2/realigner_sorted.bam -o $1$2/realigner_Recalibrate.bam -recalFile $1$2/recal.csv --default_read_group $2 --default_platform 454
	
	#Eliminar lecturas con un valor de mapeo menor de 10 despues de recalibrar. Se quita el parámetro -S, pues el input es .bam. Es muy importante aquí el -h para que siga incluyéndose la cabecera!!!
	samtools view -b -h -q 10 $1$2/realigner_Recalibrate.bam > $1$2/realigner_Recalibrate_Q10.bam
	samtools index $1$2/realigner_Recalibrate_Q10.bam
	
	#Eliminar las lecturas duplicadas
	MarkDuplicates.jar  I=$1$2/realigner_Recalibrate_Q10.bam  O=$1$2/realigner_Recalibrate_Q10_nodup.bam M=stat REMOVE_DUPLICATES=TRUE ASSUME_SORTED=FALSE
	#Es necesario ordenar el nuevo bam
	samtools sort $1$2/realigner_Recalibrate_Q10_nodup.bam $1$2/realigner_Recalibrate_Q10_nodup_sorted
	#Es necesario indexar el nuevo bam
	samtools index $1$2/realigner_Recalibrate_Q10_nodup_sorted.bam
	
	# los siguientes comentarios pertenecen a una prueba anterior
	
	# Se detecta un problema en el paso de identificación de indels, donde las lecturas no aparecen con la identificación de grupo de lectura. Se añade manualmente:
	#samtools view MID1/realigner_Recalibrate_Q10_nodup_sorted.bam > MID1/realigner_Recalibrate_Q10_nodup_sorted.sam
	#awk '{print $0,"\tRG:Z:20110221001232913"}' MID1/realigner_Recalibrate_Q10_nodup_sorted.sam > MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded.sam
	# Obtenemos la cabecera y la añadimos
	#samtools view -H MID1/realigner_Recalibrate_Q10_nodup_sorted.bam > MID1/header
	#sed 's/MDI/Sample/' MID1/header > MID1/header_std
	#cat MID1/header_std MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded.sam > MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded_headed.sam
	# Pasamos a .bam
	#samtools view -b -S -h MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded_headed.sam > MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded_headed.bam
	# se indexa nuevamente
	#samtools index MID1/realigner_Recalibrate_Q10_nodup_sorted_SMadded_headed.bam
	
	# Creación de carpeta para Indels
	if [ ! -d "$1$2/IndelsSam" ]; then
                mkdir $1$2/IndelsSam
        fi	
	
	#Identificar indels ******************este es el paso que da/daba el problema!!!!!!!!!!!!!!!!!!!!!!
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I $1$2/realigner_Recalibrate_Q10_nodup_sorted.bam -glm DINDEL -o $1$2/IndelsSam/indels.vcf
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -B:variant,VCF $1$2/IndelsSam/indels.vcf -o $1$2/IndelsSam/first_filtration_indels.vcf --filterExpression "MQ0 >=4 && ((MQ0 / (1.0 * DP ))> 0.1)" --filterName "FIRST_FILTER" --filterExpression "SB >= -1.0" --filterName "STRAND_BIAS_FILTER" --filterExpression "QUAL < 10" --filterName "QUALFILTER"
	
	# Creación de carpeta para SNPs
	if [ ! -d "$1$2/SNPsSam" ]; then
                mkdir $1$2/SNPsSam
        fi

	#Identificar SNVs
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -I $1$2/realigner_Recalibrate_Q10_nodup_sorted.bam -o $1$2/SNPsSam/snp.vcf
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar  -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -B:variant,VCF $1$2/SNPsSam/snp.vcf -o $1$2/SNPsSam/first_filtration_snps.vcf --clusterWindowSize 10 --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 *DP)) > 0.1)" --filterName "FIRST_FILTER"
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar  -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -B:variant,VCF $1$2/SNPsSam/first_filtration_snps.vcf  --filterExpression "QUAL <30.0 || QD < 5.0 || HRun > 5 " --filterName "HARD_FILTER" -o $1$2/SNPsSam/hard_filtering_snps.vcf

}

	
# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

# it will be checked if the input path finishes with "/", and the recurse function will be called in consequence.
root_dir="/data/results/Solid0065/BF18_GSJunior/"
samples=("MID1" "MID2" "MID3" "MID4" "MID5" "MID6" "MID7" "MID8")
# the input_dir is based on the input_dir, this means: root_dir/samples/input_dir

for ((k=0;k<${#samples[@]};k+=1)); do
        GATK_indels_SNPs $root_dir ${samples[k]}
done

