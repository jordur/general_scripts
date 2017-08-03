#!/bin/bash

# -------------------
# - mapping_metrics -
# -------------------

mapping_metrics ()
# Description: sequence of actions to carry out on the files for obtaining the mapping metrics
# Parameters:
# 	$1 : root directory (ending with /)
#	$2 : sample name
#	$3 : input directory (ending with /)
#	$4 : output directory sufix (ending with /)
{
	# ---------------
	# - Definitions -
	# ---------------

	# Subindexes for the different chromosomes
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y)

        # ------------------------
        # - body of the pipeline -
        # ------------------------
	
	# Creates the output directory
	if [ ! -d "$1$2/$3$4" ]; then
                mkdir $1$2/$3$4
        fi
	
	# Remove duplicates
	MarkDuplicates.jar I=$1$2/$3solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Cardio_F3_$2.csfasta.ma.bam O=$1$2/$3F3_nodup.bam M=$1$2/$3$4F3_nodup.log REMOVE_DUPLICATES=true TMP_DIR=$1$2/$3$4scratch

	# Filtering of the reads by quality values below 10
	#samtools view -b -q 10 $1$2/$3F3_nodup.bam > $1$2/$3F3_Q10_nodup.bam
	
	# Sequences (with duplicates) get piled up against the reference
	#/share/apps/samtools-0.1.9/samtools pileup -cf /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa $1$2/$3solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Cardio_F3_$2.csfasta.ma.bam > $1$2/$3F3.pileup
	
	# Sequences (without duplicates) get piled up against the reference
	/share/apps/samtools-0.1.9/samtools pileup -cf /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa $1$2/$3F3_nodup.bam > $1$2/$3$4F3_nodup.pileup
	
	# Overlapping capture intervals will get together
	for j in ${i[@]}; do
	        sg_unificar_intervalos_de_captura_solapantes.pl /data/results/Solid0065/BF11_Cardio-cancer/diseno/cardio/chr$j.bed > $1$2/$3$4chr${j}_intervalos_sondas
	done
	
	# the piled up sequences (with no duplicates) get split into the different chromosomes
	for j in ${i[@]}; do
        	awk '{if ($1=="chr'${j}'") print $0}' $1$2/$3$4F3_nodup.pileup > $1$2/$3$4chr$j.pileup
	done
	
	# the piled up sequences in range from each chromosome are obtained
	for j in ${i[@]}; do
        	sg_extraer_SNPs_fichero_tabulado.pl $1$2/$3$4chr${j}_intervalos_sondas $1$2/$3$4chr$j.pileup > $1$2/$3$4chr$j.enrango
	done

        # the piled up sequences near target from each chromosome are obtained
        for f in $1$2/$3$4chr*_intervalos_sondas; do awk '{print $1"\t"$2-100"\t"$3+100}' $f > "${f}"_near_target; done
        for j in ${i[@]}; do
                sg_extraer_SNPs_fichero_tabulado.pl $1$2/$3$4chr${j}_intervalos_sondas_near_target $1$2/$3$4chr$j.pileup > $1$2/$3$4chr$j.near_target
        done

	# Call of SNPs without duplicates and quality values over 20
	#/share/apps/samtools-0.1.9/misc/samtools.pl varFilter -D 500000 $1$2/$3$4PE_Q10_nodup.pileup | awk '$6>=20' | awk '{if ($4!="*") print $0}' > $1$2/$3$4snps_sin_duplicados_Q20
	
	# the SNPs get split into the different chromosomes
	#for j in ${i[@]}; do
        #	awk '{if ($1=="chr'${j}'") print $0}' $1$2/$3$4snps_sin_duplicados_Q20  > $1$2/$3$4chr$j.snps
	#done
	
	# the SNPs in range from each chromosome are obtained
	#for j in ${i[@]}; do
	#        sg_extraer_SNPs_fichero_tabulado.pl $1$2/$3$4chr${j}_intervalos_sondas $1$2/$3$4chr$j.snps > $1$2/$3$4chr${j}_snps.enrango
	#done
	
	# The in range SNPs from the differente chromosomes are merged into one file (todos_SNP)
	#cat $1$2/$3$4chr1_snps.enrango $1$2/$3$4chr10_snps.enrango $1$2/$3$4chr11_snps.enrango $1$2/$3$4chr12_snps.enrango $1$2/$3$4chr13_snps.enrango $1$2/$3$4chr14_snps.enrango $1$2/$3$4chr15_snps.enrango $1$2/$3$4chr16_snps.enrango $1$2/$3$4chr17_snps.enrango $1$2/$3$4chr18_snps.enrango $1$2/$3$4chr19_snps.enrango $1$2/$3$4chr2_snps.enrango $1$2/$3$4chr20_snps.enrango $1$2/$3$4chr21_snps.enrango $1$2/$3$4chr22_snps.enrango $1$2/$3$4chr3_snps.enrango $1$2/$3$4chr4_snps.enrango $1$2/$3$4chr5_snps.enrango $1$2/$3$4chr6_snps.enrango $1$2/$3$4chr7_snps.enrango $1$2/$3$4chr8_snps.enrango $1$2/$3$4chr9_snps.enrango $1$2/$3$4chrM_snps.enrango $1$2/$3$4chrX_snps.enrango $1$2/$3$4chrY_snps.enrango > $1$2/$3$4todos_SNP
	
	# Conversion from heterozygote to homozygote
	#sg_convertir_SNPs_de_heterocigotos_a_homocigotos.pl $1$2/$3$4todos_SNP > $1$2/$3$4todos_SNP_OK
	
	# The SNPs get filtered
	#variant_effect_predictor_15-2-11.pl -i $1$2/$3$4todos_SNP_OK -o $1$2/$3$4SNVs.SNVs
	
	# Annotation will be added to the SNPs file
	#sg_combinar_output_samtools_anotacion_SNP_effect.pl $1$2/$3$4todos_SNP $1$2/$3$4todos_SNP_OK
	
	# Coverage results are obtained
	for f in $1$2/$3$4chr*.pileup; do echo $f; awk '{print $8}' $f > "${f}"_coverage; sg_sumatorio_columna.pl "${f}"_coverage | grep Sumatorio; done > $1$2/$3$4results_nt_mapeables
	rm $1$2/$3$4*_coverage
	for f in $(ls ${1}${2}/${3}${4}*.enrango | grep -v snps ); do echo $f; awk '{print $8}' $f > "${f}"_coverage; sg_sumatorio_columna.pl "${f}"_coverage | grep Sumatorio; done > $1$2/$3$4results_nt_enrango
	rm $1$2/$3$4*_coverage
	for f in $(ls ${1}${2}/${3}${4}*.near_target ); do echo $f; awk '{print $8}' $f > "${f}"_coverage; sg_sumatorio_columna.pl "${f}"_coverage | grep Sumatorio; done > $1$2/$3$4results_nt_neartarget
	rm $1$2/$3$4*_coverage
	for f in $(ls ${1}${2}/${3}${4}*.enrango | grep -v snps ); do echo $f; awk '{if ( $8 >= 1) print $1}' $f | wc -l; done > $1$2/$3$4results_nt_1x
	for f in $(ls ${1}${2}/${3}${4}*.enrango | grep -v snps ); do echo $f; awk '{if ( $8 >= 10) print $1}' $f | wc -l; done > $1$2/$3$4results_nt_10x
	for f in $(ls ${1}${2}/${3}${4}*.enrango | grep -v snps ); do echo $f; awk '{if ( $8 >= 20) print $1}' $f | wc -l; done > $1$2/$3$4results_nt_20x
	
	# SNVs results are joined together
	#echo resultados de SNVs > $1$2/$3$4results_SNVs
	#echo SNVs totales: >> $1$2/$3$4results_SNVs
	#wc $1$2/$3$4todos_SNP >> $1$2/$3$4results_SNVs
	#echo SNVs nuevos: >> $1$2/$3$4results_SNVs
	#for f in $(ls ${1}${2}/${3}${4}SNVs.SNVs); do awk '{if ($9=="-") print $1"_"$2}' $f | sort -u | wc -l; done >> $1$2/$3$4results_SNVs
	#echo SNVs no-sinónimos: >> $1$2/$3$4results_SNVs
	#for f in $(ls ${1}${2}/${3}${4}SNVs.SNVs); do grep NON_SYNONYMOUS_CODING $f | awk '{print $1"_"$2}' | sort -u | wc -l; done >> $1$2/$3$4results_SNVs
	#echo SNVs zonas de splicing: >> $1$2/$3$4results_SNVs
	#for f in $(ls ${1}${2}/${3}${4}SNVs.SNVs); do grep SPLICE $f | awk '{print $1"_"$2}' | sort -u | wc -l; done >> $1$2/$3$4results_SNVs
	#echo SNVs nuevos no-sinónimos: >> $1$2/$3$4results_SNVs
	#for f in $(ls ${1}${2}/${3}${4}SNVs.SNVs); do grep NON_SYNONYMOUS_CODING $f | awk '{if ($9=="-") print $1"_"$2}' | sort -u | wc -l; done >> $1$2/$3$4results_SNVs
	#echo SNVs nuevos zonas de splicing: >> $1$2/$3$4results_SNVs
	#for f in $(ls ${1}${2}/${3}${4}SNVs.SNVs); do grep SPLICE $f | awk '{if ($9=="-") print $1"_"$2}' | sort -u | wc -l; done >> $1$2/$3$4results_SNVs
}


# -----------------------
# - SNVs_indels_calling -
# -----------------------

SNVs_indels_calling ()
# Description: function which carries out the SNVs and indels calling
# Parameters:
#       $1 : root directory (ending with /)
#       $2 : sample name
#       $3 : input directory (ending with /)
#       $4 : output directory sufix (ending with /)
{
        # ---------------
        # - Definitions -
        # ---------------

        # Subindexes for the different chromosomes
        i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 M X Y)

        # ------------------------
        # - body of the pipeline -
        # ------------------------


	samtools view -H solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Cardio_F3_BM6091.csfasta.ma.bam > header

	samtools view -F 4 -b  solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Cardio_F3_BM6091.csfasta.ma.bam  > F3-fragmentH_sort_modify.bam
	
	samtools index F3-fragmentH_sort_modify.bam
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I F3-fragmentH_sort_modify.bam  -o intervals.intervals -B:dbsnp,dbsnp /data/results/Solid0065/BF11_Cardio-cancer/datos_partida/datos_partidaJC/selected131_sort.rod
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T IndelRealigner -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I F3-fragmentH_sort_modify.bam  -targetIntervals intervals.intervals -o realigner.bam -B:dbsnp,dbsnp /data/results/Solid0065/BF11_Cardio-cancer/datos_partida/datos_partidaJC/selected131_sort.rod
	
	samtools view realigner.bam > realigner.sam
	
	perl script_ColorSpace.pl realigner.sam > realigner_withoutIndetermination.sam
	
	cat header realigner_withoutIndetermination.sam > oeoe
	
	mv oeoe realigner_withoutIndetermination.sam
	
	samtools view -b -S realigner_withoutIndetermination.sam > realigner_withoutIndetermination.bam

	samtools sort realigner_withoutIndetermination.bam  ejm
	
	mv ejm.bam realigner_withoutIndetermination.bam
	
	samtools index realigner_withoutIndetermination.bam
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T CountCovariates  -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa --DBSNP /data/results/Solid0065/BF11_Cardio-cancer/datos_partida/datos_partidaJC/selected131_sort.rod -I realigner_withoutIndetermination.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile recal.csv
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T TableRecalibration  -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -I realigner_withoutIndetermination.bam -o realigner_withoutIndetermination_Recalibrate.bam -recalFile recal.csv
	
	samtools view -b  -q 10 realigner_withoutIndetermination_Recalibrate.bam  > realigner_withoutIndetermination_Recalibrate_Q10.bam
	
	MarkDuplicates.jar  I=realigner_withoutIndetermination_Recalibrate_Q10.bam  O=realigner_withoutIndetermination_Recalibrate_Q10_nodup.bam  M=stat  REMOVE_DUPLICATES=TRUE ASSUME_SORTED=FALSE
	
	samtools index realigner_withoutIndetermination_Recalibrate_Q10_nodup.bam
	
	mkdir Indels
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -I realigner_withoutIndetermination_Recalibrate_Q10_nodup.bam  -glm DINDEL  -o Indels/indels.vcf
	
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -B:variant,VCF Indels/indels.vcf -o Indels/fist_filtration_indels.vcf --filterExpression "MQ0 >=4 && ((MQ0 / (1.0 * DP ))> 0.1)" --filterName "FIRST_FILTER" --filterExpression "SB >= -1.0" --filterName "STRAND_BIAS_FILTER" --filterExpression "QUAL < 10" --filterName "QUALFILTER"
	
	mkdir SNPs
}


# -----------------------
# - export_results -
# -----------------------

export_results ()
# Description: function to export the data for its postprocessing with other applications (excel vba)
# Parameters:
#       $1 : root directory (ending with /)
#       $2 : sample name
#       $3 : input directory (ending with /)
#       $4 : output directory sufix (ending with /)
{
        # ---------------
        # - Definitions -
        # ---------------


        # ------------------------
        # - body of the pipeline -
        # ------------------------

	if [ ! -d "$1coverage_$4" ]; then
                # Control will enter here if $DIRECTORY doesn't exist
                mkdir $1coverage_$4
        fi
        if [ ! -d "$1coverage_$4$2" ]; then
                # Control will enter here if $DIRECTORY doesn't exist
                mkdir $1coverage_$4$2
        fi

	cp $1$2/$3$4results_nt_* $1coverage_$4$2
}


# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

# it will be checked if the input path finishes with "/", and the recurse function will be called in consequence.
root_dir="/data/results/Solid0065/BF11_Cardio-cancer/cardio/"
samples=("BM3062" "BM3339" "BM3895" "BM4237" "BM5307" "BM5357" "BM6091" "BM6092" "BM6492" "BM6919" "BM8188" "BM8190")
# the input_dir is based on the input_dir, this means: root_dir/samples/input_dir
input_dir=("secondary/bioscope/maToBam/" "secondary/bioscope/maToBam/" "secondary/bioscope/maToBam/" "secondary/bioscope/maToBam/" "secondary/bioscope/maToBam/" "secondary/bioscope110214/output/F3/maToBam/" "secondary/bioscope110214/output/F3/maToBam/" "secondary/bioscope110214/output/F3/maToBam/" "secondary/bioscope20110214/output/F3/maToBam/" "secondary/bioscope110214/output/F3/maToBam/" "secondary/bioscope_110213/output/F3/maToBam/" "secondary/bioscope_110207/output/F3/maToBam/")
# the output_dir is based on the input_dir (and based dir), this means: root_dir/input_dir/samples_output_dir
output_dir="stats_F3_nodup/"

for ((k=0;k<${#samples[@]};k+=1)); do
	mapping_metrics $root_dir ${samples[k]} ${input_dir[k]} $output_dir
	export_results $root_dir ${samples[k]} ${input_dir[k]} $output_dir
done

#if [ ${1:${#1}-1} = "/" ]; then
#	recurse "${1:0:${#1}-1}"
#else
#	recurse "${1}"
#fi

