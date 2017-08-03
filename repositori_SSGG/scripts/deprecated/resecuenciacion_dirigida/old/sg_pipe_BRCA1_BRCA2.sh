#!/bin/bash

i=(BRCA1 BRCA2)
mapping ()
{
	#echo "---------------Building reference------------"
	#/share/apps/ssaha2_v2.5.3_x86_64/ssaha2Build -save ref_BRCA1_BRCA2.fa ref_BRCA1_BRCA1.fa
	echo "--------------Executing SSAHA2-mapping-----------"
        /share/apps/ssaha2_v2.5.3_x86_64/ssaha2 -454 -outfile $1.sam  -output sam  /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/ref_BRCA1_BRCA2.fa /data/results/Solid0065/BF42_BRCA_muestrasJB/primary/454Reads.$1.fastq
	echo "--------------Formatting input for GATK----------"
        cat /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/cabecera $1.sam > $1.tmp1

        awk '{print $0"\tRG:Z:20110221001322272"}' $1.tmp1  > $1.tmp2

        samtools view -S -b $1.tmp2 > $1.bam

        samtools sort $1.bam $1_sorted

        samtools index $1_sorted.bam
}

mapping_metrics ()
{
	echo "--------------Calculating metrics--------"
	/share/apps/samtools-0.1.16/samtools pileup -cf /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta $1_sorted.bam > $1.pileup
	
	for j in ${i[@]}; do
        	awk '{if ($1=="'$j'") print $0}' $1.pileup > $j\_$1.pileup
		sg_extraer_SNPs_fichero_tabulado.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$j\_v2_intervalos.tab $j\_$1.pileup > $j\_$1.enrango_tmp
		grep -v "*" $j\_$1.enrango_tmp > $j\_$1.enrango
		awk '{print $8}' $j\_$1.enrango > o
		sg_sumatorio_columna.pl o | grep sumatorio | awk '{print "coverage_en_rango\t"'$j'"\t"$0} ' >> $1\_stats
		awk '{if ( $8 >= 1) print $1}' $j\_$1.enrango | wc -l | awk '{print "nt_coverage_mas1"'$j'"\t"$0}' >> $1\_stats
		awk '{if ( $8 >= 10) print $1}' $j\_$1.enrango | wc -l | awk '{print "nt_coverage_mas10"'$j'"\t"$0}' >> $1\_stats
		awk '{if ( $8 >= 20) print $1}' $j\_$1.enrango | wc -l | awk '{print "nt_coverage_mas20"'$j'"\t"$0}' >> $1\_stats
	done
	rm o *.pileup *enrango *enrango_tmp
}


SNVs_indels_calling ()
{
	echo "-------CALLING VARIANTS-------"
	echo "----------Creating intervals to realign reads-------------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -I $1_sorted.bam  -o intervals.intervals -B:dbsnp,dbsnp /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod
	echo "----------Realigning reads-----------------"
	java -jar  /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T IndelRealigner -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -I $1_sorted.bam  -targetIntervals intervals.intervals -o stepII_realigner.bam  -B:dbsnp,dbsnp /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod
	echo "----------CountCovariates-----------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T CountCovariates  -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta --DBSNP /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod  -I stepII_realigner.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile recal.csv
	echo "-----------Recalibrating------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T TableRecalibration  -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta  -I stepII_realigner.bam -o stepII_realigner_sort_Recalibrate.bam -recalFile recal.csv
	echo "------------Variant calling--------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta  -glm SNP -nt 8  -stand_call_conf 40.0 -stand_emit_conf 20.0 -dcov 200 -I stepII_realigner_sort_Recalibrate.bam  -o snps_GATK.vcf
	echo "------------VariantFiltration--------"
	/share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 30.0" -filter "SB > 0.1" -filter "QD < 5.0" -filter "HRun >= 4"  -filter "AF < 0.05 || AF > 0.95" -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName HRUN_FILTER -filterName AFF_FILTER  -cluster 3 -window 10 -B:variant,VCF snps_GATK.vcf -o snps_GATK_filtrados.vcf
	echo "-----------IndelCalling--------------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta  -I stepII_realigner_sort_Recalibrate.bam  -glm DINDEL -nt 8 -o indels_GATK.vcf
	echo "-----------IndelFiltration----------"
	/share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 30.0" -filter "SB > -1.0" -filter "QD < 2.0" -filter "AF < 0.05 || AF > 0.95" -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName AFF_FILTER -B:variant,VCF indels_GATK.vcf -o indels_GATK_filtrados.vcf
	echo "-----------Filtering PASS SNPs and Indels---------"
	
	for x in ${i[@]}; 
	do
		awk '{if ($1=="'$x'") print $0}' indels_GATK_filtrados.vcf > $x\_$1\_indels_GATK_filtrados.vcf
		
		sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab $x\_$1\_indels_GATK_filtrados.vcf indels > $x\_$1\_indels_GATK_filtrados.vcf_tmp
		
		grep "PASS" $x\_$1\_indels_GATK_filtrados.vcf_tmp > $x\_$1\_INDELS_GATK_final.vcf_tmp
		
		awk '{if ($1=="'$x'") print $0}' snps_GATK_filtrados.vcf > $x\_$1\_snps_GATK_filtrados.vcf
		
		sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab  $x\_$1\_snps_GATK_filtrados.vcf snps >  $x\_$1\_snps_GATK_filtrados.vcf_tmp
                
		grep "PASS"  $x\_$1\_snps_GATK_filtrados.vcf > $x\_$1\_SNPs_GATK_final.vcf_tmp
		
		if [$x eq "BRCA1"]
			then
				awk '{sub($2,$2-41277500);print}' $x\_$1\_INDELS_GATK_final.vcf_tmp > $x\_$1\_INDELS_GATK_final.vcf
				awk '{sub($2,$2-41277500);print}' $x\_$1\_SNPs_GATK_final.vcf_tmp > $x\_$1\_SNPs_GATK_final.vcf
		else
			awk '{sub($2,$2+32889617);print}' $x\_$1\_INDELS_GATK_final.vcf_tmp > $x\_$1\_INDELS_GATK_final.vcf
			awk '{sub($2,$2+32889617);print}' $x\_$1\_SNPs_GATK_final.vcf_tmp > $x\_$1\_SNPs_GATK_final.vcf		
		fi
	done
	
	echo "-----------Calling SNVs and Indels with Samtools-------"
	samtools mpileup -P 454 -uf /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta   stepII_realigner_sort_Recalibrate.bam | bcftools view -bvcg - > var.raw.bcf
	
	bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > $1_SNPs_INDELS_samtools.vcf
	
	rm *_tmp *_GATK_filtrados.vcf
}



# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------


mapping $1
mapping_metrics $1
SNVs_indels_calling $1

