#!/bin/bash

##PIPELINE CANCER $1 BAM OUTPUT FROM BIOSCOPE PAIRED END $2 BAM OUTPUT FROM BIOSCOPE FRAGMENT  $3 Output Bioscope Small indel (*gff) in Paired End $4 Output Bioscope Small Indels in Fragment (*.gff) $5 Type Analysis (Exoma38 Exoma50 Cardio Cancer)

#La variable 'i' contiene el nombre de los cromosomas de todo el genoma humano. En BRCAs con 454 los cromosomas solo son el 13 y el 17, en los paneles de cardio y onco estan representados la mayoria de los cromosomas y por supuesto en el exoma completo estan todos los cromosomas excepto el mitocondrial, aunque si es cierto que hay lecturas que mapean de forma aleatoria en el genoma mitocondrial.

echo $(date)

preprocess ()
{
	#Ordenar bam
	SortSam.jar  I=$1 O=input_PicardSort.bam SO=coordinate TMP_DIR=/scratch
	#Eliminar duplicados 
	MarkDuplicates.jar  I=input_PicardSort.bam  O=input_PicardSort_nodup.bam  M=stat  TMP_DIR=/scratch REMOVE_DUPLICATES=TRUE ASSUME_SORTED=FALSE
	#Filtrar lecturas por valor de calidad de mapeo (Phred=1 --> -q 1) para eliminar las lecturas que mapean en multiples sitios
	samtools view -q 1 -b input_PicardSort_nodup.bam > input_PicardSort_nodup_Q1.bam
	#Indexar ficheros bam
	samtools index input_PicardSort_nodup.bam
	samtools index input_PicardSort_nodup_Q1.bam
	
	SortSam.jar  I=$2 O=F3_input_PicardSort.bam SO=coordinate TMP_DIR=/scratch
	MarkDuplicates.jar  I=F3_input_PicardSort.bam  O=F3_input_PicardSort_nodup.bam  M=stat  TMP_DIR=/scratch REMOVE_DUPLICATES=TRUE ASSUME_SORTED=FALSE
	samtools view -q 1 -b F3_input_PicardSort_nodup.bam > F3_input_PicardSort_nodup_Q1.bam
	samtools index F3_input_PicardSort_nodup.bam
	samtools index F3_input_PicardSort_nodup_Q1.bam
}

mapping_metrics ()
{
        echo "--------------Calculating metrics--------"
	#touch intervalos_capturas_original.bed
	#touch intervalos_capturas_original.bed.tmp

        samtools pileup -cf /share/references/genomes/human/hg19/reference/human_hg19.fa  input_PicardSort_nodup.bam  > input_PicardSort_nodup.pileup
	sed 's/chr//' input_PicardSort_nodup.pileup > input_PicardSort_nodup.pileup_tmp
        for j in ${i[@]}; do
                awk '{if ($1=="'$j'") print $0}' input_PicardSort_nodup.pileup_tmp > $j\_input_PicardSort_nodup.pileup

                sg_extraer_SNPs_fichero_tabulado.pl   $6$j$7  $j\_input_PicardSort_nodup.pileup  > $j\_input_PicardSort_nodup_enrango_tmp

                grep -v "*" $j\_input_PicardSort_nodup_enrango_tmp > $j\_input_PicardSort_nodup_enrango
                awk '{print $8}' $j\_input_PicardSort_nodup_enrango_tmp > o
                sg_sumatorio_columna.pl o | grep Sumatorio | awk '{print "coverage_en_rango\t"'$j'"\t"$0} ' >> mapping_metric_stats
                awk '{if ( $8 >= 1) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas1_chr"'$j'"\t"$0}' >> mapping_metric_stats
                awk '{if ( $8 >= 10) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas10_chr"'$j'"\t"$0}' >> mapping_metric_stats
                awk '{if ( $8 >= 20) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas20_chr"'$j'"\t"$0}' >> mapping_metric_stats

		cat intervalos_capturas_original.bed.tmp $6$j$7 >> intervalos_capturas_original.bed.tmp2
		mv intervalos_capturas_original.bed.tmp2 intervalos_capturas_original.bed.tmp
        done
	mv intervalos_capturas_original.bed.tmp intervalos_capturas_original.bed
        rm o *.pileup *enrango *enrango_tmp
}

SNVs_indels_calling ()
{
	        echo "-------CALLING VARIANTS-------"
#       echo "----------Creating intervals to realign reads-------------"
#       java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -I $1_sorted.bam  -o intervals.intervals -B:dbsnp,dbsnp /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod
#       echo "----------Realigning reads-----------------"
#       java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T IndelRealigner -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta -I $1_sorted.bam  -targetIntervals intervals.intervals -o stepII_realigner.bam  -B:dbsnp,dbsnp /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod
#       echo "----------CountCovariates-----------"
#       java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T CountCovariates  -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta --DBSNP /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/rod.rod  -I stepII_realigner.bam -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile recal.csv
#       echo "-----------Recalibrating------"
#       java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T TableRecalibration  -R /data/results/Solid0065/BF42_BRCA_muestrasJB/datos_partida/referencia/refOK.fasta  -I stepII_realigner.bam -o stepII_realigner_sort_Recalibrate.bam -recalFile recal.csv





        echo "------------Variant calling Paired End--------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -glm SNP -nt 8  -stand_call_conf 20.0 -stand_emit_conf 20.0 -dcov 200 -I input_PicardSort_nodup_Q1.bam  -o snps_GATK.vcf

	echo "------------Variant calling Fragment--------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -glm SNP -nt 8  -stand_call_conf 20.0 -stand_emit_conf 20.0 -dcov 200 -I F3_input_PicardSort_nodup_Q1.bam  -o F3_snps_GATK.vcf

        echo "------------VariantFiltration Paired End--------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T VariantFiltration -R /share/references/genomes/human/hg19/reference/human_hg19.fa -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -filterName HRUN_FILTER  -B:variant,VCF snps_GATK.vcf -o snps_GATK_filtrados.vcf
	
	echo "------------VariantFiltration Fragment--------"

	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian/GenomeAnalysisTK.jar -T VariantFiltration -R /share/references/genomes/human/hg19/reference/human_hg19.fa -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0"  -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -filterName HRUN_FILTER  -B:variant,VCF F3_snps_GATK.vcf -o F3_snps_GATK_filtrados.vcf

        echo "-----------IndelCalling Paired End--------------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian_old/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -I input_PicardSort_nodup_Q1.bam  -glm DINDEL -nt 8 -o indels_GATK.vcf

	echo "-----------IndelCalling Fragment--------------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian_old/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -I F3_input_PicardSort_nodup.bam  -glm DINDEL -nt 8 -o F3_indels_GATK.vcf

        echo "-----------IndelFiltration Paired End----------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian_old/GenomeAnalysisTK.jar -T VariantFiltration -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -filterName HRUN_FILTER  -cluster 3 -window 10 -B:variant,VCF indels_GATK.vcf -o indels_GATK_filtrados.vcf

	echo "-----------IndelFiltration Fragment----------"
	java -Djava.io.tmpdir=/scratch -jar /share/apps/local/GATK_whiterussian_old/GenomeAnalysisTK.jar -T VariantFiltration -R /share/references/genomes/human/hg19/reference/human_hg19.fa  -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName HRUN_FILTER  -cluster 3 -window 10  -B:variant,VCF F3_indels_GATK.vcf -o F3_indels_GATK_filtrados.vcf

	grep "#" indels_GATK_filtrados.vcf > header
 
	echo "-----------Calling SNVs and Indels with Samtools Paired End-------"
       samtools mpileup -g -uf /share/references/genomes/human/hg19/reference/human_hg19.fa   input_PicardSort_nodup.bam  | bcftools view -bvcg - > var.raw.bcf
       bcftools view var.raw.bcf | vcfutils.pl varFilter -D2000 > SNPs_INDELS_samtools.vcf

	echo "-----------Calling SNVs and Indels with Samtools Fragment-------"
	samtools mpileup -g -uf /share/references/genomes/human/hg19/reference/human_hg19.fa F3_input_PicardSort_nodup.bam | bcftools view -bvcg - > F3_var.raw.bcf
	bcftools view F3_var.raw.bcf | vcfutils.pl varFilter -D2000 > F3_SNPs_INDELS_samtools.vcf

	echo "chromosomes = "$i
	

	echo "----------Bioscope Small Indels Preprocess ------------"

	sg_convertir_GFF_en_VCF.sh $4 bioscope_INDELS.vcf 
	sg_convertir_GFF_en_VCF.sh $5 F3_bioscope_INDELS.vcf 

	echo "-----------Filtering PASS SNPs and Indels---------"
        for x in ${i[@]};
        do


                awk '{if ($1=="chr'$x'") print $0}' indels_GATK_filtrados.vcf > $x\_indels_GATK_filtrados.vcf
                awk '{if ($1=="chr'$x'") print $0}' snps_GATK_filtrados.vcf > $x\_snps_GATK_filtrados.vcf
                awk '{if ($1=="chr'$x'") print $0}' SNPs_INDELS_samtools.vcf > $x\_SNPs_INDELS_samtools.vcf_tmp
		awk '{if ($1=="chr'$x'") print $0}' bioscope_INDELS.vcf > $x\_bioscope_INDELS.vcf

		awk '{if ($1=="chr'$x'") print $0}' F3_indels_GATK_filtrados.vcf > $x\_F3_indels_GATK_filtrados.vcf
		awk '{if ($1=="chr'$x'") print $0}' F3_snps_GATK_filtrados.vcf > $x\_F3_snps_GATK_filtrados.vcf
		awk '{if ($1=="chr'$x'") print $0}' F3_SNPs_INDELS_samtools.vcf > $x\_F3_SNPs_INDELS_samtools.vcf_tmp
		awk '{if ($1=="chr'$x'") print $0}' F3_bioscope_INDELS.vcf > $x\_F3_bioscope_INDELS.vcf

                grep "INDEL;" $x\_SNPs_INDELS_samtools.vcf_tmp > $x\_INDELS_samtools.vcf_tmp
		grep "INDEL;" $x\_F3_SNPs_INDELS_samtools.vcf_tmp > $x\_F3_INDELS_samtools.vcf_tmp

                grep -v "INDEL;" $x\_SNPs_INDELS_samtools.vcf_tmp > $x\_SNPs_samtools.vcf_tmpx
		grep -v "INDEL;" $x\_F3_SNPs_INDELS_samtools.vcf_tmp > $x\_F3_SNPs_samtools.vcf_tmpx
                grep -v "#" $x\_SNPs_samtools.vcf_tmpx > $x\_SNPs_samtools.vcf_tmp
		grep -v "#" $x\_F3_SNPs_samtools.vcf_tmpx > $x\_F3_SNPs_samtools.vcf_tmp

                echo "---INDELS_samtools.vcf_tmp2--"
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7 $x\_INDELS_samtools.vcf_tmp indels > $x\_INDELS_samtools.vcf_tmp2
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7 $x\_F3_INDELS_samtools.vcf_tmp indels > $x\_F3_INDELS_samtools.vcf_tmp2

                echo "---SNPs_samtools.vcf_tmp2--"
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_SNPs_samtools.vcf_tmp snps > $x\_SNPs_samtools.vcf_tmp2
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_F3_SNPs_samtools.vcf_tmp snps > $x\_F3_SNPs_samtools.vcf_tmp2

                echo "----indel_GATK.vcf_tmp2--"
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_indels_GATK_filtrados.vcf indels > $x\_indels_GATK_filtrados.vcf_tmp2
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_F3_indels_GATK_filtrados.vcf indels > $x\_F3_indels_GATK_filtrados.vcf_tmp2

                echo "----SNPs_GATK.vcf_tmp2--"
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_snps_GATK_filtrados.vcf snps >  $x\_snps_GATK_filtrados.vcf_tmp2
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_F3_snps_GATK_filtrados.vcf snps >  $x\_F3_snps_GATK_filtrados.vcf_tmp2

		echo "----indel_bioscope_tmp2--"
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_bioscope_INDELS.vcf indels > $x\_bioscope_INDELS.vcf_tmp2
		sg_extraer_SNPs_fichero_tabulado28032011.pl $6$x$7  $x\_F3_bioscope_INDELS.vcf indels > $x\_F3_bioscope_INDELS.vcf_tmp2


                grep "PASS" $x\_indels_GATK_filtrados.vcf_tmp2 > $x\_INDELS_GATK_final.vcf_tmp3
		grep "PASS" $x\_F3_indels_GATK_filtrados.vcf_tmp2 > $x\_F3_INDELS_GATK_final.vcf_tmp3

                grep "PASS" $x\_snps_GATK_filtrados.vcf_tmp2 >  $x\_SNPs_GATK_final.vcf_tmp3
		grep "PASS" $x\_F3_snps_GATK_filtrados.vcf_tmp2 >  $x\_F3_SNPs_GATK_final.vcf_tmp3

		mv $x\_INDELS_samtools.vcf_tmp2 $x\_INDELS_samtools_final.vcf_tmp3
		mv $x\_SNPs_samtools.vcf_tmp2 $x\_SNPs_samtools_final.vcf_tmp3
		mv $x\_F3_INDELS_samtools.vcf_tmp2 $x\_F3_INDELS_samtools_final.vcf_tmp3
		mv $x\_F3_SNPs_samtools.vcf_tmp2 $x\_F3_SNPs_samtools_final.vcf_tmp3


		cat header $x\_INDELS_GATK_final.vcf_tmp3 | vcf-sort > $x\_INDELS_GATK_tmp5
		cat header $x\_F3_INDELS_GATK_final.vcf_tmp3 | vcf-sort > $x\_F3_INDELS_GATK_tmp5

                cat header $x\_SNPs_GATK_final.vcf_tmp3 | vcf-sort > $x\_SNPs_GATK_tmp5
		cat header $x\_F3_SNPs_GATK_final.vcf_tmp3 | vcf-sort > $x\_F3_SNPs_GATK_tmp5

                cat header $x\_INDELS_samtools_final.vcf_tmp3 | vcf-sort > $x\_INDELS_samtools_tmp5
		cat header $x\_F3_INDELS_samtools_final.vcf_tmp3 | vcf-sort > $x\_F3_INDELS_samtools_tmp5

                cat header $x\_SNPs_samtools_final.vcf_tmp3 | vcf-sort > $x\_SNPs_samtools_tmp5
		cat header $x\_F3_SNPs_samtools_final.vcf_tmp3 | vcf-sort > $x\_F3_SNPs_samtools_tmp5

		cat header $x\_bioscope_INDELS.vcf_tmp2 | vcf-sort > $x\_bioscope_INDELS_tmp5
		cat header $x\_F3_bioscope_INDELS.vcf_tmp2 | vcf-sort > $x\_F3_bioscope_INDELS_tmp5


                bgzip $x\_INDELS_GATK_tmp5
 		bgzip $x\_F3_INDELS_GATK_tmp5

		bgzip $x\_SNPs_GATK_tmp5
		bgzip $x\_F3_SNPs_GATK_tmp5

		bgzip $x\_INDELS_samtools_tmp5
		bgzip $x\_F3_INDELS_samtools_tmp5

		bgzip $x\_SNPs_samtools_tmp5
		bgzip $x\_F3_SNPs_samtools_tmp5

		bgzip $x\_bioscope_INDELS_tmp5
		bgzip $x\_F3_bioscope_INDELS_tmp5


                tabix -f -p vcf $x\_INDELS_GATK_tmp5.gz
		tabix -f -p vcf $x\_F3_INDELS_GATK_tmp5.gz
		
                tabix -f -p vcf $x\_SNPs_GATK_tmp5.gz
		tabix -f -p vcf $x\_F3_SNPs_GATK_tmp5.gz

                tabix -f -p vcf $x\_INDELS_samtools_tmp5.gz
		tabix -f -p vcf $x\_F3_INDELS_samtools_tmp5.gz

                tabix -f -p vcf $x\_SNPs_samtools_tmp5.gz
		tabix -f -p vcf $x\_F3_SNPs_samtools_tmp5.gz

		tabix -f -p vcf $x\_bioscope_INDELS_tmp5.gz 
		tabix -f -p vcf $x\_F3_bioscope_INDELS_tmp5.gz


		# If one of the SNPS/indels files is empty, then its name shouldn't be passed as argument to vcf-isec:
                if [ -s "${x}_INDELS_GATK_final.vcf_tmp3" ]; then
                        file1="${x}_INDELS_GATK_tmp5.gz"
                else
                        file1=""
                fi
                if [ -s "${x}_INDELS_samtools_final.vcf_tmp3" ]; then
                        file2="${x}_INDELS_samtools_tmp5.gz"
                else
                        file2=""
                fi
                if [ -s "${x}_bioscope_INDELS.vcf_tmp2" ]; then
                        file3="${x}_bioscope_INDELS_tmp5.gz"
                else
                        file3=""
                fi
                if [ -s "${x}_F3_INDELS_GATK_final.vcf_tmp3" ]; then
                        file4="${x}_F3_INDELS_GATK_tmp5.gz"
                else
                        file4=""
                fi
                if [ -s "${x}_F3_INDELS_samtools_final.vcf_tmp3" ]; then
                        file5="${x}_F3_INDELS_samtools_tmp5.gz"
                else
                        file5=""
                fi
                if [ -s "${x}_F3_bioscope_INDELS.vcf_tmp2" ]; then
                        file6="${x}_F3_bioscope_INDELS_tmp5.gz"
                else
                        file6=""
                fi

		echo "vcf-isec -n +1 $file1 $file2 $file3 $file4 $file5 $file6 > ${x}_INDELS_END.vcf"
		vcf-isec -n +1 $file1 $file2 $file3 $file4 $file5 $file6 > ${x}_INDELS_END.vcf

                if [ -s "${x}_SNPs_GATK_final.vcf_tmp3" ]; then
                        file1="${x}_SNPs_GATK_tmp5.gz"
                else
                        file1=""
                fi
                if [ -s "${x}_SNPs_samtools_final.vcf_tmp2" ]; then
                        file2="${x}_SNPs_samtools_tmp5.gz"
                else
                        file2=""
                fi
                if [ -s "${x}_F3_SNPs_GATK_final.vcf_tmp3" ]; then
                        file3="${x}_F3_SNPs_GATK_tmp5.gz"
                else
                        file3=""
                fi
                if [ -s "${x}_F3_SNPs_samtools_final.vcf_tmp2" ]; then
                        file4="${x}_F3_SNPs_samtools_tmp5.gz"
                else
                        file4=""
                fi

                vcf-isec -n +1 $file1 $file2 $file3 $file4 > ${x}_SNPs_END.vcf

	done
	
	# for test purposes, the _tmp and _fitrados.vcf files won't be deleted!!!!!!!!!!!
	# rm *_tmp* *_filtrados.vcf*

        cat *INDELS_END.vcf | grep -v "#" > cambios_INDELS.vcf
	cat *SNPs_END.vcf | grep -v "#" > cambios_SNPs.vcf


	# This step is only included here for test purposes. In the future should be included in the "anotacion" function


	sg_integration_vcfFotmat.pl cambios_INDELS.vcf > cambios_INDELS_integrados.vcf ##GENERAR UN FORMATO UNIFICADO DE VCF PARA INDELS

	sg_integration_vcfFotmat.pl cambios_SNPs.vcf  > cambios_SNPs_integrados.vcf ##GENERAR UN FORMATO UNIFICADO DE VCF PARA SNPs
}

anotacion ()
{
        echo "--------------Variant annotation--------"
	echo "SNP Annotation   "$(date)
       	variant_effect_predictor_Ensembl62_111013.pl --hgvs --condel=b -i cambios_SNPs_integrados.vcf  -o variant_effect_output.txt_SNPs_tmp1 -format vcf -b 1000
	echo "Indel Annotation   "$(date)
	variant_effect_predictor_Ensembl62_111013.pl --hgvs --condel=b -i cambios_INDELS_integrados.vcf -o variant_effect_output.txt_INDELS_tmp1 -format vcf -b 1000
	

	sg_sustituir_porespacio_annotacion.pl variant_effect_output.txt_SNPs_tmp1 > variant_effect_output.txt_SNPs_tmp
	sg_sustituir_porespacio_annotacion.pl variant_effect_output.txt_INDELS_tmp1 > variant_effect_output.txt_INDELS_tmp
        
	for c in ${i[@]}; do
                awk '{if ($1=="chr'$c'") print $0}' cambios_SNPs_integrados.vcf  | sort -k2n > $c\_SNPs.vcf
		awk '{if ($1=="chr'$c'") print $0}' cambios_INDELS_integrados.vcf  | sort -k2n > $c\_INDELS.vcf

                awk '{if ($1=="'$c'") print $0}' variant_effect_output.txt_SNPs_tmp | sort -k2n > $c\_SNPs.annot
		awk '{if ($1=="'$c'") print $0}' variant_effect_output.txt_INDELS_tmp | sort -k2n > $c\_INDELS.annot

                echo "-----------Combining annotations-------------"
                sg_combinar_vcf_anotacion_APIs.pl $c\_SNPs.vcf $c\_SNPs.annot > $c\_SNPs.comb
		sg_combinar_vcf_anotacion_APIs.pl $c\_INDELS.vcf $c\_INDELS.annot > $c\_INDELS.comb

                echo "-----------Biobase annotation----------------"
                sg_extraerBiobaseInformation.pl $c\_SNPs.comb /share/references/biobase/bioBaseSortUniques_chr$c.txt > $c\_SNPs.annot.biobase
		sg_extraerBiobaseInformation.pl $c\_INDELS.comb /share/references/biobase/bioBaseSortUniques_chr$c.txt > $c\_INDELS.annot.biobase
        done

        cat *SNPs.annot.biobase | grep -v \"#\" > cambios_anotados_SNPs_tmp > cambios_anotados_tmp_SNPs.xls
	cat *INDELS.annot.biobase | grep -v \"#\" > cambios_anotados_INDELS_tmp > cambios_anotados_tmp_INDELS.xls

	grep -v -P "Gene_name\t" cambios_anotados_tmp_SNPs.xls > cambios_anotados_tmp2_SNPs.xls
	grep -v -P "Gene_name\t" cambios_anotados_tmp_INDELS.xls > cambios_anotados_tmp2_INDELS.xls

	rm *.comb
	
	mkdir tmp

	mv *.bam* tmp/
	mv cambios_INDELS_integrados.vcf tmp/
	mv cambios_SNPs_integrados.vcf tmp/
	mv variant_effect_output.txt_SNPs_tmp1 tmp/
	mv variant_effect_output.txt_INDELS_tmp1 tmp/

	rm *.vcf
	rm *.idx


	 printf "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tBiobase_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tCDS_position\taa_position\taa_change\tGERP_conservarion_score\tHGVs_ID\tCondel_Prediction\tCoverage\tGenotype_quality\tMapping_quality\tInterpro_ID\tInterpro_Description\tDisease\tHGVS_ID\tPubMed_ID\n" > header

	if [ $6 == "cardio" ]; then
               cat header  cambios_anotados_tmp2_SNPs.xls >  cardio_SNPs.xls
               cat header  cambios_anotados_tmp2_INDELS.xls >  cardio_INDELS.xls
	fi

        if [ $6 == "cancer" ]; then
		 cat header  cambios_anotados_tmp2_SNPs.xls >  cancer_SNPs.xls
		 cat header  cambios_anotados_tmp2_INDELS.xls >  cancer_INDELS.xls
        fi

	if [ $6 == "Exoma38_hg18" ]; then
                cat header  cambios_anotados_tmp2_SNPs.xls >  Exoma38_hg18_SNPs.xls
		cat header  cambios_anotados_tmp2_INDELS.xls >  Exoma38_hg18_INDELS.xls
        fi

        if [ $6 == "Exoma38_hg19" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls >  Exoma38_hg19_SNPs.xls
		cat header cambios_anotados_tmp2_INDELS.xls >  Exoma38_hg19_INDELS.xls
        fi

	if [ $6 == "Exoma50" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls >  Exoma50_SNPs.xls
		cat header cambios_anotados_tmp2_INDELS.xls >  Exoma50_INDELS.xls
        fi

        if [ $6 == "Exoma51" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls > Exoma51_SNPs.xls
                cat header cambios_anotados_tmp2_INDELS.xls > Exoma51_INDELS.xls
        fi

        if [ $6 == "Exoma71" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls > Exoma71_SNPs.xls
                cat header cambios_anotados_tmp2_INDELS.xls > Exoma71_INDELS.xls
        fi

        if [ $6 == "hipoacusias" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls >  Hipoacusias_SNPs.xls
                cat header cambios_anotados_tmp2_INDELS.xls >  Hipoacusias_INDELS.xls
        fi

        if [ $6 == "cardio2" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls >  Cardio2_SNPs.xls
                cat header cambios_anotados_tmp2_INDELS.xls >  Cardio2_INDELS.xls
        fi
	
        if [ $6 == "nf2" ]; then
                cat header cambios_anotados_tmp2_SNPs.xls >  NF2_SNPs.xls
                cat header cambios_anotados_tmp2_INDELS.xls >  NF2_INDELS.xls
        fi


        rm *_tmp* *biobase *annot *gz
}


# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

EXPECTED_ARGS=5 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


##PIPELINE CANCER $1 BAM OUTPUT FROM BIOSCOPE PAIRED END $2 BAM OUTPUT FROM BIOSCOPE FRAGMENT  $3 Output Bioscope Small indel (*gff) in Paired End $4 Output Bioscope Small Indels in Fragment (*.gff) $5 Type Analysis (Exoma38 Exoma50 Cardio Cancer)

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Name: sg_GeneralPipe_SNPs_Indels.sh'
        echo 'Description: General pipeline for SNPs & Indels calling'
        echo 'Mandatory parameters:'
        echo '       $1 : .bam output PE'
        echo '       $2 : .bam output F3'
	echo '       $3 : .gff output Bioscope Small Indels PE'
	echo '       $4 : .gff output Bioscope Small Indels F3'
        echo '       $5 : application type (Exoma38_hg18, Exoma38_hg19, Exoma50, Exoma51, Exoma71, cardio, cancer, hipoacusias, cardio2, nf2)'
        exit $E_BADARGS
fi


if [ $5 == "Exoma38_hg19" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4 
#	mapping_metrics $i $1 $2 $3 $4       "/share/references/intervalos_captura/exoma38_hg19/hg19_unicos_chr" "_intervalos.bed"
        SNVs_indels_calling $i $1 $2 $3 $4  "/share/references/intervalos_captura/exoma38_hg19/hg19_unicos_chr" "_intervalos.bed" 
	anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "Exoma38_hg18" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4
#	mapping_metrics $i $1 $2 $3 $4      "/share/references/intervalos_captura/exoma38_hg18/chr" "_intervalos_sondas"
        SNVs_indels_calling $i $1 $2 $3 $4  "/share/references/intervalos_captura/exoma38_hg18/chr" "_intervalos_sondas"
	anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "Exoma50" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4
#	mapping_metrics $i $1 $2 $3 $4     "/share/references/intervalos_captura/exoma50_hg19/exoma50_chr" "_intervalos_sondas"
	SNVs_indels_calling $i $1 $2 $3 $4 "/share/references/intervalos_captura/exoma50_hg19/exoma50_chr" "_intervalos_sondas"
	anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "Exoma51" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
        preprocess $1 $2 $3 $4
#       mapping_metrics $i $1 $2 $3 $4     "/share/references/intervalos_captura/exoma51_hg19/exoma51_chr" "_intervalos_sondas"
        SNVs_indels_calling $i $1 $2 $3 $4     "/share/references/intervalos_captura/exoma51_hg19/exoma51_chr" "_intervalos_sondas"
        anotacion $i $1 $2 $3 $4 $5 $6
fi

if [ $5 == "Exoma71" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
        preprocess $1 $2 $3 $4
#       mapping_metrics $i $1 $2 $3 $4     "/share/references/intervalos_captura/exoma71_hg19/nuevos-intervalos_" ""
        SNVs_indels_calling $i $1 $2 $3 $4     "/share/references/intervalos_captura/exoma71_hg19/nuevos-intervalos_" ""
        anotacion $i $1 $2 $3 $4 $5 $6
fi

if [ $5 == "cardio" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4
#	mapping_metrics $i $1 $2 $3 $4     "/share/references/intervalos_captura/cardio/chr" "_intervalos_sondas"
	SNVs_indels_calling $i $1 $2 $3 $4 "/share/references/intervalos_captura/cardio/chr" "_intervalos_sondas"
	anotacion $i $1 $2 $3 $4 $5
fi

if [ $5 == "cancer" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4
#	mapping_metrics $i $1 $2 $3 $4      "/share/references/intervalos_captura/cancer/chr" "_intervalos_sondas"
 	SNVs_indels_calling $i $1 $2 $3 $4  "/share/references/intervalos_captura/cancer/chr" "_intervalos_sondas"
	anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "hipoacusias" ]; then
        i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)
	preprocess $1 $2 $3 $4
#       mapping_metrics $i $1 $2 $3 $4      "/share/references/intervalos_captura/hipoacusias/chr" "_intervalos_sondas"
        SNVs_indels_calling $i $1 $2 $3 $4  "/share/references/intervalos_captura/hipoacusias/chr" "_intervalos_sondas"
        anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "cardio2" ]; then
	i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y MT)   	
	preprocess $1 $2 $3 $4
#       mapping_metrics $i $1 $2 $3 $4     "/share/references/intervalos_captura/cardio2/chr" "_intervalos_sondas"
        SNVs_indels_calling $i $1 $2 $3 $4 "/share/references/intervalos_captura/cardio2/chr" "_intervalos_sondas"
        anotacion $i $1 $2 $3 $4 $5 
fi

if [ $5 == "nf2" ]; then
	i=(22)
        preprocess $1 $2 $3 $4
#       mappiing_metrics $i $1 $2 $3 $4     "/data/results/Solid0065/OnGoing/DB11_long_pcr/NF2_chr" "_intervalos_sondas"
        SNVs_indels_calling $i $1 $2 $3 $4 "/data/results/Solid0065/OnGoing/DB11_long_pcr/NF2_chr" "_intervalos_sondas"
        anotacion $i $1 $2 $3 $4 $5
fi

echo $(date)



