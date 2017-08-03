
#!/bin/bash

##PIPELINE CANCER $1 BAM OUTPUT FROM BIOSCOPE PAIRED END $2 BAM OUTPUT FROM BIOSCOPE FRAGMENT  $3 Output Bioscope Small indel (*gff) in Paired End $BED Output Bioscope Small Indels in Fragment (*.gff) $SUFFIX Type Analysis (Exoma38 Exoma50 Cardio Cancer)

#La variable 'i' contiene el nombre de los cromosomas de todo el genoma humano. En BRCAs con 454 los cromosomas solo son el 13 y el 17, en los paneles de cardio y onco estan representados la mayoria de los cromosomas y por supuesto en el exoma completo estan todos los cromosomas excepto el mitocondrial, aunque si es cierto que hay lecturas que mapean de forma aleatoria en el genoma mitocondrial.

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

echo $(date)

preprocess ()
{
	#Ordenar bam
	/share/apps/picard-tools-1.31/SortSam.jar  I=$1 O=input_PicardSort.bam SO=coordinate TMP_DIR=/data/results/solid0065/jm_trash/
	#Eliminar duplicados 
	/share/apps/picard-tools-1.31/MarkDuplicates.jar  I=input_PicardSort.bam  O=input_PicardSort_nodup.bam  M=stat  TMP_DIR=/data/results/solid0065/jm_trash REMOVE_DUPLICATES=TRUE ASSUME_SORTED=FALSE
	#Filtrar lecturas por valor de calidad de mapeo (Phred=1 --> -q 1) para eliminar las lecturas que mapean en multiples sitios
	/share/apps/samtools-0.1.16/samtools view -q 1 -b input_PicardSort_nodup.bam > input_PicardSort_nodup_Q1.bam
	#Indexar ficheros bam
	/share/apps/samtools-0.1.16/samtools index input_PicardSort_nodup.bam
	/share/apps/samtools-0.1.16/samtools index input_PicardSort_nodup_Q1.bam
}

mapping_metrics ()
{
	echo "--------------Calculating metrics--------"
	#touch intervalos_capturas_original.bed
	rm intervalos_capturas_original.bed.tmp
	touch intervalos_capturas_original.bed.tmp

	/share/apps/samtools-0.1.16/samtools pileup -cf /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  input_PicardSort_nodup.bam  > input_PicardSort_nodup.pileup
	
	awk '$3 == "*"' input_PicardSort_nodup.pileup > input_PicardSort_nodup.samtools.indels
	
	sed 's/chr//' input_PicardSort_nodup.pileup > input_PicardSort_nodup.pileup_tmp
	for j in ${i[@]}; do
		awk '{if ($1=="'$j'") print $0}' input_PicardSort_nodup.pileup_tmp > $j\_input_PicardSort_nodup.pileup
		echo $BED$j$SUFFIX
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado.pl $BED$j$SUFFIX $j\_input_PicardSort_nodup.pileup  > $j\_input_PicardSort_nodup_enrango_tmp

		grep -v "*" $j\_input_PicardSort_nodup_enrango_tmp > $j\_input_PicardSort_nodup_enrango
		awk '{print $8}' $j\_input_PicardSort_nodup_enrango_tmp > o
		/share/apps/scripts/sg_sumatorio_columna.pl o | grep Sumatorio | awk '{print "coverage_en_rango\t"'$j'"\t"$0} ' >> mapping_metric_stats
		awk '{if ( $8 >= 1) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas1_chr"'$j'"\t"$0}' >> mapping_metric_stats
		awk '{if ( $8 >= 10) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas10_chr"'$j'"\t"$0}' >> mapping_metric_stats
		awk '{if ( $8 >= 20) print $1}' $j\_input_PicardSort_nodup_enrango | wc -l | awk '{print "nt_coverage_mas20_chr"'$j'"\t"$0}' >> mapping_metric_stats

		cat intervalos_capturas_original.bed.tmp $BED$j$SUFFIX >> intervalos_capturas_original.bed.tmp2
		mv intervalos_capturas_original.bed.tmp2 intervalos_capturas_original.bed.tmp
	done
	#mv intervalos_capturas_original.bed.tmp intervalos_capturas_original.bed
	#rm o *.pileup *enrango *enrango_tmp
}

SNVs_indels_calling ()
{
	echo "------------Variant calling Fragment--------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -glm SNP -nt 8  -stand_call_conf 20.0 -stand_emit_conf 20.0 -dcov 1000 -I input_PicardSort_nodup_Q1.bam  -o snps_GATK.vcf

	echo "------------VariantFiltration Fragment--------"
	/share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -filterName HRUN_FILTER  -B:variant,VCF snps_GATK.vcf -o snps_GATK_filtrados.vcf

	echo "-----------IndelCalling Fragment--------------"
	java -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -I input_PicardSort_nodup.bam  -glm DINDEL -nt 8 -o indels_GATK.vcf

	echo "-----------IndelFiltration Fragment----------"
	/share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T VariantFiltration -R /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa  -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 20.0" -filter "SB > 0.1" -filter "QD < 2.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -filterName HRUN_FILTER  -cluster 3 -window 10 -B:variant,VCF indels_GATK.vcf -o indels_GATK_filtrados.vcf

	grep "#" indels_GATK_filtrados.vcf > header

	echo "-----------Calling SNVs and Indels with Samtools Fragment-------"
	/share/apps/samtools-0.1.16/samtools mpileup -g -uf /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa   input_PicardSort_nodup.bam  | /share/apps/samtools-0.1.16/bcftools/bcftools view -bvcg - > var.raw.bcf
	
	/share/apps/samtools-0.1.16/bcftools/bcftools view var.raw.bcf | /share/apps/samtools-0.1.16/bcftools/vcfutils.pl varFilter -D 2000 > SNPs_INDELS_samtools.vcf_raw

	sg_parsing_samtools_indels.pl -p input_PicardSort_nodup.samtools.indels -v SNPs_INDELS_samtools.vcf_raw > SNPs_INDELS_samtools.vcf
	
	#rm input_PicardSort_nodup.samtools.indels
	
	echo "----------Bioscope Small Indels Preprocess ------------"

	/share/apps/scripts/resecuenciacion_dirigida/sg_convertir_GFF_en_VCF.sh $2 bioscope_INDELS.vcf 

	echo "-----------Filtering PASS SNPs and Indels---------"
	for x in ${i[@]}; 
	do

		awk '{if ($1=="chr'$x'") print $0}' indels_GATK_filtrados.vcf > $x\_indels_GATK_filtrados.vcf
		awk '{if ($1=="chr'$x'") print $0}' snps_GATK_filtrados.vcf > $x\_snps_GATK_filtrados.vcf
		awk '{if ($1=="chr'$x'") print $0}' SNPs_INDELS_samtools.vcf > $x\_SNPs_INDELS_samtools.vcf_tmp
		awk '{if ($1=="chr'$x'") print $0}' bioscope_INDELS.vcf > $x\_bioscope_INDELS.vcf

		grep "INDEL;" $x\_SNPs_INDELS_samtools.vcf_tmp > $x\_INDELS_samtools.vcf_tmp

		grep -v "INDEL;" $x\_SNPs_INDELS_samtools.vcf_tmp > $x\_SNPs_samtools.vcf_tmpx
		grep -v "#" $x\_SNPs_samtools.vcf_tmpx > $x\_SNPs_samtools.vcf_tmp

		echo "---INDELS_samtools.vcf_tmp2--"
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado28032011.pl $BED$x$SUFFIX $x\_INDELS_samtools.vcf_tmp indels > $x\_INDELS_samtools.vcf_tmp2

		echo "---SNPs_samtools.vcf_tmp2--"
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado28032011.pl $BED$x$SUFFIX  $x\_SNPs_samtools.vcf_tmp snps > $x\_SNPs_samtools.vcf_tmp2

		echo "----indel_GATK.vcf_tmp2--"
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado28032011.pl $BED$x$SUFFIX  $x\_indels_GATK_filtrados.vcf indels > $x\_indels_GATK_filtrados.vcf_tmp2

		echo "----SNPs_GATK.vcf_tmp2--"
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado28032011.pl $BED$x$SUFFIX  $x\_snps_GATK_filtrados.vcf snps >  $x\_snps_GATK_filtrados.vcf_tmp2

		echo "----indel_bioscope_tmp2--"
		/share/apps/scripts/sg_extraer_SNPs_fichero_tabulado28032011.pl $BED$x$SUFFIX  $x\_bioscope_INDELS.vcf indels > $x\_bioscope_INDELS.vcf_tmp2

		grep "PASS" $x\_indels_GATK_filtrados.vcf_tmp2 > $x\_INDELS_GATK_final.vcf_tmp3

		grep "PASS" $x\_snps_GATK_filtrados.vcf_tmp2 >  $x\_SNPs_GATK_final.vcf_tmp3

		mv $x\_INDELS_samtools.vcf_tmp2 $x\_INDELS_samtools_final.vcf_tmp3
		mv $x\_SNPs_samtools.vcf_tmp2 $x\_SNPs_samtools_final.vcf_tmp3

		cat header $x\_INDELS_GATK_final.vcf_tmp3 | /share/apps/vcftools_0.1.4a/perl/vcf-sort > $x\_INDELS_GATK_tmp5

		cat header $x\_SNPs_GATK_final.vcf_tmp3 | /share/apps/vcftools_0.1.4a/perl/vcf-sort > $x\_SNPs_GATK_tmp5

		cat header $x\_INDELS_samtools_final.vcf_tmp3 | /share/apps/vcftools_0.1.4a/perl/vcf-sort > $x\_INDELS_samtools_tmp5

		cat header $x\_SNPs_samtools_final.vcf_tmp3 | /share/apps/vcftools_0.1.4a/perl/vcf-sort > $x\_SNPs_samtools_tmp5

		cat header $x\_bioscope_INDELS.vcf_tmp2 | /share/apps/vcftools_0.1.4a/perl/vcf-sort > $x\_bioscope_INDELS_tmp5
	

		/data/results/Solid0065/jtrivino/tabix-0.2.5/bgzip $x\_INDELS_GATK_tmp5

		/data/results/Solid0065/jtrivino/tabix-0.2.5/bgzip $x\_SNPs_GATK_tmp5

		/data/results/Solid0065/jtrivino/tabix-0.2.5/bgzip $x\_INDELS_samtools_tmp5

		/data/results/Solid0065/jtrivino/tabix-0.2.5/bgzip $x\_SNPs_samtools_tmp5

		/data/results/Solid0065/jtrivino/tabix-0.2.5/bgzip $x\_bioscope_INDELS_tmp5

		/data/results/Solid0065/jtrivino/tabix-0.2.5/tabix -f -p vcf $x\_INDELS_GATK_tmp5.gz

		/data/results/Solid0065/jtrivino/tabix-0.2.5/tabix -f -p vcf $x\_SNPs_GATK_tmp5.gz

		/data/results/Solid0065/jtrivino/tabix-0.2.5/tabix -f -p vcf $x\_INDELS_samtools_tmp5.gz

		/data/results/Solid0065/jtrivino/tabix-0.2.5/tabix -f -p vcf $x\_SNPs_samtools_tmp5.gz

		/data/results/Solid0065/jtrivino/tabix-0.2.5/tabix -f -p vcf $x\_bioscope_INDELS_tmp5.gz 

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

		/share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file1 $file2 $file3 > ${x}_INDELS_END.vcf

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
		
		#Modificar#
		/share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file1 $file2 > ${x}_SNPs_END.vcf

	done

	# for test purposes, the _tmp and _fitrados.vcf files won't be deleted!!!!!!!!!!!
	#rm *_tmp* *_filtrados.vcf*

	cat *INDELS_END.vcf *SNPs_END.vcf | grep -v "#" > cambios.vcf
}

realignment () {

	echo "--------------Indel relalignment--------"
	sg_integration_vcfFotmat.pl cambios.vcf > cambios_integrados.vcf_raw 

	java -Djava.io.tmpdir=/data/results/solid0065/jm_trash/ -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T UnifiedGenotyper -R  /data/results/Solid0065/referencia_genoma_humano/GRch37.58/hg19.fa -nt 4 -I input_PicardSort_nodup_Q1.bam  -o allbases.vcf_raw --output_mode EMIT_ALL_SITES -dcov 1000
	
	sg_realign_indel_vcf.pl -v cambios_integrados.vcf_raw -p allbases.vcf_raw -j > cambios_integrados.vcf
	
#	echo "--------------Pileup per base--------"
#	sg_GATK_pileup_junior.pl -i $1\_allbases.vcf -o $1\_GATK
	
	if [ -d pileup ]; then 
		mv *_GATK.pileup pileup/
		mv *_allbases.vcf pileup/
  	else  
		mkdir pileup
		mv *_GATK.pileup pileup/
		mv *_allbases.vcf pileup/
	fi
	rm *raw*
}
	
anotacion ()
{

	echo "--------------Variant annotation--------"
	/share/apps/scripts/variant_effect_predictor_Ensembl62_111013.pl --hgvs --condel=b -i cambios_integrados.vcf -o variant_effect_output.txt_tmp1 -format vcf -b 1000

	/share/apps/scripts/sg_sustituir_porespacio_annotacion.pl variant_effect_output.txt_tmp1 > variant_effect_output.txt_tmp

	for c in ${i[@]}; do
		awk '{if ($1=="chr'$c'") print $0}' cambios_integrados.vcf  | sort -k2n > $c.vcf

		awk '{if ($1=="'$c'") print $0}' variant_effect_output.txt_tmp | sort -k2n > $c.annot

		echo "-----------Combining annotations-------------"
		/share/apps/scripts/sg_combinar_vcf_anotacion_APIs.pl $c.vcf $c.annot > $c.comb

		echo "-----------Biobase annotation----------------"
		/share/apps/scripts/sg_extraerBiobaseInformation.pl $c.comb /data/results/Solid0065/info_resecuenciacion/databases/HGMD-Biobase/bioBaseSortUniques_chr$c.txt > $c.annot.biobase
	done

	cat *.annot.biobase | grep -v \"#\" > cambios_anotados_tmp.xls

	grep -v -P "Gene_name\t" cambios_anotados_tmp.xls > cambios_anotados_tmp2.xls

	rm *.comb

	mkdir tmp

	mv *.bam* tmp/
	mv cambios_integrados.vcf tmp/
	mv variant_effect_output.txt_tmp1 tmp/

	rm *.vcf
	rm *.idx


	printf "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tCDS_position\taa_position\taa_change\tGERP_conservarion_score\tHGVs_ID\tCondel_Prediction\tCoverage\tGenotype_quality\tMapping_quality\tInterprot_ID\tInterprot_Description\tDisease\tHGVS_ID\tPubMed_ID\n" > header

	#Ensembl v59

	if [ $3 == "cardio" ]; then
		awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tHGVs_ID\tCondel_Score\tCoverage\tGenotype_quality\tMapping_quality\tInterpro_ID\tInterpro_Description\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000070846") || ($11=="ENST00000155840") || ($11=="ENST00000200639") || ($11=="ENST00000211998") || ($11=="ENST00000217381") || ($11=="ENST00000223528") || ($11=="ENST00000228841") || ($11=="ENST00000232975") || ($11=="ENST00000238682") || ($11=="ENST00000243457") || ($11=="ENST00000252321") || ($11=="ENST00000261201") || ($11=="ENST00000261349") || ($11=="ENST00000261448") || ($11=="ENST00000261590") || ($11=="ENST00000262186") || ($11=="ENST00000262464") || ($11=="ENST00000262631") || ($11=="ENST00000265968") || ($11=="ENST00000266732") || ($11=="ENST00000269881") || ($11=="ENST00000271348") || ($11=="ENST00000280904") || ($11=="ENST00000281456") || ($11=="ENST00000282541") || ($11=="ENST00000287878") || ($11=="ENST00000290310") || ($11=="ENST00000290378") || ($11=="ENST00000292327") || ($11=="ENST00000295754") || ($11=="ENST00000299328") || ($11=="ENST00000299333") || ($11=="ENST00000306077") || ($11=="ENST00000307128") || ($11=="ENST00000309889") || ($11=="ENST00000310128") || ($11=="ENST00000316623") || ($11=="ENST00000324501") || ($11=="ENST00000324727") || ($11=="ENST00000330010") || ($11=="ENST00000333535") || ($11=="ENST00000334785") || ($11=="ENST00000337385") || ($11=="ENST00000337851") || ($11=="ENST00000342992") || ($11=="ENST00000343849") || ($11=="ENST00000344887") || ($11=="ENST00000348997") || ($11=="ENST00000354410") || ($11=="ENST00000355349") || ($11=="ENST00000356239") || ($11=="ENST00000356287") || ($11=="ENST00000357077") || ($11=="ENST00000357525") || ($11=="ENST00000357998") || ($11=="ENST00000361308") || ($11=="ENST00000361594") || ($11=="ENST00000366574") || ($11=="ENST00000366578") || ($11=="ENST00000366783") || ($11=="ENST00000367317") || ($11=="ENST00000367319") || ($11=="ENST00000367895") || ($11=="ENST00000369519") || ($11=="ENST00000371372") || ($11=="ENST00000372066") || ($11=="ENST00000372980") || ($11=="ENST00000373960") || ($11=="ENST00000375985") || ($11=="ENST00000375994") || ($11=="ENST00000376480") || ($11=="ENST00000377329") || ($11=="ENST00000379802") || ($11=="ENST00000392770") || ($11=="ENST00000393931") || ($11=="ENST00000395869") || ($11=="ENST00000396576") || ($11=="ENST00000399249") || ($11=="ENST00000399655") || ($11=="ENST00000403994") || ($11=="ENST00000405093") || ($11=="ENST00000432085") || ($11=="ENST00000432168") || ($11=="ENST00000433631") || ($11=="ENST00000452339") || ($11=="ENST00000508053") || ($11=="ENST00000525550") || ($11=="ENST00000545968")) print $0}' cambios_anotados_tmp2.xls > cardio.xls
	fi
	
	if [ $3 == "cancer" ]; then
	#Ensembl v59
		awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tHGVs_ID\tCondel_Score\tCoverage\tGenotype_quality\tMapping_quality\tInterpro_ID\tInterpro_Description\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000257430") || ($11=="ENST00000278616") || ($11=="ENST00000380152") || ($11=="ENST00000234420") || ($11=="ENST00000231790") || ($11=="ENST00000233146") || ($11=="ENST00000371953") || ($11=="ENST00000326873") || ($11=="ENST00000260947") || ($11=="ENST00000441310") || ($11=="ENST00000265849") || ($11=="ENST00000269305") || ($11=="ENST00000372115") || ($11=="ENST00000265433") || ($11=="ENST00000261769") || ($11=="ENST00000323929") || ($11=="ENST00000265335") || ($11=="ENST00000328354") || ($11=="ENST00000309486") || ($11=="ENST00000261584") || ($11=="ENST00000259008") || ($11=="ENST00000337432") || ($11=="ENST00000357654") || ($11=="ENST00000380152")  || ($11=="ENST00000404276")  || ($11=="ENST00000261584")  || ($11=="ENST00000259008")  || ($11=="ENST00000445888")  || ($11=="ENST00000371953")  || ($11=="ENST00000326873")  || ($11=="ENST00000261769")  || ($11=="ENST00000278616")  || ($11=="ENST00000260947")  || ($11=="ENST00000231790")  || ($11=="ENST00000323929")  || ($11=="ENST00000233146")  || ($11=="ENST00000234420")  || ($11=="ENST00000372115")  || ($11=="ENST00000265433")  || ($11=="ENST00000441310")  || ($11=="ENST00000265849")  || ($11=="ENST00000265335")  || ($11=="ENST00000337432")  || ($11=="ENST00000257430")) print $0}' cambios_anotados_tmp2.xls > cancer.xls
	fi
	
	if [ $3 == "Exoma38_hg18" ]; then
		cat header  cambios_anotados_tmp2.xls >  Exoma38_hg18.xls
	fi

	if [ $3 == "Exoma38_hg19" ]; then
		cat header cambios_anotados_tmp2.xls >  Exoma38_hg19.xls
	fi
	if [ $3 == "Exoma50" ]; then
		cat header cambios_anotados_tmp2.xls >  Exoma50.xls
	fi

	rm *_tmp* *biobase *annot *gz
}


# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

EXPECTED_ARGS=3 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


##PIPELINE CANCER $1 BAM OUTPUT FROM BIOSCOPE PAIRED END $2 BAM OUTPUT FROM BIOSCOPE FRAGMENT  $3 Output Bioscope Small indel (*gff) in Paired End $BED Output Bioscope Small Indels in Fragment (*.gff) $SUFFIX Type Analysis (Exoma38 Exoma50 Cardio Cancer)

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo 'Name: sg_GeneralPipe_SNPs_Indels.sh'
	echo 'Description: General pipeline for SNPs & Indels calling'
	echo 'Mandatory parameters:'
	echo '       $1 : .bam output F3'
	echo '       $2 : .gff output Bioscope Small Indels F3'
	echo '       $3 : application type (Exoma38_hg18, Exoma38_hg19, Exoma50, cardio,cancer)'
	exit $E_BADARGS
fi


if [ $3 == "Exoma38_hg19" ]; then
	BED="/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma38_hg19/hg19_unicos_chr"
	SUFFIX="_intervalos.bed"
fi

if [ $3 == "Exoma38_hg18" ]; then
	BED="/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma38_hg18/chr"
	SUFFIX="_intervalos_sondas"
fi

if [ $3 == "Exoma50" ]; then
	BED="/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma50_hg19/exoma50_chr" 
	SUFFIX="_intervalos_sondas"
fi

if [ $3 == "cardio" ]; then
	BED="/data/results/Solid0065/info_resecuenciacion/intervalos_captura/cardio/chr"
	SUFFIX="_intervalos_sondas"
fi

if [ $3 == "cancer" ]; then
	BED="/data/results/Solid0065/info_resecuenciacion/intervalos_captura/cancer/chr"
	SUFFIX="_intervalos_sondas"
fi

echo $1 $2 $3 $BED $SUFFIX
preprocess $1 $2  
mapping_metrics $1 $2 $BED $SUFFIX
SNVs_indels_calling  $1 $2 $BED $SUFFIX
realignment
anotacion $1 $2 $3 $BED

echo $(date)



