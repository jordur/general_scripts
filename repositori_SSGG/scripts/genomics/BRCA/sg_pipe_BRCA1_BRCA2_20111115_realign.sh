#!/bin/bash

echo Starting BRCA pipeline on
echo $(date)

# Global parameters definition 
i=(BRCA1 BRCA2)

mapping ()
{
	#echo "---------------Building reference------------"
	#/share/apps/ssaha2_v2.5.3_x86_64/ssaha2Build -save ref_BRCA1_BRCA2.fa ref_BRCA1_BRCA1.fa
	echo "--------------Executing Smalt-mapping-----------"
	/share/apps/smalt-0.5.7/smalt_x86_64 map -f sam  -o $sample_name.sam  $smalt_index  $fastq
	echo "--------------Formatting input for GATK----------"
	cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera $sample_name.sam > $sample_name.tmp1

	awk '{print $0"\tRG:Z:20110221001322272"}' $sample_name.tmp1  > $sample_name.tmp2

	samtools view -S -b $sample_name.tmp2 > $sample_name.bam

	samtools sort $sample_name.bam ${sample_name}_sorted

	samtools index ${sample_name}_sorted.bam

	samtools view -b -q 1 ${sample_name}_sorted.bam > ${sample_name}_sorted_Q1.bam
	
	samtools index ${sample_name}_sorted_Q1.bam
}

mapping_metrics ()
{
	echo "--------------Calculating metrics--------"
	/share/apps/samtools-0.1.16/samtools pileup -cf $reference  ${sample_name}_sorted_Q1.bam > $sample_name.pileup
	awk '$3 == "*"' $sample_name.pileup > $sample_name.samtools.indels
	
	for j in ${i[@]}; do
		awk '{if ($1=="'$j'") print $0}' $sample_name.pileup > $j\_$sample_name.pileup
		sg_extraer_SNPs_fichero_tabulado.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$j\_v2.1_intervalos.tab $j\_$1.pileup > $j\_$1.enrango_tmp
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
	java -Djava.io.tmpdir=/scratch/ -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T UnifiedGenotyper -R  $reference  -glm SNP -nt 4  -stand_call_conf 40.0 -stand_emit_conf 20.0 -dcov 200 -I ${sample_name}_sorted_Q1.bam  -o snps_GATK.vcf
	echo "------------VariantFiltration--------"
	java -Djava.io.tmpdir=/scratch/ -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T VariantFiltration -R $reference  -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 30.0" -filter "SB > 0.1" -filter "QD < 5.0" -filter "HRun >= 4"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName HRUN_FILTER  -cluster 3 -window 10 -B:variant,VCF snps_GATK.vcf -o snps_GATK_filtrados.vcf
	echo "-----------IndelCalling--------------"
	java -Djava.io.tmpdir=/scratch/  -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T UnifiedGenotyper -R $reference -I ${sample_name}_sorted_Q1.bam  -glm DINDEL -nt 4 -o indels_GATK.vcf
	echo "-----------IndelFiltration----------"
	java -Djava.io.tmpdir=/scratch/  -jar /share/apps/GenomeAnalysisTK-1.0.5083/GenomeAnalysisTK.jar -T VariantFiltration -R $reference -filter "MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)" -filter "QUAL < 30.0" -filter "SB > -1.0" -filter "QD < 2.0"  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -B:variant,VCF indels_GATK.vcf -o indels_GATK_filtrados.vcf
	echo "-----------Calling SNVs and Indels with Samtools-------"
	samtools mpileup -g -uf $reference  ${sample_name}_sorted_Q1.bam | bcftools view -bvcg - > var.raw.bcf
	bcftools view var.raw.bcf | vcfutils.pl varFilter -D 1000 > $1_SNPs_INDELS_samtools.vcf_raw
	sg_parsing_samtools_indels.pl -p $sample_name.samtools.indels -v $1_SNPs_INDELS_samtools.vcf_raw > $1_SNPs_INDELS_samtools.vcf
	
	rm $sample_name.samtools.indels
	
	
	echo "-----------Filtering PASS SNPs and Indels---------"
	for x in ${i[@]}; 
	do
		
		awk '{if ($1=="'$x'") print $0}' indels_GATK_filtrados.vcf > $x\_$1\_indels_GATK_filtrados.vcf
		awk '{if ($1=="'$x'") print $0}' snps_GATK_filtrados.vcf > $x\_$1\_snps_GATK_filtrados.vcf
		awk '{if ($1=="'$x'") print $0}' $1_SNPs_INDELS_samtools.vcf > $x\_$1\_SNPs_INDELS_samtools.vcf_tmp
		grep "INDEL;" $x\_$1\_SNPs_INDELS_samtools.vcf_tmp > $x\_$1\_INDELS_samtools.vcf_tmp
		grep -v "INDEL;" $x\_$1\_SNPs_INDELS_samtools.vcf_tmp > $x\_$1\_SNPs_samtools.vcf_tmpx
		grep -v "#" $x\_$1\_SNPs_samtools.vcf_tmpx > $x\_$1\_SNPs_samtools.vcf_tmp 
		echo "---INDELS_samtools.vcf_tmp2--"
		#sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab $x\_$1\_INDELS_samtools.vcf_tmp indels > $x\_$1\_INDELS_samtools.vcf_tmp2
		echo "---SNPs_samtools.vcf_tmp2--"
		#sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab $x\_$1\_SNPs_samtools.vcf_tmp snps > $x\_$1\_SNPs_samtools.vcf_tmp2
		echo "----indel_GATK.vcf_tmp2--"
		#sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab $x\_$1\_indels_GATK_filtrados.vcf indels > $x\_$1\_indels_GATK_filtrados.vcf_tmp2
		echo "----SNPs_GATK.vcf_tmp2--"
		#sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_resecuenciacion/intervalos_captura/BRCA/$x\_v2_intervalos.tab $x\_$1\_snps_GATK_filtrados.vcf snps >  $x\_$1\_snps_GATK_filtrados.vcf_tmp2
		mv $x\_$1\_INDELS_samtools.vcf_tmp $x\_$1\_INDELS_samtools.vcf_tmp2
		mv $x\_$1\_SNPs_samtools.vcf_tmp $x\_$1\_SNPs_samtools.vcf_tmp2
		mv $x\_$1\_indels_GATK_filtrados.vcf $x\_$1\_indels_GATK_filtrados.vcf_tmp2
		mv $x\_$1\_snps_GATK_filtrados.vcf $x\_$1\_snps_GATK_filtrados.vcf_tmp2
		grep "PASS" $x\_$1\_indels_GATK_filtrados.vcf_tmp2 > $x\_$1\_INDELS_GATK_final.vcf_tmp3
		grep "PASS" $x\_$1\_snps_GATK_filtrados.vcf_tmp2 >  $x\_$1\_SNPs_GATK_final.vcf_tmp3

		if [ $x = 'BRCA1' ];
			then
				echo "AWK de $x"
				#awk '{sub($2,41277501-$2);sub($1,"chr17");print}' $x\_$1\_INDELS_GATK_final.vcf_tmp3 > $x\_$1\_INDELS_GATK_final.vcf_tmp4
				#awk '{sub($2,41277501-$2);sub($1,"chr17");print}' $x\_$1\_SNPs_GATK_final.vcf_tmp3 > $x\_$1\_SNPs_GATK_final.vcf_tmp4
				#awk '{sub($2,41277501-$2);sub($1,"chr17");print}' $x\_$1\_INDELS_samtools.vcf_tmp2 > $x\_$1\_INDELS_samtools_final.vcf_tmp4
				#awk '{sub($2,41277501-$2);sub($1,"chr17");print}' $x\_$1\_SNPs_samtools.vcf_tmp2 > $x\_$1\_SNPs_samtools_final.vcf_tmp4
				awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $x\_$1\_INDELS_GATK_final.vcf_tmp3 > $x\_$1\_INDELS_GATK_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $x\_$1\_SNPs_GATK_final.vcf_tmp3 > $x\_$1\_SNPs_GATK_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $x\_$1\_INDELS_samtools.vcf_tmp2 > $x\_$1\_INDELS_samtools_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $x\_$1\_SNPs_samtools.vcf_tmp2 > $x\_$1\_SNPs_samtools_final.vcf_tmp4
				#mv  $x\_$1\_INDELS_GATK_final.vcf_tmp3 $x\_$1\_INDELS_GATK_final.vcf_tmp4
				#mv $x\_$1\_SNPs_GATK_final.vcf_tmp3 $x\_$1\_SNPs_GATK_final.vcf_tmp4
				#mv $x\_$1\_INDELS_samtools.vcf_tmp2 $x\_$1\_INDELS_samtools_final.vcf_tmp4
				#mv $x\_$1\_SNPs_samtools.vcf_tmp2 $x\_$1\_SNPs_samtools_final.vcf_tmp4
				echo "anado cabecera y ordeno con vcf $x\_$1"
				cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_INDELS_GATK_final.vcf_tmp4 | vcf-sort > $x\_$1\_INDELS_GATK_tmp5
		                cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_SNPs_GATK_final.vcf_tmp4 | vcf-sort > $x\_$1\_SNPs_GATK_tmp5
		                cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_INDELS_samtools_final.vcf_tmp4 | vcf-sort > $x\_$1\_INDELS_samtools_tmp5
                		cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_SNPs_samtools_final.vcf_tmp4 | vcf-sort > $x\_$1\_SNPs_samtools_tmp5
				echo "BGZIP"
				bgzip $x\_$1\_INDELS_GATK_tmp5
		                bgzip $x\_$1\_SNPs_GATK_tmp5
		                bgzip $x\_$1\_INDELS_samtools_tmp5
                		bgzip $x\_$1\_SNPs_samtools_tmp5
				echo "TABIX"
		                tabix -f -p vcf $x\_$1\_INDELS_GATK_tmp5.gz
                		tabix -f -p vcf $x\_$1\_SNPs_GATK_tmp5.gz
		                tabix -f -p vcf $x\_$1\_INDELS_samtools_tmp5.gz
                		tabix -f -p vcf $x\_$1\_SNPs_samtools_tmp5.gz
		                echo "entro IF $x\_$1\_INDELS_GATK_final.vcf_tmp4, $x\_$1\_SNPs_GATK_final.vcf_tmp4, $x\_$1\_INDELS_samtools_final.vcf_tmp4, $x\_$1\_SNPs_samtools_final.vcf_tmp4"
				if [ -s "${x}_${1}_INDELS_GATK_final.vcf_tmp4" ]; then
                		        file1="${x}_${1}_INDELS_GATK_tmp5.gz"
		                else
                		        file1=""
		                fi
                		if [ -s "${x}_${1}_SNPs_GATK_final.vcf_tmp4" ]; then
		                        file2="${x}_${1}_SNPs_GATK_tmp5.gz"
                		else
		                        file2=""
                		fi
		                if [ -s "${x}_${1}_INDELS_samtools_final.vcf_tmp4" ]; then
                		        file3="${x}_${1}_INDELS_samtools_tmp5.gz"
		                else
                		        file3=""
		                fi
                		if [ -s "${x}_${1}_SNPs_samtools_final.vcf_tmp4" ]; then
		                        file4="${x}_${1}_SNPs_samtools_tmp5.gz"
                		else
		                        file4=""
                		fi


                                /share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file1 $file3 > $x\_$1_INDELS_END.vcf;
                                /share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file2 $file4 > $x\_$1_SNPs_END.vcf;
		else
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $x\_$1\_INDELS_GATK_final.vcf_tmp3 > $x\_$1\_INDELS_GATK_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $x\_$1\_SNPs_GATK_final.vcf_tmp3 > $x\_$1\_SNPs_GATK_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $x\_$1\_INDELS_samtools.vcf_tmp2 > $x\_$1\_INDELS_samtools_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $x\_$1\_SNPs_samtools.vcf_tmp2 > $x\_$1\_SNPs_samtools_final.vcf_tmp4
#			mv  $x\_$1\_INDELS_GATK_final.vcf_tmp3 $x\_$1\_INDELS_GATK_final.vcf_tmp4
 #                       mv $x\_$1\_SNPs_GATK_final.vcf_tmp3 $x\_$1\_SNPs_GATK_final.vcf_tmp4
  #                      mv $x\_$1\_INDELS_samtools.vcf_tmp2 $x\_$1\_INDELS_samtools_final.vcf_tmp4
   #                     mv $x\_$1\_SNPs_samtools.vcf_tmp2 $x\_$1\_SNPs_samtools_final.vcf_tmp4

			cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_INDELS_GATK_final.vcf_tmp4 | vcf-sort > $x\_$1\_INDELS_GATK_tmp5
	                cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_SNPs_GATK_final.vcf_tmp4 | vcf-sort > $x\_$1\_SNPs_GATK_tmp5
        	        cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_INDELS_samtools_final.vcf_tmp4 | vcf-sort > $x\_$1\_INDELS_samtools_tmp5
                	cat /share/apps/scripts/resecuenciacion_dirigida/Junior/cabecera_vcf $x\_$1\_SNPs_samtools_final.vcf_tmp4 | vcf-sort > $x\_$1\_SNPs_samtools_tmp5
			bgzip $x\_$1\_INDELS_GATK_tmp5
		        bgzip $x\_$1\_SNPs_GATK_tmp5
		      	bgzip $x\_$1\_INDELS_samtools_tmp5
              		bgzip $x\_$1\_SNPs_samtools_tmp5
	        	tabix -f -p vcf $x\_$1\_INDELS_GATK_tmp5.gz
		       	tabix -f -p vcf $x\_$1\_SNPs_GATK_tmp5.gz
		        tabix -f -p vcf $x\_$1\_INDELS_samtools_tmp5.gz
        		tabix -f -p vcf $x\_$1\_SNPs_samtools_tmp5.gz
			if [ -s "${x}_${1}_INDELS_GATK_final.vcf_tmp4" ]; then
        	                file1="${x}_${1}_INDELS_GATK_tmp5.gz"
                	else
                        	file1=""
			fi
			if [ -s "${x}_${1}_SNPs_GATK_final.vcf_tmp4" ]; then
        	                file2="${x}_${1}_SNPs_GATK_tmp5.gz"
                	else
                        	file2=""
	                fi
			if [ -s "${x}_${1}_INDELS_samtools_final.vcf_tmp4" ]; then
                	        file3="${x}_${1}_INDELS_samtools_tmp5.gz"
	                else
        	                file3=""
                	fi
			if [ -s "${x}_${1}_SNPs_samtools_final.vcf_tmp4" ]; then
        	                file4="${x}_${1}_SNPs_samtools_tmp5.gz"
                	else
                        	file4=""
	                fi
			/share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file1 $file3 > $x\_$1_INDELS_END.vcf;
			/share/apps/vcftools_0.1.4a/perl/vcf-isec -n +1 $file2 $file4 > $x\_$1_SNPs_END.vcf;
		fi	
	
	done

	cat  *_$1_INDELS_END.vcf *_$1_SNPs_END.vcf | grep -v "#" > $1\_cambios.vcf
	
	rm *tmp* BRCA1_* BRCA2_*
    rm *_filtrados.vcf*  var.raw.bcf

}

realignment () 
{
	echo "--------------Indel relalignment--------"
	sg_integration_vcfFotmat.pl $1\_cambios.vcf > $1\_cambios_integrados.vcf_raw ##GENERAR UN FORMATO UNIFICADO DE VCF

	java -Djava.io.tmpdir=/scratch/ -jar /share/apps/GenomeAnalysisTK-1.0.5777/GenomeAnalysisTK.jar -T UnifiedGenotyper -R  $reference -nt 4 -I $1\_sorted_Q1.bam  -o $1\_allbases.vcf_raw --output_mode EMIT_ALL_SITES -dcov 1000
	
	sg_realign_indel_vcf.pl -v $1\_cambios_integrados.vcf_raw -p $1\_allbases.vcf_raw -j > $1\_cambios_integrados.vcf
	
	echo "--------------Pileup per base--------"
	sg_GATK_pileup_junior.pl -i $1\_allbases.vcf -o $1\_GATK
	
	if [ -d pileup ]; then 
		mv $1\_GATK.pileup pileup/
		mv $1\_allbases.vcf pileup/
  	else  
		mkdir pileup
		mv $1\_GATK.pileup pileup/
		mv $1\_allbases.vcf pileup/
	fi
	rm *raw*
}

anotacion ()
{
	echo "--------------Variant annotation--------"

	/share/apps/scripts/variant_effect_predictor_Ensembl62_111013.pl --hgvs --condel=b -i $1\_cambios_integrados.vcf  -o $1\_variant_effect_output.txt -format vcf -b 1000
	
	awk '{if (($4=="ENST00000380152") || ($4=="ENST00000544455") || ($4=="ENST00000412061") || ($4=="ENST00000309486"))print $0}' $1\_variant_effect_output.txt > $1\_variant_effect_output.txt_tmp1
	sg_sustituir_porespacio_annotacion.pl $1\_variant_effect_output.txt_tmp1 > $1\_variant_effect_output.txt_tmp2
	f=(13 17)
        for c in ${f[@]}; do
	        awk '{if ($1=="chr'$c'") print $0}' $1\_cambios_integrados.vcf  | sort -k2n > $c.vcf
		sed 's/_/\t/gi' $1\_variant_effect_output.txt_tmp2 > $1\_variant_effect_output.txt_tmp
                awk '{if ($1=="'$c'") print $0}' $1\_variant_effect_output.txt_tmp | sort -k2n > $c.annot
		echo "-----------Combining annotations-------------"
		sg_combinar_vcf_anotacion_APIs.pl $c.vcf $c.annot > $c.comb
		echo "-----------Biobase annotation----------------"
		sg_extraerBiobaseInformation.pl $c.comb /data/results/Solid0065/info_resecuenciacion/databases/HGMD-Biobase/bioBaseSortUniques_chr$c.txt > $c.annot.biobase	
        done

	cat *.annot.biobase | grep -v "Gene_name" | sort -u |  awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservation_score\tHGVS_ID\tCondel_Score\tCoverage\tGenotype_quality\tMapping_quality\tDisease\tHGVS_ID\tPubMed_ID";} {print $0}'   > $1\_cambios_anotados.xls
	
	if [ -d tmp ]; then 
		mv *bam* tmp/
       	mv *.sam tmp/
        mv $1\_cambios_integrados.vcf tmp/
        mv $1\_variant_effect_output.txt tmp/
		mv $1\_allbases.vcf tmp/
  	else  
		mkdir tmp
		mv *bam* tmp/
		mv *.sam tmp/
		mv $1\_cambios_integrados.vcf tmp/
		mv $1\_variant_effect_output.txt tmp/
		mv $1\_allbases.vcf tmp/
	fi

	rm *.vcf *_tmp* *anno* *.comb *.tbi *.idx *.gz *.bcf *.biobase *.tmp* *.txt

#	sort -u $1\_cambios_anotados_tmp.xls | sort -k1n -k2n > $1\_cambios_anotados.xls ##ELIMINAR UNA DE LAS DOS CABECRAS QUE APARECEN
}

# ---------------------------------------------------------------
# ---------------- sg_pipe_BRCA1_BRCA2.sh ----------------
# ---------------------------------------------------------------

# Name: sg_pipe_BRCA1_BRCA2.sh
# Description: BRCA1 and BRCA2 pipeline for SNPs & Indels calling
# Mandatory parameters:
#       $1 : MID*

# ---------------
# - Definitions -
# ---------------

EXPECTED_ARGS=4 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo ''
	echo ''
	echo 'Name: sg_pipe_BRCA1_BRCA2.sh'
	echo 'Description: BRCA pipeline for SNPs & Indels calling'
	echo 'Mandatory parameters:'
	echo ''
	echo '$1 : Nombre de la muestra '
	echo '$2 : Reference *.fa'
	echo '$3 : Fastq input'
	echo '$4 : Smalt reference Index (default is /share/apps/script/resecuenciacion/Junior/reference)'
	echo ''
	echo ''
	exit $E_BADARGS
fi

# environment definitions


# variables definitions
sample_name=$1
reference=$2
fastq=$3
smalt_index=$4

# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

mapping $1 $2 $3 $4
mapping_metrics $1
SNVs_indels_calling $1 $2
realignment $1
anotacion $1 $2

echo Finishing pipeline on
echo $(date)
