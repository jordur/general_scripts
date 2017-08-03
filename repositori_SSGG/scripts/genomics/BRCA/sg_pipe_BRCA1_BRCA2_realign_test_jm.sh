#!/bin/bash

echo Starting BRCA pipeline on
echo $(date)

# Global parameters definition 
i=(BRCA1 BRCA2)


primary_quality ()
{
	echo "---------Primary Quality------------------"
	sg_get_metrics.sh $qual_input $threads_number $junior_cod
}


mapping ()
{
	source whiterussian
	#echo "---------------Building reference------------"
	echo "--------------Executing Smalt-mapping-----------"
	smalt_x86_64 map -f sam  -o $sample_name.sam  $smalt_index  $fastq
	echo "--------------Formatting input for GATK----------"
	awk '{print $0"\tRG:Z:20110221001322272"}' $sample_name.sam  > $sample_name.sam.tmp

	cat /share/references/genomes/human/BRCA1_BRCA2_junior/cabecera $sample_name.sam.tmp | samtools view -S -b - > $sample_name.bam

	samtools sort $sample_name.bam ${sample_name}_sorted

	samtools index ${sample_name}_sorted.bam

	samtools view -b -q 1 ${sample_name}_sorted.bam > ${sample_name}_sorted_Q1.bam
	
	samtools index ${sample_name}_sorted_Q1.bam
	
	echo "--------------Calculating pileup--------"
	samtools pileup -cf $reference  ${sample_name}_sorted_Q1.bam > $sample_name.pileup
	awk '$3 == "*"' $sample_name.pileup > $sample_name.samtools.indels
	
	samtools view -h ${sample_name}_sorted_Q1.bam | sed 's/BRCA1/chr17/g' | sed 's/BRCA2/chr13/g' | awk '{if ($3 == "chr17") $4 += 41196311; else if ($3 == "chr13") $4 += 32889616; else if ($3 == "LN:81189") $3 = "LN:50000000"; else if ($3 == "LN:84193") $3 = "LN:50000000"; print $0}' | sed 's/ /\t/g'  | samtools view -Shb - > ${sample_name}_sorted_Q1_rebuilt.bam
	samtools index ${sample_name}_sorted_Q1_rebuilt.bam
}

SNVs_indels_calling ()
{
	source whiterussian
	echo "-------CALLING VARIANTS-------"
	gatk_toolkit "-T UnifiedGenotyper -R  $reference  -glm SNP -nt 2  -stand_call_conf 40.0 -stand_emit_conf 20.0 -dcov 200 -I ${sample_name}_sorted_Q1.bam  -o snps_GATK.vcf"
	echo "------------VariantFiltration--------"
	gatk_toolkit "-T VariantFiltration -R $reference  -filter 'MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)' -filter 'QUAL < 30.0' -filter 'SB > 0.1' -filter 'QD < 5.0' -filter 'HRun >= 4'  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName HRUN_FILTER  -cluster 3 -window 10 -B:variant,VCF snps_GATK.vcf -o snps_GATK_filtrados.vcf"

	# Even an older version of GATK is UNFORTUNATELY needed...
	source whiterussian_old

	echo "-----------IndelCalling--------------"
	gatk_toolkit "-T UnifiedGenotyper -R $reference -I ${sample_name}_sorted_Q1.bam  -glm DINDEL -nt 2 -o indels_GATK.vcf"
	echo "-----------IndelFiltration----------"
	gatk_toolkit "-T VariantFiltration -R $reference -filter 'MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)' -filter 'QUAL < 30.0' -filter 'SB > -1.0' -filter 'QD < 2.0'  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -B:variant,VCF indels_GATK.vcf -o indels_GATK_filtrados.vcf"
	# And back to not that old version...
	source whiterussian

	echo "-----------Calling SNVs and Indels with Samtools-------"
	samtools mpileup -g -uf $reference  ${sample_name}_sorted_Q1.bam | bcftools view -bvcg - > var.raw.bcf
	bcftools view var.raw.bcf | vcfutils.pl varFilter -D 1000 > ${1}_SNPs_INDELS_samtools.vcf_raw
	sg_parsing_samtools_indels.pl -p $sample_name.samtools.indels -v ${1}_SNPs_INDELS_samtools.vcf_raw > ${1}_SNPs_INDELS_samtools.vcf
	
	rm $sample_name.samtools.indels
	
	
	echo "-----------Filtering PASS SNPs and Indels---------"
	for x in ${i[@]}; 
	do
		
		awk '{if ($1=="'$x'") print $0}' indels_GATK_filtrados.vcf > ${x}_${1}_indels_GATK_filtrados.vcf
		awk '{if ($1=="'$x'") print $0}' snps_GATK_filtrados.vcf > ${x}_${1}_snps_GATK_filtrados.vcf
		awk '{if ($1=="'$x'") print $0}' ${1}_SNPs_INDELS_samtools.vcf > ${x}_${1}_SNPs_INDELS_samtools.vcf_tmp
		grep "INDEL;" ${x}_${1}_SNPs_INDELS_samtools.vcf_tmp > ${x}_${1}_INDELS_samtools.vcf_tmp
		grep -v "INDEL;" ${x}_${1}_SNPs_INDELS_samtools.vcf_tmp > ${x}_${1}_SNPs_samtools.vcf_tmpx
		grep -v "#" ${x}_${1}_SNPs_samtools.vcf_tmpx > ${x}_${1}_SNPs_samtools.vcf_tmp 
		
		mv ${x}_${1}_INDELS_samtools.vcf_tmp ${x}_${1}_INDELS_samtools.vcf_tmp2
		mv ${x}_${1}_SNPs_samtools.vcf_tmp ${x}_${1}_SNPs_samtools.vcf_tmp2
		mv ${x}_${1}_indels_GATK_filtrados.vcf ${x}_${1}_indels_GATK_filtrados.vcf_tmp2
		mv ${x}_${1}_snps_GATK_filtrados.vcf ${x}_${1}_snps_GATK_filtrados.vcf_tmp2
		grep "PASS" ${x}_${1}_indels_GATK_filtrados.vcf_tmp2 > ${x}_${1}_INDELS_GATK_final.vcf_tmp3
		grep "PASS" ${x}_${1}_snps_GATK_filtrados.vcf_tmp2 >  ${x}_${1}_SNPs_GATK_final.vcf_tmp3

		if [ $x = 'BRCA1' ];
			then
				echo "AWK de $x"
				awk '{sub($2,$2+41196311);sub($1,"chr17");print}' ${x}_${1}_INDELS_GATK_final.vcf_tmp3 > ${x}_${1}_INDELS_GATK_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' ${x}_${1}_SNPs_GATK_final.vcf_tmp3 > ${x}_${1}_SNPs_GATK_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' ${x}_${1}_INDELS_samtools.vcf_tmp2 > ${x}_${1}_INDELS_samtools_final.vcf_tmp4
                                awk '{sub($2,$2+41196311);sub($1,"chr17");print}' ${x}_${1}_SNPs_samtools.vcf_tmp2 > ${x}_${1}_SNPs_samtools_final.vcf_tmp4
				echo "anado cabecera y ordeno con vcf ${x}_$1"
				cat /share/references/BRCA_junior_454/old/cabecera_vcf ${x}_${1}_INDELS_GATK_final.vcf_tmp4 | vcf-sort > ${x}_${1}_INDELS_GATK_tmp5
		                cat /share/references/BRCA_junior_454/old/cabecera_vcf ${x}_${1}_SNPs_GATK_final.vcf_tmp4 | vcf-sort > ${x}_${1}_SNPs_GATK_tmp5
		                cat /share/references/BRCA_junior_454/old/cabecera_vcf ${x}_${1}_INDELS_samtools_final.vcf_tmp4 | vcf-sort > ${x}_${1}_INDELS_samtools_tmp5
                		cat /share/references/BRCA_junior_454/old/cabecera_vcf ${x}_${1}_SNPs_samtools_final.vcf_tmp4 | vcf-sort > ${x}_${1}_SNPs_samtools_tmp5
				echo "BGZIP"
				bgzip ${x}_${1}_INDELS_GATK_tmp5
		                bgzip ${x}_${1}_SNPs_GATK_tmp5
		                bgzip ${x}_${1}_INDELS_samtools_tmp5
                		bgzip ${x}_${1}_SNPs_samtools_tmp5
				echo "TABIX"
		                tabix -f -p vcf ${x}_${1}_INDELS_GATK_tmp5.gz
                		tabix -f -p vcf ${x}_${1}_SNPs_GATK_tmp5.gz
		                tabix -f -p vcf ${x}_${1}_INDELS_samtools_tmp5.gz
                		tabix -f -p vcf ${x}_${1}_SNPs_samtools_tmp5.gz
		                echo "entro IF ${x}_${1}_INDELS_GATK_final.vcf_tmp4, ${x}_${1}_SNPs_GATK_final.vcf_tmp4, ${x}_${1}_INDELS_samtools_final.vcf_tmp4, ${x}_${1}_SNPs_samtools_final.vcf_tmp4"
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


                                vcf-isec -n +1 $file1 $file3 > ${x}_${1}_INDELS_END.vcf;
                                vcf-isec -n +1 $file2 $file4 > ${x}_${1}_SNPs_END.vcf;
		else
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' ${x}_${1}_INDELS_GATK_final.vcf_tmp3 > ${x}_${1}_INDELS_GATK_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' ${x}_${1}_SNPs_GATK_final.vcf_tmp3 > ${x}_${1}_SNPs_GATK_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' ${x}_${1}_INDELS_samtools.vcf_tmp2 > ${x}_${1}_INDELS_samtools_final.vcf_tmp4
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' ${x}_${1}_SNPs_samtools.vcf_tmp2 > ${x}_${1}_SNPs_samtools_final.vcf_tmp4

			cat /share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf ${x}_${1}_INDELS_GATK_final.vcf_tmp4 | vcf-sort > ${x}_${1}_INDELS_GATK_tmp5
	        cat /share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf ${x}_${1}_SNPs_GATK_final.vcf_tmp4 | vcf-sort > ${x}_${1}_SNPs_GATK_tmp5
        	cat /share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf ${x}_${1}_INDELS_samtools_final.vcf_tmp4 | vcf-sort > ${x}_${1}_INDELS_samtools_tmp5
            cat /share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf ${x}_${1}_SNPs_samtools_final.vcf_tmp4 | vcf-sort > ${x}_${1}_SNPs_samtools_tmp5
			bgzip ${x}_${1}_INDELS_GATK_tmp5
	        bgzip ${x}_${1}_SNPs_GATK_tmp5
	      	bgzip ${x}_${1}_INDELS_samtools_tmp5
      		bgzip ${x}_${1}_SNPs_samtools_tmp5
        	tabix -f -p vcf ${x}_${1}_INDELS_GATK_tmp5.gz
	       	tabix -f -p vcf ${x}_${1}_SNPs_GATK_tmp5.gz
	        tabix -f -p vcf ${x}_${1}_INDELS_samtools_tmp5.gz
    		tabix -f -p vcf ${x}_${1}_SNPs_samtools_tmp5.gz

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
			vcf-isec -n +1 $file1 $file3 > ${x}_${1}_INDELS_END.vcf;
			vcf-isec -n +1 $file2 $file4 > ${x}_${1}_SNPs_END.vcf;
		fi	
	
	done

	cat  *_${1}_INDELS_END.vcf *_${1}_SNPs_END.vcf | grep -v "#" > ${1}_cambios.vcf
	rm *tmp* BRCA1_* BRCA2_*
	rm *_filtrados.vcf*  var.raw.bcf
}

realignment () 
{
	source whiterussian
	echo "--------------Indel relalignment--------"

	sort -k1 -k2 ${1}_cambios.vcf > ${1}_cambios_integrados.vcf_raw

	gatk_toolkit "-T UnifiedGenotyper -R  $reference -nt 2 -I ${1}_sorted_Q1.bam  -o ${1}_allbases.vcf_raw --output_mode EMIT_ALL_SITES -dcov 1000"

	sg_realign_indel_vcf.pl -v ${1}_cambios_integrados.vcf_raw -p ${1}_allbases.vcf_raw -j > ${1}_cambios_integrados.vcf
	
	echo "--------------Pileup per base--------"
	sg_GATK_pileup_junior.pl -i ${1}_allbases.vcf -o ${1}_GATK
	
	if [ -d pileup ]; then 
		mv ${1}_GATK.pileup pileup/
		mv ${1}_allbases.vcf pileup/
  	else  
		mkdir pileup
		mv ${1}_GATK.pileup pileup/
		mv ${1}_allbases.vcf pileup/
	fi
	rm *raw*
}

anotacion ()
{
	echo "--------------Variant annotation--------"
	
	sg_parsing_vcf_indels.pl -v ${1}_cambios_integrados.vcf > ${1}_cambios_integrados_parsed.vcf

	sg_variant_effect_predictor_Ensembl64_test.pl -i ${1}_cambios_integrados.vcf --parsed_file ${1}_cambios_integrados_parsed.vcf -o ${1}_cambios_anotados.xls --host blackrussian --user bioinfo --password A29bcd1234# --port 3306 --bed_file /share/references/genomes/human/BRCA1_BRCA2_junior/BRCA_nomenclature.bed
	
	 if [ -d tmp ]; then 
		 mv *bam* tmp/
	  	 mv *.sam tmp/
		 mv ${1}_cambios_integrados.vcf tmp/
		 mv ensembl64.tmp tmp/${1}_ensembl64.tmp
		 mv ensembl64_tmp.log tmp/${1}_ensembl64_tmp.log
 	 else  
		 mkdir tmp
		 mv *bam* tmp/
		 mv *.sam tmp/
		 mv ${1}_cambios_integrados.vcf tmp/
	     mv ensembl64.tmp tmp/${1}_ensembl64.tmp
		 mv ensembl64_tmp.log tmp/${1}_ensembl64_tmp.log
 
	 fi

	rm *.vcf *_tmp* *anno* *.comb *.tbi *.idx *.gz *.bcf *.biobase *.tmp* *.txt

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

EXPECTED_ARGS=7 # Number of arguments expected in the script call
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
	echo '$5 : MIDX.qual (primary Stats)'
	echo '$6 : Threads number'
	echo '$7 : Junior, yes/no'
	echo ''
	echo ''
	exit $E_BADARGS
fi


# variables definitions
sample_name=$1
reference=$2
fastq=$3
smalt_index=$4
qual_input=$5
threads_number=$6
junior_cod=$7

# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

source whiterussian

#primary_quality $5 $6 $7
mapping $1 $2 $3 $4
SNVs_indels_calling $1 $2
realignment $1
anotacion $1 $2

echo Finishing pipeline on
echo $(date)
