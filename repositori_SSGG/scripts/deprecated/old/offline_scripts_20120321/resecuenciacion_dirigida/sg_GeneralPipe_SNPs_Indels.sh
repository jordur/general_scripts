#!/bin/bash

# --------------------------
# ------- indels_SNPs ------
# --------------------------

indels_SNPs ()
# Description: sequence of actions to carry out on the files for obtaining the indels & SNPs metrics
# Parameters:
#       $1 : .vcf file containing the SNVs
#       $2 : .vcf file containing the indels
#	$3 : name prefix of the capture intervals
#	$4 : name sufix of the capture intervals
#	$5 : application type
#	$6 : sample name
{
        # ---------------
        # - Definitions -
        # ---------------
	
	# The chromosomes to consider in the analysis are application type dependant!
	if [ $5 == "BRCA1-BRCA2" ]; then
		chrom=(13 17)
	fi

	if [ $5 == "Exoma38_hg19" ]; then
		chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
	fi

	if [ $5 == "Exoma38_hg18" ]; then
                chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
        fi

	if [ $5 == "Exoma50" ]; then
		chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
	fi

	if [ $5 == "Panel_cardio" ]; then
		chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
	fi

	if [ $5 == "Panel_cancer" ]; then
		chrom=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
	fi


        # ------------------------
        # - body of the pipeline -
        # ------------------------

	awk '{if ($7=="PASS") print $0}' $1 > OK.SNVs
	awk '{if ($7=="PASS") print $0}' $2 > OK.indels
	sed 's/chr//' OK.SNVs > o
	sed 's/chr//' OK.indels > oo
	for c in ${chrom[@]}; do
        	awk '{if ($1=="'$c'") print $0}' o | sort -k2n > chr$c.SNV.vcf
	        awk '{if ($1=="'$c'") print $0}' oo | sort -k2n > chr$c.Indel.vcf
	done
	
	rm o oo

	for c in ${chrom[@]}; do
	        sg_extraer_SNPs_fichero_tabulado28032011.pl $3$c$4 chr$c.SNV.vcf  snps > chr$c\_SNVs_N_rango
        	sg_extraer_SNPs_fichero_tabulado28032011.pl $3$c$4 chr$c.Indel.vcf  indels > chr$c\_Indels_N_rango
	done

	cat chr*_N_rango > variantes_N_rango.vcf

	sg_variant_effect_predictor_Ensembl_190411.pl -i variantes_N_rango.vcf  -o variant_effect_output.txt -format vcf -b 1000

	rm chr*.SNV.vcf chr*.Indel.vcf chr*_SNVs_N_rango chr*_Indels_N_rango

	for c in ${chrom[@]}; do
	        awk '{if ($1=="'$c'") print $0}' variantes_N_rango.vcf | sort -k2n > chr$c.vcf
        	sed 's/_/\t/' variant_effect_output.txt > ook
	        sed 's/_/\t/' ook > ok
	        awk '{if ($1=="'$c'") print $0}' ok | sort -k2n > chr$c.annot
	done

	rm ok ook
	for c in ${chrom[@]}; do
	        sg_combinar_vcf_anotacion_APIs.pl chr$c.vcf chr$c.annot > chr$c.comb
	done

	for c in ${chrom[@]}; do
	        sg_extraerBiobaseInformation.pl chr$c.comb /data/results/Solid0065/info_resecuenciacion/databases/HGMD-Biobase/bioBaseSortUniques_chr$c.txt > chr$c.annot.biobase
	done
	
	cat chr*.annot.biobase  > variantes_anotadas_combinadas.tab_tmp
	rm chr*.comb chr*.annot chr*.vcf chr*.annot.biobase

	grep -v -P "Gene_name\t" variantes_anotadas_combinadas.tab_tmp > variantes_anotadas_combinadas.tab #Nombre_muestra.csv

	# The last step depends on the application type:

	if [ $5 == "BRCA1-BRCA2" ]; then
		awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tCoverage\tGenotype_quality\tMapping_quality\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000380152") || ($11=="ENST00000412061") || ($11=="ENST00000309486"))print $0}' variantes_anotadas_combinadas.tab > variants_$5_$6.xls
	fi

	if [ $5 == "Exoma38" ]; then
		cp variantes_anotadas_combinadas.tab variants_$5_$6.xls
	fi

	if [ $5 == "Exoma50" ]; then
		cp variantes_anotadas_combinadas.tab variants_$5_$6.xls
	fi

	if [ $5 == "Panel_cardio" ]; then
	        awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tCoverage\tGenotype_quality\tMapping_quality\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000070846") || ($11=="ENST00000155840") || ($11=="ENST00000200639") || ($11=="ENST00000211998") || ($11=="ENST00000217381") || ($11=="ENST00000223528") || ($11=="ENST00000228841") || ($11=="ENST00000232975") || ($11=="ENST00000238682") || ($11=="ENST00000243457") || ($11=="ENST00000252321") || ($11=="ENST00000261201") || ($11=="ENST00000261349") || ($11=="ENST00000261448") || ($11=="ENST00000261590") || ($11=="ENST00000262186") || ($11=="ENST00000262464") || ($11=="ENST00000262631") || ($11=="ENST00000265968") || ($11=="ENST00000266732") || ($11=="ENST00000269881") || ($11=="ENST00000271348") || ($11=="ENST00000280904") || ($11=="ENST00000281456") || ($11=="ENST00000282541") || ($11=="ENST00000287878") || ($11=="ENST00000290310") || ($11=="ENST00000290378") || ($11=="ENST00000292327") || ($11=="ENST00000295754") || ($11=="ENST00000299328") || ($11=="ENST00000299333") || ($11=="ENST00000306077") || ($11=="ENST00000307128") || ($11=="ENST00000309889") || ($11=="ENST00000310128") || ($11=="ENST00000316623") || ($11=="ENST00000324501") || ($11=="ENST00000324727") || ($11=="ENST00000330010") || ($11=="ENST00000333535") || ($11=="ENST00000334785") || ($11=="ENST00000337385") || ($11=="ENST00000337851") || ($11=="ENST00000342992") || ($11=="ENST00000343849") || ($11=="ENST00000344887") || ($11=="ENST00000348997") || ($11=="ENST00000354410") || ($11=="ENST00000355349") || ($11=="ENST00000356239") || ($11=="ENST00000356287") || ($11=="ENST00000357077") || ($11=="ENST00000357525") || ($11=="ENST00000357998") || ($11=="ENST00000361308") || ($11=="ENST00000361594") || ($11=="ENST00000366574") || ($11=="ENST00000366578") || ($11=="ENST00000366783") || ($11=="ENST00000367317") || ($11=="ENST00000367319") || ($11=="ENST00000367895") || ($11=="ENST00000369519") || ($11=="ENST00000371372") || ($11=="ENST00000372066") || ($11=="ENST00000372980") || ($11=="ENST00000373960") || ($11=="ENST00000375985") || ($11=="ENST00000375994") || ($11=="ENST00000376480") || ($11=="ENST00000377329") || ($11=="ENST00000379802") || ($11=="ENST00000392770") || ($11=="ENST00000393931") || ($11=="ENST00000395869") || ($11=="ENST00000396576") || ($11=="ENST00000399249") || ($11=="ENST00000399655") || ($11=="ENST00000403994") || ($11=="ENST00000405093") || ($11=="ENST00000432085") || ($11=="ENST00000432168") || ($11=="ENST00000433631") || ($11=="ENST00000452339") || ($11=="ENST00000508053") || ($11=="ENST00000525550") || ($11=="ENST00000545968")) print $0}' variantes_anotadas_combinadas.tab > variants_$5_$6.xls
	fi

	if [ $5 == "Panel_cancer" ]; then
		#Ensembl v59
	        awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tCoverage\tGenotype_quality\tMapping_quality\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000257430") || ($11=="ENST00000278616") || ($11=="ENST00000380152") || ($11=="ENST00000234420") || ($11=="ENST00000231790") || ($11=="ENST00000233146") || ($11=="ENST00000371953") || ($11=="ENST00000326873") || ($11=="ENST00000260947") || ($11=="ENST00000441310") || ($11=="ENST00000265849") || ($11=="ENST00000269305") || ($11=="ENST00000372115") || ($11=="ENST00000265433") || ($11=="ENST00000261769") || ($11=="ENST00000323929") || ($11=="ENST00000265335") || ($11=="ENST00000328354") || ($11=="ENST00000309486") || ($11=="ENST00000261584") || ($11=="ENST00000259008") || ($11=="ENST00000337432") || ($11=="ENST00000357654") || ($11=="ENST00000380152")  || ($11=="ENST00000404276")  || ($11=="ENST00000261584")  || ($11=="ENST00000259008")  || ($11=="ENST00000445888")  || ($11=="ENST00000371953")  || ($11=="ENST00000326873")  || ($11=="ENST00000261769")  || ($11=="ENST00000278616")  || ($11=="ENST00000260947")  || ($11=="ENST00000231790")  || ($11=="ENST00000323929")  || ($11=="ENST00000233146")  || ($11=="ENST00000234420")  || ($11=="ENST00000372115")  || ($11=="ENST00000265433")  || ($11=="ENST00000441310")  || ($11=="ENST00000265849")  || ($11=="ENST00000265335")  || ($11=="ENST00000337432")  || ($11=="ENST00000257430")) print $0}' variantes_anotadas_combinadas.tab > variants_$5_$6.xls
	fi

	rm variantes_anotadas_combinadas.tab
}



# ---------------------------------------------------------------
# ---------------- sg_GeneralPipe_SNPs_Indels.sh ----------------
# ---------------------------------------------------------------

# Name: sg_GeneralPipe_SNPs_Indels.sh
# Description: General pipeline for SNPs & Indels calling
# Mandatory parameters:
#       $1 : .vcf file containing the SNVs
#       $2 : .vcf file containing the indels
#	$3 : application type (BRCA1-BRCA2, Exoma38, Exoma50, Panel_cardio, Panel_cancer)
#	$4 : sample name

# ---------------
# - Definitions -
# ---------------

EXPECTED_ARGS=4	# Number of arguments expected in the script call
E_BADARGS=65	# Error code in case of bad arguments


# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo 'Name: sg_GeneralPipe_SNPs_Indels.sh'
	echo 'Description: General pipeline for SNPs & Indels calling'
	echo 'Mandatory parameters:'
	echo '       $1 : .vcf file containing the SNVs'
	echo '       $2 : .vcf file containing the indels'
	echo '       $3 : application type (BRCA1-BRCA2, Exoma38, Exoma50, Panel_cardio, Panel_cancer)'
	echo '       $4 : sample name'
	exit $E_BADARGS
fi

if [ $3 == "BRCA1-BRCA2" ]; then
        indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/BF18_GSJunior/chr" "_intervalos_sondas" $3 $4
fi

if [ $3 == "Exoma38_hg19" ]; then
	indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma38_hg19/hg19_unicos_chr" "_intervalos.bed" $3 $4
fi

if [ $3 == "Exoma38_hg18" ]; then
        indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma38_hg18/chr" "_intervalos_sondas" $3 $4 
fi

if [ $3 == "Exoma50" ]; then
        indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/exoma38_hg19/hg19_unicos_chr" "_intervalos.bed" $3 $4
fi

if [ $3 == "Panel_cardio" ]; then
        indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/cardio/chr" "_intervalos_sondas" $3 $4
fi

if [ $3 == "Panel_cancer" ]; then
        indels_SNPs $1 $2 "/data/results/Solid0065/info_resecuenciacion/intervalos_captura/cancer/chr" "_intervalos_sondas" $3 $4
fi

