#!/bin/bash

awk '{if ($7=="PASS") print $0}' $1 > OK.SNVs
awk '{if ($7=="PASS") print $0}' $2 > OK.indels 
b=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
c=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
sed 's/chr//' OK.SNVs > o
sed 's/chr//' OK.indels > oo
for b in ${b[@]}
do
	awk '{if ($1=="'$b'") print $0}' o | sort -k2n > chr$b.SNV.vcf
        awk '{if ($1=="'$b'") print $0}' oo | sort -k2n > chr$b.Indel.vcf
done


for c in ${c[@]}
do

        sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_paneles/cancer/chr$c\_intervalos_sondas chr$c.SNV.vcf  snps > chr$c\_SNVs_N_rango
        sg_extraer_SNPs_fichero_tabulado28032011.pl /data/results/Solid0065/info_paneles/cancer/chr$c\_intervalos_sondas chr$c.Indel.vcf  indels > chr$c\_Indels_N_rango

done

cat chr*_N_rango > variantes_N_rango.vcf

sg_variant_effect_predictor_Ensembl_190411.pl -i variantes_N_rango.vcf  -o variant_effect_output.txt -format vcf -b 1000

rm chr*.SNV.vcf chr*.Indel.vcf chr*_SNVs_N_rango chr*_Indels_N_rango

u=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
a=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)
for i in ${i[@]}  
do
	awk '{if ($1=="'$i'") print $0}' variantes_N_rango.vcf | sort -k2n > chr$i.vcf
	sed 's/_/\t/' variant_effect_output.txt > ook
	sed 's/_/\t/' ook > ok
	awk '{if ($1=="'$i'") print $0}' ok | sort -k2n > chr$i.annot
done

rm ok ook
for a in ${a[@]}
do
	sg_combinar_vcf_anotacion_APIs.pl chr$a.vcf chr$a.annot > chr$a.comb
	
done



for u in ${u[@]}
do
	sg_extraerBiobaseInformation.pl chr$u.comb /data/results/Solid0065/info_paneles/HGMD-Biobase/bioBaseSortUniques_chr$u.txt > chr$u.annot.biobase
done
cat chr*.annot.biobase  > variantes_anotadas_combinadas.tab_tmp
rm chr*.comb chr*.annot chr*.vcf chr*.annot.biobase


grep -v -P "Gene_name\t" variantes_anotadas_combinadas.tab_tmp > variantes_anotadas_combinadas.tab
awk 'BEGIN {print "Gene_name\tChr\tPosition\tReference_genotype\tSample_genotype\tVariation_type\tVariation_length\tZygosity\tVariant_ID\tEnsembl_Gene_ID\tEnsembl_Transcript_ID\tGene_description\tConsequence_on_transcripts\tcDNA_position\taa_position\taa_change\tGERP_conservarion_score\tCoverage\tGenotype_quality\tMapping_quality\tDisease\tHGVS_ID\tPubMed_ID";} {if (($11=="ENST00000357654") || ($11=="ENST00000544455")  || ($11=="ENST00000404276")  || ($11=="ENST00000261584")  || ($11=="ENST00000259008")  || ($11=="ENST00000445888")  || ($11=="ENST00000371953")  || ($11=="ENST00000326873")  || ($11=="ENST00000261769")  || ($11=="ENST00000278616")  || ($11=="ENST00000260947")  || ($11=="ENST00000231790")  || ($11=="ENST00000323929")  || ($11=="ENST00000233146")  || ($11=="ENST00000234420")  || ($11=="ENST00000372115")  || ($11=="ENST00000265433")  || ($11=="ENST00000441310")  || ($11=="ENST00000265849")  || ($11=="ENST00000265335")  || ($11=="ENST00000337432")  || ($11=="ENST00000257430")) print $0}' variantes_anotadas_combinadas.tab > variantes_cancer.xls 

rm variantes_anotadas_combinadas.tab 
