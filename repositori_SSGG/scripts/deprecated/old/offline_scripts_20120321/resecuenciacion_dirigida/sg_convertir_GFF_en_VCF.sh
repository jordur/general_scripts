#!/bin/bash

#echo $(date)

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

grep -v "#" $1 | sed 's/chrchr/chr/' | sed 's/ID=.*;del_len=//' | sed 's/ID=.*;ins_len=//' | sed 's/;clustered-.*allele-call=/\t/' | sed 's/;allele-pos=.*alleles=/\t/' | sed 's/;allele-counts=/\t/' | sed 's/;tight_chrom_pos=.*;no_nonred_reads=/\t/' | sed 's/;coverage_ratio=.*;experimental-zygosity=/\t/' | sed 's/;experimental-zygosity-score.*$//' | sed 's/ //gi ' | sed 's/tight//' | sed 's/loose//' > GFF_VCF_tmp1

for j in ${i[@]}
do
	#echo "CHR$j"
	awk '{if ($1=="chr'$j'") print $0}' GFF_VCF_tmp1 >  $j\_GFF_VCF_tmp1
	awk '{if ($3=="insertion_site") print $1"_+_"$4"_"$5"\t"$4"\t"$4"\t+"; else print $1"_+_"$4"_"$5"\t"$4-1"\t"$5"\t+"}' $j\_GFF_VCF_tmp1 > $j\_GFF_VCF_coords_tmp2
	#echo "extracting sequences in CHR$j"
	/usr/local/bin/scripts/sg_extract_seq.pl $j\_GFF_VCF_coords_tmp2 /data/results/solid0065/tomato_Solanum_lycopersicum_2.40/genome/chr$j | grep -v ">" > $j\_GFF_VCF_bases_tmp3
	paste $j\_GFF_VCF_tmp1 $j\_GFF_VCF_bases_tmp3 > $j\_GFF_VCF_tmp4
	#echo "formatting output in CHR$j"
	sg_formatear_fichero_pseudo_GFF.pl $j\_GFF_VCF_tmp4 > $j\_GFF_VCF_tmp5	
done

cat *_GFF_VCF_tmp5 > indels_bioscope.vcf
mv indels_bioscope.vcf $2

#echo "cleaning"
#rm *tmp*
#echo $(date)
