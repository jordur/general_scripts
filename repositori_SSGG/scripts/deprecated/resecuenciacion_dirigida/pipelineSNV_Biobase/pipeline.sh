#!/bin/bash

#REQUIREMENT:
	# sg_coordenadas_no_cubiertas_250411.py GUILLE SCRIPT
	# bed2pileup.pl BED <-> PILEUP
	# sg_noRepeat_intervals.pl TRANSFORM INTERVALS TO NO_REPEAT/LOW COMLPEXITY INTERVALS
	# sg_extraer_SNPs_fichero_Biobase.pl #OBTAIN SNV OF BIOBASE IN INTERVALS

##$1 = PATH WHERE ARE cardioCoordenate_chr*.bed
##$2 = PATH WHERE ARE chr*_intervalos_sondas
##$3 = PATH WHERE ARE MASK FASTA GENOME

mkdir tmp

rm ALL_SNV_Biobase_noCubiertos
touch ALL_SNV_Biobase_noCubiertos

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

for i in ${i[@]}
do
	sort -k2n -k3n $1cardioCoordenate_chr$i.bed > tmp/cardioCoordenate_chr$i.ord.bed #OK DONE!
	sg_unificar_intervalos_de_captura_solapantes.pl tmp/cardioCoordenate_chr$i.ord.bed > tmp/cardioCoordenate_chr$i.ord.noSolap.bed
	
	cat $2chr$i"_intervalos_sondas" tmp/cardioCoordenate_chr$i.ord.noSolap.bed >  tmp/combinate_chr$i
	sort -k2n -k3n tmp/combinate_chr$i > tmp/oeoe
	mv tmp/oeoe tmp/combinate_chr$i
	sg_unificar_intervalos_de_captura_solapantes.pl tmp/combinate_chr$i > tmp/combinate_chr$i.noSolap
	perl bed2pileup.pl $2chr$i"_intervalos_sondas" > tmp/chr$i.pileup
	sg_coordenadas_no_cubiertas_250411.py tmp/chr$i.pileup tmp/combinate_chr$i.noSolap
	sg_miRNA_extract_seq.pl chr$i"_coordenadas_no_cubiertas" $3Homo_sapiens.GRCh37.62.dna_rm.chromosome.$i.fa > tmp/seq_interval_noCubiertos_chr$i.fa
	perl sg_noRepeat_intervals.pl   tmp/seq_interval_noCubiertos_chr$i.fa > tmp/seq_noCubiertos_Nopeat_NolowComplexy_chr$i.fa
	grep ">" tmp/seq_noCubiertos_Nopeat_NolowComplexy_chr$i.fa | sed 's/>/chr'$i'\t/' > tmp/interval_noCubiertos_Nopeat_NolowComplexy.$i.bed
#	perl sg_extraer_SNPs_fichero_Biobase.pl tmp/interval_noCubiertos_Nopeat_NolowComplexy.$i.bed /data/results/Solid0065/info_paneles/HGMD-Biobase/bioBaseSortUniques_chr$i.txt > chr$i.SNV_Biobase_noCubiertos

	perl sg_extraer_SNPs_fichero_Biobase.pl tmp/interval_noCubiertos_Nopeat_NolowComplexy.$i.bed /data/results/Solid0065/info_paneles/HGMD-Biobase/bioBaseSortUniques_chr$i.txt >> ALL_SNV_Biobase_noCubiertos
	#rm chr$i"_coordenadas_no_cubiertas" 
done

