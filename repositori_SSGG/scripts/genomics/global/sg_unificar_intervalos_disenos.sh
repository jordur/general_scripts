#/bin/bash

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


i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

sed 's/\t0*/\t/gi' $1 | sed 's/chr//' > spl_tmp_0

for j in ${i[@]}
do
	
	awk '{if ($1=="'$j'") print $0}' spl_tmp_0 > spl_tmp1_$j
	sg_unificar-intervalos-parte1.pl spl_tmp1_$j | sort -u -k2n  > spl_tmp2_$j
	awk '{print "chr"$1"\t"$2"\t"$8}' spl_tmp2_$j > spl_tmp3_$j
	sg_intervalos_comunes_071111.pl spl_tmp3_$j "covered"   > spl_tmp4_$j
	 
done
cat spl_tmp4_* > $2.bed
rm spl_tmp*

