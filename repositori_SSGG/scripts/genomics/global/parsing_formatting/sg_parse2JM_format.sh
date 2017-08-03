#!/bin/bash

#####################################################
# Copyright 2012 - Sistemas Genómicos S.L.
# @Desc: Script for collecting BioScope indels and transform them into the standard format, so that its information can be added to the standard bioinformatics analysis output files
# @Author: Arbol
# @Contributors: Sheila Zuniga, JM Rosa
#####################################################

# ---------------
# -- functions --
# ---------------

function usage ()
# Displays the script usage
{
	# Parameters:
	scriptname=$1

	echo "$scriptname is a script for collecting BioScope indels and transform them into the standard format, so that its information can be added to the standard bioinformatics analysis output files."
	echo "Usage:"
	echo "  $scriptname <options> input1 input2"
	echo
	echo "* Input (mandatory parameters): *"
	echo " input1     gff file containing BioScope indels"
	echo
	echo "* Options: *"
	echo " -s sample_name   Sample name to include in the standard format results (default is \"sample\")"
	echo " -o output_name   Output file name (default is )"
	echo
	echo "* Examples of running: *"
	echo " \$$scriptname -s BM04832 -o BM04832/tmp/indels/indels_PE_bioscope.vcf BioScope.gff"
	exit
}

# -----------------------
# -- global parameters --
# -----------------------
OPTERROR=65
cwd=`pwd`

# ------------------------
# -- default parameters --
# ------------------------
sample="sample"
output="$cwd/indels_bioscope.vcf"

# -----------------------------
# -- Name, version and usage --
# -----------------------------
scriptname=`basename $0`
if [ $# -eq 0 ]; then
	usage $scriptname
	exit $OPTERROR
fi

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------
while getopts ":o:s:" Option; do
	case $Option in
		o) output=`readlink -f $OPTARG`;;
		s) sample="$OPTARG";;
		*) echo "Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
	esac
done

args=($@)
let index='OPTIND-1'
gff=`readlink -f ${args[index]}`

if [ "$gff" == "" ]; then
	echo "Not enough parameters in call to $scriptname!!"
	echo "Call $scriptname with no parameters for help."
	exit
fi

# ---------------
# -- main body -- 
# ---------------

#echo $(date)

i=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M)

RAND1=$(date +%N)
mkdir temp_$RAND1

cd temp_$RAND1
#grep -v "#" $gff | sed 's/chrchr/chr/' | sed 's/ID=.*;del_len=//' | sed 's/ID=.*;ins_len=//' | sed 's/;clustered-.*allele-call=/\t/' | sed 's/;allele-pos=.*alleles=/\t/' | sed 's/;allele-counts=/\t/' | sed 's/;tight_chrom_pos=.*;no_nonred_reads=/\t/' | sed 's/;coverage_ratio=.*;experimental-zygosity=/\t/' | sed 's/;experimental-zygosity-score.*$//' | sed 's/ //gi ' | sed 's/tight//' | sed 's/loose//' > GFF_VCF_tmp1
grep -v "#" $gff | sed 's/chrchr/chr/' | sed 's/ID=.*;del_len=//' | sed 's/ID=.*;ins_len=//' | sed 's/;clustered-.*allele-call=/\t/' | sed 's/;allele-pos=.*alleles=/\t/' | sed 's/;allele-counts=/\t/' | sed 's/;tight_chrom_pos=.*;no_nonred_reads=/\t/' | sed 's/;coverage_ratio=.*;experimental-zygosity=/\t/' | sed 's/;experimental-zygosity-score=/\t/' |  sed 's/;run_names=.*;strands=/\t/' | sed 's/;tags=.*$//' | sed 's/ /_/gi' > GFF_VCF_tmp1

for j in ${i[@]}
do
	#echo "CHR$j"
	awk '{if ($1=="chr'$j'") print $0}' GFF_VCF_tmp1 >  $j\_GFF_VCF_tmp1
	awk '{if ($3=="insertion_site") print $1"_+_"$4"_"$5"\t"$4"\t"$4"\t+"; else print $1"_+_"$4"_"$5"\t"$4-1"\t"$5"\t+"}' $j\_GFF_VCF_tmp1 > $j\_GFF_VCF_coords_tmp2
	#echo "extracting sequences in CHR$j"
	sg_extract_seq.pl $j\_GFF_VCF_coords_tmp2 /share/references/genomes/human/hg19/reference/chr$j | grep -v ">" > $j\_GFF_VCF_bases_tmp3
	paste $j\_GFF_VCF_tmp1 $j\_GFF_VCF_bases_tmp3 > $j\_GFF_VCF_tmp4
	#echo "formatting output in CHR$j"
	sg_parse2JM_format.pl $j\_GFF_VCF_tmp4 $sample > $j\_GFF_VCF_tmp5

	if [ "$j" == "1" ]; then
		cat $j\_GFF_VCF_tmp5 > "$output"
	else
		grep -v "#Chr" $j\_GFF_VCF_tmp5 >> "$output"
	fi
done

#echo "cleaning"
rm *tmp*
rmdir $cwd/temp_$RAND1
#echo $(date)
