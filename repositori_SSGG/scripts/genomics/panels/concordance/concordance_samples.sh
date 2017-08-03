#!/bin/bash

#####################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to perform concordance studies between samples
# @Author: Arbol, based on Sheila's script
# @Date: June 2013
# @Version: v0.1 
#####################################################

# ---------------
# -- functions --
# ---------------


usage ()
{
	scriptname=$1
	echo ""
	echo "$scriptname is a script to perform concordance studies between samples"
	echo ""
	echo "Usage:"
	echo "  $scriptname -v <pseudovcf_file> -f <first_sample> -s <second_sample>"
	echo
	echo "If belonging to an routine analysis, be sure to launch this script from the analysis path (where analysis, logs, jobs and trash folders are created)"
	echo "* Mandatory parameters: *"
	echo " -v	pseudoVCF file containing the samples"
	echo " -f	First sample. Ex: NA12144_replicate1"
	echo " -s	Second sample"
	echo
	exit
}

pseudoVCFToVCF ()
{
	allsamples=$1
	# Sample directory must be given by json template
	vcf=$2
	output=$3
	samplenumber=0
	#To know the number of samples we count the commas in the list of samples provided from the command line
	countcomas=`grep -o "," <<<"$allsamples" | wc -l`
	countcomas=0
	sepbycomas=`sed 's/,/ /gi' <<< "$allsamples"` 
	sg_parse_vcf.pl $vcf > $tmp_dir/${allsamples}_tmp.vcf 
	for i in $sepbycomas; do
		#A check to verify sample presence in vcf file should be made
		findsample=`grep "#" $vcffile | grep $i`
		if [ -z "$findsample" ]; then
			echo "ERROR: Sample $i does not exist in the VCF file"
			usage
			exit
		fi
		samplenumber=$(($samplenumber+1))
        if [ "$samplenumber" == "$(($countcomas+1))" ]; then
            samplelist+="$i"
        else
            samplelist+="$i	"
        fi
    done
	# Once file validation has been completed, the script generates a real VCF file from the pseudoVCF file
	echo "##fileformat=VCFv4.1" > $tmp_dir/$output
	export exec_time=`date`
	echo "##fileDate=$exec_time" >> $tmp_dir/$output
	echo "##source=$vcffile" >> $tmp_dir/$output
	echo "##reference=$refhg19" >> $tmp_dir/$output
	echo "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >> $tmp_dir/$output
	echo "##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Genotype Quality">" >> $tmp_dir/$output
	echo "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">" >> $tmp_dir/$output
	echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$allsamples" >> $tmp_dir/$output
	sg_fromPseudoVCFToVCF.pl $tmp_dir/${allsamples}_tmp.vcf $tmp_dir/$output $allsamples
}

compare_samples_genotypes ()
{
	sample1=$1
	sample2=$2
	bgzip $tmp_dir/$sample1.vcf
	tabix -p vcf -f $tmp_dir/$sample1.vcf.gz
	bgzip $tmp_dir/$sample2.vcf
	tabix -p vcf -f $tmp_dir/$sample2.vcf.gz
	vcf-merge $tmp_dir/$sample2.vcf.gz $tmp_dir/$sample1.vcf.gz > $tmp_dir/sample2_vs_sample1.vcf
	echo "Ignore error \"Broken VCF: empty columns\", results look fine."
	sg_compare_samples_Vs_HapMap.pl -i $tmp_dir/sample2_vs_sample1.vcf -o $tmp_dir/
	sepbycomas=`sed 's/,/ /gi' <<< "$sample2"`
    echo "SAMPLE_NAME	All $sample2 variants	All reference-$sample1 variants	TP (true positives)	DTP (discordant true positives)	CTP (concordant true positives)	FN (false negatives)	TN (true negatives)	FP (false positives)" > $outputpath/concordance_stats.tsv
	for i in $sepbycomas; do
		# Positions not interrogated by arrays in HapMap data that vary in control samples are initially classified as "false positives" but in order to be more precise the annotation of these variants is launched and from that results the number of potential 'false positive' (which of course includes 'novel' variants) will be determined.
		all2=`grep -v "#" $tmp_dir/sample2_vs_sample1.vcf | awk '{print $10}' | grep -v "0\/0" | wc -l`
		all1=`grep -v "#" $tmp_dir/sample2_vs_sample1.vcf | awk '{print $11}' | grep -v "0\/0" | wc -l`
		dtp=`wc -l $tmp_dir/DTP_$i.vcf | awk '{print $1}'`
		ctp=`wc -l $tmp_dir/CTP_$i.vcf | awk '{print $1}'`
		tp=$(( $dtp + $ctp ))
		fn=`wc -l $tmp_dir/FN_$i.vcf | awk '{print $1}'`
		tn=`wc -l $tmp_dir/TN_$i.vcf | awk '{print $1}'`
		fp=`wc -l $tmp_dir/FP_$i.vcf | awk '{print $1}'`
		tn=`wc -l $tmp_dir/TN_$i.vcf | awk '{print $1}'`
		echo $i"	"$all2"	"$all1"	"$tp"	"$dtp"	"$ctp"	"$fn"	"$tn"	"$fp>> $outputpath/concordance_stats.tsv
	done
}


# **************************************************************************************************
# **************************************************************************************************

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------

while getopts ":v:f:s:" Option; do
    case $Option in
        v) vcffile=$OPTARG;;
		f) sample1=$OPTARG;;
		s) sample2=$OPTARG;;
        *) echo "Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
    esac
done

# -----------------------------
# -- Name, version and usage --
# -----------------------------
scriptname=`basename $0`
version="0.1"
if [ $# -eq 0 ]; then
    usage $scriptname $version
    exit $OPTERROR
elif [ -z "$sample1" ]; then
    echo "Sample1 is needed"
elif [ -z "$vcffile" ]; then
    echo "An VCF file is needed"
elif [ -z "$sample2" ]; then
    echo "Sample2 is needed"
fi

# -----------------------
# -- global parameters --
# -----------------------
OPTERROR=65


# ------------------------
# -- default parameters --
# ------------------------

cwd=`pwd`
outputpath=$cwd/analysis/exploit/concordance
tmp_dir=$cwd/trash/exploit/concordance
# Creating output and temporal folders
if [ ! -d $tmp_dir ]; then
    mkdir -p $tmp_dir
fi 

if [ ! -d $outputpath ]; then
    mkdir -p $outputpath
fi

# ---------------
# -- main body -- 
# ---------------

export exec_time=`date`
echo $scriptname Script started at $exec_time
echo "Command launched: $scriptname -f $sample1 -s $sample2 -p $panel -v $vcffile"

pseudoVCFToVCF $sample1 $vcffile $sample1.vcf
pseudoVCFToVCF $sample2 $vcffile $sample2.vcf
compare_samples_genotypes $sample1 $sample2

export exec_time=`date`
echo Script finished at $exec_time