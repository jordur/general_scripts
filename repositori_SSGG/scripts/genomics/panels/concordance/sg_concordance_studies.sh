#!/bin/bash

#####################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script to perform concordance studies in HapMap samples
# @Author: Sheila
# @Contributors: Arbol
# @Date: April 2013
# @Version: v0.1 
#####################################################

# ---------------
# -- functions --
# ---------------


usage ()
{
	scriptname=$1
	echo ""
	echo "$scriptname is a script to perform concordance studies between HapMap samples and data from the HapMap and 1000Genomes projects"
	echo ""
	echo "Usage:"
	echo "  $scriptname -v <pseudovcf_file> -s <sample_list> -d <production panel/exome name> "
	echo
	echo "* Mandatory parameters: *"
	echo " -v	pseudoVCF file with samples"
	echo " -s	A list of comma-separated samples on which to perform the analysis. Ex: NA12144_replicate1,NA12144_replicate2"
	echo " -d	Name of the production panel or exome version. Currently, the script supports: Cardio1,Cardio2,Oto1,Oto2,CardioMio1,OsteoNeuro1,OncoChile and Mama2"
	echo " -h	Name of the HapMap cell lines. Ex: NA12144"
	echo " -o	Output folder (default: cwd/concordance)"
	echo " -t	Temporary folder (default: cwd/trash/concordance)"
	echo
	exit
}

pseudoVCFToVCF ()
{
	allsamples=$1
	# Sample directory must be given by json template
	vcf=$2
	samplenumber=0
	#To know the number of samples we count the comas in the list of samples provided from the command line
	countcomas=`grep -o "," <<<"$allsamples" | wc -l`
	sepbycomas=`sed 's/,/ /gi' <<< "$allsamples"` 
	sg_parse_vcf.pl $vcf > $tmp_dir/pseudoVCF_tmp.vcf 
	for i in $sepbycomas
    do
		#A check to verify sample presence in vcf file should be made
		findsample=`grep "#" $vcffile | grep $i`
		if [ -z "$findsample" ]
		then
			echo "Sample $i does not exist in the VCF file"
			usage
			exit
		fi
		samplenumber=$(($samplenumber+1))
        if [ "$samplenumber" == "$(($countcomas+1))" ]
        then
            samplelist+="$i"
        else
            samplelist+="$i	"
        fi
    done
	# Once file validation has been completed, the script generates a real VCF file from the pseudoVCF file
	echo "##fileformat=VCFv4.1" > $tmp_dir/pseudoVCFToVCF.vcf
	export exec_time=`date`
	echo "##fileDate=$exec_time" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "##source=$vcffile" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "##reference=$refhg19" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Genotype Quality">" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">" >> $tmp_dir/pseudoVCFToVCF.vcf
	echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$samplelist" >> $tmp_dir/pseudoVCFToVCF.vcf
	sg_fromPseudoVCFToVCF.pl $tmp_dir/pseudoVCF_tmp.vcf $tmp_dir/pseudoVCFToVCF.vcf $allsamples
}

get_HapMap_data ()
{
	hapmapsample=$1
	colcounter=0
	
	# Get header (samplecolumn gets he hapmadata header, that is first line)
	samplecolumn=`zcat $hapmapdata | head -1`
	
	# Look for the column containing the hapmap data
	IFS=' ' 
	for a in $samplecolumn; do
		colcounter=$(($colcounter+1))
		if [ "$a" == "$hapmapsample" ]; then
			samplecol=$colcounter
			break
		fi
	done
	zcat $hapmapdata | awk -v initial_line=1 '{if (NR > initial_line) print $3"\t"$1"\texon\t"$4"\t"$4"\t0\t+\t.\t"$'$samplecol'}' > $tmp_dir/$hapmapsample\_gff_tmp

	for i in ${chrs[@]}; do
  	 	awk '{if ($1=="chr'$i'") print $0}' $tmp_dir/$hapmapsample\_gff_tmp > $tmp_dir/$i\_$hapmapsample\_gff_tmp
    	awk '{print $1"_"$4"\t"$4"\t"$4"\t+"}' $tmp_dir/$i\_$hapmapsample\_gff_tmp > $tmp_dir/$i\_$hapmapsample\_gff_tmp2
  		sg_extract_seq.pl $tmp_dir/$i\_$hapmapsample\_gff_tmp2 $refhg18path/chr$i | grep -v ">" | tr "[:lower:]" "[:upper:]" > $tmp_dir/$i\_allele_in_reference_tmp
    	paste $tmp_dir/$i\_$hapmapsample\_gff_tmp $tmp_dir/$i\_allele_in_reference_tmp > $tmp_dir/$i\_$hapmapsample\_gff_tmp3
	done
	cat $tmp_dir/*_$hapmapsample\_gff_tmp3 >  $tmp_dir/$hapmapsample\_gff_all_chrs_tmp
	liftOver -gff $tmp_dir/$hapmapsample\_gff_all_chrs_tmp $liftOverhg18Tohg19 $tmp_dir/$hapmapsample\_gff_all_chrs_in_hg19_tmp $tmp_dir/$hapmapsample\_gff_all_chrs_no_match_hg19_tmp
	echo "## Ignore Warning \"WARNING: -gff is not recommended\". It comes from the liftOver tool - It warns you not to use a gff file to convert genomic coordinates but results are fine and this is the best way to convert HapMap data genomics positions from hg18 to hg19."
	sg_comprobar_genotipo.pl $tmp_dir/$hapmapsample\_gff_all_chrs_in_hg19_tmp > $tmp_dir/$hapmapsample.vcf
	intersectBed -a $tmp_dir/$hapmapsample.vcf -b $target > $tmp_dir/$hapmapsample\_$prefix.vcf_tmp0
	## Some positions changed from hg18 to hg19 in the reference and so ref allele must be changed.
	for j in ${chrs[@]}
	do
    	awk '{if ($1=="chr'$j'") print $1"\t"$2"\t"$2"\t+"}' $tmp_dir/$hapmapsample\_$prefix.vcf_tmp0 > $tmp_dir/$j\_tmp1
    	#awk '{if ($1=="chr'$j'") print $0}' $1  > $j\_tmp1
    	sg_extract_seq.pl $tmp_dir/$j\_tmp1 $refhg19path/chr$j.fa | grep -v ">" > $tmp_dir/$j\_tmp2
    	awk '{if ($1=="chr'$j'") print $0}' $tmp_dir/$hapmapsample\_$prefix.vcf_tmp0  > $tmp_dir/$j\_tmp3
    	paste $tmp_dir/$j\_tmp3 $tmp_dir/$j\_tmp2 | sed 's/\t\t/\t/' > $tmp_dir/$j\_tmp4
    	awk '{print $1"\t"$2"\t"$3"\t"$11"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' $tmp_dir/$j\_tmp4 > $tmp_dir/$j\_tmp5
	done
	cat $tmp_dir/*_tmp5 > $tmp_dir/$hapmapsample\_$prefix.vcf_tmp
	echo "##fileformat=VCFv4.1" > $tmp_dir/$hapmapsample\_$prefix.vcf
    export exec_time=`date`
    echo "##fileDate=$exec_time" >> $tmp_dir/$hapmapsample\_$prefix.vcf
    echo "##source=$tmp_dir/$hapmapsample.vcf" >> $tmp_dir/$hapmapsample\_$prefix.vcf
    echo "##reference=$refhg19" >> $tmp_dir/$hapmapsample\_$prefix.vcf
    echo "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >> $tmp_dir/$hapmapsample\_$prefix.vcf
    echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$hapmapsample" >> $tmp_dir/$hapmapsample\_$prefix.vcf
	vcf-sort $tmp_dir/$hapmapsample\_$prefix.vcf_tmp >> $tmp_dir/$hapmapsample\_$prefix.vcf
#	rm $tmp_dir/*_tmp*  $tmp_dir/$hapmapsample.vcf
}

compare_HapMap_against_samples ()
{
	hapmapsample=$1
	allsamples=$2
	bgzip $tmp_dir/$hapmapsample\_$prefix.vcf
	tabix -p vcf -f $tmp_dir/$hapmapsample\_$prefix.vcf.gz
	bgzip $tmp_dir/pseudoVCFToVCF.vcf
	tabix -p vcf -f $tmp_dir/pseudoVCFToVCF.vcf.gz
	vcf-merge $tmp_dir/pseudoVCFToVCF.vcf.gz $tmp_dir/$hapmapsample\_$prefix.vcf.gz > $tmp_dir/all_samples_Vs_$hapmapsample.vcf
	echo "Ignore error \"Broken VCF: empty columns\", results look fine."
	sg_compare_samples_Vs_HapMap.pl -i $tmp_dir/all_samples_Vs_$hapmapsample.vcf -o $tmp_dir/
	sepbycomas=`sed 's/,/ /gi' <<< "$allsamples"`
    echo "Sample name	DTP	CTP	FN	TN	FP	FP+NV	TP-ANNOT" > $outputpath/concordance_stats.tsv
	for i in $sepbycomas
    do
		# Positions not interrogated by arrays in HapMap data that vary in control samples are innitially classified as "false positives" but in order to be more precise the annotation of these variants is launched and from that results the number of potential 'false positive' (which of course includes 'novel' variants) will be determined.
		variant_effect_predictor_64.pl -i $tmp_dir/FP_$i.vcf -format vcf --check_existing --host blackrussian --user bioinfo --password A29bcd1234# --port 3306 -o $tmp_dir/FP_$i.annot
		dtp=`wc -l $tmp_dir/DTP_$i.vcf | awk '{print $1}'`
		ctp=`wc -l $tmp_dir/CTP_$i.vcf | awk '{print $1}'`
		fn=`wc -l $tmp_dir/FN_$i.vcf | awk '{print $1}'`
		tn=`wc -l $tmp_dir/TN_$i.vcf | awk '{print $1}'`
		fp=`wc -l $tmp_dir/FP_$i.vcf | awk '{print $1}'`
		fpnv=`grep -v "#" $tmp_dir/FP_$i.annot | awk '{if ($13=="-") print $1}' |  sort -u | wc -l `
		annot=`grep -v "#" $tmp_dir/FP_$i.annot | awk '{if ($13!="-") print $1}' |  sort -u | wc -l`
		tn=`wc -l $tmp_dir/TN_$i.vcf | awk '{print $1}'`
		echo $i"	"$dtp"	"$ctp"	"$fn"	"$tn"	"$fp"	"$fpnv"	"$annot >> $outputpath/concordance_stats.tsv
	done
		echo "DTP=True positives but discordant in genotype against control data (HapMap)" >> $outputpath/concordance_stats.tsv
		echo "CTP=True positives and concordant in genotype against control data (HapMap)" >> $outputpath/concordance_stats.tsv
		echo "FN=False negative variants" >> $outputpath/concordance_stats.tsv
		echo "TN=True negative variants" >> $outputpath/concordance_stats.tsv
		echo "FP="False positive variants", i.e., variants that were not genotyped in HapMap sample so this category includes known and true variants + novel variants + false positive variants " >> $outputpath/concordance_stats.tsv
		echo "FP+NV=False positive + novel variants. From FP classified variants those that after annotating against ENSEMBL do not get an ID" >> $outputpath/concordance_stats.tsv
		echo "TP-ANNOT=True positive variants. From FP classified variants those with a known ID against ENSEMBL" >> $outputpath/concordance_stats.tsv
}


# **************************************************************************************************
# **************************************************************************************************

# ------------------------
# -- default parameters --
# ------------------------

cwd=`pwd`
outputpath=$cwd/concordance
tmp_dir=$cwd/trash/concordance

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------

while getopts ":v:s:h:d:o:t:" Option; do
    case $Option in
        v) vcffile=$OPTARG;;
		s) samples=$OPTARG;;
		h) hapmap=$OPTARG;;
        d) panel=$OPTARG;;
        o) outputpath=$OPTARG;;
        t) tmp_dir=$OPTARG;;
        *) echo "ERROR: Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
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
elif [ -z "$samples" ]
then
    echo "A sample list is needed"
elif [ -z "$panel" ]
then
    echo "A panel name or exome version is needed"
	echo "Currently, this script supports Cardio1,Cardio2,Oto1,Oto2,CardioMio1,OsteoNeuro1 and Mama2"
elif [ -z "$vcffile" ]
then
    echo "An VCF file is needed"
elif [ -z "$hapmap" ]
then
    echo "A valid HapMap name is needed. Name should be something like ex: NA12144"
fi

# -----------------------
# -- global parameters --
# -----------------------
OPTERROR=65
# Parameters
refhg18="/share/references/genomes/human/hg18/hg18.fa"
refhg19="/share/references/genomes/human/hg19/reference/human_hg19.fa"
refhg18path="/share/references/genomes/human/hg18/"
refhg19path="/share/references/genomes/human/hg19/reference/"
liftOverhg18Tohg19="/share/references/genomes/liftOver/hg18ToHg19.over.chain"
hapmapdata="/share/dbs/HapMap_Genotype_Data/CEU.hmap.gz"
if [[ $panel == "Cardio1" ]]
then
    target="/share/references/target_seq/cardio/PanelCardio_hg19.bed"
    capture="/share/references/target_seq/cardio/PanelCardio_hg19.bed"
	chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Cardio1"
elif [[ $panel == "Cardio2" ]]
then
    target="/share/references/target_seq/Cardio2/target_regions.bed"
    capture="/share/references/target_seq/Cardio2/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Cardio2"
elif [[ $panel == "Onco1" ]]
then
    target="/share/references/target_seq/Onco1/target_regions.bed"
    capture="/share/references/target_seq/Onco1/capture_regions.bed"
	chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Onco1"
elif [[ $panel == "Onco2" ]]
then
    target="/share/references/target_seq/Onco2/target_regions.bed"
    capture="/share/references/target_seq/Onco2/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Onco2"
elif [[ $panel == "Mama2" ]]
then
    target="/share/references/target_seq/BRCA/Multiplicom_v2.1/target_regions.bed"
    capture="/share/references/target_seq/BRCA/Multiplicom_v2.1/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Mama2"
elif [[ $panel == "Oto1" ]]
then
    target="/share/references/target_seq/Oto1/target_regions.bed"
    capture="/share/references/target_seq/Oto1/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Oto1"
elif [[ $panel == "Oto2" ]]
then
    target="/share/references/target_seq/Oto2/target_regions.bed"
    capture="/share/references/target_seq/Oto2/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="Oto2"
elif [[ $panel == "CardioMio1" ]]
then
    target="/share/references/target_seq/CardioMio1/target_regions.bed"
    capture="/share/references/target_seq/CardioMio1/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="CardioMio1"
elif [[ $panel == "OsteoNeuro1" ]]
then
    target="/share/references/target_seq/OsteoNeuro1/target_regions.bed"
    capture="/share/references/target_seq/OsteoNeuro1/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
    prefix="OsteoNeuro1"
elif [[ $panel == "canada" ]]; then
    target="/share/gluster/OnGoing/BF133_panel_Canada/target-intervals/target.bed"
    capture="/share/gluster/OnGoing/BF133_panel_Canada/target-intervals/capture.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
	prefix="CND"
elif [[ $panel == "OncoChile" ]]; then
	target="/share/references/target_seq/OncoChile/target_regions.bed"
    capture="/share/references/target_seq/OncoChile/capture_regions.bed"
    chrs=`awk '{print $1}' $target | sed 's/chr//' | sort -u | awk '{printf("%s ", $0)}'`
	prefix="OncoChile"
else
	echo 
	echo "ERROR - HEY! $panel is not a valid argument for panel name or exome version"
	echo 
	usage
	exit
fi

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
echo "Command launched: $scriptname -s $samples -p $panel -h $hapmap -v $vcffile"

pseudoVCFToVCF $samples $vcffile
get_HapMap_data $hapmap
compare_HapMap_against_samples $hapmap $samples

export exec_time=`date`
echo Script finished at $exec_time