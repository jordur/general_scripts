#!/bin/bash
#####################################################
### Bioinformatics Department                     ###
### Sistemas Genomicos April 2013                 ###
### BRCA pipeline for GSJunior + Multiplicom      ###
### Version: 0.1								  ###
### Authors: Sheila                               ###
#####################################################


EXPECTED_ARGS=3 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Name: sg_pipe_GSJunior_BRCA.sh'
        echo 'Description: General pipeline for GSJunior-BRCA with Multiplicom kit. This pipeline must be executed from the project root directory. The sff file must be copied into this directory.'
        echo 'Mandatory parameters:'
        echo '	$1 : .sff file'
        echo '	$2 : capture kit version, either 2.0 or 2.1'
		echo '	$3 : list of MID numbers (comma-separated). Ex. 1,5,6'
        exit $E_BADARGS
fi


# Defining input parameters and files
# Raw data file 
if [[ $1 == "" ]]; then
    echo "SFF file needed" 
	exit
else
    sff=$1
fi

# Multimplicom BRCA capture Kit version definition
if [[ $2 == "2.0" ]]
then    
    target="/share/references/target_seq/BRCA/Multiplicom_v2.0/target_regions.bed"
	probes="/share/references/target_seq/BRCA/Multiplicom_v2.0/probes.bed"
	reference="/share/references/genomes/human/BRCA1_BRCA2_junior/referencia.fasta"
	index="/share/references/genomes/human/BRCA1_BRCA2_junior/smalt-index"
	cabecera="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera-multiplicom"
	cabvcf="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf"
	refseqtranscripts="/share/references/genomes/human/BRCA1_BRCA2_junior/BRCA_nomenclature.bed"
	genes=(BRCA1 BRCA2)

elif [[ $2 == "2.1" ]]
then
    target="/share/references/target_seq/BRCA/Multiplicom_v2.1/target_regions.bed"
	probes="/share/references/target_seq/BRCA/Multiplicom_v2.1/probes.bed"
	reference="/share/references/genomes/human/BRCA1_BRCA2_junior/referencia.fasta"
	index="/share/references/genomes/human/BRCA1_BRCA2_junior/smalt-index"
	cabecera="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera-multiplicom"
	cabvcf="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf"
	refseqtranscripts="/share/references/genomes/human/BRCA1_BRCA2_junior/BRCA_nomenclature.bed"
	probespergene="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_gene_v2.1.bed"
	probespergenename="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_gene_name_v2.1.bed"
	probespermultiplex="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_multiplex_v2.1.bed"
	probespermultiplexname="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_multiplex_name_v2.1.bed"
	genes=(BRCA1 BRCA2)
	primersfiveprime=(AAGACTCGGCAGCATCTCCA GCGATCGTCACTGTTCTCCA)
	primersthreeprime=(TGGAGATGCTGCCGAGTCTT TGGAGAACAGTGACGATCGC)
elif [[ $2 == "" ]]
then
	echo "Capture kit version needed"
	exit
fi

# Number of the samples to be analyzed 
if [[ $3 == "" ]]; then
	echo "A comma-separated list of sample numbers are needed. Ex: 2,4,5,6 - those numbers must match to MID numbers"
	exit
fi

# Defining working directory
cwd=`pwd`

# Creating directory for input files

if [ ! -d $cwd/rawdata ]
then
        mkdir -p $cwd/rawdata/tmp  
fi

# Creating output folders
if [ ! -d $cwd/analysis ]; then
	mkdir -p $cwd/analysis/stats/coverage/
	midsamples=$3
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}	
	do
		mkdir -p $cwd/analysis/stats/primary/454Reads.MID$s
		mkdir -p $cwd/analysis/annotation/454Reads.MID$s
		mkdir -p $cwd/analysis/variants/454Reads.MID$s
		mkdir -p $cwd/analysis/mapping/bam/454Reads.MID$s	
	done
fi  

if [ ! -d $cwd/trash ]; then
    midsamples=$3
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
		mkdir -p $cwd/trash/mapping/454Reads.MID$s
		mkdir -p $cwd/trash/variants/454Reads.MID$s
		mkdir -p $cwd/trash/annotation/454Reads.MID$s
	done
fi
	

decompressing ()
{
	mv $sff $cwd/rawdata
	cd $cwd/rawdata/tmp
	sfffile -c 200 -r -s -xlr $cwd/rawdata/$sff
	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
	for s in ${sepbycomas[@]}
    do
        sffinfo -s $cwd/rawdata/tmp/454Reads.MID$s.sff > $cwd/rawdata/tmp/454Reads.MID$s.fasta
		sffinfo -q $cwd/rawdata/tmp/454Reads.MID$s.sff > $cwd/rawdata/tmp/454Reads.MID$s.qual
		sg_fastaQual2fastq.pl 454Reads.MID$s.fasta 454Reads.MID$s.qual
    done
}

primary_quality_stats ()
{
    echo "--> Raw data quality statistics calculculation starts <--"
	date
	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
	do
		for p in ${primersfiveprime[@]}
        do            
			primerListfiveprime+="-g $p "
		done

		for prim in ${primersthreeprime[@]}
        do
            primerListthreeprime+="-a $prim "
        done
#		cutadapt -b AAGACTCGGCAGCATCTCCA $cwd/rawdata/tmp/454Reads.MID$s.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp1.fastq && cutadapt -b GCGATCGTCACTGTTCTCCA $cwd/rawdata/tmp/454Reads.MID$s\_tmp1.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp2.fastq && cutadapt -b TGGAGATGCTGCCGAGTCTT $cwd/rawdata/tmp/454Reads.MID$s\_tmp2.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp3.fastq && cutadapt -b TGGAGAACAGTGACGATCGC $cwd/rawdata/tmp/454Reads.MID$s\_tmp3.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp4.fastq
		cutadapt -O 5 $primerListfiveprime $cwd/rawdata/tmp/454Reads.MID$s.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_trim1.fastq
		cutadapt -O 5 $primerListthreeprime $cwd/rawdata/tmp/454Reads.MID$s\_trim1.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_trim.fastq
		prinseq-lite.pl -trim_qual_right 30 -max_len 400 -fastq $cwd/rawdata/tmp/454Reads.MID$s\_trim.fastq -out_good $cwd/rawdata/tmp/454Reads.MID$s\_clean -out_bad $cwd/rawdata/tmp/454Reads.MID$s\_bad
		fastqc $cwd/rawdata/tmp/454Reads.MID$s\_clean.fastq -o $cwd/analysis/stats/primary/454Reads.MID$s/
		echo "454Reads.MID$s"
		line=`wc -l $cwd/rawdata/tmp/454Reads.MID$s\_clean.fastq | awk '{print $1}'`
		reads=$(($line/4))
		echo "454Reads.MID$s" $reads >> $cwd/analysis/stats/primary/total_reads.txt
		cd $cwd
	done
	echo "--> Raw data quality statistics calculculation ends <--"
	date
}

mapping ()
{
	source whiterussian
    echo "--> Mapping starts <--"
	date

	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
        cd $cwd/trash/mapping/454Reads.MID$s
    	smalt_x86_64 map -f sam  -o $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam  $index  $cwd/rawdata/tmp/454Reads.MID$s\_clean.fastq
    
		# Needed format for Samtools
		awk '{if ($1 !~ /#/) print $0"\tRG:Z:20110221001322272"; else print $0}' $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam  > $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam.tmp

		samtools view -S -b $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam.tmp > $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.bam.tmp
		samtools reheader $cabecera $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.bam.tmp > $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.bam.tmp2
		samtools view -b -q 1 $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.bam.tmp2 > $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1.bam

    	samtools sort $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1.bam $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted

	    samtools index $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam

	    samtools pileup -cf $reference  $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam > $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.pileup
	    samtools view -h $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam | sed 's/BRCA1/chr17/g' | sed 's/BRCA2/chr13/g' | awk '{if ($3 == "chr17") $4 += 41196311; else if ($3 == "chr13") $4 += 32889616; else if ($3 == "LN:81189") $3 = "LN:81195210"; else if ($3 == "LN:84193") $3 = "LN:115169878"; print $0}' | sed 's/ /\t/g'  | samtools view -Shb - > $cwd/analysis/mapping/bam/454Reads.MID$s/454Reads.MID$s.bam
    	samtools index $cwd/analysis/mapping/bam/454Reads.MID$s/454Reads.MID$s.bam
	done
	
	echo "--> Mapping ends <--"
	date
}

variant_calling ()
{
	echo "--> Variant calling starts <--"
	date
	
	source whiterussian
	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
        cd $cwd/trash/variants/454Reads.MID$s
   		gatk_toolkit "-T UnifiedGenotyper -R  $reference  -glm SNP -nt 2  -stand_call_conf 40.0 -stand_emit_conf 20.0 -dcov 200 -I $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam -o $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK.vcf"
    	gatk_toolkit "-T VariantFiltration -R $reference  -filter 'MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)' -filter 'QUAL < 30.0' -filter 'SB > 0.1' -filter 'QD < 5.0' -filter 'HRun >= 4'  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER -filterName HRUN_FILTER  -cluster 3 -window 10 -B:variant,VCF $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK.vcf -o $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados.vcf"

	    source whiterussian_old
		cd $cwd/trash/variants/454Reads.MID$s
    	gatk_toolkit "-T UnifiedGenotyper -R $reference -I $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam  -glm DINDEL -nt 4 -o $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK.vcf"
	    gatk_toolkit "-T VariantFiltration -R $reference -filter 'MQ0 >= 4 && ((MQ0 /(1.0*DP))> 0.1)' -filter 'QUAL < 30.0' -filter 'SB > -1.0' -filter 'QD < 2.0'  -filterName HARD -filterName QUAL_FILTER -filterName SB_FILTER -filterName QD_FILTER  -B:variant,VCF $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK.vcf -o $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados.vcf"
    
		source whiterussian
		cd $cwd/trash/variants/454Reads.MID$s
		samtools mpileup -g -uf $reference  $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam | bcftools view -bvcg - > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_var.raw.bcf
    	bcftools view $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_var.raw.bcf | vcfutils.pl varFilter -D 1000 > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_INDELS_samtools.vcf_raw
		grep "INDEL;" $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_INDELS_samtools.vcf_raw > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools.vcf_tmp0
		awk '$3 == "*"' $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.pileup > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s.samtools.indels
		sg_parsing_samtools_indels-modificado.pl -p $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s.samtools.indels -v $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools.vcf_tmp0 > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools.vcf
        grep -v "INDEL;" $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_INDELS_samtools.vcf_raw |  grep -v "#"  > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools.vcf

		for x in ${genes[@]};
    	do
        awk '{if ($1=="'$x'") print $0}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados.vcf | grep "PASS" > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp
        awk '{if ($1=="'$x'") print $0}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados.vcf  | grep "PASS" > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp
        awk '{if ($1=="'$x'") print $0}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools.vcf > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp
		awk '{if ($1=="'$x'") print $0}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools.vcf > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp

		if [ $x = 'BRCA1' ]
		then
			awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp2
            awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp2
            awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp2
            awk '{sub($2,$2+41196311);sub($1,"chr17");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp2
		else
			awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp2 
            awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp2
            awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp2
            awk '{sub($2,$2+32889616);sub($1,"chr13");print}' $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp2
		fi
            cat $cabvcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp2 | vcf-sort > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp3
            cat $cabvcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp2 | vcf-sort > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp3 
            cat $cabvcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp2 | vcf-sort > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp3
			cat $cabvcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp2 | vcf-sort > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp3
                
			bgzip $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp3 &&  bgzip $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp3 && bgzip $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp3 && bgzip $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp3 
                       
			tabix -f -p vcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp3.gz &&  tabix -f -p vcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp3.gz &&  tabix -f -p vcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp3.gz &&  tabix -f -p vcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp3.gz

			if [ -s $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp2 ]
			then
				file1=$cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_indels_GATK_filtrados_$x.vcf_tmp3.gz
            else
				file1=""
			fi
            if [ -s $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp2 ]
			then
            	file2=$cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_GATK_filtrados_$x.vcf_tmp3.gz
			else
				file2=""
			fi
			if [ -s $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp2 ]
			then
                file3=$cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_INDELS_samtools_$x.vcf_tmp3.gz
            else
                file3=""
            fi
            if [ -s $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp2 ]
			then
				file4=$cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_SNVs_samtools_$x.vcf_tmp3.gz
            else
                file4=""
            fi

			vcf-isec -n +1 $file1 $file3 > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_$x\_FINAL_INDELS.vcf
			vcf-isec -n +1 $file2 $file4 > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_$x\_FINAL_SNVs.vcf
		done	
		cat $cwd/trash/variants/454Reads.MID$s/*_FINAL_INDELS.vcf $cwd/trash/variants/454Reads.MID$s/*_FINAL_SNVs.vcf | grep -v "#" | sort -u | vcf-sort >  $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_ALL_VARIANTS.vcf
    done

	echo "--> Variant calling ends <--"
	date
}

realignment ()
{
    # Sometimes works sometimes it does not work...
	source whiterussian
    echo "--> Indel realignment starts <--"
	date
    midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
		gatk_toolkit "-T UnifiedGenotyper -R  $reference -nt 2 -I $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_Q1_sorted.bam  -o $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_allbases.vcf_raw --output_mode EMIT_ALL_SITES -dcov 1000"
    	sg_realign_indel_vcf.pl -v $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_ALL_VARIANTS.vcf -p $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_allbases.vcf_raw -j > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_ALL_VARIANTS_realigned.vcf
	done

	# Para continuar con el análisis por base se tendría que coger el resultado de GATK del pileup y meterlo en el programa sg_pipe_coordenadas_comunes_151111.sh como pileup.  
    #sg_GATK_pileup_junior.pl -i ${1}_allbases.vcf -o ${1}_GATK ---> El fichero 454Reads.MID$s\_allbases.vcf sale del comando anterior

    #if [ -d pileup ]; then
     #   mv ${1}_GATK.pileup pileup/
      #  mv ${1}_allbases.vcf pileup/
    #else
     #   mkdir pileup
      #  mv ${1}_GATK.pileup pileup/
       # mv ${1}_allbases.vcf pileup/
    #fi
    #rm *raw*


	echo "--> Indel realignment ends <--"
    date
}

coverage_stats()
{
    echo "--> Generation of files for coverage statistics starts <--"
	date
    midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
	for s in ${sepbycomas[@]}
    do
	#This is the way to tell bash to read a file line by line. File name is given after the done statement of the while loop. $probespergene in this case.
	while read line
	do
    	echo -e "$line" > $cwd/trash/mapping/reg_tmp
    	# This line extracts PCR name
		exname=`echo -e "$line" | awk '{print $4}'`
		# This line obtains PCR size
    	size=`echo -e "$line" | awk '{print $3-$2}'`
    	samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/analysis/mapping/bam/454Reads.MID$s/454Reads.MID$s.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
   		sum=`sg_suma_columna.pl $cwd/trash/mapping/count_tmp`
    	if [ $sum -eq "0" ]
    	then
        	meandepth=0
    	else
        	meandepth=$(($sum/$size))
    	fi
    	echo $meandepth >> $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_depth_by_probe_gene.tab
	done <$probespergene
	
	while read line
    do
        echo -e "$line" > $cwd/trash/mapping/reg_tmp
        exname=`echo -e "$line" | awk '{print $4}'`
        size=`echo -e "$line" | awk '{print $3-$2}'`
        samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/analysis/mapping/bam/454Reads.MID$s/454Reads.MID$s.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
        sum=`sg_suma_columna.pl $cwd/trash/mapping/count_tmp`
        if [ $sum -eq "0" ]
        then
            meandepth=0
        else
            meandepth=$(($sum/$size))
        fi
        echo $meandepth >> $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s\_depth_by_probe_multiplex.tab
    done <$probespermultiplex
	
		sampleslistgene+=" $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s""_depth_by_probe_gene.tab"
		sampleslistmultiplex+=" $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s""_depth_by_probe_multiplex.tab"
		allsamples+="MID$s"
	done
	maximumdepth=`cat $cwd/trash/mapping/454Reads.MID*/454Reads.MID*depth_by_probe_multiplex.tab | sort -k1n | tail -n 1`
	# This constructs the header for the file that will be the input for the script in R to generate final graphics
	echo "AMPLICON$allsamples" | sed 's/MID/\tMID/'gi > $cwd/trash/mapping/samples_tmp
	
	# This will generate the file with the proper header and the values for all samples ordered by gene
	paste $probespergenename $sampleslistgene > $cwd/trash/mapping/gene_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/gene_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab
	
	# This will generate the file with the proper header and the values for all samples ordered by multiplex
	paste $probespermultiplexname $sampleslistmultiplex > $cwd/trash/mapping/multiplex_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/multiplex_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab
	
	# Next lines will normalize the two matrix that have previously been created, the one that summarizes depth of coverage per PCR ordered by gene and the one that summarizes mean depth of coverage per PCR ordered by multiplex. 

	# This script will normalize samples per line.
	sg_normalize_samples_by_row_454.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_gene -mid $midsamples
	sg_normalize_samples_by_row_454.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex -mid $midsamples

	# The following lines will normalize the sample by columns.
	matrixfiles=($cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_row.tab $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_row.tab)
	for file in ${matrixfiles[@]}
	do
		grep -v "AMPLICON" $file > $file\_tmp
		num_col=$(head -n 1 $file | awk '{ print NF}')
		for ((a=2; a<=$num_col; a++))
		do
    		valor="0"
		    awk -v temp=$a '{print $temp}' $file\_tmp > $file\_col_tmp
    		median=$(sg_calcular_mediana_columna.pl $file\_col_tmp)
    		awk -v temp=$a '{if ($median==0) {print $valor} else {$temp=$temp/'$median'; print}}' $file\_tmp > $file\_normalized_col_tmp
			awk '{print $1/'$median'}' $file\_col_tmp > $file\_col_tmp_$a
		    mv $file\_normalized_col_tmp $file\_tmp
		done
		sed 's/\t/ /gi' $cwd/trash/mapping/samples_tmp > $cwd/trash/mapping/samples_tmp1
		cat $cwd/trash/mapping/samples_tmp1 $file\_tmp > $file\_normalized_by_row_and_column
	done
	maximumdepthnorm=`cat $cwd/trash/mapping/*_col_tmp_* | sort -k1n | tail -n 1`

	sg_pipe_BRCA_coverage_plot.R $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab $cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_row.tab_normalized_by_row_and_column $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_row.tab_normalized_by_row_and_column $maximumdepth $maximumdepthnorm $cwd/analysis/stats/coverage
	echo "--> Generation of files for coverage statistics ends <--"
	date
}

annotation ()
{
	echo "--> Variant annotation starts <--"
	date
	bash
	cd $cwd/trash/annotation/
	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
		# Two different files are joined in this step, one comes from indel realignment and the other one from the original .vcf file without realignment around indes. This is due to an increase in false negatives when realignment is performed with our 'in house' tools.
		cat $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_ALL_VARIANTS_realigned.vcf $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_ALL_VARIANTS.vcf | sort -u | vcf-sort > $cwd/analysis/variants/454Reads.MID$s/454Reads.MID$s.vcf
		sg_parsing_vcf_indels.pl -v $cwd/analysis/variants/454Reads.MID$s/454Reads.MID$s.vcf > $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_parsed.vcf
		sg_variant_effect_predictor_Ensembl64_test.pl -i $cwd/analysis/variants/454Reads.MID$s/454Reads.MID$s.vcf --parsed_file $cwd/trash/variants/454Reads.MID$s/454Reads.MID$s\_parsed.vcf -o $cwd/analysis/annotation/454Reads.MID$s/454Reads.MID$s\_annotated.csv --host blackrussian --user bioinfo --password A29bcd1234# --port 3306 --bed_file $refseqtranscripts
	done

	echo "--> Variant annotation ends <--"
	date
}

cleaning ()
{
	rm -rf $cwd/trash $cwd/rawdata/tmp
}

##### MAIN BODY #####

echo "Launching GSJUnior-BRCAs Multiplicom pipeline for: SFF_file=$1"
echo "Kit_version=$2"
echo "MID number to be analyzed=$3"

decompressing $3
primary_quality_stats $3
mapping $3
variant_calling $3
realignment $3
coverage_stats $3
annotation $3
#cleaning
