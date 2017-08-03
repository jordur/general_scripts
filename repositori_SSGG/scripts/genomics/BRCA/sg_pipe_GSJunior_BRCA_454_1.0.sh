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
	#reference="/share/references/genomes/human/BRCA1_BRCA2_junior/referencia.fasta"
	reference="/share/references/genomes/human/hg19/reference/human_hg19.fa"
#	index="/share/references/genomes/human/BRCA1_BRCA2_junior/smalt-index"
	index="/share/references/genomes/human/hg19/reference/bowtie2/hg19"
	cabecera="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera-multiplicom"
	cabvcf="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf"
	refseqtranscripts="/share/references/genomes/human/BRCA1_BRCA2_junior/BRCA_nomenclature.bed"
	genes=(BRCA1 BRCA2)
	prefix="BRCA-454"

elif [[ $2 == "2.1" ]]
then
    target="/share/references/target_seq/BRCA/Multiplicom_v2.1/target_regions.bed"
	probes="/share/references/target_seq/BRCA/Multiplicom_v2.1/probes.bed"
	#reference="/share/references/genomes/human/BRCA1_BRCA2_junior/referencia.fasta"
	reference="/share/references/genomes/human/hg19/reference/human_hg19.fa"
	#index="/share/references/genomes/human/BRCA1_BRCA2_junior/smalt-index"
	index="/share/references/genomes/human/hg19/reference/bowtie2/hg19"
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
	prefix="BRCA-454"

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

raw_data_preprocessing_and_stats ()
{
    echo "--> Raw data preprocessing and stats calculation start <--"
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
#       cutadapt -b AAGACTCGGCAGCATCTCCA $cwd/rawdata/tmp/454Reads.MID$s.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp1.fastq && cutadapt -b GCGATCGTCACTGTTCTCCA $cwd/rawdata/tmp/454Reads.MID$s\_tmp1.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp2.fastq && cutadapt -b TGGAGATGCTGCCGAGTCTT $cwd/rawdata/tmp/454Reads.MID$s\_tmp2.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp3.fastq && cutadapt -b TGGAGAACAGTGACGATCGC $cwd/rawdata/tmp/454Reads.MID$s\_tmp3.fastq > $cwd/rawdata/tmp/454Reads.MID$s\_tmp4.fastq
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
	echo "-->  Raw data preprocessing and stats calculation start <--"
	date
}

mapping ()
{
	#source whiterussian
    echo "--> Mapping starts <--"
	date

	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
    for s in ${sepbycomas[@]}
    do
		bowtie2 -p 5 --very-sensitive -x $index -U $cwd/rawdata/tmp/454Reads.MID$s\_clean.fastq -S $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam
 		AddOrReplaceReadGroups.jar I=$cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s.sam O=$cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s-RG.sam SO=coordinate LB=MID$s PL=454 PU=MID$s SM=MID$s   
		samtools view -h -b -S $cwd/trash/mapping/454Reads.MID$s/454Reads.MID$s-RG.sam >  $cwd/analysis/mapping/454Reads.MID$s-RG.bam
		samtools index $cwd/analysis/mapping/454Reads.MID$s-RG.bam
	done
	echo "--> Mapping ends <--"
	date
}

variant_calling ()
{
	echo "--> Variant calling starts <--"
    date
	midsamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$midsamples"`
	for s in ${sepbycomas[@]}
    do
		bams+="$cwd/analysis/mapping/454Reads.MID$s-RG.bam "
		bamsgatk+="-I $cwd/analysis/mapping/454Reads.MID$s-RG.bam "
		mids+="MID$s,"
	done
	samtools mpileup -l $target -B -A -D -S -m 3 -F 0.001 -f $reference $bams > $cwd/trash/mapping/454Reads.mpileup
	VarScan mpileup2cns $cwd/trash/mapping/454Reads.mpileup --min-coverage 3 --min-freq-for-hom 0.85 --min-var-freq 0 --p-value 0.99 --strand-filter 0 --variants --min-reads2 3 > $cwd/trash/variants/varscan.snvs
	sg_parsing_varscan.pl -i $cwd/trash/variants/varscan.snvs -o $cwd/trash/variants/parse-varscan --sample $mids
	gatk_toolkit -T UnifiedGenotyper -R $reference -nt 8 $bamsgatk -o $cwd/trash/variants/gatk.vcf -glm BOTH -dcov 10000 -L $target --min_base_quality_score 0 -minIndelCnt 3 -stand_call_conf 0 -stand_emit_conf 0 -minIndelFrac 0
	sg_parsing_gatk.pl -v $cwd/trash/variants/gatk.vcf -s $mids -o $cwd/trash/variants/parse-gatk
	sg_collecting_variants.pl -vp $cwd/trash/variants/parse-varscan.snvs -gp $cwd/trash/variants/parse-gatk.snvs -o $cwd/trash/variants/vars
	cd $cwd/trash/variants/
	sg_filtering_variants.pl -not_annot -freq_indels 0.1 -gatk 1e+100 -samtools_up 1e+100 -samtools_low 0  -c $cwd/trash/variants/vars.collected.vcf
}

realignment ()
{
    # Sometimes works sometimes it does not work...
	#source whiterussian
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
    	samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/analysis/mapping/454Reads.MID$s-RG.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
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
        samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/analysis/mapping/454Reads.MID$s-RG.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
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
	sg_normalize_samples_by_row.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_gene -mid $midsamples
	sg_normalize_samples_by_row.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex -mid $midsamples

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
	source ensembl71
	#sg_variant_effect_predictor_Ensembl64_parallel.pl -i $cwd/trash/variants/filtered_not_annotated.vcf -o $cwd/trash/variants/all_samples_with_indels_annotated --samples $1 --check_existing --threads 4 --host blackrussian --user bioinfo --password A29bcd1234# --freq 0.05 --port 3306
	#sg_parsing_freqs_pops.pl -i $cwd/trash/variants/all_samples_with_indels_annotated -o $cwd/analysis/annotation/all_samples_annotated.psv $chrs	
	vep_71_launcher.pl -i $cwd/trash/variants/filtered_not_annotated.vcf -o $cwd/analysis/annotation/${prefix}_annotated.tsv -t 8 -b 5000
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
raw_data_preprocessing_and_stats $3
mapping $3
variant_calling $3
#realignment $3
coverage_stats $3
annotation $3
#cleaning
