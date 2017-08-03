#!/bin/bash
#####################################################
### Bioinformatics Department                     ###
### Sistemas Genomicos June  2013                 ###
### BRCA pipeline for Illumina + Multiplicom      ###
### Version: 0.1								  ###
### Authors: Sheila                               ###
#####################################################


EXPECTED_ARGS=2 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
       echo 'Name: sg_pipe_BRCA_Multiplicom_MiSeq_0.1.sh'
        echo 'Description: General pipeline for BRCA with Multiplicom kit in MiSeq. This pipeline must be executed from the project root directory. fastq.gz files must be copied into a directory called rawdata'
        echo 'Mandatory parameters:'
        echo '	$1 : list of comma-separated samples. Ex: sample1_S1,sample2_S2'
        echo '	$2 : capture kit version, either 2.0 or 2.1'
        exit $E_BADARGS
fi


# Defining input parameters and files

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
	reference="/share/references/genomes/human/hg19/reference/human_hg19.fa"
	index="/share/references/genomes/human/hg19/reference/bowtie2/hg19"
	cabecera="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera-multiplicom"
	cabvcf="/share/references/genomes/human/BRCA1_BRCA2_junior/cabecera_vcf"
	refseqtranscripts="/share/references/genomes/human/BRCA1_BRCA2_junior/BRCA_nomenclature.bed"
	probespergene="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_gene_v2.1.bed"
	probespergenename="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_gene_name_v2.1.bed"
	probespermultiplex="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_multiplex_v2.1.bed"
	probespermultiplexname="/share/references/genomes/human/BRCA1_BRCA2_junior/probes_per_multiplex_name_v2.1.bed"
	genes=(BRCA1 BRCA2)
	#primersfiveprime=(AAGACTCGGCAGCATCTCCA GCGATCGTCACTGTTCTCCA)
	primersfiveprime="/share/references/target_seq/BRCA/Multiplicom_v2.1/primers_5prime.fasta"
	primersthreeprime="/share/references/target_seq/BRCA/Multiplicom_v2.1/primers_3prime.fasta"
#    primersthreeprime=(TGGAGATGCTGCCGAGTCTT TGGAGAACAGTGACGATCGC)
	prefix="BRCA"

elif [[ $2 == "" ]]
then
	echo "Capture kit version needed"
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
if [ ! -d $cwd/Analysis ]; then
	mkdir -p $cwd/Analysis/stats/coverage
	mkdir -p $cwd/Analysis/stats/mapping
	mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
    for s in ${sepbycomas[@]}	
	do
		mkdir -p $cwd/Analysis/stats/primary/$s
		mkdir -p $cwd/Analysis/annotation
		#mkdir -p $cwd/Analysis/variants/$s
		mkdir -p $cwd/Analysis/mapping/bam	
	done
fi  

if [ ! -d $cwd/trash ]; then
    mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
    for s in ${sepbycomas[@]}
    do
		mkdir -p $cwd/trash/rawdata-preprocessing
		mkdir -p $cwd/trash/mapping/$s
		mkdir -p $cwd/trash/variants/$s
		#mkdir -p $cwd/trash/annotation/$s
	done
fi
	

raw_data_preprocessing_and_stats ()
{
    echo "--> Raw data preprocessing and stats calculation start <--"
	date
	mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
    
	for s in ${sepbycomas[@]}
	do
		
		#samplefileR1=`ls $cwd/rawdata/$s\_L001_R1_001.fastq.gz`
		echo "-->  Sample $s <--"	
		echo $cwd/rawdata/$s\_L001_R1_001.fastq.gz
		if [ -s $cwd/rawdata/$s\_L001_R1_001.fastq.gz ]
# -s $cwd/rawdata/$s\_L001_R2_001.fastq.gz] ]
		then
			if [ ! -d $cwd/trash/rawdata-preprocessing/$s ]
			then
				mkdir -p $cwd/trash/rawdata-preprocessing/$s
			fi
		
			for p in ${primersfiveprime[@]}
        	do
            	primerListfiveprime+="-g $p "
        	done

        	for prim in ${primersthreeprime[@]}
        	do
            	primerListthreeprime+="-a $prim "
        	done
        	

		echo $cwd/rawdata/$s\_L001_R1_001.fastq.gz
		java org.usadellab.trimmomatic.TrimmomaticPE -threads 5 -phred33 -trimlog $cwd/trash/rawdata-preprocessing/$s/trimlog_primersfiveprime $cwd/rawdata/$s\_L001_R1_001.fastq.gz $cwd/rawdata/$s\_L001_R2_001.fastq.gz $cwd/trash/rawdata-preprocessing/$s/$s\_1-paired.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_unpaired_R1_5prime.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_2-paired.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_unpaired_R2_5prime.fastq ILLUMINACLIP:$primersfiveprime:2:5:5 MINLEN:30
		
		java org.usadellab.trimmomatic.TrimmomaticPE -threads 5 -phred33 -trimlog $cwd/trash/rawdata-preprocessing/$s/trimlog_primersthreeprime $cwd/trash/rawdata-preprocessing/$s/$s\_1-paired.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_2-paired.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_R1.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_unpaired_R1_3prime.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_R2.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_unpaired_R2_3prime.fastq ILLUMINACLIP:$primersthreeprime:2:5:5 MINLEN:30

		fastqc $cwd/trash/rawdata-preprocessing/$s/$s\_R1.fastq $cwd/trash/rawdata-preprocessing/$s/$s\_R2.fastq  -o $cwd/Analysis/stats/primary/$s/
			line=`zcat $cwd/rawdata/$s\_L001_R1_001.fastq.gz | wc -l | awk '{print $1}'`
			reads=$(($line/2))
			echo "$s" $reads >> $cwd/Analysis/stats/primary/total_reads.txt
		else
			echo "Raw data for sample $s is not correct. Please check whether data for this samples exists in the rawdata dir under the root directory with the name you're giving as input in the command line"
			exit
		fi
	done
	echo "-->  Raw data preprocessing and stats calculation start <--"
	date
}

mapping ()
{
    echo "--> Mapping starts <--"
	date

	mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
    for s in ${sepbycomas[@]}
    do
		bowtie2 -p 5 --very-sensitive -x $index -1 $cwd/trash/rawdata-preprocessing/$s/$s\_R1.fastq -2 $cwd/trash/rawdata-preprocessing/$s/$s\_R2.fastq  -S $cwd/trash/mapping/$s/$s.sam
 		AddOrReplaceReadGroups.jar I=$cwd/trash/mapping/$s/$s.sam O=$cwd/trash/mapping/$s/$s-RG.sam SO=coordinate LB=$s PL=Illumina PU=$s SM=$s ID=$s 
		samtools view -h -b -S $cwd/trash/mapping/$s/$s-RG.sam >  $cwd/trash/mapping/$s/$s-RG.bam
		samtools index $cwd/trash/mapping/$s/$s-RG.bam
	done
	echo "--> Mapping ends <--"
	date
}

variant_calling ()
{
	echo "--> Variant calling starts <--"
    date
	mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
	for s in ${sepbycomas[@]}
    do
		gatk_toolkit -T RealignerTargetCreator -I $cwd/trash/mapping/$s/$s-RG.bam -R $reference -known:dbsnp,vcf /share/references/realign_recalibrate_hapmap/common_all.vcf -o $cwd/trash/mapping/$s/$s.intervals -L $target
		gatk_toolkit -T IndelRealigner -targetIntervals $cwd/trash/mapping/$s/$s.intervals -R $reference -maxReads 1000000 -filterMBQ -I $cwd/trash/mapping/$s/$s-RG.bam -o $cwd/trash/mapping/$s/$s\_realign.bam
		gatk_toolkit -T BaseRecalibrator -rf BadCigar -I $cwd/trash/mapping/$s/$s\_realign.bam -R $reference -knownSites /share/references/realign_recalibrate_hapmap/common_all.vcf --disable_indel_quals -o $cwd/trash/mapping/$s/$s\_realign_recaldata.grp -L $target
		gatk_toolkit -T PrintReads -R $reference -I $cwd/trash/mapping/$s/$s\_realign.bam -o $cwd/Analysis/mapping/$s.bam -BQSR $cwd/trash/mapping/$s/$s\_realign_recaldata.grp
		gatk_toolkit -T UnifiedGenotyper -R $reference -nt 6 -I $cwd/Analysis/mapping/$s.bam -o $cwd/trash/variants/$s/$s\_gatk.vcf -glm BOTH -dcov 10000 -L $target --min_base_quality_score 0 -minIndelCnt 3 -stand_call_conf 0 -stand_emit_conf 0 -minIndelFrac 0
        bgzip $cwd/trash/variants/$s/$s\_gatk.vcf
        tabix -f -p vcf $cwd/trash/variants/$s/$s\_gatk.vcf.gz
        bams+="$cwd/Analysis/mapping/$s.bam "
        vcfs+="$cwd/trash/variants/$s/$s""_gatk.vcf.gz "
		samp+="$s,"
	done
	## Varscan
	samtools mpileup -l $target -d 2000 -A -D -S -m 1 -F 0.001 -f $reference $bams > $cwd/trash/mapping/allsamples.mpileup
	VarScan mpileup2cns $cwd/trash/mapping/allsamples.mpileup --min-coverage 3 --min-freq-for-hom 0.85 --min-var-freq 0.01 --p-value 0.99 --strand-filter 0 --variants --min-reads2 3 > $cwd/trash/variants/varscan.snvs
	sg_parsing_varscan.pl -i $cwd/trash/variants/varscan.snvs -o $cwd/trash/variants/parse-varscan --sample $samp
	
	vcf-merge $vcfs > $cwd/trash/variants/gatk.vcf
	sg_parsing_gatk.pl -v $cwd/trash/variants/gatk.vcf -s $samp -o $cwd/trash/variants/parse-gatk
	sg_collecting_variants.pl -vp $cwd/trash/variants/parse-varscan.snvs -gp $cwd/trash/variants/parse-gatk.snvs -o $cwd/trash/variants/vars
	cd $cwd/trash/variants/
	sg_filtering_variants.pl -not_annot -c $cwd/trash/variants/vars.collected.vcf
}

coverage_stats()
{
    echo "--> Generation of files for coverage statistics starts <--"
	date
    mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`
	for s in ${sepbycomas[@]}
    do
	#This is the way to tell bash to read a file line by line. File name is given after the done statement of the while loop. $probespergene in this case.
	totreads=`grep $s $cwd/Analysis/stats/primary/total_reads.txt | awk '{print $2}'`

	while read line
	do
    	echo -e "$line" > $cwd/trash/mapping/reg_tmp
    	# This line extracts PCR name
		exname=`echo -e "$line" | awk '{print $4}'`
		# This line obtains PCR size
    	size=`echo -e "$line" | awk '{print $3-$2}'`
    	samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/Analysis/mapping/$s.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
   		sum=`sg_suma_columna.pl $cwd/trash/mapping/count_tmp`
    	if [ $sum -eq "0" ]
    	then
        	meandepth=0
    	else
			# Mean depth by PCR
			meandepth=`awk 'BEGIN{printf("%0.5f",('$sum' / '$size'))}'`
			# Mean depth by PCR normalized by total reads
			normByreads=`awk 'BEGIN{printf("%0.5f", (('$sum' / '$size') / '$totreads' ))}'`
    	fi
    	echo $meandepth >> $cwd/trash/mapping/$s/$s\_depth_by_probe_gene.tab
		echo $normByreads >> $cwd/trash/mapping/$s/$s\_depth_by_probe_gene_normalized_by_total_reads.tab 
	done <$probespergene
	
	while read line
    do
        echo -e "$line" > $cwd/trash/mapping/reg_tmp
        exname=`echo -e "$line" | awk '{print $4}'`
        size=`echo -e "$line" | awk '{print $3-$2}'`
        samtools depth -b $cwd/trash/mapping/reg_tmp $cwd/Analysis/mapping/$s.bam | awk '{print $3}' > $cwd/trash/mapping/count_tmp
        sum=`sg_suma_columna.pl $cwd/trash/mapping/count_tmp`
        if [ $sum -eq "0" ]
        then
            meandepth=0
        else
			
			meandepth=`awk 'BEGIN{printf("%0.5f",('$sum' / '$size'))}'`
			normByreads=`awk 'BEGIN{printf("%0.5f", (('$sum' / '$size') / '$totreads' ))}'`
        fi
        echo $meandepth >> $cwd/trash/mapping/$s/$s\_depth_by_probe_multiplex.tab
		echo $normByreads >> $cwd/trash/mapping/$s/$s\_depth_by_probe_multiplex_normalized_by_total_reads.tab
    done <$probespermultiplex
	
		sampleslistgene+=" $cwd/trash/mapping/$s/$s""_depth_by_probe_gene.tab"
		sampleslistmultiplex+=" $cwd/trash/mapping/$s/$s""_depth_by_probe_multiplex.tab"
		sampleslistgenenorm+=" $cwd/trash/mapping/$s/$s""_depth_by_probe_gene_normalized_by_total_reads.tab"
		sampleslistmultiplexnorm+=" $cwd/trash/mapping/$s/$s""_depth_by_probe_multiplex_normalized_by_total_reads.tab"
		allsamples+="	$s"
	done
		
	maximumdepth=`cat $cwd/trash/mapping/*/*depth_by_probe_multiplex.tab | sort -k1n | tail -n 1`
	# This constructs the header for the file that will be the input for the script in R to generate final graphics
	echo "AMPLICON$allsamples"  > $cwd/trash/mapping/samples_tmp
	
	# This will generate the file with the proper header and the values for all samples ordered by gene
	paste $probespergenename $sampleslistgene > $cwd/trash/mapping/gene_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/gene_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab
	paste $probespergenename $sampleslistgenenorm > $cwd/trash/mapping/gene_norm_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/gene_norm_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_total_reads.tab
	
	# This will generate the file with the proper header and the values for all samples ordered by multiplex
	paste $probespermultiplexname $sampleslistmultiplex > $cwd/trash/mapping/multiplex_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/multiplex_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab
	paste $probespermultiplexname $sampleslistmultiplexnorm > $cwd/trash/mapping/multiplex_norm_tmp
	cat $cwd/trash/mapping/samples_tmp $cwd/trash/mapping/multiplex_norm_tmp > $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_total_reads.tab
	
	# Next lines will normalize the two matrix that have previously been created, the one that summarizes depth of coverage per PCR ordered by gene and the one that summarizes mean depth of coverage per PCR ordered by multiplex. 

	# This script will normalize samples per line.
	sg_normalize_samples_by_row_illum.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_total_reads.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_gene -s $mysamples
	sg_normalize_samples_by_row_illum.pl -i $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_total_reads.tab -o $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex -s $mysamples

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
			#The following line will divide matrix values between median and total number of reads per sample
    		awk -v temp=$a '{if ($median==0) {print $valor} else {$temp=$temp/'$median'; print}}' $file\_tmp > $file\_normalized_col_tmp
			awk '{print $1/'$median'}' $file\_col_tmp > $file\_col_tmp_$a
		    mv $file\_normalized_col_tmp $file\_tmp
		done
		sed 's/\t/ /gi' $cwd/trash/mapping/samples_tmp > $cwd/trash/mapping/samples_tmp1
		cat $cwd/trash/mapping/samples_tmp1 $file\_tmp > $file\_normalized_by_row_and_column
	done
	maximumdepthnorm=`cat $cwd/trash/mapping/*_col_tmp_* | sort -k1n | tail -n 1`

	sg_pipe_BRCA_coverage_plot.R $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab $cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_row.tab_normalized_by_row_and_column $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_row.tab_normalized_by_row_and_column $maximumdepth $maximumdepthnorm $cwd/Analysis/stats/coverage
	mv $cwd/trash/mapping/all_samples_per_probe_depth_by_gene.tab $cwd/Analysis/stats/coverage/depth_by_gene.tab
	mv $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex.tab $cwd/Analysis/stats/coverage/depth_by_multiplex.tab
	mv $cwd/trash/mapping/all_samples_per_probe_depth_by_gene_normalized_by_row.tab_normalized_by_row_and_column $cwd/Analysis/stats/coverage/depth_by_gene_normalized.tab
	mv $cwd/trash/mapping/all_samples_per_probe_depth_by_multiplex_normalized_by_row.tab_normalized_by_row_and_column $cwd/Analysis/stats/coverage/depth_by_multiplex_normalized.tab
	echo "--> Generation of files for coverage statistics ends <--"
	date
}

mapping_stats ()
{
    echo "--> Mapping stats calculation start <--"
    date
    mysamples=$1
    sepbycomas=`sed 's/,/ /gi' <<< "$mysamples"`

    echo "Sample_name   Mapped_reads    Mapped_properly_paired" >> $cwd/Analysis/stats/mapping/mapping_stats.tsv
    for s in ${sepbycomas[@]}
    do
        mapped=`samtools flagstat $cwd/trash/mapping/$s/$s-RG.bam | grep " mapped (" | awk '{print $1}'`
        proppaired=`samtools flagstat $cwd/trash/mapping/$s/$s-RG.bam | grep " properly paired " | awk '{print $1}'`
        echo "$s    $mapped $proppaired" >> $cwd/Analysis/stats/mapping/mapping_stats.tsv
        bams+="$cwd/trash/mapping/$s/$s-RG.bam "
    done
    echo $sepbycomas >> $cwd/trash/mapping/correlation.tsv
    samtools depth -b $target $bams | awk '{for (i=3; i <= NF; i++) printf FS$i; print NL }' | sed 's/^ //' >> $cwd/trash/mapping/correlation.tsv
    sg_calculate_correlation_values.R $cwd/trash/mapping/correlation.tsv $cwd/Analysis/stats/mapping
    echo "--> Mapping stats calculation ends <--"
    date
}

annotation ()
{
	echo "--> Variant annotation starts <--"
	date
	source ensembl71
	vep_71_launcher.pl -i $cwd/trash/variants/filtered_not_annotated.vcf -o $cwd/Analysis/annotation/${prefix}_annotated.tsv -t 8 -b 5000
	echo "--> Variant annotation ends <--"
	source genetonic
	date
}

cleaning ()
{
	echo "--> House cleaning <--"
	date
	rm -rf $cwd/trash
	echo "--> House cleaning ends <--"
    date
}

##### MAIN BODY #####

echo "Launching BRCAs Multiplicom pipeline for samples: $1"
echo "Kit_version=$2"

raw_data_preprocessing_and_stats $1
mapping $1
variant_calling $1
coverage_stats $1
mapping_stats $1
annotation $1
cleaning

