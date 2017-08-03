!#/bin/bash

echo 'Description: Trimming and Assembly of Paired Ends Illumina reads'
echo 'Starting Transcriptome Paired End reads Assembly pipeline on'
echo $(date) 



## Variable definitions ##                                                                                                              
##########################
dir='/media/bec2-jcalvete/Elements/rokyta_adamanteus'
samples=$1

### Directories structure ###
#############################
create_dir_check_raw_files()
{ 	
### Logs directory    
    if [ ! -d $dir/logs ]
	then
	mkdir -p $dir/logs
    fi
    
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do
	
	echo "---> Creating Directories structure and Rawdata Quality Control <---"
	
### Directory results/stats/preprocessing
	if [ ! -d $dir/${s}/results/stats/preprocessing ]
	    then
	    mkdir -p $dir/${s}/results/stats/preprocessing && echo "Directory ${s}/results/stats/preprocessing created." || echo "ERROR: Failed to create ${s}/results/stats/preprocessing directory." 
	else
	    echo "WARNING: ${s}/results/stats/preprocessing directory exists!" 
	fi

### Directory mapping
	if [ ! -d $dir/${s}/mapping ]; then
	    mkdir -p $dir/${s}/mapping && echo "Directory ${s}/mapping created." || echo "ERROR: Failed to create ${s}/mapping directory."
	else
    echo "WARNING: ${s}/mapping directory exists!"
	fi
	
### Directory results/stats/rawdata
	if [ ! -d $dir/${s}/results/stats/rawdata ]
	    then
	    mkdir -p $dir/${s}/results/stats/rawdata > /dev/null 2>&1 && echo "Directory ${s}/results/stats/rawdata created." ||  echo "ERROR: Failed to create ${s}/results/stats/rawdata directory."
	else
	    echo "WARNING: ${s}/results/stats/rawdata directory exists!" 
	fi
	
### Directory results/preprocessing
	if [ ! -d $dir/${s}/results/preprocessing ]
	    then
	mkdir -p $dir/${s}/results/preprocessing > /dev/null 2>&1 && echo "Directory ${s}/results/preprocessing created." ||  echo "ERROR: Failed to create ${s}/results/preprocessing directory." 
	else
	    echo "WARNING: ${s}/results/preprocessing directory exists!" 
	fi

### Directory results/assembly
	if [ ! -d $dir/${s}/results/assembly ]
	    then
	    mkdir -p $dir/${s}/results/assembly > /dev/null 2>&1 && echo "Directory ${s}/results/assembly created." ||  echo "ERROR: Failed to create ${s}/results/assembly directory."
	else
	    echo "WARNING: ${s}/results/assembly directory exists!"
	fi
	
	echo "--> Rawdata Quality Control for Sample $s <---"
	#fastqc $dir/rawdata/${s}/${s}_R1*.fastq.gz $dir/rawdata/${s}/${s}_R2*.fastq.gz -o $dir/$s/results/stats/rawdata
	fastqc $dir/rawdata/${s}/${s}_R1*.fastq $dir/rawdata/${s}/${s}_R2*.fastq -o $dir/$s/results/stats/rawdata
	echo "---> Rawdata Quality Control finished for Sample $s <---"
echo "done"

done
}

### Error correction & Filtering & Trimming ### 
###############################################

filtering_and_trimming()
{
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do
    	
### RCorrector ###
echo "---> Executing Rcorrector <---"
echo $(date)
/home/bec2-jcalvete/Feina_Jordi/programes/Rcorrector-master/run_rcorrector.pl -t 4 -1 $dir/rawdata/${s}/${s}_R1.fastq -2 $dir/rawdata/${s}/${s}_R2.fastq -od $dir/$s/results/preprocessing
echo "---> Rcorrector done <---"
echo $(date)

echo "--> Rawdata Quality Control for Sample $s <---"
fastqc $dir/${s}/results/preprocessing/${s}_R1.cor.fq $dir/${s}/results/preprocessing/${s}_R2.cor.fq -o $dir/$s/results/stats/preprocessing
echo "---> Rawdata Quality Control finished for Sample $s <---"


### Trimmomatic
	echo "---> Executing Trimmomatic in sample $s <---"
	date
	java -jar /home/bec2-jcalvete/Feina_Jordi/programes/Trimmomatic-0.35/trimmomatic-0.35.jar PE -threads 4 -phred33 $dir/$s/results/preprocessing/${s}_R1.cor.fq $dir/$s/results/preprocessing/${s}_R2.cor.fq $dir/$s/results/preprocessing/${s}_R1-paired.fastq $dir/$s/results/preprocessing/${s}_R1-unpaired.fastq $dir/$s/results/preprocessing/${s}_R2-paired.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired.fastq ILLUMINACLIP:/home/bec2-jcalvete/Feina_Jordi/programes/Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:40:15 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:90 HEADCROP:9
	fastqc $dir/$s/results/preprocessing/${s}_R1-paired.fastq $dir/$s/results/preprocessing/${s}_R1-unpaired.fastq $dir/$s/results/preprocessing/${s}_R2-paired.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired.fastq -o $dir/$s/results/stats/preprocessing 
	echo "---> Trimmomatic DONE in sample $s <---"
	date
 
### Prinseq sobre les unpaired ###
	echo "---> Executing Prinseq in sample $s <---"
	date
	perl /home/bec2-jcalvete/Feina_Jordi/programes/prinseq-lite-0.20.4/prinseq-lite.pl -min_len 90 -fastq $dir/$s/results/preprocessing/${s}_R1-unpaired.fastq
	perl /home/bec2-jcalvete/Feina_Jordi/programes/prinseq-lite-0.20.4/prinseq-lite.pl -min_len 90 -fastq $dir/$s/results/preprocessing/${s}_R2-unpaired.fastq
	fastqc $dir/$s/results/preprocessing/${s}_R1-unpaired_prinseq_good*.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired_prinseq_good*.fastq $dir/$s/results/preprocessing/${s}_R1-unpaired_prinseq_bad*.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired_prinseq_bad*.fastq -o $dir/$s/results/stats/preprocessing 
	echo "---> Prinseq DONE in sample $s <---"
	date

### Resincronitzar R1 i R2 ###
	echo "----> Synchronizing R1 and R2 sample $s <---"
	date
	/home/bec2-jcalvete/Feina_Jordi/scripts/scripts_normandeau_data/fastqCombinePairedEnd.py  $dir/$s/results/preprocessing/${s}_R1-paired.fastq $dir/$s/results/preprocessing/${s}_R2-paired.fastq
	echo "---> Synchronization DONE in sample $s <---"
	date
	
### Opcional ###
	echo "---> Executing FASTQC on preprocessing data in sample $s <---"
	date
	fastqc  $dir/$s/results/preprocessing/${s}_R1-paired*pairs_R1.fastq  $dir/$s/results/preprocessing/${s}_R2-paired*pairs_R2.fastq $dir/$s/results/preprocessing/${s}_R1-unpaired*good*.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired*good*.fastq -o $dir/$s/results/stats/preprocessing
	echo "---> FASTQC done on preprocessing in sample $s <---"  
	date  
	
done

}


### Reducing complexity ###                                                                                                          
###########################

reducing_complexity()
{
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do
	
	### FLASH ###
	echo "---> Executing FLASH in sample $s <---" 
	date 
	/home/bec2-jcalvete/Feina_Jordi/programes/FLASH-1.2.11/flash -t 5 $dir/${s}/results/preprocessing/${s}_R1-paired.fastq $dir/$s/results/preprocessing/${s}_R2-paired.fastq  -r 125 -f 200 -s 20 -d $dir/$s/results/preprocessing
	echo "---> FLASH done in sample $s <---" 
	date 
	
	echo "---> Concatenem lectures singles <---"
	date
	cat $dir/$s/results/preprocessing/${s}_R1-unpaired_prinseq_good*.fastq $dir/$s/results/preprocessing/${s}_R2-unpaired_prinseq_good*.fastq $dir/$s/results/preprocessing/out.extendedFrags.fastq > $dir/$s/results/preprocessing/${s}_singles.fastq
	echo "---> Singletons fets <---"
	date
	
	### Colapsar-lecturas
#	echo "---> Executing sg_colapsar-lecturas-paired-end.sh in sample $s <---"  
#	date 
#	sg_colapsar-lecturas-paired-end.sh $dir/$s/results/preprocessing/${s}_1-paired-mask.fastq $dir//$s/results/preprocessing/${s}_2-paired-mask.fastq $dir//$s/results/preprocessing
#	echo "---> sg_colapsar-lecturas-paired-end.sh done in sample $s <---" 
#	date 
	
	### FASTX-COLLAPSER de singletons
#	echo "---> Executing FastX Collapser in sample $s <---" 
#	fastx_collapser -Q33 -i $dir/$s/results/preprocessing/${s}_singles.fastq -o $dir/$s/results/preprocessing/${s}_singles_collapse.fasta
#	echo "---> FastX Collapser done in sample $s <---" 
#	date 
	
done
    echo "---> Preprocessing done <---" 
}


### Assembly Oases ###
######################

assembly_oases()
{
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do	
	date 
	echo "---> Executing OASES Assembly in sample $s <---"  
	date 
	
	oases_pipeline.py -m 23 -M 71 -o $dir/$s/results/assembly/oases -d "-strand_specific -fastq -shortPaired -separate $dir/$s/results/preprocessing/out.notCombined_1.fastq $dir/$s/results/preprocessing/out.notCombined_2.fastq -fastq -short $dir/$s/results/preprocessing/${s}_singles.fastq" -p "-scaffolding yes -cov_cutoff auto "
	
	echo "---> OASES Assembly done in sample $s <---" 
	date 
	
done
}


#### Assembly Trinity ####
##########################

assembly_trinity()
{
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do
	echo "---> Executing Trinity assembler in sample $s <---" 
	date 
	Trinity.pl --seqType fq --JM 50G --left $dir/$s/rawdata/Trinity_R1.fq.gz--right $dir/$s/rawdata/Trinity_R2.fq.gz --SS_lib_type RF --output $dir/$s/results/assembly/Trinity --CPU 5 --min_contig_length 100 --full_cleanup
	echo "---> Trinity assembly done in sample $s <---" 
	date 
done
}

### Oases and Trinity stats ###
###############################

stats_oases_trinity()
{
sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
for s in ${sepbycomas[@]}; do
	echo "---> Executing Oases and Trinity stats in sample $s <---"
	date
	
	out=$dir/${s}/results/stats/assembly
	echo -e "K-mer,Total_contigs,Contigs>200nt,Contigs>500nt,Contigs>1000nt,Contigs>2000nt,Contigs>5000n,Contigs>8000nt,Max_length,Total_nt,Avrge_length,N50" > $out/header	 
	
	for i in $dir/${s}/results/assembly/Trinity/Trinity.Trinity.fasta;do
	    resultado=`grep -c ">" $i`
	    /home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/length_fasta.pl $i | awk '{print $2}' | sort -nr > $out/t
	    long=`awk '{if ($1>200) print}' $out/t | wc -l`
	    long2=`awk '{if ($1>500) print}' $out/t | wc -l`
	    long3=`awk '{if ($1>1000) print}' $out/t | wc -l`
	    long4=`awk '{if ($1>2000) print}' $out/t | wc -l`
	    long6=`awk '{if ($1>5000) print}' $out/t | wc -l`
	    long7=`awk '{if ($1>8000) print}' $out/t | wc -l`
	    max=`head -n 1 $out/t | awk '{print}'`
	    res2=`/home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/sg_suma_columna.pl $out/t`
	    avrge=`echo $res2/$resultado | bc`
	    n50=`/home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/quality_stats/sg_calcularN50.pl $i | head -n 1`
	    echo -e "Trinity\t$resultado\t$long\t$long2\t$long3\t$long4\t$long6\t$long7\t$max\t$res2\t$avrge\t$n50" > $out/t
		cat $out/header $out/t > $out/trinity_stats.tsv
		
	    /home/bec2-jcalvete/Feina_Jordi/programes/trinity/util/TrinityStats.pl $dir/${s}/results/assembly/Trinity/Trinity.Trinity.fasta > $out/Trinity_stats.txt
	    rm $out/t
	done

	for i in `seq 25 2 71`;
		do
		cd $dir/${s}/results/assembly/Oases/;
		name=`echo $i | awk -F "/" '{print $NF}'`
#		resultado=`grep -c ">" ${i}/${i}_transcripts.fa`
		resultado=`grep -c ">" oases_${i}/transcripts.fa`
 		/home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/length_fasta.pl oases_${i}/transcripts.fa | awk '{print $2}'| sort -nr > $out/t
 		long=`awk '{if ($1>200) print}' $out/t | wc -l`
  		long2=`awk '{if ($1>500) print}' $out/t | wc -l`
  		long3=`awk '{if ($1>1000) print}' $out/t | wc -l`
  		long4=`awk '{if ($1>2000) print}' $out/t | wc -l`
  		long6=`awk '{if ($1>5000) print}' $out/t | wc -l`
  		long7=`awk '{if ($1>8000) print}' $out/t | wc -l`
  		max=`head -n 1 $out/t | awk '{print}'`
 		res2=`/home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/sg_suma_columna.pl $out/t`
		avrge=`echo $res2/$resultado | bc`
 		n50=`/home/bec2-jcalvete/Feina_Jordi/scripts/repositori_SSGG/scripts/deprecated/quality_stats/sg_calcularN50.pl oases_${i}/transcripts.fa | head -n 1` 		echo -e "$name,$resultado,$long,$long2,$long3,$long4,$long6,$long7,$max,$res2,$avrge,$n50" > $out/stats_${name}
 		rm $out/t 
		done
		 
	cat $out/header $out/stats_* > $out/oases_final_stats.csv
	rm $out/stats_*
	
	echo "---> Oases and Trinity stats done in sample $s <---"
	date
    done
}


### K-mer selection ###
#######################

kmer_selection()
{
   sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
   for s in ${sepbycomas[@]}; do
	echo "---> Executing k-mer selection in sample $s <---" 
	date 
	echo "Please, introduce k-mer selection values space-separated"
	read intervalo
	for j in  ${intervalo[@]};
	  do
	  /share/gluster/tests/jdurban/scripts/sg_seleccionar_lecturas_conf1.pl $dir//$s/results/assembly/oases-51-71_$j/transcripts.fa > $dir//$s/results/assembly/oases-51-71_$j/$j-filtered_transcripts.fa
	  transcritos="$transcritos $dir//$s/results/assembly/oases-51-71_$j/$j-filtered_transcripts.fa"
	done
	cat $transcritos > $dir//$s/results/assembly/oases-transcripts.fa
	
	echo "---> K-mer selection done in sample $s <---"  
	date 
   done

}

### Assembly CAP3 ###
#####################


assembly_cap3()
{
   sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
   for s in ${sepbycomas[@]}; do
	if [ ! -d "$dir//$s/results/assembly/cap3/" ]
	    then
	    mkdir -p $dir//$s/results/assembly/cap3/  > /dev/null 2>&1 && echo "Directory cap3 created." ||  echo "ERROR: Failed to create logs directory." 
       else
	    echo "WARNING: cap3 directory exists!"
       fi
	echo "---> Performing CAP3 Assembly in sample $s <---" 
	date 
	
	cat $dir//$s/results/assembly/oases-transcripts.fa $dir//$s/results/assembly/Trinity/Trinity.fasta > $dir//$s/results/assembly/cap3/all.contigs.fa	
	cap3 $dir//$s/results/assembly/cap3/all.contigs.fa -o 70 -p 98 -z 3 -f 2
	
	echo "---> CAP3 assembly done in sample $s <---" 
   done
}


### Assembly stats ###
######################


assembly_stats()
{
    sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
    for s in ${sepbycomas[@]}; do
	python /media/bec2-jcalvete/Disc_Dur/Feina_Jordi/programes/bin/quast.py --min-contig 100 -o $dir/$s/results/stats/assembly/Assembly_transfuse_stats -t 4 $dir/${s}/results/assembly/transfuse/transfuse_17
    done
    echo "Your Quast statistics for samples $samples have been done"
    
}


### Mapping reads back ###


map_back()
{
    sepbycomas=`sed 's/,/ /gi' <<< "$samples"`
    for s in ${sepbycomas[@]}; do
	bowtie2-build -f $dir/${s}/results/assembly/${s}-rename-contigs.fa $dir/tmp/${s}    
	index="$dir/tmp/${s}"

	bowtie2 -p 5 --very-sensitive -x $index -1 $dir/rawdata//Sample_${s}/${s}_1.fastq.gz -2 $dir/rawdata//Sample_${s}/${s}_2.fastq.gz -S $dir//${s}/mapping/${s}.sam

	samtools view -S -b -q 1 $dir//${s}/mapping/${s}.sam > $dir//${s}/mapping/${s}_Q1.bam

	SortSam.jar I=$dir//${s}/mapping/${s}_Q1.bam O=$dir//${s}/mapping/${s}_Q1_Sort.bam SO=coordinate TMP_DIR=$dir/tmp

	java -Xms2048m -Xmx12288m  -jar /share/apps/local/picard-tools/MarkDuplicates.jar I=$dir//${s}/mapping/${s}_Q1_Sort.bam  O=$dir//${s}/mapping/${s}_Q1_Sort_NoDup.bam M=$dir/${s}/mapping/stats_${s} TMP_DIR=$dir/tmp  REMOVE_DUPLICATES=TRUE CREATE_INDEX=TRUE ASSUME_SORTED=FALSE MAX_RECORDS_IN_RAM=40000000 MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
	
    done
       
}





# ----------------------------------------------------
# ------- Human Assembly pipeline --------------------
# -----------------------------------------------------

EXPECTED_ARGS=1 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments  

##It will be checked for the proper number of files ##

if [ ${#} -ne $EXPECTED_ARGS ]
    then
    echo ''
    echo 'ERROR : Check input parameters'
    echo 'Name: pipeline_ensamblaje_humano.sh'
    echo 'Usage: sh pipeline_ensamblaje_humano.sh sample1,sample2'
    echo ''
    exit $E_BADARGS
fi


# -----------------------------------------
# ------------ main body ------------------
# -----------------------------------------
#create_dir_check_raw_files ${@}
#filtering_and_trimming ${@}
#reducing_complexity ${@}
#assembly_oases ${@}
#assembly_trinity ${@}
stats_oases_trinity ${@}
#kmer_selection ${@}
#assembly_cap3 ${@}
#assembly_stats ${@}
#map_back ${@}
