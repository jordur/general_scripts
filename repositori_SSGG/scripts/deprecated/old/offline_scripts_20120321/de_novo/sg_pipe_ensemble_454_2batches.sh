#!/bin/bash


# -----------------
# - extract_path --
# -----------------

function extract_path_ext_filename ()
# Description: function to extract the path, extension and filename (whithout extension) from a whole filename
# Its results are passed as global variables ($filename $extension $pathname)
# Parameters
#       $1 : whole filename
{
	# ---------------
	# - Definitions -
	# ---------------


	# ------------------------
	# - body of the pipeline -
	# ------------------------

	filename=$(basename $1)
	extension=${filename##*.}
	pathname=`echo $1 | sed "s/$filename//" | sed "s/\/*$//g"`
	filename=$(basename $1 .$extension)
}


# -------------------
# - dff_to_fasta ----
# -------------------

dff_to_fasta ()
# Description: sequence of actions to transform the dff files into fasta format
# Parameters
#       $1 : root directory
#       $2 : first sequence directory sufix
#       $3 : second sequence directory sufix
#       $4 : sample name
#       $5 : output directory sufix
{
	# ---------------
	# - Definitions -
	# ---------------


	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	for j in $1/$2/$4/*.sff; do
		/opt/454/bin/sffinfo -s $j > $j.fasta
		/opt/454/bin/sffinfo -q $j > $j.qual
		/opt/454/bin/sffinfo -m $j > $j.mft
	done
}


# --------------------
# - assemble_454 -----
# --------------------

assemble_454 ()
# Description: sequence of actions to carry out the assembly with the 454 software
# Parameters
#	$1 : root directory
#	$2 : first sequence directory sufix
#	$3 : second sequence directory sufix
#	$4 : sample name
#	$5 : output directory sufix (this parameter also indicates the type of assembly which is going to be carried out)
#	$6 : intermediate directory sufix
#	$7 : original investigator name
{
	# ---------------
	# - Definitions -
	# ---------------
	
	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	#Create the output directories for the new assemblies
	if [ ! -d "$1/$6" ]; then
		mkdir $1/$6
	fi

	if [ ! -d "$1/$6/$4" ]; then
		mkdir $1/$6/$4
	fi

	# Create a new directory structure for the assembly
	newAssembly -force $1/$6/$4/$5 #$1/$5/$4
	
	# The different sequences reads will be concatenated
	# It is considered if the fasta and qual files from both batches are needed or not
	if [ "$5" != "2nd_batch" ]; then
		if [ -e $1/$6/$4/$7_batch1.fasta ]; then
			rm $1/$6/$4/$7_batch1.fasta
		fi
		for j in $1/$2/$4/*_sequence.fna; do
			if [ ! -e $1/$6/$4/$7_batch1.fasta ]; then
				cat $j > $1/$6/$4/$7_batch1.fasta
			else
				cat $j >> $1/$6/$4/$7_batch1.fasta
			fi
		done
		if [ -e $1/$6/$4/$7_batch1.qual ]; then
			rm $1/$6/$4/$7_batch1.qual
		fi
		for j in $1/$2/$4/*_sequence.qual; do
			if [ ! -e $1/$6/$4/$7_batch1.qual ]; then
				cat $j > $1/$6/$4/$7_batch1.qual
			else
				cat $j >> $1/$6/$4/$7_batch1.qual
			fi
		done
	fi

	if [ "$5" != "1st_batch" ]; then
		if [ -e $1/$6/$4/$7_batch2.fasta ]; then
			rm $1/$6/$4/$7_batch2.fasta
		fi
		for j in $1/$3/$4/*_sequence.fna; do
			if [ ! -e $1/$6/$4/$7_batch2.fasta ]; then
				cat $j > $1/$6/$4/$7_batch2.fasta
			else
				cat $j >> $1/$6/$4/$7_batch2.fasta
			fi
		done
		if [ -e $1/$6/$4/$7_batch2.qual ]; then
			rm $1/$6/$4/$7_batch2.qual
		fi
		for j in $1/$3/$4/*_sequence.qual; do
			if [ ! -e $1/$6/$4/$7_batch2.qual ]; then
				cat $j > $1/$6/$4/$7_batch2.qual
			else
				cat $j >> $1/$6/$4/$7_batch2.qual
			fi
		done
	fi

	#It's needed to change to the project directory:
	cd "$1/$6/$4/$5"

	# 2E+1E strategy: first the 2nd batch is added and ensembled, then 1st batch is added and reensembled
	if [ "$5" = "2E+1E" ]; then
	
		#Anyadir el fichero fasta con las lecturas. Usar la extension .fna para identificar este fichero. El fichero de calidad asociado a la muestra se tiene que llamar igual que el fichero fasta que contiene las lecturas. Ejemplo: NG-5250_4_in.454.fna NG-5250_4_in.454.qual Si no se llaman igual el programa no va a reconocerlos!!
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch2.fasta"

		#Se ejecuta el ensamblaje de esta muestra. El numero de cpu se puede cambiar dependiendo del nodo que estemos usando y su disponibilidad.
		# Please notice that the maximum allowed value for -minlen is 45!!!!!
		runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip
	
		#Para anadir mas lecturas se ejecuta de nuevo el segundo paso
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch1.fasta"

		#Para unir estas nuevas lecturas al ensamblaje anterior se ejecuta el paso 3
		runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip

		#Si nos equivocamos de fichero .fna y queremos eliminarlo del proyecto
		#removeRun nombre_directorio_proyecto fichero_fasta_a_eliminar.fna
	fi
	
	# 1+2 strategy: first of all, batches first and second are added, then get ensembled
	if [ "$5" = "1+2" ]; then
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch2.fasta"
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch1.fasta"
		runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip
	fi
	
	# 2nd_batch strategy: only 2nd batch is added and ensembled
	if [ "$5" = "2nd_batch" ]; then
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch2.fasta"
		runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip
	fi

	# 1st_batch: only 1st batch is added and ensembled
	if [ "$5" = "1st_batch" ]; then
		addRun "$1/$6/$4/$5"  "$1/$6/$4/$7_batch1.fasta"
		runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip
	fi
}


# --------------------
# - assemble_results -
# --------------------

assemble_results ()
# Description: sequence of actions to create the results files of the assembly
# Parameters
#	$1 : root directory
#	$2 : first sequence directory sufix
#	$3 : second sequence directory sufix
#	$4 : sample name
#	$5 : output directory sufix (this parameter also indicates the type of assembly which is going to be carried out)
#	$6 : intermediate directory sufix
#	$7 : original investigator name
{

    # ---------------
	# - Definitions -
	# ---------------

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# Creation of the histogram data files:
	# Read lengths histogram
	awk '{if ( $2 == "Assembled" ) {if ( $5 == "-" ) {print $4-$7} else {print $7-$4}}}' "$1/$6/$4/$5/assembly/454ReadStatus.txt" > "$1/$6/$4/$5/assembly/histogram_$4.csv"

	# Read lengths histogram including "PartiallyAssembled"
	awk '{if ( $2 == "Assembled" || $2 =="PartiallyAssembled" ) {if ( $5 == "-" ) {print $4-$7} else {print $7-$4}}}' "$1/$6/$4/$5/assembly/454ReadStatus.txt" > "$1/$6/$4/$5/assembly/histogram_with_partially_$4.csv"

	# Creation of the file with contig numbers and their sizes
	/share/apps/scripts/length_fasta.pl "$1/$6/$4/$5/assembly/454AllContigs.fna" | sort -k2nr > "$1/$6/$4/$5/assembly/contigs_$4.csv"

	# Creation of the results file:
	# -----------------------------
	echo "Assembly type:        $5" > "$1/$6/$4/$5/assembly/results_$4.csv"
	echo "Ensemble metrics & results for sample       $4" >> "$1/$6/$4/$5/assembly/results_$4.csv"
	echo "--------------------------------------------------------------------------------" >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Total number of reads
	if [ "$5" != "2nd_batch" ]; then
		echo "Total number of reads in source file $7_batch1.fasta:" >> "$1/$6/$4/$5/assembly/results_$4.csv"
		grep -c ">" "$1/$6/$4/$7_batch1.fasta" >> "$1/$6/$4/$5/assembly/results_$4.csv"
	fi
	
	if [ "$5" != "1st_batch" ]; then
		echo "Total number of reads in source file $7_batch2.fasta:" >> "$1/$6/$4/$5/assembly/results_$4.csv"
		grep -c ">" "$1/$6/$4/$7_batch2.fasta" >> "$1/$6/$4/$5/assembly/results_$4.csv"
	fi

	# Obtain the N50 value
	echo "N50 value (from 454NewblerMetrics):" >> "$1/$6/$4/$5/assembly/results_$4.csv"
	grep "N50ContigSize" "$1/$6/$4/$5/assembly/454NewblerMetrics.txt" | awk '{print $3}' >> "$1/$6/$4/$5/assembly/results_$4.csv"
	echo "N50 value (from sg_calcularN50):" >> "$1/$6/$4/$5/assembly/results_$4.csv"
	/share/apps/scripts/sg_calcularN50.pl "$1/$6/$4/$5/assembly/454AllContigs.fna" | grep N50 >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Total number of contigs
	echo Total number of contigs: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	grep "" "$1/$6/$4/$5/assembly/contigs_$4.csv" | wc -l >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Number of contigs above 1kbp
	echo Contigs above 1kbp: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	awk '{if ( $2 >= 1000 ) print $1}' "$1/$6/$4/$5/assembly/contigs_$4.csv" | wc -l >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Size of the largest contig
	echo Size of the largest contig: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	head -1 "$1/$6/$4/$5/assembly/contigs_$4.csv" | awk '{print $2}' >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Number of assembled reads
	echo Number of assembled reads: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	grep -c "" "$1/$6/$4/$5/assembly/histogram_with_partially_$4.csv" >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Number of assembled bases
	echo Number of assembled bases: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	grep "numAlignedBases" "$1/$6/$4/$5/assembly/454NewblerMetrics.txt" | awk '{print $3}' | tail -1 >> "$1/$6/$4/$5/assembly/results_$4.csv"
	#WRONG TRY!!!!!: awk '{SUM+=$1} END {print SUM}' "$1/$6/$4/$5/assembly/histogram_with_partially_$4.csv" >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Estimated genome size
	echo Estimated genome length: >> "$1/$6/$4/$5/assembly/results_$4.csv"
	awk '{SUM+=$2} END {print SUM}' "$1/$6/$4/$5/assembly/contigs_$4.csv" >> "$1/$6/$4/$5/assembly/results_$4.csv"
}


# --------------------
# - assemble_report --
# --------------------

assemble_report ()
# Description: creation of a single report for each type of assembly with the results of all the samples
# Parameters
#	$1 : root directory
#	$2 : output directory sufix (this parameter also indicates the type of assembly which is going to be carried out)
#	$3 : intermediate directory sufix
#	$4 : sample number
#	$5 : output directory
#	$6 : intermediate directory
#	$7 : original investigator name
{

	# ---------------
	# - definitions -
	# ---------------

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	echo "Sample: $4 (original name: $7 )" >> $1/$6/$4/report_$4.xls
	echo ------------------------------- >> $1/$6/$4/report_$4.xls

	cat "$1/$6/$4/$5/assembly/results_$4.csv" >> $1/$6/$4/report_$4.xls
	echo ------------------------------- >> $1/$6/$4/report_$4.xls
	echo                                 >> $1/$6/$4/report_$4.xls
}


# ---------------------------
# - create_assemble_report --
# ---------------------------

create_assemble_report ()
# Description: creation of a single report for each type of assembly with the results of all the samples
# Parameters:
#	the function uses the parameters defined in main (global and sample parameters)
{
	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	for ((i=0;i<${#output_dir[@]};i+=1)); do
		for ((k=0;k<${#samples[@]};k+=1)); do
			# delete the existent assemble reports
			if [ -e "$root_dir/$interm_dir/${samples[k]}/report_${samples[k]}.xls" ] && [ $i -eq 0  ]; then
				rm "$root_dir/$interm_dir/${samples[k]}/report_${samples[k]}.xls"
			fi
			assemble_report $root_dir $firstseq $secondseq ${samples[k]} ${output_dir[i]} $interm_dir ${original_names[k]}
		done
	done
}


# -------------------
# - ORFs_prediction -
# -------------------

ORFs_prediction ()
# Description: Carries out the prediction and annotation of open reading frames (ORFs)
# Parameters
#	See definitions section
{
	# ---------------
	# - definitions -
	# ---------------

	local sample=$1
	local reference=$2
	local gene_coordinate=$3
	local original_name=$4
	extract_path_ext_filename $reference
	local reference_dir=$pathname
	local reference_ext=$extension
	local reference_name=$filename
		
	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# --------------
	# - Prediction -
	# --------------

	# Cambiar el nombre de las secuencias. Entregaremos los ficheros modificados
	sed 's/  length.*$//' $root_dir/$interm_dir/$sample/assembly/454AllContigs.fna > $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna
	sed 's/  length.*$//' $root_dir/$interm_dir/$sample/assembly/454AllContigs.qual > $root_dir/$interm_dir/$sample/secondary/454AllContigs.qual
	# Extract known genes sequences in fasta format
	if [ ! -e "${reference_dir}/${reference_name}_knowngenes.${reference_ext}" ]; then
		local contigs_number=`grep -c ">" $reference`
		if [ $contigs_number -eq 1 ]; then
			# If only one contig is present, then "extract" can be employed:
			extract $reference $gene_coordinate > ${reference_dir}/${reference_name}_knowngenes.${reference_ext}
		else
			# If more than one contig is present, then "multi-extract" is required. Multi-extract needs a different coordinates file which has to be created:
			awk '{sum=sum+1; print sum"\t"$0}' $gene_coordinate > ${gene_coordinate}_2
			multi-extract -w $reference ${gene_coordinate}_2 > ${reference_dir}/${reference_name}_knowngenes.${reference_ext}
		fi
	fi
	
	# Obtain likelihood model of coding sequences (ICM)
	if [ ! -e "${reference_dir}/${reference_name}.icm" ]; then
		build-icm -r ${reference_dir}/${reference_name}.icm < ${reference_dir}/${reference_name}_knowngenes.${reference_ext}
	fi

	# The upstream sequences of the known genes are obtained, so that the ribosome binding site can be calculated by means of a probability matrix
	if [ ! -e "${reference_dir}/${reference_name}.upstream" ]; then
		# The "extract" program was supposed to require non-inverted values when those are in the negative strand (for this reason the "awk" was added initially and now is commented):
		#upstream-coords.awk 25 0 $gene_coordinate | awk '{ if ($2<$3) print $1"\t"$2"\t"$3; else print $1"\t"$3"\t"$2 }' | extract $reference - > ${reference_dir}/${reference_name}.upstream
		upstream-coords.awk 25 0 $gene_coordinate | extract $reference - > ${reference_dir}/${reference_name}.upstream
	fi
	if [ ! -e "${reference_dir}/${reference_name}.motif" ]; then
		elph ${reference_dir}/${reference_name}.upstream LEN=6 | get-motif-counts.awk > ${reference_dir}/${reference_name}.motif
	fi

	# Calculate the frecuency of the 3 possible start codons (this parameter will be passed to the next step)
	local params=`start-codon-distrib -3 ${reference_dir}/${reference_name}_knowngenes.${reference_ext} $gene_coordinate`

	# Utilizar los ficheros motifs/-P/icm/ensamblaje que corresponda (-P abrir *.start-codon-distrib)
	local contigs_number=`grep -c ">" $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna`
	if [ $contigs_number -eq 1 ]; then
		glimmer3 -o50 -g110 -t30 -b ${reference_dir}/${reference_name}.motif -P $params $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna ${reference_dir}/${reference_name}.icm $root_dir/$interm_dir/$sample/assembly/orfs
	else
		# Althought bacterial (prokaryotic) genome is circular, if it's not complete (several contigs) then it has to be processed without wraparound (as linear)
		glimmer3 -l -o50 -g110 -t30 -b ${reference_dir}/${reference_name}.motif -P $params $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna ${reference_dir}/${reference_name}.icm $root_dir/$interm_dir/$sample/assembly/orfs
	fi

	# --------------
	# - Annotation -
	# --------------
	
	#Ejecutar un parser para la salida de Glimmer
	sg_glimmer_sacar_coordenadas_ORFs_predichos.pl $root_dir/$interm_dir/$sample/assembly/orfs.predict > $root_dir/$interm_dir/$sample/assembly/ORFs_extraer

	#Extraer las secuencias de los contigs
	multi-extract -d -w $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna $root_dir/$interm_dir/$sample/assembly/ORFs_extraer > $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs.fasta

	#Final gene predictions in tab separated format (xls):
	sg_reorderFastaFraccGenome.pl $root_dir/$interm_dir/$sample/assembly/orfs.predict > $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORF_prediction.xls

	# In order to make the annotation process faster, the _ORFs.fasta file will be splitted if greater than 1000 seqs:
	local export nseqs=`grep -c ">" $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs.fasta`
	if [ $nseqs -gt 1000 ]; then
		cd $root_dir/$interm_dir/$sample/assembly/
		/share/apps/scripts/sg_split_fasta_file_en_varios_ficheros.pl $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs.fasta 1000
	else
		# if the _ORFs.fasta file is less than 1000 seqs long, then it is copied to the assembly directory wiht "resultado1" name, so that the programm recognises the file to process
		cp $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs.fasta $root_dir/$interm_dir/$sample/assembly/resultado1
	fi

	# the annotation processes will be launched by means of "sg_commands_nodes_sequencer.sh". This script needs as argument a file containing the commands to run, and will therefore be created
    if [ -e $root_dir/$interm_dir/$sample/assembly/commands_${original_name}.txt ]; then
		rm $root_dir/$interm_dir/$sample/assembly/commands_${original_name}.txt
	fi

	# Iteration over all the "resultado*" files created by sg_split_fasta_file_en_varios_ficheros.pl
	for j in $root_dir/$interm_dir/$sample/assembly/resultado*; do
		extract_path_ext_filename $j
		local file_name=$filename
		# XML annotations and normal annotations
		echo blastall -p blastx -d $blastall_db -i $j -o $root_dir/$interm_dir/$sample/assembly/ORFs_${original_name}_${file_name}.fasta.blastx_XML -m 7 -a 4 -e 1e-20 -K 1 >> $root_dir/$interm_dir/$sample/assembly/commands_${original_name}.txt
		echo blastall -p blastx -d $blastall_db -i $j -o $root_dir/$interm_dir/$sample/assembly/ORFs_${original_name}_${file_name}.fasta.blastx -a 4 -e 1e-20 -K 1 >> $root_dir/$interm_dir/$sample/assembly/commands_${original_name}.txt
	done
	sg_commands_nodes_sequencer.sh 2 $root_dir/$interm_dir/$sample/assembly/commands_${original_name}.txt

	# The "resultado" files will be deleted:
	for j in $root_dir/$interm_dir/$sample/assembly/resultado*; do
		rm $j
	done
		
	# if the ORFs prediction blastx file already exists, it will be deleted
	if [ -e $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls ]; then
		rm $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls
	fi

	# Blast result will be parsed:
	for j in $root_dir/$interm_dir/$sample/assembly/ORFs_${original_name}_resultado*.fasta.blastx; do
		extract_path_ext_filename $j
		local file_name=$filename
		sg_parsear_resultados_BLASTX.pl $root_dir/$interm_dir/$sample/assembly/${file_name}.blastx > ${j}_parseado
		cat ${j}_parseado >> $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls
		rm $j
	done
		
	# Elimination of redundant hits:
	# Not working yet!!!!!!
	sg_blastx.parseado2blastx.parseado.filtrado_Interactive.pl $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls 70 > $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls2
	mv $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls2 $root_dir/$interm_dir/$sample/tertiary/${original_name}_ORFs_blastx.xls
}


# --------------------
# - tRNAs_prediction -
# --------------------

tRNAs_prediction ()
# Description: Carries out the prediction and annotation of mRNAs
# Parameters
#       $1 : current sample
{
	# ---------------
	# - Definitions -
	# ---------------
	local sample=$1
	local original_name=$2

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	cd /share/apps/tRNAscan-SE-1.3/
	./tRNAscan-SE -B -o $root_dir/$interm_dir/$sample/tertiary/${original_name}_tRNA_prediction.xls $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna
	awk '{print $1"_"$3"_"$4"\t"$1"\t"$3"\t"$4}' $root_dir/$interm_dir/$sample/tertiary/${original_name}_tRNA_prediction.xls > $root_dir/$interm_dir/$sample/tertiary/ex
	multi-extract -w secondary/454AllContigs.fna $root_dir/$interm_dir/$sample/tertiary/ex > $root_dir/$interm_dir/$sample/tertiary/${original_name}_tRNAs.fasta
	rm $root_dir/$interm_dir/$sample/tertiary/ex
}


# -------------------
# - enrichment test -
# -------------------

enrichment_test ()
# Description: Carries out the comparisonal enrichment test between the reference organism and the sequenced organism
# Parameters
#       $1 : current sample
{
	# ---------------
	# - definitions -
	# ---------------
	local sample=$1
	local original_name=$2

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	/share/apps/scripts/transcriptoma/sg_binomial_test_enrichment.pl gene_ontology_annotation8.txt ../../../mira/a_baumani/functionalAnnotation/a_baumanniReference_GOannotation.txt 0.05
}


# -----------------
# - fastqc_report -
# -----------------

fastqc_report ()
# Description: Creates a single fastq file from all the existent fasta and qual files from the different batches, and then it analyzes the fastq file and creates the corresponding metrics directory (in the primary directory) 
# Parameters
#       $1 : current sample
#       $2 : original name of the sample given by the investigator
{
	# ---------------
	# - definitions -
	# ---------------
	local sample=$1
	local original_name=$2

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	# The following lines are not anymore needed, since the fasta_qual file names include already the original investigator name!!!

	# First the fasta_qual file names will be changed, including the original name of the investigator at its beginning
#	if [ ! -e $root_dir/$interm_dir/$sample/primary/names_updated.flag ]; then
#		echo "Names of the files in the directory already changed!!" > $root_dir/$interm_dir/$sample/primary/names_updated.flag
#		for file in $root_dir/$interm_dir/$sample/primary/${sample}_*.fasta; do
#			extract_path_ext_filename $file
#			local file_dir=$pathname
#			local file_ext=$extension
#		    local file_name=$filename
#			mv $file $file_dir/${original_name}_${file_name}.$extension
#			mv $file_dir/${file_name}.qual $file_dir/${original_name}_${file_name}.qual
#		done
#	fi

	# The fastq files will be created (from the fasta & qual files) and will be analyzed. The metrics get output to the primary directory
	local fastq_files=""
	cd $root_dir/$interm_dir/$sample/primary/

	for file in $root_dir/$interm_dir/$sample/primary/${original_name}_*.fasta; do
		/share/apps/scripts/sg_fastaQual2fastq.pl $file
		extract_path_ext_filename $file
		local file_dir=$pathname
		local file_ext=$extension
		local file_name=$filename
		export fastq_files="$fastq_files ${file_dir}/${file_name}.fastq"
	done
	/share/apps/FastQC/fastqc $fastq_files
}


# -------------------------------
# - genes_prediction_annotation -
# -------------------------------

genes_prediction_annotation ()
# Description: predicts the genes ORFs& and carries out its annotation
# Parameters
#	see definitions section
{
	# ---------------
	# - definitions -
	# ---------------
	
	local sample=$1
	local chosen_assembly=$2
	local reference=$3 
	local gene_coordinate=$4 
	local original_name=$5 
	local interm_dir=$6 

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# Preliminaries, always important, keep this in mind

	# It will be checked if the different assembly directories still exist. If so, then the chosen assembly files will be moved to $sample/assembly
	#and the existent assembly directories will be deleted
	if [ -d "$root_dir/$interm_dir/$sample/$chosen_assembly" ]; then
		# Directory to store all the scratch-existent files is created ("to_delete_date")
		local date_suffix=`date +"%Y%m%d%H%M%S"`

		if [ ! -d "$root_dir/$interm_dir/$sample/to_delete_${date_suffix}" ]; then
			mkdir $root_dir/$interm_dir/$sample/to_delete_${date_suffix}
		fi
		# if there is an existent assembly directory, it is moved to "to_delete"
		if [ -d "$root_dir/$interm_dir/$sample/assembly" ]; then
			mv $root_dir/$interm_dir/$sample/assembly $root_dir/$interm_dir/$sample/to_delete_${date_suffix}
		fi
	
		mv $root_dir/$interm_dir/$sample/$chosen_assembly/assembly $root_dir/$interm_dir/$sample/

		# Construction of the directories structure (primary, secondary, etc.)
		if [ ! -d "$root_dir/$interm_dir/$sample/primary" ]; then
			mkdir $root_dir/$interm_dir/$sample/primary
		fi
		mv $root_dir/$interm_dir/$sample/${original_name}_*.fasta $root_dir/$interm_dir/$sample/primary
		mv $root_dir/$interm_dir/$sample/${original_name}_*.qual $root_dir/$interm_dir/$sample/primary

		if [ ! -d "$root_dir/$interm_dir/$sample/secondary" ]; then
			mkdir $root_dir/$interm_dir/$sample/secondary
		fi
		#cp $root_dir/$interm_dir/$sample/assembly/454AllContigs.* $root_dir/$interm_dir/$sample/secondary
		
		if [ ! -d "$root_dir/$interm_dir/$sample/tertiary" ]; then
			mkdir $root_dir/$interm_dir/$sample/tertiary
		fi
		
		# move existent assembly directories to "to_delete"
		for ((i=0;i<${#output_dir[@]};i+=1)); do
			mv $root_dir/$interm_dir/$sample/${output_dir[i]} $root_dir/$interm_dir/$sample/to_delete_${date_suffix}
		done
	else
		echo "Chosen assembly directory non existent. It is assumed that the chosen assembly directory is the current assembly dir!"
	fi

	# Create fastqc report (quality of incoming data)
	#fastqc_report $sample $original_name

	# ORFs (open reading frames) prediction
	ORFs_prediction $sample $reference $gene_coordinate $original_name

	# tRNAs prediction
	#tRNAs_prediction $sample $original_name

	# Enrichment test
	#enrichment_test
	
	# Scaffolds generation from reference sequence:
	cd $root_dir/$interm_dir/$sample/assembly/
	promer $reference $root_dir/$interm_dir/$sample/secondary/454AllContigs.fna
	mv $root_dir/$interm_dir/$sample/assembly/out.delta $root_dir/$interm_dir/$sample/assembly/${original_name}_out.delta
	
	# Pseudomolecule creation and scaffolding results:
	show-tiling -c -p $root_dir/$interm_dir/$sample/tertiary/${original_name}_pseudomolecule.fasta $root_dir/$interm_dir/$sample/assembly/${original_name}_out.delta > $root_dir/$interm_dir/$sample/assembly/${original_name}_promer_output
	sg_reorderFastaFraccGenome.pl $root_dir/$interm_dir/$sample/assembly/${original_name}_promer_output > $root_dir/$interm_dir/$sample/tertiary/${original_name}_scaffolding.xls
	#rm $root_dir/$interm_dir/$sample/secondary/out.delta $root_dir/$interm_dir/$sample/secondary/promer_output
}


# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

# GLOBAL VARIABLES:
# -----------------
# Environment variables (those which depend on the directory structure, for example, locations of references, annotations, etc.)
blastall_db="/data/results/Solid0065/blast_db/refseq_protein"

# Sample variables (those dependant on the sample directories)
root_dir="/data/results/Solid0065/BF5_Trocar/GATC_454"
firstseq="secuencias_primer_batch"
secondseq="secuencias_segundo_batch"

samples=("NG-5250_7_sp" "NG-5250_8_junii")
original_names=("41-2" "209")
chosen_assemblies=("1+2" "1+2")
references=("/data/results/solid0065/BF5_Trocar/GATC_454/mira/a_sp/a_sp.fna" "/data/results/Solid0065/BF5_Trocar/GATC_454/mira/a_junii/a_junii.fna")
genes_coordinates=("/data/results/solid0065/BF5_Trocar/GATC_454/mira/a_sp/a_sp.train" "/data/results/Solid0065/BF5_Trocar/GATC_454/mira/a_junii/a_junii.train")

output_dir=("1st_batch" "2nd_batch" "2E+1E" "1+2")
interm_dir="NUEVOS_ENSAMBLAJES"

# -------------------
# - FUNCTION CALLS: -
# -------------------
# Functions/procedures calls
for ((i=0;i<${#output_dir[@]};i+=1)); do
	for ((k=0;k<${#samples[@]};k+=1)); do
		#dff_to_fasta $root_dir $firstseq $secondseq ${samples[k]} $output_dir #please review dff_to_fasta function before its use, and notice that usually you get the fasta and qual files from the 454 system
		export exec_time=`date`
		echo Starting assembly of ${samples[k]} - ${output_dir[i]} on $exec_time
		#assemble_454 $root_dir $firstseq $secondseq ${samples[k]} ${output_dir[i]} $interm_dir ${original_names[k]}
		
		export exec_time=`date`
		echo Starting creation of result files of ${samples[k]} - ${output_dir[i]} on $exec_time
		#assemble_results $root_dir $firstseq $secondseq ${samples[k]} ${output_dir[i]} $interm_dir ${original_names[k]}
	done
done

# call of the genes prediction & annotation function
for ((k=0;k<${#samples[@]};k+=1)); do
	export exec_time=`date`
	echo "Starting genes prediction & annotation of ${samples[k]} - ${chosen_assemblies[k]} creation of assemble report on $exec_time"
	genes_prediction_annotation ${samples[k]} ${chosen_assemblies[k]} ${references[k]} ${genes_coordinates[k]} ${original_names[k]} $interm_dir
done

# create the assemble reports:
export exec_time=`date`
echo Starting creation of assemble report on $exec_time
#create_assemble_report

export exec_time=`date`
echo Pipeline finished on $exec_time
