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
#       $1 : sample name
{
	# ---------------
	# - Definitions -
	# ---------------

    local sample=$1

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# Create directories structure
	if [ ! -d $root_dir/$sample/primary ]; then
		mkdir -p $root_dir/$sample/primary/sff
	fi

	if [ ! -d $root_dir/$sample/secondary ]; then
		mkdir -p $root_dir/$sample/secondary
	fi

	if [ ! -d $root_dir/$sample/tertiary ]; then
		mkdir -p $root_dir/$sample/tertiary
	fi
	
	if [ -e $root_dir/$sample/*.sff ]; then
		mv $root_dir/$sample/*.sff $root_dir/$sample/primary/sff
	fi

	# Extract fasta and qual files from sff:
	for j in $root_dir/$sample/primary/sff/*.sff; do
		extract_path_ext_filename $j
		local sff_name=$filename

		sffinfo -s $j > $root_dir/$sample/primary/$sff_name.fasta
		sffinfo -q $j > $root_dir/$sample/primary/$sff_name.qual
		sffinfo -m $j > $root_dir/$sample/primary/$sff_name.mft
	done
}


# --------------------
# - assemble_454 -----
# --------------------

assemble_454 ()
# Description: sequence of actions to carry out the assembly with the 454 software
# Parameters
#	$1 : sample name
{
	# ---------------
	# - Definitions -
	# ---------------
	
	local sample=$1

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	# Create a new directory structure for the assembly
	newAssembly -force $root_dir/$sample/assembly_454
	
	#It's needed to change to the project directory:
	cd "$root_dir/$sample/assembly_454"

	# Add reads file (fasta) to project
	addRun "$root_dir/$sample/assembly_454" "$root_dir/$sample/primary/$sample.fasta"

	# Run sample assembly
	# Please notice that the maximum allowed value for -minlen is 45!!!!!
	runProject -cpu 4 -minlen 45 -at -pair -l 1000 -rip	

	# The assembly file is saved to the secondary folder
	#mv $root_dir/$sample/assembly_454/assembly/454AllContigs.fna $root_dir/$sample/secondary/$sample.fasta
	#ln -sf $root_dir/$sample/secondary/$sample.fasta $root_dir/$sample/assembly_454/assembly/454AllContigs.fna

	# The sequences names get changed, so that these information won't be supplied to customers
	sed 's/  length.*$//' $root_dir/$sample/assembly_454/assembly/454AllContigs.fna > $root_dir/$sample/secondary/${sample}_assembly.fasta
	sed 's/  length.*$//' $root_dir/$sample/assembly_454/assembly/454AllContigs.qual > $root_dir/$sample/secondary/${sample}_assembly.qual
}


# --------------------
# - assemble_results -
# --------------------

assemble_results ()
# Description: sequence of actions to create the results files of the assembly
# Parameters
#	$1 : sample name
{

    # ---------------
	# - Definitions -
	# ---------------

	local sample=$1

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# Creation of the histogram data files:
	# Read lengths histogram
	awk '{if ( $2 == "Assembled" ) {if ( $5 == "-" ) {print $4-$7} else {print $7-$4}}}' "$root_dir/$sample/assembly_454/assembly/454ReadStatus.txt" > "$root_dir/$sample/assembly_454/assembly/histogram_$sample.csv"

	# Read lengths histogram including "PartiallyAssembled"
	awk '{if ( $2 == "Assembled" || $2 =="PartiallyAssembled" ) {if ( $5 == "-" ) {print $4-$7} else {print $7-$4}}}' "$root_dir/$sample/assembly_454/assembly/454ReadStatus.txt" > "$root_dir/$sample/assembly_454/assembly/histogram_with_partially_$sample.csv"

	# Creation of the file with contig numbers and their sizes
	length_fasta.pl "$root_dir/$sample/assembly_454/assembly/454AllContigs.fna" | sort -k2nr > "$root_dir/$sample/assembly_454/assembly/contigs_$sample.csv"

	# Creation of the results file:
	# -----------------------------
	echo "Ensemble metrics & results for sample       $sample" > "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo "Assembly by means of 454 Roche software" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo "--------------------------------------------------------------------------------" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Total number of reads
	echo "Total number of reads in source file $sample.fasta:" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	grep -c ">" "$root_dir/$sample/primary/${sample}.fasta" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo "Total number of reads in assembly file $sample.fasta:" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	grep -c ">" "$root_dir/$sample/secondary/${sample}_assembly.fasta" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	
	# Obtain the N50 value
	echo "N50 value (from 454NewblerMetrics):" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	grep "N50ContigSize" "$root_dir/$sample/assembly_454/assembly/454NewblerMetrics.txt" | awk '{print $3}' >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo "N50 value (from sg_calcularN50):" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	/share/apps/scripts/sg_calcularN50.pl "$root_dir/$sample/assembly_454/assembly/454AllContigs.fna" | grep N50 >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Total number of contigs
	echo Total number of contigs: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	cat "$root_dir/$sample/assembly_454/assembly/contigs_$sample.csv" | wc -l >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Number of contigs above 1kbp
	echo Contigs above 1kbp: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	awk '{if ( $2 >= 1000 ) print $1}' "$root_dir/$sample/assembly_454/assembly/contigs_$sample.csv" | wc -l >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Size of the largest contig
	echo Size of the largest contig: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	head -1 "$root_dir/$sample/assembly_454/assembly/contigs_$sample.csv" | awk '{print $2}' >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Number of assembled reads
	echo Number of assembled reads: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	cat "$root_dir/$sample/assembly_454/assembly/histogram_with_partially_$sample.csv" | wc -l >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Number of assembled bases
	echo Number of assembled bases: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	grep "numAlignedBases" "$root_dir/$sample/assembly_454/assembly/454NewblerMetrics.txt" | awk '{print $3}' | tail -1 >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	#WRONG TRY!!!!!: awk '{SUM+=$1} END {print SUM}' "$1/$6/$4/$5/assembly/histogram_with_partially_$4.csv" >> "$1/$6/$4/$5/assembly/results_$4.csv"

	# Estimated genome size
	echo Estimated genome length: >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	awk '{SUM+=$2} END {print SUM}' "$root_dir/$sample/assembly_454/assembly/contigs_$sample.csv" >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"
	echo >> "$root_dir/$sample/assembly_454/assembly/results_$sample.csv"

	# Print assembly results to general results file:
	cat "$root_dir/$sample/assembly_454/assembly/results_$sample.csv" >> $root_dir/report_TROCAR_plasmids.xls
}


tryyo ()
{

	local sample=$1
    local reference=$2
    local gene_coordinate=$3



	if [ -e $root_dir/$sample/secondary/${sample}_unique_contig.fasta ]; then
        fastas="$root_dir/$sample/secondary/${sample}_unique_contig.fasta"
        names="${sample}_unique_contig"
    else
        fastas="$root_dir/$sample/secondary/${sample}.fasta"
        names="$sample"
    fi
    if [ -e $root_dir/$sample/secondary/${sample}_rest_contigs.fasta ]; then
        fastas="$fastas $root_dir/$sample/secondary/${sample}_rest_contigs.fasta"
        names="$names ${sample}_rest_contigs"
    fi


	for fasta in $fastas; do
		echo "$fasta"
        extract_path_ext_filename $fasta
        local name=$filename
		echo "$name"
	done

	if [ "$reference" == "" ]; then
		echo "PAAAASSOOOOOOOOO"
	fi
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
	
	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# --------------
	# - Prediction -
	# --------------

	# First of all, it will be searched for the contigs .fasta files, since it is possible to have preprocessed them before predicting genes. 
	# It will be searched for fasta files in the secondary directory containing _unique_contig.fasta and _rest_contigs.fasta:
	if [ -e $root_dir/$sample/secondary/${sample}_unique_contig.fasta ]; then
		fastas="$root_dir/$sample/secondary/${sample}_unique_contig.fasta"
	else
		fastas="$root_dir/$sample/secondary/${sample}_assembly.fasta"
	fi
	if [ -e $root_dir/$sample/secondary/${sample}_rest_contigs.fasta ]; then 
		fastas="$fastas $root_dir/$sample/secondary/${sample}_rest_contigs.fasta"
	fi

	# Now the genes prediction for the fasta files is carried out:
	for fasta in $fastas; do
		extract_path_ext_filename $fasta
		local name=$filename

		# If a reference is provided, then the prediction software will be trained:
		if [ ! "$reference" == "" ]; then 

			extract_path_ext_filename $reference
			local reference_dir=$pathname
			local reference_ext=$extension
			local reference_name=$filename

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
			local contigs_number=`grep -c ">" $fasta`
			if [ $contigs_number -eq 1 ]; then
				glimmer3 -o50 -g110 -t30 -b ${reference_dir}/${reference_name}.motif -P $params $fasta ${reference_dir}/${reference_name}.icm $root_dir/$sample/assembly_454/assembly/${name}.orfs
			else
				# Althought bacterial (prokaryotic) genome is circular, if it's not complete (several contigs) then it has to be processed without wraparound (as linear)
				glimmer3 -l -o50 -g110 -t30 -b ${reference_dir}/${reference_name}.motif -P $params $fasta ${reference_dir}/${reference_name}.icm $root_dir/$sample/assembly_454/assembly/${name}.orfs
			fi

		# If a reference is not provided:
		else
		
			local contigs_number=`grep -c ">" $fasta`
			if [ $contigs_number -eq 1 ]; then
				long-orfs -n -t 1.15 $fasta $root_dir/$sample/assembly_454/assembly/${name}.longorfs
				# If only one contig is present, then "extract" can be employed:
				extract -t $fasta $root_dir/$sample/assembly_454/assembly/${name}.longorfs > $root_dir/$sample/assembly_454/assembly/${name}.train
			else
				mkdir -p $root_dir/$sample/assembly_454/assembly/contigs_$sample
				cd $root_dir/$sample/assembly_454/assembly/contigs_$sample
				all2many.pl $fasta 1
				if [ -e $root_dir/$sample/assembly_454/assembly/$name.train ]; then
					rm $root_dir/$sample/assembly_454/assembly/$name.train
				fi
				for contig in *; do
					long-orfs -l -n -t 1.15 $contig ${contig}_longorfs
					extract -t $contig ${contig}_longorfs >> $root_dir/$sample/assembly_454/assembly/$name.train
				done
				cd $root_dir
				rm -rf $root_dir/$sample/assembly_454/assembly/contigs_$sample
			fi

			build-icm -r $root_dir/$sample/assembly_454/assembly/$name.icm < $root_dir/$sample/assembly_454/assembly/$name.train

			if [ $contigs_number -eq 1 ]; then
				glimmer3 -o50 -g110 -t30 $fasta $root_dir/$sample/assembly_454/assembly/$name.icm $root_dir/$sample/assembly_454/assembly/${name}_orfs
			else
				# Althought bacterial (prokaryotic) genome is circular, if it's not complete (several contigs) then it has to be processed without wraparound (as linear)
				glimmer3 -l -o50 -g110 -t30 $fasta $root_dir/$sample/assembly_454/assembly/$name.icm $root_dir/$sample/assembly_454/assembly/${name}_orfs
			fi
		fi

		# --------------
		# - Annotation -
		# --------------
	
		# Run parser on Glimmer's output:
		sg_glimmer_sacar_coordenadas_ORFs_predichos.pl $root_dir/$sample/assembly_454/assembly/${name}_orfs.predict > $root_dir/$sample/assembly_454/assembly/${name}_ORFs_extraer

		# Extract contigs sequences:
		multi-extract -d -w $fasta $root_dir/$sample/assembly_454/assembly/${name}_ORFs_extraer > $root_dir/$sample/tertiary/${name}_ORFs.fasta

		#Final gene predictions in tab separated format (xls):
		sg_reorderFastaFraccGenome.pl $root_dir/$sample/assembly_454/assembly/${name}_orfs.predict > $root_dir/$sample/tertiary/${name}_ORF_prediction.xls

		# In order to make the annotation process faster, the _ORFs.fasta file will be splitted if greater than 1000 seqs:
		local export nseqs=`grep -c ">" $root_dir/$sample/tertiary/${name}_ORFs.fasta`
		if [ $nseqs -gt 1000 ]; then
			cd $root_dir/$sample/assembly_454/assembly/
			/share/apps/scripts/sg_split_fasta_file_en_varios_ficheros.pl $root_dir/$sample/tertiary/${name}_ORFs.fasta 1000
		else
			# if the _ORFs.fasta file is less than 1000 seqs long, then it is copied to the assembly directory wiht "resultado1" name, so that the programm recognises the file to process
			cp $root_dir/$sample/tertiary/${name}_ORFs.fasta $root_dir/$sample/assembly_454/assembly/resultado1
		fi

		# the annotation processes will be launched by means of "sg_commands_nodes_sequencer.sh". This script needs as argument a file containing the commands to run, and will therefore be created
		if [ -e $root_dir/$sample/assembly_454/assembly/commands_${name}.txt ]; then
			rm $root_dir/$sample/assembly_454/assembly/commands_${name}.txt
		fi

		# Iteration over all the "resultado*" files created by sg_split_fasta_file_en_varios_ficheros.pl
		for j in $root_dir/$sample/assembly_454/assembly/resultado*; do
			extract_path_ext_filename $j
			local file_name=$filename
			# XML annotations and normal annotations
			echo blastall -p blastx -d $blastall_db -i $j -o $root_dir/$sample/assembly_454/assembly/ORFs_${name}_${file_name}.fasta.blastx_XML -m 7 -a 4 -e 1e-20 -K 1 >> $root_dir/$sample/assembly_454/assembly/commands_${name}.txt
			echo blastall -p blastx -d $blastall_db -i $j -o $root_dir/$sample/assembly_454/assembly/ORFs_${name}_${file_name}.fasta.blastx -a 4 -e 1e-20 -K 1 >> $root_dir/$sample/assembly_454/assembly/commands_${name}.txt
		done
		sg_commands_nodes_sequencer.sh 2 $root_dir/$sample/assembly_454/assembly/commands_${name}.txt

		# The "resultado" files will be deleted:
		for j in $root_dir/$sample/assembly_454/assembly/resultado*; do
			rm $j
		done
		
		# if the ORFs prediction blastx file already exists, it will be deleted
		if [ -e $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls ]; then
			rm $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls
		fi

		# Blast result will be parsed:
		for j in $root_dir/$sample/assembly_454/assembly/ORFs_${name}_resultado*.fasta.blastx; do
			extract_path_ext_filename $j
			local file_name=$filename
			sg_parsear_resultados_BLASTX.pl $root_dir/$sample/assembly_454/assembly/${file_name}.blastx > ${j}_parseado
			cat ${j}_parseado >> $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls
			rm $j
		done
		
		# Elimination of redundant hits:
		# Not working yet!!!!!!
		sg_blastx.parseado2blastx.parseado.filtrado_Interactive.pl $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls 70 > $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls2
		mv $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls2 $root_dir/$sample/tertiary/${name}_ORFs_blastx.xls
	done
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

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	cd /share/apps/tRNAscan-SE-1.3/
	./tRNAscan-SE -B -o $root_dir/$sample/tertiary/${sample}_tRNA_prediction.xls $root_dir/$sample/secondary/${sample}_assembly.fasta
	awk '{print $1"_"$3"_"$4"\t"$1"\t"$3"\t"$4}' $root_dir/$sample/tertiary/${sample}_tRNA_prediction.xls > $root_dir/$sample/tertiary/ex
	multi-extract -w $root_dir/$sample/secondary/${sample}_assembly.fasta $root_dir/$sample/tertiary/ex > $root_dir/$sample/tertiary/${sample}_tRNAs.fasta
	rm $root_dir/$sample/tertiary/ex
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
{
	# ---------------
	# - definitions -
	# ---------------
	local sample=$1

	# ------------------------
	# - body of the pipeline -
	# ------------------------
	
	# The fastq files will be created (from the fasta & qual files) and will be analyzed. The metrics get output to the primary directory
	local fastq_files=""
	cd $root_dir/$sample/primary/

	for file in $root_dir/$sample/primary/${sample}*.fasta; do
		sg_fastaQual2fastq.pl $file
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
	local reference=$2
	local gene_coordinate=$3

	# ------------------------
	# - body of the pipeline -
	# ------------------------

	# Preliminaries, always important, keep this in mind

	# ORFs (open reading frames) prediction
	ORFs_prediction $sample $reference $gene_coordinate

	# tRNAs prediction
	tRNAs_prediction $sample

	# Enrichment test
	#enrichment_test
	
	# If a reference is provided, then scaffolds are generated:
	if [ ! "$reference" == "" ]; then
	
		# Scaffolds generation from reference sequence:
		cd $root_dir/$sample/assembly_454/assembly/
		promer $reference $root_dir/$sample/secondary/${sample}_assembly.fasta
		mv $root_dir/$sample/assembly_454/assembly/out.delta $root_dir/$sample/assembly_454/assembly/${sample}_out.delta
	
		# Pseudomolecule creation and scaffolding results:
		show-tiling -c -p $root_dir/$sample/tertiary/${sample}_pseudomolecule.fasta $root_dir/$sample/assembly_454/assembly/${sample}_out.delta > $root_dir/$sample/assembly_454/assembly/${sample}_promer_output
		sg_reorderFastaFraccGenome.pl $root_dir/$sample/assembly_454/assembly/${sample}_promer_output > $root_dir/$sample/tertiary/${sample}_scaffolding.xls
		#rm $root_dir/$interm_dir/$sample/secondary/out.delta $root_dir/$interm_dir/$sample/secondary/promer_output
	fi
}


# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

# GLOBAL VARIABLES:
# -----------------
# Environment variables (those which depend on the directory structure, for example, locations of references, annotations, etc.)
blastall_db="/data/results/Solid0065/blast_db/refseq_protein"

# Sample variables (those dependant on the sample directories)
root_dir="/data/results/Solid0065/BF5_Trocar/plasmidos"
firstseq="secuencias_primer_batch"
secondseq="secuencias_segundo_batch"

samples=("CTXM1" "H67" "HEC9" "pA1" "p1433" "pY41")
references=("" "" "" "" "" "")
genes_coordinates=("" "" "" "" "" "")

# -------------------
# - FUNCTION CALLS: -
# -------------------
# Functions/procedures calls
echo "Assembly results:" > $root_dir/report_TROCAR_plasmids.xls
echo "-----------------" >> $root_dir/report_TROCAR_plasmids.xls
echo >> $root_dir/report_TROCAR_plasmids.xls
for ((k=0;k<${#samples[@]};k+=1)); do
	
	export exec_time=`date`
	echo Starting conversion to fasta/qual from dff of ${samples[k]} on $exec_time
	#dff_to_fasta ${samples[k]} 
	
	export exec_time=`date`
	echo Starting fastqc report of ${samples[k]} on $exec_time
	#fastqc_report ${samples[k]}
	
	export exec_time=`date`
	echo Starting assembly of ${samples[k]} on $exec_time
	#assemble_454 ${samples[k]}
	
	export exec_time=`date`
	echo Starting creation of result files of ${samples[k]} on $exec_time
	#assemble_results ${samples[k]}

	# call of the genes prediction & annotation function
	export exec_time=`date`
	echo "Starting genes prediction & annotation of ${samples[k]} on $exec_time"
	genes_prediction_annotation ${samples[k]} ${references[k]} ${genes_coordinates[k]}
done

export exec_time=`date`
echo Pipeline finished on $exec_time
