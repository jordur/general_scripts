#!/bin/bash

# -------------------------------------------
# ---------------- main body ----------------
# -------------------------------------------

# Function to extract the desired reads (reads) from the target csfasta/fasta and qual files (samples)
# The user only has to redefine the following parameters:


echo 'Function name: sg_extract_reads_fasta_qual.sh'
echo 'Usage: the script sg_extract_reads_fasta_qual.sh has to be edited in order to add the root path, csfasta/fasta & qual files, the desired reads to be extracted and the output file names.'
echo 'Dont worry, just go to the parameters section of the script and edit them to your convenience.'
echo 'By the way, if the output files already exist, you will be asked to delete them. If you dont do so, the new reads will be added to the existing output files, so BE AWARE!!'


# -------------------------------------------
# ---------------- parameters ----------------
# -------------------------------------------


# Root directory where the files are stored in
root_dir="./"
# Names of the csfasta/fasta and qual files to work on
samples=("solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_06S213.csfasta" "solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F3_QV_06S213.qual" "solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F5-BC_06S213.csfasta" "solid0065_20110110_PE_BC_FC2_Cardio_MamaColon_Mama_Colon_F5-BC_QV_06S213.qual")
# Names of the reads to be searched in the csfasta/fasta and/or qual files
reads=("217_3_1311" "713_4_419" "538_3_1779" "433_4_577" "73_2_1515")
# Names of the output files to be created
output_files=("F3.csfasta" "F3.qual" "F5.csfasta" "F5.qual")

for ((k=0;k<${#samples[@]};k+=1)); do
	if [ -e $root_dir${output_files[k]} ]; then
		echo "Do you want to remove the $root_dir${output_files[k]} file?? (y/n) followed by [ENTER]:"
		read answer
		if [ "$answer" == "y" ] || [ "$answer" == "Y" ]; then
			rm $root_dir${output_files[k]}
		fi
	fi
	for ((j=0;j<${#reads[@]};j+=1)); do
                echo ${reads[$j]} $1${samples[k]} $1${output_files[k]}
		grep -A 1 ${reads[$j]} $root_dir${samples[k]} >> $root_dir${output_files[k]}
                #grep -A 1 ${reads[$j]} $1$2 >> $1$3
        done
done
