#!/bin/bash
#####################################################
# Copyright 2012 - Sistemas Genomicos S.L.
# @Desc: Template script for developing bash scripts
# @Author: Arbol
# @Contributors:
#####################################################

# ---------------
# -- functions --
# ---------------

function usage ()
# Displays the script usage
{
	# Parameters:
	scriptname=$1

	echo "$scriptname is a template script for developing bash scripts."
	echo "Usage:"
	echo "  $scriptname <options> input1 input2"
	echo
	echo "* Input (mandatory parameters): *"
	echo " input1     csfasta file"
	echo " input2     qual file"
	echo
	echo "* Options: *"
	echo " -p prefix_			Adds the given string as prefix for the ouput results"
	echo " -o output_folder		Creates the output files in the given output folder, creating it if needed"
	echo
	echo "* Examples of running: *"
	echo " \$$scriptname -prefix ohio_ file1.csfasta file2_QV.qual"
	exit
}


# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
fooGlobalVar=`echo $PATH`
OPTERROR=65

# ------------------------
# -- default parameters --
# ------------------------
prefix=""
outdir=`pwd`
csfasta=""
qual=""

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
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag 
while getopts ":o:p:m:v:s:t" Option; do
	case $Option in
		o) outdir="$OPTARG";;
		p) prefix="$OPTARG";;
		*) echo "ERROR: Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
	esac
done
args=($@)
let index='OPTIND-1'
csfasta=${args[index]}
qual=${args[index+1]}

if [ "$csfasta" == "" ] || [ "$qual" == "" ]; then
	echo "ERROR: Not enough parameters in call to $scriptname!!"
	echo "Call $scriptname with no parameters for help."
	exit
fi


# ---------------
# -- main body -- 
# ---------------

export exec_time=`date`
echo Script started at $exec_time

# Add here the script code | | | | | |
#                          V V V V V V
echo "Hi there! Command launched: $0 prefix=$prefix outdir=$outdir csfasta=$csfasta qual=$qual"

export exec_time=`date`
echo Script finished at $exec_time
