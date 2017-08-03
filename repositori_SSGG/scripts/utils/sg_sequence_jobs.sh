#!/bin/bash
#####################################################
# Copyright 2012 - Sistemas Genómicos S.L.
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

        echo "$scriptname is a script for launching a sequence of commands over the queue/resources manager."
        echo "Usage:"
        echo "  $scriptname <options> input1"
        echo
        echo "* Input (mandatory parameters): *"
        echo " input1     file containing the list of commands/jobs/tasks to launch in sequence (one after another)"
        echo
        echo "* Options: *"
	echo " -t interpreter_type	type of interpreter to be used. At the moment one of the following: bash, perl or python"
	echo " -w wait_time		time to be waited between checks for jobs to be finished"
        echo
        echo "* Examples of running: *"
        echo " \$$scriptname -n 5 tasks.txt"

        exit
}


# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
EXPECTED_ARGS=2 # Number of arguments expected in the script call
OPTERROR=65

# ------------------------
# -- default parameters --
# ------------------------

interpreter="bash" # type of interpreter to be used
wait_time=1 # number of seconds to wait between processes status checks
tasks_file=""

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
while getopts "t:w:N:" Option; do
	case $Option in
		t)
			interpreter="$OPTARG";; 
		w)
			wait_time="$OPTARG";;
		N)
			jobname="$OPTARG";;
		*)
			echo "Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
	esac
done
args=($@)
let index='OPTIND-1'
task_file=${args[index]}

if [ "$task_file" == "" ]; then
        echo "Not enough parameters in call to $scriptname!!"
        echo "Call $scriptname with no parameters for help."
        exit
fi

if [ ! -e $task_file ]; then
	echo "Not existing file $task_file!!"
	exit
fi

if [ "$interpreter" != "bash" ] && [ "$interpreter" != "perl" ] && [ "$interpreter" != "perl" ]; then
	echo "Invalid value for parameter -t in call to $scriptname!!"
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

# Array for processes IDs:
declare -a proc_commands # array which contains the different commands to launch

number_commands=`wc -l $task_file | awk '{print $1}'` # number of commands to launch
wait_time=1 # number of seconds to wait between processes status checks

# the file with the command lines will be read. First the field separator environment variable is set to "new line".
old_IFS=$IFS
IFS=$'\n'
proc_commands=( $( < $task_file ) ) 
IFS=$old_IFS


# ----------------------------
# --- body of the function ---
# ----------------------------

export exec_time=`date`
echo Starting launching the tasks in file $task_file on $exec_time

# Iteration over all the command lines included in the arguments file
options=""
for ((k=0;k<${number_commands};k+=1)); do
	
	export exec_time=`date`

	# Store the command in a file that is gonna be submitted
	if [ ! -e bash_job ]; then
        echo "${proc_commands[ ( $k ) ]}" > bash_job
    else
        echo "Problem with existent file bash_job!! Please rename the existent file"
        exit
    fi

	# Run the command and save the new job identification in the $job_info variable
	if [ "$interpreter" == "bash" ]; then
		echo "Launching command:      ${proc_commands[ ( $k ) ]} with submit_bash_job"
		job_info=`submit_bash_job $options -N command$k bash_job | awk '{print $3}'`
	else 
		if [ "$interpreter" == "perl" ]; then
			echo "Launching command:      submit_perl_job ${proc_commands[ ( $k ) ]}"
			job_info=`submit_perl_job $options -N command$k bash_job | awk '{print $3}'`
		else
			echo "Launching command:      submit_python_job ${proc_commands[ ( $k ) ]}"
			job_info=`submit_python_job $options -N command$k bash_job | awk '{print $3}'`
		fi
	fi
	echo "on $exec_time"
	
	rm bash_job

	# the new job identification must be added to the $options variable
	if [ $job_info != "" ]; then
		options=`echo -hold_jid $job_info`
	else
		echo "Problem with command ${proc_commands[ ( $k ) ]}. The rest of commands won't continue executing."
		exit
	fi
done

export exec_time=`date`
echo Script finished at $exec_time

