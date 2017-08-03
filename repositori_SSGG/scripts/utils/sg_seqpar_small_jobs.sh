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

        echo "$scriptname is a script for launching in parallel a list of commands over the queue/resources manager."
        echo "Usage:"
        echo "  $scriptname <options> input1"
        echo
        echo "* Input (mandatory parameters): *"
        echo " input1     file containing the list of commands/jobs/tasks to launch in parallel"
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
while getopts ":n" Option; do
	case $Option in
		t) interpreter="$OPTARG";; 
		w) wait_time="$OPTARG";;
		*) echo "Bad option: $Option! Call $scriptname with no parameters for help."; exit;;
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
	if [ ! -e bash_sequence_job ]; then
        	echo "${proc_commands[ ( $k ) ]}" > bash_sequence_job
	else
		echo "Problem with existent file bash_sequence_job.sh!! Please rename the existent file"
		exit
	fi

	# Run the sequence commands and save the last job identification in the $job_info variable
	echo "Launching sequence commands: ${proc_commands[($k)]}"
	job_info=`../sg_sequence_small_jobs.sh ${proc_commands[ ( $k ) ]} | awk '{if($1=="Last_job_ID:") print $2}'`
	echo "on $exec_time"
	
	rm bash_sequence_job
	
	# the new job identification must be added to the $options variable
	echo $job_info
	if [ $job_info != "" ]; then
		options=`echo $options -hold_jid $job_info`
	else
		echo "Problem with commands in ${proc_commands[ ( $k ) ]}. The rest of commands will continue executing."
	fi
done

# A watcher-job (sleep 0) will be submitted. The watcher-job will finish only when all the parallel jobs are finished:
echo sleep 0 > watcher-job
job_id=`submit_small_bash_job $options -N watcher-job watcher-job | awk '{print $3}'`

# Waiting until watcher-job finishes:
while [ "`qstat | awk -v val=$job_id '{ if ($1==val) print 1 }'`" == "1" ]; do
	sleep $wait_time
done

rm watcher-job

export exec_time=`date`
echo Script finished at $exec_time

