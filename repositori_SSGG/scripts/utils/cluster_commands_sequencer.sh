#!/bin/bash

# ---------------------
# ---- definitions ----
# ---------------------

EXPECTED_ARGS=3 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo 'Name: cluster_commands_sequencer.sh'
	echo 'Description: script for launching simultaneously a list of commands, considering the maximum amount of jobs which can stay (either running or waiting) simultaneously in a cluster queue'
	echo 'Mandatory parameters:'
	echo '  $1 : maximum number of jobs that can simultaneously stay on the cluster queue'
	echo '  $2 : cluster queue to be checked'
	echo '  $3 : file containing the list of commands to run. This file should include ONLY the commands (including submit_job if necessary), that is, without comments and other scrap!'
	exit $E_BADARGS
fi

max_threads=$1 # maximum number of threads allowed to run simultaneously in one node
queue=$2
number_commands=`wc -l $3 | awk '{print $1}'` # number of commands to launch
wait_time=5 # number of seconds to wait between processes status checks

# the file with the command lines will be read. First the field separator environment variable is set to "new line".
old_IFS=$IFS
IFS=$'\n'
proc_commands=( $( < $3 ) ) 
IFS=$old_IFS


# ----------------------------
# --- body of the function ---
# ----------------------------

export exec_time=`date`
echo Starting launching the line commands in file $2 on $exec_time

jobs=0
k=1

# Iteration over all the command lines included in the arguments file
while [ $jobs -lt $max_threads ] || [ $k -le $number_commands ]; do

	# Check amount of jobs of current user after some seconds
	user=`whoami`
	jobs=`qstat -u "$user" -q $queue | wc -l`
	jobs=$(( $jobs -2 ))

	# a new command is launched
	# -------------------------
	if [ $jobs -lt $max_threads ]; then

		export exec_time=`date`
		echo "Launching command:     ${proc_commands[ ( $k - 1 ) ]} & on $exec_time"

		# Command to run:
		${proc_commands[ ( $k - 1 ) ]} &
		
		# Go on with next command
		k=$(( $k + 1 ))
	fi
	
	# Wait $wait_time seconds until launching next job
	sleep $wait_time

	# Check if all jobs were processed
	if [ $k -gt $number_commands ]; then
		jobs=$max_threads
	fi
done
