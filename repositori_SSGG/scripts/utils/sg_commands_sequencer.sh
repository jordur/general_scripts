#!/bin/bash

# ---------------------
# ---- definitions ----
# ---------------------

EXPECTED_ARGS=2 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments

# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
	echo 'Name: sg_commands_sequencer.sh'
	echo 'Description: script for launching simultaneously a list of commands, considering the maximum amount of processes which can run simultaneously in the machine'
	echo 'Mandatory parameters:'
	echo '	$1 : maximum number of processes to run simultaneously on the machine'
	echo '	$2 : file containing the list of commands to run. This file should include ONLY the commands, that is, without comments and other scrap!'
	exit $E_BADARGS
fi

# Array for processes IDs:
declare -a proc_IDs # identification of the processes currently running
declare -a proc_nodes # identification of the node where a process is currently running on
declare -a proc_commands # array which contains the different commands to launch

max_threads=$1 # maximum number of threads allowed to run simultaneously in one node
number_nodes=1 # number of nodes in the cluster (the first version will work only on the director node, since the passwordless login is not yet active in the different nodes
number_commands=`wc -l $2 | awk '{print $1}'` # number of commands to launch
wait_time=1 # number of seconds to wait between processes status checks

# the file with the command lines will be read. First the field separator environment variable is set to "new line".
old_IFS=$IFS
IFS=$'\n'
proc_commands=( $( < $2 ) ) 
IFS=$old_IFS


# ----------------------------
# --- body of the function ---
# ----------------------------

export exec_time=`date`
echo Starting launching the line commands in file $2 on $exec_time

# Each "max_threads" threads of bfast match will be launched simultaneously over the "splitted" files, and then it will be waited for launching a new thread until one of the threads is finished
num_threads=0 #number of threads currently running
threads_to_wait=0 #number of threads to wait for until the loop-for-waiting-for-processes-to-be-finished is abandoned

# initialisation of proc_IDs and proc_nodes vector:
for ((k=0;k<${max_threads};k+=1)); do
	proc_IDs[$k]=0;
done

for ((k=0;k<${number_nodes};k+=1)); do
	proc_nodes[$k]=0;
done

# Iteration over all the command lines included in the arguments file

for ((k=1;k<${number_commands}+1;k+=1)); do
	
	# a new command is launched
	# -------------------------

	# First it has to be decided the node where the process will run on
	#i=0
	#while [ ( ! ${proc_nodes[i]} -eq 0 ) && ( $i -lt $number_nodes ) ]; do
	#	let i+=1
	#done
	#if [ $i -lt $number_nodes ]; then
	#	ssh compute-0-${i+1}
	#	echo ${proc_commands[ ( $k - 1 ) ]} &
	#	${proc_commands[ ( $k - 1 ) ]} &
	#	proc_IDs[$num_threads]=$!
	#	let num_threads+=1
	#fi

	export exec_time=`date`
	echo "Launching command:     ${proc_commands[ ( $k - 1 ) ]} &"
	echo "on $exec_time"

	# Command to run:
	${proc_commands[ ( $k - 1 ) ]} &

	# the new process identification must be saved in the proc_IDs vector in the first gap found (position of the vector containing 0)
	aux=$!
	m=0
	while [ ${proc_IDs[${m}]} -ne 0 ] && [ $m -lt $max_threads ]; do
		let m+=1
	done
	if [ $m -lt $max_threads ]; then
		proc_IDs[$m]=$aux
	else
		echo "Problem by launching of new process: no gap found in proc_IDs vector!!"
		echo "The last process ($aux) will be killed and later relaunched!"
		kill -9 $aux
		let k=$k-1
	fi
	let num_threads+=1
	
	# Here it's decided for how many threads the script should wait to be finished
	if [ $k -eq $number_commands ]; then
		# when the last command line is being processed, then the script has to wait for continueing until all the active processes are finished
		threads_to_wait=$num_threads 
	else
		if [ $num_threads -lt $max_threads ]; then
			# if the maximal number of threads hasn't been reached, then new processes can be launched
			threads_to_wait=0
		else
			# if the maximal number of threads has been reached, then it is needed to wait for one process to be finished
			threads_to_wait=1
		fi
	fi
	echo "New process: ${proc_IDs[m]} num_treads: $num_threads max_threads: $max_threads threads_to_wait: $threads_to_wait"
	
	# Now the script decides if it has to wait:
	# -------------------------------------------------------
	# ---- loop-for-waiting-for-processes-to-be-finished ----
	# -------------------------------------------------------

	while [ $threads_to_wait -gt 0 ]; do
		#echo Loop: num_treads: $num_threads max_threads: $max_threads threads_to_wait: $threads_to_wait
		# Each "wait_time" seconds the programm will check the processes status. If a process is finished, then a new process can be launched
		sleep $wait_time

		#echo "Waiting till end of one of the processes with IDs:"

		num_threads=0
		for ((l=0;l<${max_threads};l+=1)); do
			if [ ${proc_IDs[l]} -ne 0 ]; then
				#echo process: ${proc_IDs[l]} \t;
				
				# the kill -0 command just tests the status of the process, not killing it. Its error message is redirected to /dev/null, so that no message will be output in case that the process doesn't exist anymore (2>&1 would redirect error to stdout)
				kill -0 ${proc_IDs[l]} 2>/dev/null
				if [ $? -eq 1 ]; then
					echo "Process finished: ${proc_IDs[l]}"
					proc_IDs[$l]=0
				else
					let num_threads+=1
				fi
			fi
		done
		if [ $k -eq $number_commands ]; then
			# when the last command line is being processed, then the script has to wait for continueing until all the active processes are finished
			threads_to_wait=$num_threads
		else
			if [ $num_threads -lt $max_threads ]; then
				# if the maximal number of threads hasn't been reached, then new processes can be launched
				threads_to_wait=0
			else
				# if the maximal number of threads has been reached, then it is needed to wait for one process to be finished
				threads_to_wait=1
			fi
		fi
	done
done

