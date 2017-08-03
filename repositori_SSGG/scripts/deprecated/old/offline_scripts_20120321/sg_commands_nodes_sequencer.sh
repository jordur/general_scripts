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
	echo '	$1 : maximum number of processes allowed to run simultaneously on one node of the machine'
	echo '	$2 : file containing the list of commands to run. This file should include ONLY the commands, that is, without comments and other scrap!'
	exit $E_BADARGS
fi


commands_file=$2

# Array for processes IDs:
declare -a proc_IDs # identification of the processes currently running
declare -a proc_nodes # identification of the node where a process is currently running on
declare -a proc_nodes_load # array which contains the number of processes currently running/waiting on each node (node load)
declare -a proc_nodes_procs # array which contains the number of commands currently running on each node
declare -a proc_commands # array which contains the different commands to launch

max_threads=$1 # maximum number of threads allowed to run simultaneously in one node
number_nodes=5 # number of nodes in the cluster, frontend node excluded (the first version will work only on the director node, since the passwordless login is not yet active in the different nodes
number_cores=8 # number of processing cores (CPUs) per cluster node
let total_max_threads=${max_threads}*${number_nodes}
number_commands=`wc -l $commands_file | awk '{print $1}'` # number of commands to launch
wait_time=1 # number of seconds to wait between processes status checks
wait_time_load=60 # number of seconds to wait between nodes status checks when they were overload


# the file with the command lines will be read. First the field separator environment variable is set to "new line".
old_IFS=$IFS
IFS=$'\n'
proc_commands=( $( < $commands_file ) ) 
IFS=$old_IFS


# ----------------------------
# --- body of the function ---
# ----------------------------

export exec_time=`date`
echo Starting launching the line commands in file $commands_file on $exec_time

# Each "max_threads" threads of bfast match will be launched simultaneously over the "splitted" files, and then it will be waited for launching a new thread until one of the threads is finished
num_threads=0 #number of threads currently running
threads_to_wait=0 #number of threads to wait for until the loop-for-waiting-for-processes-to-be-finished is abandoned

# initialisation of proc_IDs and proc_nodes vector:
for ((k=0;k<${total_max_threads};k+=1)); do
	proc_IDs[$k]=0
	proc_nodes[$k]=0
done

for ((k=0;k<${number_nodes};k+=1)); do
	proc_nodes_procs[$k]=0
done

# Iteration over all the command lines included in the arguments file

#for ((k=1;k<${number_commands}+1;k+=1)); do
k=1
let aux_var=${number_commands}+1
while [ $k -lt $aux_var ]; do
	
	# a new command is launched
	# -------------------------

	# First it has to be decided the node where the process will run on. The node load will be checked.
	# If all the nodes are overload, then the overload flag will be set (and keep) equal to 1
	overload=1
	for ((l=0;l<${number_nodes};l+=1)); do
		# the "1 minute load" will be consider
		loadavg=`ssh compute-0-$(($l+1)) uptime | awk -F',' '{print $4}' | awk -F' ' '{print $3}'`
		proc_nodes_load[$l]=`echo $loadavg | awk -F \. '{print $1}'`
		echo "node: compute-0-$(($l + 1)) - load: ${proc_nodes_load[$l]}"
		if [ ${proc_nodes_load[l]} -lt $number_cores ] && [ ${proc_nodes_procs[l]} -lt $max_threads ]; then
			overload=0
			#echo 'The system is underload (overload=0)'
		fi
	done
	
	i=0
	while [[ ${proc_nodes_load[$i]} -ge $number_cores || ${proc_nodes_procs[$i]} -ge $max_threads && $i -lt $number_nodes ]]; do
		let i+=1
	done

	if [ $i -lt $number_nodes ]; then
	
		export exec_time=`date`
		let u=$i+1
		echo "On node compute-0-$u - Launching command:     ssh compute-0-$u ${proc_commands[ ( $k - 1 ) ]} &"
		echo "on $exec_time - node load: ${proc_nodes_load[i]}"

		# Command to run:
		export path=`pwd`
		path_existent=0
		#ssh compute-0-$u '[[ ! -d $path ]] && path_existent=1'
		if [ $path_existent -eq 1 ]; then
			echo "Directory not found in node compute-0-$u!!!"
			echo "The last process will be relaunched!"
		else
			ssh compute-0-$u "cd $path; ${proc_commands[ ( $k - 1 ) ]}" &

			# the new process identification must be saved in the proc_IDs vector in the first gap found (position of the vector containing 0)
			aux=$!
			m=0
			while [ ${proc_IDs[m]} -ne 0 ] && [ $m -lt $total_max_threads ]; do
				let m+=1
			done
			if [ $m -lt $total_max_threads ]; then
				proc_IDs[$m]=$aux
				proc_nodes[$m]=$((i+1))
				let proc_nodes_procs[$i]+=1
				#echo "proc_ID: ${proc_IDs[m]} - stored in position $m - running on node ${proc_nodes[m]} with ${proc_nodes_procs[i]} processes in total"
			else
				echo "Problem by launching of new process: no gap found in proc_IDs vector!!"
				echo "The last process ($aux) will be killed and later relaunched!"
				ssh compute-0-$u kill -9 $aux
				let k=$k-1
			fi
			let num_threads+=1

			#for ((z=0;z<$number_nodes;z+=1)); do
			#	echo Node "$(($z + 1)), processes: ${proc_nodes_procs[z]}"
			#done
		
			# Here it's decided for how many threads the script should wait to be finished
			if [ $k -eq $number_commands ]; then
				# when the last command line is being processed, then the script has to wait for continueing until all the active processes are finished
				threads_to_wait=$num_threads 
			else
				if [ $num_threads -lt $total_max_threads ]; then
					# if the maximal number of threads hasn't been reached, then new processes can be launched
					threads_to_wait=0
				else
					# if the maximal number of threads has been reached, then it is needed to wait for one process to be finished
					threads_to_wait=1
				fi
			fi
			echo "Node compute-0-$(($i + 1)) new process: ${proc_IDs[m]} processes running in node: ${proc_nodes_procs[i]} num_treads: $num_threads max_threads_per_node: $max_threads threads_to_wait: $threads_to_wait"
		fi
	fi
	
	# Now the script decides if it has to wait:
	# -------------------------------------------------------
	# ---- loop-for-waiting-for-processes-to-be-finished ----
	# -------------------------------------------------------

	while [[ $threads_to_wait -gt 0 ]]; do
		#echo Loop: num_treads: $num_threads max_threads: $max_threads threads_to_wait: $threads_to_wait
		# Each "wait_time" seconds the programm will check the processes status. If a process is finished, then a new process can be launched
		sleep $wait_time

		#echo "Waiting till end of one of the processes with IDs:"

		num_threads=0
		for ((l=0;l<${total_max_threads};l+=1)); do
			if [ ${proc_IDs[l]} -ne 0 ]; then
				#echo process: ${proc_IDs[l]} \t;
				
				# the kill -0 command just tests the status of the process, not killing it. Its error message is redirected to /dev/null, so that no message will be output in case that the process doesn't exist anymore (2>&1 would redirect error to stdout)
				# Please note that there is a parent process for each process launched in the nodes, so that it isn't needed to check the processes in each node (ssh compute-0-${proc_nodes[l]} kill -0 ${proc_IDs[l]} 2>/dev/null), and it is enough to check in the frontend node:
				kill -0 ${proc_IDs[l]} 2>/dev/null
				if [ $? -eq 1 ]; then
					export exec_time=`date`
					echo "On node compute-0-${proc_nodes[l]} - Process finished: ${proc_IDs[l]} on $exec_time"
					proc_IDs[$l]=0
					u=${proc_nodes[$l]}-1
					let proc_nodes_procs[$u]-=1
					proc_nodes[$l]=0
				else
					let num_threads+=1
				fi
			fi
		done
		if [ $k -eq $number_commands ]; then
			# when the last command line is being processed, then the script has to wait for continueing until all the active processes are finished
			threads_to_wait=$num_threads
		else
			if [ $num_threads -lt $total_max_threads ]; then
				# if the maximal number of threads hasn't been reached, then new processes can be launched
				threads_to_wait=0
			else
				# if the maximal number of threads has been reached, then it is needed to wait for one process to be finished
				threads_to_wait=1
			fi
		fi
	done

	# If all the nodes are overload (overload=1) then it will be waited "wait_time_load" seconds. Otherwise, a new process can be launched
	if [ $overload -ne 0 ]; then
		echo "System overload. Waiting before launching new processes"
		sleep $wait_time_load

		# Obtain the number of processes currently running
		num_threads=0
		for ((l=0;l<${total_max_threads};l+=1)); do
			if [ ${proc_IDs[l]} -ne 0 ]; then
				#echo process: ${proc_IDs[l]} \t;

				# the kill -0 command just tests the status of the process, not killing it. Its error message is redirected to /dev/null, so that no message will be output in case that the process doesn't exist anymore (2>&1 would redirect error to stdout)
				# Please note that there is a parent process for each process launched in the nodes, so that it isn't needed to check the processes in each node (ssh compute-0-${proc_nodes[l]} kill -0 ${proc_IDs[l]} 2>/dev/null), and it is enough to check in the frontend node:
				kill -0 ${proc_IDs[l]} 2>/dev/null
				if [ $? -eq 1 ]; then
					export exec_time=`date`
					echo "On node compute-0-${proc_nodes[l]} - Process finished: ${proc_IDs[l]} on $exec_time"
					proc_IDs[$l]=0
					u=${proc_nodes[$l]}-1
					let proc_nodes_procs[$u]-=1
					proc_nodes[$l]=0
				else
					let num_threads+=1
				fi
			fi
		done
		if [ $k -eq $number_commands ]; then
			# when the last command line is being processed, then the script has to wait for continueing until all the active processes are finished
			threads_to_wait=$num_threads
		else
			if [ $num_threads -lt $total_max_threads ]; then
				# if the maximal number of threads hasn't been reached, then new processes can be launched
				threads_to_wait=0
			else
				# if the maximal number of threads has been reached, then it is needed to wait for one process to be finished
				threads_to_wait=1
			fi
		fi
	else
		let k+=1
	fi
done

