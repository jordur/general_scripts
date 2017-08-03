#!/bin/bash
folder=$1
if [[ "$folder" != "" ]]; then
	job_list=`qstat | awk '{print $1}' | grep -v job-ID | grep -v '\-\-\-\-'`
	csv_list=""
	echo "List of jobs that are going to be deleted:"
	for id in $job_list; do
		script_file=`qstat -j $id | grep script_file | awk '{print $2}'`
		if [[ ${script_file:0:1} == "/" ]]; then
			parameter="script_file"
		else
			parameter="sge_o_workdir"
		fi
		count=`qstat -j $id | grep $parameter | grep -c $folder`;
		if [[ $count > 0 ]]; then
			echo qdel $id;
			if [[ "$csv_list" == "" ]]; then
				csv_list="$id"
			else
				csv_list="$csv_list,$id"
			fi
		fi
	done
	while true; do
		read -p "Are you sure you want to delete the jobs? " yn
		case $yn in
			[Yy]* ) echo qdel $csv_list; qdel $csv_list; break;;
			[Nn]* ) exit;;
			* ) echo "Please answer yes or no.";;
		esac
	done
fi
