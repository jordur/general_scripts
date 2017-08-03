#!/bin/bash
#####################################################
# Copyright 2013 - Sistemas Genomicos S.L.
# @Desc: Script for updating repositories in the production cluster
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

	echo "$scriptname is a script for updating repositories in the production cluster, in both testing or production environments."
	echo "Usage:"
	echo "  $scriptname environment repository"
	echo
	echo "* Input (mandatory parameters): *"
	echo " environment  name of the environment to update (testing, production, both)"
	echo " repository   name of the repository to update (scripts, pipeta, biolib, panels, all)"
	echo
	echo "* Examples of running: *"
	echo " \$$scriptname testing biolib"
	exit
}

function UpdateRepository ()
#
{
	environment=$1
	repository=$2
	# User name is needed for having different testing environments for each user
	read user <<< $USER
	if [ "$environment" == "both" ]; then
		for env in ${environments[@]}; do
			UpdateRepository $env $repository
		done
	else
		if [ "$repository" == "all" ]; then
			for rep in ${repositories[@]}; do
				UpdateRepository $environment $rep
			done
		else
			# validate environment
			env_ok=0
			for env in ${environments[@]}; do
				if [ "$environment" == "$env" ]; then
					env_ok=1
				fi
			done
			if [ $env_ok -eq 0 ]; then
				echo "ERROR: Wrong environment name. Call $scriptname with no parameters for help."
				exit
			fi
			if [ "$environment" == "testing" ]; then
				envdir="testing_$user"
				branch="dev-$user"
			else
				envdir="scripts"
				branch="master"
			fi
			
			# validate repository
			rep_ok=0
			for rep in ${repositories[@]}; do
				if [ "$repository" == "$rep" ]; then
					rep_ok=1
				fi
			done
			if [ $rep_ok -eq 0 ]; then
				echo "ERROR: Wrong repository name. Call $scriptname with no parameters for help."
				exit
			fi
			if [ "$repository" == "scripts" ]; then
				path="/share/apps/$envdir"
				
				# If the testing folder doesn't yet exist, it will be created
				if [ ! -d $path ]; then
					git clone git@beefeater.sgnet.local:bioinfo/scripts.git $path
					cd $path
					git submodule init
					git submodule update
				fi
			fi
			if [ "$repository" == "biolib" ]; then
				path="/share/apps/$envdir/biolib"
			fi
			if [ "$repository" == "panels" ]; then
				path="/share/apps/$envdir/genomics/panels/panel_design"
			fi
			if [ "$repository" == "pipeta" ]; then
				path="/share/apps/$envdir/utils/PipeTA"
			fi
			
			# Update repository
			cd $path
	
			# Fetch repository information if not yet found 
			if [ "$environment" == "testing" ]; then
				found=`git branch | grep dev-$user`
				if [ "$found" == "" ]; then
					git fetch
				fi
			fi
			
			git checkout $branch
			git stash
			git pull origin $branch
			git fetch --tags
		fi
	fi
}

# **************************************************************************************************
# **************************************************************************************************

# -----------------------
# -- global parameters --
# -----------------------
OPTERROR=65
environments=( production testing )
repositories=( scripts pipeta biolib panels)


# ------------------------
# -- default parameters --
# ------------------------
environment=""
repository=""


# -----------------------------
# -- Name, version and usage --
# -----------------------------
scriptname=`basename $0`
if [ $# -lt 2 ]; then
	usage $scriptname
	exit $OPTERROR
fi

# ------------------------------------------
# -- Setting input options and parameters --
# ------------------------------------------
# When you want getopts to expect an argument for an option, just place a : (colon) after the proper option flag 
args=($@)
let index=0 #'OPTIND-1'
environment=${args[index]}
repository=${args[index+1]}

if [ "$environment" == "" ] || [ "$repository" == "" ]; then
	echo "ERROR: Not enough parameters in call to $scriptname!!"
	echo "Call $scriptname with no parameters for help."
	exit
fi


# ---------------
# -- main body -- 
# ---------------

export exec_time=`date`
echo Script started at $exec_time

UpdateRepository $environment $repository

export exec_time=`date`
echo Script finished at $exec_time
