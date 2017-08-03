#!/bin/bash
if [ "$1" == "" ]; then
	directory="./"
else
	if [[ "$1" =~ ^[0-9]+$ ]] ; then
		directory="$PROJECTS/BF$1*/*/logs"
	else
		directory=$1
	fi
fi
echo "Checking for errors in $directory..."
echo
grep -Ri error $directory
grep -Ri fatal $directory
grep -Ri fault $directory
