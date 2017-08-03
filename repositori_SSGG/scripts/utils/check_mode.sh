#!/bin/bash
echo "For changing working mode, use \"source testing\" or \"source production\""
echo
repo=`which pipeta.pl | awk -F"/" '{print $4}'`
if [ $repo == "scripts" ]; then 
	echo "INFO: currently on PRODUCTION MODE";
else 
	echo "INFO: currently on TESTING MODE"; 
fi
