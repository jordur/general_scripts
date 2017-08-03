#!/bin/bash

pipeline ()
# sequence of actions to carry out on the files
{
	echo $2 $1
	$2 $1
	return
}

recurse ()
# function to recursively access all the files in a directory path
{
	#echo $2
	for file in ${1}/*; do
		if [ -d "$file" ]; then
				pipeline "$file" "$2"
        		recurse "$file" "$2"
		else
			if [ -f "$file" ]; then	
	        		pipeline "$file" "$2"
			fi
		fi
	done
	return
}

# main body

# ---------------
# - Definitions -
# ---------------

EXPECTED_ARGS=2 # Number of arguments expected in the script call
E_BADARGS=65    # Error code in case of bad arguments


# It will be checked for the proper number of arguments
if [ ${#} -ne $EXPECTED_ARGS ]; then
        echo 'Name: sg_recurse.sh'
        echo 'Description: Pipeline for carrying out an action recursively over all the files in a given directory.'
	echo '		Thought it may seem an unnecessary script, since many actions can be carried out recursively with "awk" and "ls -R", this is not the case for some actions. Just imagine how would you carry out a "chmod 775" recursively over all the files in a directory, since you cannot use "ls -R directory | chmod 775"'
	echo 'Mandatory parameters:'
        echo '       $1 : target directory'
        echo '       $2 : action to carry out over the files in the target directory. Please be sure to include the action between qutotation marks (simple or double) when it contains more than one word. It can also be a sequences of actions inside an executable script'
        exit $E_BADARGS
fi

#if [ ${*:${#*}-1} = "/" ]; then
#	recurse ${*:0:${#*}-2}
#else
	recurse $1 "$2"
#fi
