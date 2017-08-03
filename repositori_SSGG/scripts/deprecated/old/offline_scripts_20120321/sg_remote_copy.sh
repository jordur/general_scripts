#!/bin/bash

#########################################
## Bash Script to find files on remote ##
## server and copy them to your local  ##
## machine, using find and scp         ##
## Licenced under the GNU GPLv3 licence##
## Author: Guillermo Garron            ##
## Date: 2009-04-10                    ##
#########################################

usage()
{
echo "usage copy-remote.sh [user@server] [path/to/look/in] ['file-pattern-to-look-for'] [output-directory]"
echo "Remeber to use the ('') in the pattern"
}

if [ $# != 4 ] 
then
    usage
    exit 1
fi
ssh $1 "find $2 -iname '$3'" > ./file-list.txt
while read filename; do
	if [ ! -d $4 ]; then
		mkdir -p $4
	fi
        scp $1:$filename $4/;
done < ./file-list.txt
