#!/bin/bash

#in case of SGNET user (e.g. SGNET\arbol), get only username (e.g. arbol):
read user <<< $USER
user=`echo $user | sed 's/SGNET//g'`

if [ "$1" != "" ]; then
	branch=`git branch | grep "*" | awk '{print $2}'`
	if [ "$branch" != "dev-$user" ]; then
		echo "WARNING: Currently not working in develop branch (dev-$user)!!!!"
		echo "Move to the develop branch (git checkout dev-$user) and be sure that you already saved your changes before submitting them"
		exit
	fi
	chmod -R 775 *
	git pull origin dev-$user
	git add .
	git commit -m "$@"
	git push origin dev-$user
	git checkout master
	git pull origin master
	git merge dev-$user
	git add .
	git commit -m "$@"
	git push origin master
	git checkout dev-$user
	# Finally get other master changes in the development branch
	git pull origin master
else
	echo "ERROR: No changes committed to repository!! Please be sure to add a descriptive comment for your commitment"
fi
