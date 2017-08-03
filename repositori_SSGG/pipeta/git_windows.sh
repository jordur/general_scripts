#!/bin/bash

#in case of SGNET user (e.g. SGNET\arbol), get only username (e.g. arbol):
read user <<< $USER
user=`echo $user | sed 's/SGNET//g'`

if [ "$1" != "" ]; then
	branch=`git branch | grep "*" | awk '{print $2}'`
	if [ "$branch" != "dev-$user" ]; then
		read -p "WARNING: Your develop branch (dev-$user) doesn't yet exist. Do you want me to create it and include all your current branch changes in it (y/n)? Keep in mind that all not yet committed changes will be ONLY included in your develop branch, BUT they can later be committed to the production branch. So, yes (y) or no (n)??" dev
		case $dev in
			[Nn]* ) exit;;
			[Yy]* ) git branch dev-$user; git checkout dev-$user; git add *; git commit -m "Creation of branch dev-$user"; git push origin dev-$user; break;;
			* ) echo "Please answer yes or no"
		esac
	fi
	chmod -R 775 *
	git pull origin dev-$user
	git add *
	git commit -m "$@"
	git push origin dev-$user
else
	echo "ERROR: No changes committed to repository!! Please be sure to add a descriptive comment for your commitment"
fi
