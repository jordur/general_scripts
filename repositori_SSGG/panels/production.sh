#!/bin/bash
if [ "$1" != "" ]; then
	git checkout master
	git merge develop
	git add *
	git commit -m "$@"
	git pull origin master
	git push origin master
	git checkout develop
else
	echo "ERROR: No changes committed to repository!! Please be sure to add a descriptive comment for your commitment"
fi