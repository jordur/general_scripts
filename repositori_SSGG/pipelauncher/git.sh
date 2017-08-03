#!/bin/bash
if [ "$1" != "" ]; then
	git checkout develop
	chmod -R 775
	git add *
	git commit -m "$@"
	git pull origin develop
	git push origin develop
else
	echo "ERROR: No changes committed to repository!! Please be sure to add a descriptive comment for your commitment"
fi
