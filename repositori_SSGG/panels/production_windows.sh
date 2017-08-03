git checkout master
git merge develop
git add *
git commit -m "$@"
git pull origin master
git push origin master
git checkout develop
