git checkout master
git merge develop
git add *
git commit -m "$@"
git pull bitbucket master
git push bitbucket master
git checkout develop
