git checkout develop
git add *
git commit -m "$@"
git pull bitbucket develop
git push bitbucket develop
