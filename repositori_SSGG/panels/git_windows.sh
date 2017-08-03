git checkout develop
git add *
git commit -m "$@"
git pull origin develop
git push origin develop
