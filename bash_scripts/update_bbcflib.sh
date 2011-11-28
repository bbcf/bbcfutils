#!/bin/sh

################################################################
# Update all bbcflib projects (git pull) at once.
# You must have cloned them once, and they must exist 
# in the same root directory.
# bbcflib, bein, track, bbcfutils
################################################################


if [ -z $1 ]; then
    echo 'git directory where all bbcf git projects are located is missing'
    echo 'USAGE : sh update_bbcflib.sh /path/to/bbcf/projects/home'
    exit 1
fi
CUR=$PWD
BBCFLIB_HOME=$1

# depending on your configuration you can have other names
GIT_REMOTE='origin'
GIT_BRANCH='master'

cd $BBCFLIB_HOME

# update bbcflib
echo '[x] bbcflib [x]'
cd bbcflib;git pull $GIT_REMOTE $GIT_BRANCH

# update bein
echo '[x] bein [x]'
cd ../bein;git pull $GIT_REMOTE $GIT_BRANCH

# update track
echo '[x] track [x]'
cd ../track;git pull $GIT_REMOTE $GIT_BRANCH

# update bbcfutils
echo '[x] bbcfutils [x]'
cd ../bbcfutils;git pull $GIT_REMOTE $GIT_BRANCH


# end 
cd $CUR
exit 0