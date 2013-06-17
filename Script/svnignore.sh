#!/bin/sh

# cd ~/Workspace/PathControl/
ignore_file=./Script/svnignore
svn propset svn:ignore -F $ignore_file . -R

svn propset svn:ignore '*' ./Build/Release
svn propset svn:ignore '*' ./Build/Debug

svn propset svn:ignore '*' ./Lib/Release
svn propset svn:ignore '*' ./Lib/Debug

svn propset svn:ignore '*' ./Bin/Release
svn propset svn:ignore '*' ./Bin/Debug