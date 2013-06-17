#! /usr/bin/env python

###############################################################################
# Add all of the new files to repository, and remove all of the missed files.
# And if parameter is provided, the parameter will be use as svn's message 
# and a commit command will be performed.
# usage: 
# 1) ./Script/svn_add_del_commit.py
# 2) ./Script/svn_add_del_commit.py "message for commiting is here"
###############################################################################

import os
import sys

svn_status_cmd = "svn status "
svn_add_cmd = "svn add "
svn_del_cmd = "svn rm "
svn_commit_cmd = "svn commit -m "

svn_status = os.popen(svn_status_cmd)
svn_status = svn_status.read()

print "the original svn status***********************"
print svn_status

svn_status = svn_status.split("\n");
print "add new files*********************************"
for oneline in svn_status:
    if len(oneline) > 2 and oneline[0] == '?':
            os.system(svn_add_cmd + oneline[1:])

print "remove missing files**************************"
for oneline in svn_status:
    if len(oneline) > 2 and oneline[0] == '!':
            os.system(svn_del_cmd + oneline[1:])

if len(sys.argv)==2:
    os.system( svn_commit_cmd + "\""+ sys.argv[1]+ "\"" )
