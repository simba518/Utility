#! /usr/bin/env python

import os
import sys

parameter_is_correct = True

rootdir = "./Bin/"
if len(sys.argv)==2:
    if sys.argv[1] == "-d":
        rootdir = rootdir + "Debug"
        print "Run Debug Test Cases\n"
    elif sys.argv[1] == "-r":
        rootdir = rootdir + "Release"
        print "Run Release Test Cases\n"
    else :
        print "error: the parameters for this python script is uncorrect"
        print "usage: runalltest [-r|-d]\n"
        parameter_is_correct = False

if parameter_is_correct:
    fileList = []
    logifle = "./Doc/test.log"
    for root, subFolders, files in os.walk(rootdir):
        for file in files:
            fileList.append(os.path.join(root,file))
    for files in fileList:
        if (files[-4:]).lower() == 'test' :
            print "TestCase Name: " +  files
            os.system(files)
