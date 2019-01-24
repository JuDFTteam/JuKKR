#!/usr/bin/env python

import os, sys

# default path
path = '../source/common/'

# read in path from first given command line argument
if len(sys.argv)>1:
    path = sys.argv[1]
if path[-1]!='/': path+='/'

# if a second command line argument is given then the occurences of 'use' without only statements are printed
if len(sys.argv)>2:
    print_occurences = True
else:
    print_occurences = False

missing_all = 0
# walk through all subdirectories
for dirpath, subdirs, files in os.walk(path):
    if files != []: # do only something if a file was found
        files = [f for f in files if '.f' in f or '.F' in f] # sort out only fortran source files
        for f in files:
            with open(dirpath+'/'+f) as ifile:
                txt = ifile.readlines() # read file 
                i_noonly = 0
                for line in txt: # loop over all lines to find 'use' statements
                    tmpline = line.lower() # remove case sensitivity
                    if ('use' in tmpline and 'only' not in tmpline # look for use without only
                        and tmpline.split()[0]=='use' # make sure use is the first statement (otherwise it appears as a comment) 
                        and 'mpi' not in tmpline # exclude 'use mpi'
                        and 'omp_lib' not in tmpline # exclude 'use omp_lib'
                       ):
                        i_noonly+=1 # count number of missing only statements
                        if print_occurences: print line # print line if option is chosen
                # print result for each (sub) directory
                if i_noonly>0:
                    print (dirpath+'/'+f).replace('//','/'), i_noonly
                    missing_all+=i_noonly

print path, 'total', missing_all

