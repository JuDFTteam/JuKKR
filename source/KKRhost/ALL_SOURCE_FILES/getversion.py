#! /usr/local/bin/python
from subprocess import check_output
version=check_output('git describe',shell=True).strip()
head = ["      character(len=*), parameter :: version='%s'"%version]
file = open('ALL_SOURCE_FILES/version.f90','w')
file.writelines(head)
file.close()
print 'version file written for',version
