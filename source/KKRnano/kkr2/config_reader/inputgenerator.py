#!/usr/bin/python

"""
Reads file defs.txt and generates Fortran code for reading an inputfile.
Author: Elias Rabel, 2013
"""

import sys

# Some constants
CONFIG_READER_DICT_VALUE_LENGTH = 192
FILE_UNIT = 67

configname = sys.argv[1]
deffilename = sys.argv[2]
deffile = open(deffilename, 'r')


typedict = {'i':'integer','d':'double precision',
            's':'character(len='+str(CONFIG_READER_DICT_VALUE_LENGTH)+')',
            'l':'logical', 'dv':'double precision', 'iv':'integer'}

getterNamesdict = {'i':'getValueInteger','d':'getValueDouble',
                   's':'getValueString','l':'getValueLogical',
                   'dv':'getValueDoubleVector', 'iv':'getValueIntVector'}
                   
vectors = ('dv', 'iv')


HEADER = ("""
!------------------------------------------------------------------------------
! Automatically generated source file. Do not edit by hand.
! To add/remove/modify input parameters:
""" +
"! Edit " +  deffilename + " and run \n! 'inputgenerator.py " + configname + " " + deffilename + " > source.f90'" +
"""
! to generate source code.
!------------------------------------------------------------------------------

""")

print HEADER

print 'module ' + configname + '_mod'

#------------Generate type declaration -------------------
print 'type ' + configname

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    typename = typedict[splitted_line[0]]
    print '  ' + typename + " :: " + splitted_line[1],

    if (splitted_line[0] in vectors):
      print '(' + str(splitted_line[2]) + ')'
    else:
      print

print 'end type ' + configname

print

deffile.close()

deffile = open(deffilename, 'r')

print 'CONTAINS'

#----------- Generate code for retrieving config values -
print '!'+'-'*79
print 'integer function get' + configname + 'Values(filename, confvalues) result(ierror)'

print '  use Config_Reader'
print '  implicit none'
print
print '  character(len=*) :: filename'
print '  type (ConfigReader) :: conf'
print '  type (' + configname + '), intent(inout) :: confvalues'
print
print """  ierror = 0
  call createConfigReader(conf)
  call parseFile(conf, filename, ierror)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile", filename
    call destroyConfigReader(conf)
    return
  end if

"""


for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    gettername = getterNamesdict[splitted_line[0]]

    if (splitted_line[0] in vectors):
      print '  call ' + gettername + \
            '(conf, "' + splitted_line[1] + '", confvalues%' + splitted_line[1] + \
            ', ' + splitted_line[2] + ', ierror)'      
    else:
      print '  call ' + gettername + \
            '(conf, "' + splitted_line[1] + '", confvalues%' + splitted_line[1] + \
            ', ierror)'

    print '  if (ierror /= 0) then'
    print '    write(*,*) "Bad/no value given for ' + splitted_line[1] + '."'
    print '    call destroyConfigReader(conf)'
    print '    return'
    print '  end if'

deffile.close()

print '  call destroyConfigReader(conf)'
print 'end function'
print

#----------- Generate code for reading config values from unformatted file -
print '!'+'-'*79
print 'integer function read' + configname + 'FromFile(filename, confvalues) result(ierror)'
print '  implicit none'
print '  character(len=*), intent(in) :: filename'
print '  type (' + configname + '), intent(inout) :: confvalues'
print
print '  integer, parameter :: FILEHANDLE = ' + str(FILE_UNIT)
print
print '  ierror = 0'
print '  open(FILEHANDLE, file=filename, form="unformatted")'

deffile = open(deffilename, 'r')

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    print '  read(FILEHANDLE) confvalues%' + splitted_line[1]

deffile.close()

print '  close(FILEHANDLE)'
print 'end function'
print

#----------- Generate code for reading config values from unformatted file -
print '!'+'-'*79
print 'integer function write' + configname + 'ToFile(filename, confvalues) result(ierror)'
print '  implicit none'
print '  character(len=*), intent(in) :: filename'
print '  type (' + configname + '), intent(inout) :: confvalues'
print
print '  integer, parameter :: FILEHANDLE = ' + str(FILE_UNIT)
print
print '  ierror = 0'
print '  open(FILEHANDLE, file=filename, form="unformatted")'

deffile = open(deffilename, 'r')

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    print '  write(FILEHANDLE) confvalues%' + splitted_line[1]

deffile.close()

print '  close(FILEHANDLE)'
print 'end function'
print

print 'end module'

#------------ Assignment statements ------------------------------
#deffile = open(deffilename, 'r')
#
#for line in deffile:
#    if line[0] == '#': continue  #skip comment
#    splitted_line = line.split()
#
#    print '! '+ splitted_line[1] + ' = cf%' + splitted_line[1]
#
#deffile.close()
