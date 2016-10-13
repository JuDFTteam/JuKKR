#!/usr/bin/python

"""
Reads file defs.txt and generates Fortran code for reading an inputfile.
Author: Elias Rabel, 2013
"""

import sys

# Some constants
CONFIG_READER_DICT_VALUE_LENGTH = 96
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
! Automatically generated source file. Do not edit manually!
! To add/remove/modify input parameters:
""" + "! Edit " + deffilename + " and run \n! 'inputgenerator.py " + configname + " " + deffilename + " > " + configname + "_mod.F90'" + """
! to generate source code.
!------------------------------------------------------------------------------

""")

print HEADER

print 'module ' + configname + '_mod'
print '  implicit none'
print '  private'
print '  public :: ' + configname + ', load, store, getValues'
print
print

#------------Generate type declaration -------------------
print '  type ' + configname

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    typename = typedict[splitted_line[0]]
    print '    ' + typename + " :: " + splitted_line[1],

    if (splitted_line[0] in vectors):
      print '(' + str(splitted_line[2]) + ')'
    else:
      print

print '  endtype ! ' + configname
print

deffile.close()

deffile = open(deffilename, 'r')

print
print '  interface load'
print '    module procedure read' + configname + 'FromFile'
print '  endinterface'
print
print '  interface store'
print '    module procedure write' + configname + 'ToFile'
print '  endinterface'
print
print '  contains'

#----------- Generate code for retrieving config values -
print '!'+'-'*79
print 'integer function getValues(filename, values) result(ierror)'
print '  use ConfigReader_mod, only: ConfigReader, create, destroy'
print '  use ConfigReader_mod, only: not_found => CONFIG_READER_ERR_VAR_NOT_FOUND'
print '  use ConfigReader_mod, only: use_default => CONFIG_READER_USE_DEFAULT_VALUE'
print '  use ConfigReader_mod, only: getValue, parseFile'

print
print '  character(len=*), intent(in) :: filename'
print
print '  type(' + configname + '), intent(inout) :: values'
print '  type(ConfigReader) :: cr'
print
print """  ierror = 0
  write(*,*) "Reading information from input.conf..."
  call create(cr) ! createConfigReader
#define destroy_and_return   call destroy(cr) ; return
  ierror = parseFile(cr, filename)
  if (ierror /= 0) then
    write(*,*) "Error reading configfile ", trim(filename)
    destroy_and_return
  endif

"""


for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    gettername = 'getValue'
    default_value = None

    if (splitted_line[0] in vectors):
      if len(splitted_line) > 3:
          default_value = splitted_line[3]    
    else:
      if len(splitted_line) > 2:
          default_value = splitted_line[2]             

    if default_value is not None:
        print '  ierror = ' + gettername + '(cr, "' + splitted_line[1] + '", values%' + splitted_line[1],
        print ', def=' + default_value + ')'
        print '  if (ierror == use_default) then'
        print '    write(*,*) "WARNING: Bad/no value given for ' + splitted_line[1] + '. Set ' + splitted_line[1] + ' to ' + default_value + '"' 
        print '    ierror = 0 ! ok, no error'
        print '  elseif (ierror /= 0) then'
    else:
        print '  ierror = ' + gettername + '(cr, "' + splitted_line[1] + '", values%' + splitted_line[1] + ')'
        print '  if (ierror /= 0) then'
    print '    write(*,*) "Bad/no value given for ' + splitted_line[1] + '."'
    print '    destroy_and_return' # program will die
    print '  endif'
    print

deffile.close()

print '  write(*,*) "Finished reading information from input.conf"'
print '  destroy_and_return'
print '#undef destroy_and_return'
print 'endfunction ! get'
print

#----------- Generate code for reading config values from unformatted file -
print '!'+'-'*79
print 'integer function read' + configname + 'FromFile(values, filename) result(ios)'
print '  type(' + configname + '), intent(inout) :: values'
print '  character(len=*), intent(in) :: filename'
print
print '  integer, parameter :: fu = ' + str(FILE_UNIT)
print
print '  open(fu, file=filename, form="unformatted", action="read", status="old", iostat=ios)'
print '  if (ios /= 0) return'

deffile = open(deffilename, 'r')

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    print '  read(fu) values%' + splitted_line[1]

deffile.close()

print '  close(fu)'
print 'endfunction ! load'
print

#----------- Generate code for reading config values from unformatted file -
print '!'+'-'*79
print 'subroutine write' + configname + 'ToFile(values, filename)'
print '  type(' + configname + '), intent(inout) :: values'
print '  character(len=*), intent(in) :: filename'
print
print '  integer, parameter :: fu = ' + str(FILE_UNIT)
print '  integer :: ios'
print
print '  open(fu, file=filename, form="unformatted", action="write", iostat=ios)'
print '  if (ios /= 0) return'

deffile = open(deffilename, 'r')

for line in deffile:
    if line[0] == '#': continue  #skip comment
    splitted_line = line.split()

    print '  write(fu) values%' + splitted_line[1]

deffile.close()

print '  close(fu)'
print 'endsubroutine'
print

print 'endmodule !' + configname

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
