#!/usr/bin/env python3

####################################################
### Installation script for the Juelich KKR code ###
####################################################

# import modules
import os
import sys
import getopt

##########################################################################
 
def greeting():
   """ Prints greeting. """
   print("***********************************************************")
   print("Welcome to the Installation script of the Juelich KKR code.")
   print("\nYou can find useful additional information in our KKR wiki: \033[4mhttps://kkr.iff.kfa-juelich.de/doku.php\033[0m")
   print("and in our code's documentation: \033[4mhttps://kkr.iffgit.fz-juelich.de/kkrjm/index.html\033[0m")
   print("***********************************************************\n")

##########################################################################
 
def goodbye():
   """ Prints final remarks before exiting. """
   print("\n***********************************************************")
   print("Congratulations you have now set up the build directory.")
   print("To compile the code please go to the build directory and type")
   print("    'make'")
   print(" or")
   print("    'make -jn'")
   print(" where n is the number of tasks you want to use in the parallel compilation process.")
   print("***********************************************************")

##########################################################################

def usage():
   """ Prints usage info and exists. """
   print("Acceptable arguments:")
   print("  -i or --interactive        Use interactive installation script asking the user for input (default if none are given).")
   print("  -d or --debug              Set up debug mode of the code.")
   print("  -h or --help               Print this message and exit.")
   print("  --machine=name             Use a predefined set of settings for a specific machine where 'name' is one of 'iff', 'claix', 'jureca'.")
   print("  --compiler=name            Use a specific compiler, 'name' could for example be 'mpiifort' or 'gfortran'.")
   print("  --parallelization=scheme   Use either MPI, OpenMP or both (hybrid) parallelization: 'scheme should be one of 'mpi', 'openmp', 'hybrid'.")
   print("  --flags=flag1,flag2        Add cmake flags manually (can be combined with either -m or the settings with -c and -p).")
   sys.exit()

##########################################################################

def read_interactive(flags):
   """ Interactively ask user for input and return the neede options. """

   print("Please input the compiler name (e.g. 'gfortran', 'ifort', 'mpiifort'). Empty input will try the system's default.")
   compiler = input()
   if compiler=='': compiler = None

   print("Please input the parallelization scheme ('serial', 'mpi', 'openmp', 'hybrid')")
   while True:
      parallelization = input()
      if parallelization not in ['serial', 'mpi', 'openmp', 'hybrid']:
         print("your input is not a valid parallelization mode. Please chose a valid scheme.")
      else:
         break

   print("Please input additional flags for cmake (empty input skips adding additional flags)")
   while True:
      print("Next flag? (empty line exists read-in)")
      tmpinput = input()
      if tmpinput=='':
         break
      else:
         if 'enable_mpi' in tmpinput.lower() or 'enable_omp' in tmpinput.lower():
            print("Skipping flag {} since it messes with the parallelization settings.")
         else:
            flags.append(tmpinput)

   print("Use debug flags for compilation? (leave empty for 'Release' build)")
   debug = input()
   if debug != '':
      flags.append("CMAKE_BUILD_TYPE=Debug")

   print("Summary of inputs:\n----------\n Compiler: {}\n Parallelization scheme: {}\n Cmake flags: {}\n".format(compiler, parallelization, flags))

   return compiler, parallelization, flags

##########################################################################

def read_machine(flags, mname):
   """ Read options from machine file. """
   print("read_machine: Not implemented yet.")
   sys.exit()
   return compiler, parallelization, flags

##########################################################################

def check_dependencies():
   """ Check if all necessary dependencies are available (e.g. cmake). """
   print("check_dependecies: Not implemented yet.")
   pass

##########################################################################

def create_build_dir():
   """ Create a build directory and store eventually previously existing ones. """
   import shutil
   # check for existence of build dir
   if 'build' in os.listdir('.'):
      i = 1
      oldbuild = 'build_{}'.format(i)
      while oldbuild in os.listdir('.'):
         i += 1
         oldbuild = 'build_{}'.format(i)
      print("Found old build directory. Moving this to {}".format(oldbuild))
      shutil.move('build', oldbuild)

   # create build dir
   os.makedirs('build')

##########################################################################

def run_cmake(compiler, parallelization, flags, verbose):
   """ Runs cmake step to create makefile etc in build directory. """
   from subprocess import call

   # start creating cmake task
   task = "cd build && "

   # add compiler if given
   if compiler is not None:
      task += "FC=$(which {}) ".format(compiler)

   # add cmake command
   task += "cmake  "

   # add parallelization flags
   if parallelization=='hybrid':
      task += "-DENABLE_MPI=ON  "
      task += "-DENABLE_OMP=ON  "
   elif parallelization=='openmp':
      task += "-DENABLE_MPI=OFF  "
      task += "-DENABLE_OMP=ON  "
   elif parallelization=='mpi':
      task += "-DENABLE_MPI=ON  "
      task += "-DENABLE_OMP=OFF  "
   else:
      task += "-DENABLE_MPI=OFF  "
      task += "-DENABLE_OMP=OFF  "

   # add additional flags if given
   for flag in flags:
      task += "-D{} ".format(flag)

   # finalize task and run
   task += '..'

   print("\nNow run cmake command: '{}'".format(task))
   print("\n***********************************************************\n")

   call(task, shell=True)

##########################################################################

def main(argv):
   # print greeting lines
   greeting()
   # process script options
   if 1:#try:
      if len(argv)==0: usage()
      opts, args = getopt.getopt(argv, "ivhdmcpf:", ["interactive", "verbose", "help", "debug", "machine=", "compiler=", "parallelization=", "flags="])
      print(argv, len(argv), opts,args)
   else:#except getopt.GetoptError:
      usage()

   # define defaults
   compiler = None
   parallelization = None
   flags = []
   verbose = False

   # now deal with different input options:
   for opt, arg in opts:
      if opt in ("-v", "--verbose"):
         verbose = True

      if opt in ("-h", "--help"):
         usage()

      elif opt in ("-i", "--interactive"):
         compiler, parallelization, flags = read_interactive(flags)

      elif opt in ("-m", "--machine"):
         compiler, parallelization, flags = read_machine(flags, args)

      else: # read in options from compiler and parallelization info
         if opt in ("-c", "--compiler"):
            compiler = arg
            if compiler=='': compiler = None
         
         if opt in ("-p", "--parallelization"):
            parallelization = arg
            if parallelization=='': parallelization = None
         
      # finally allow to add flags additionally
      if opt in ("-f", "--flags") and "-i" not in opts and "--interactive" not in opts:
         for iarg in arg.split(':'):
             flags.append(iarg)

      if opt in ("-d", "--debug"):
          flags.append("CMAKE_BUILD_TYPE=Debug")

   # check dependencies
   create_build_dir()

   # create build dir
   check_dependencies()

   # run cmake with previously read-in options
   run_cmake(compiler, parallelization, flags, verbose)

   # print final messages
   goodbye()

##########################################################################

# run script
if __name__ == "__main__":
   main(sys.argv[1:])
