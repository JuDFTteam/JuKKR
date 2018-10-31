#!/usr/bin/env python3

####################################################
### Installation script for the Juelich KKR code ###
####################################################

# import modules
import os
import sys
import getopt
import shutil

##########################################################################
#global settings:

codeslist = ['kkrhost', 'kkrimp', 'voronoi'] # still to add 'kkrsusc', pkkprime

##########################################################################
 
def greeting():
   """ Prints greeting. """
   print("***********************************************************")
   print("Welcome to the Installation script of the Juelich KKR codes.")
   print("\nYou can find useful additional information in our KKR wiki: \033[4mhttps://kkr.iff.kfa-juelich.de/doku.php\033[0m")
   print("and in our code's documentation: \033[4mhttps://kkr.iffgit.fz-juelich.de/jukkr/index.html\033[0m")
   print("\nCurrently this script only supports individual builds of one of these codes:")
   for code in codeslist:
      print("  * {}".format(code))
   print("\nTo build a different code please go to the source/code subdirectory and build the code as in a previous version (see wiki).\n")
   print("***********************************************************\n")

##########################################################################
 
def goodbye():
   """ Prints final remarks before exiting. """
   print("\n***********************************************************")
   print("Congratulations you have now set up the build directory.")
   print("To compile the code please go to the build directory and type")
   print("    'make'")
   print(" or")
   print("    'make -j n'")
   print(" where n is the number of tasks you want to use in the parallel compilation process.")
   print("***********************************************************")

##########################################################################

def usage():
   """ Prints usage info and exists. """
   print("Acceptable arguments:")
   print("  -h or --help               Print this message and exit.")
   print("  -i or --interactive        Use interactive installation script asking the user for input.")
   print("  -v or --verbose            Verbose mode.")
   print("  -d or --debug              Set up debug mode of the code.")
   print("  --program=name             Set up build of this program. Available are: {}".format(codeslist))
   print("  --machine=name             Use a predefined set of settings for a specific machine where 'name' is one of 'iff', 'claix', 'jureca'.")
   print("  --compiler=name            Use a specific compiler, 'name' could for example be 'mpiifort' or 'gfortran'.")
   print("  --parallelization=scheme   Use either MPI, OpenMP or both (hybrid) parallelization: 'scheme should be one of 'mpi', 'openmp', 'hybrid'.")
   print("  --flags=flag1,flag2        Add cmake flags manually (can be combined with either -m or the settings with -c and -p).")
   sys.exit()

##########################################################################

def read_interactive(flags, codeslist):
   """ Interactively ask user for input and return the neede options. """

   print("Please input the the codename of the code you want to compile (one of {}).".format(codeslist))
   while True:
      code = input()
      if code not in codeslist:
         print("your input is not a valid code name. Please chose a valid scheme.")
      else:
         break

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

   return compiler, parallelization, flags, code

##########################################################################

def read_machine(flags, mname):
   """ Read options from machine file. """
   print("read_machine: Not implemented yet.")
   sys.exit()
   return compiler, parallelization, flags

##########################################################################

def check_dependencies(verbose):
   """ Check if all necessary dependencies are available (e.g. cmake). """
   print("check-dependencies...")
   # check make
   print("make availabale?")
   if shutil.which('make') is None:
      print("Command 'make' not found. Maybe you have to import some modules?")
      sys.exit()
   print("OK")
   # check cmake
   print("cmake availabale?")
   if shutil.which('cmake') is None:
      print("Command 'cmake' not found. Maybe you have to import some modules?")
      sys.exit()
   print("OK")

##########################################################################

def create_build_dir(verbose):
   """ Create a build directory and store eventually previously existing ones. """
   # check for existence of build dir
   if 'build' in os.listdir('.'):
      i = 1
      oldbuild = 'build_{}'.format(i)
      while oldbuild in os.listdir('.'):
         i += 1
         oldbuild = 'build_{}'.format(i)
      print("Found old build directory. Moving this to {}".format(oldbuild))
      if verbose:
         print("Continue (will rename old build dir)? [Y/n]")
         answer = input()
         if 'n' in answer.lower(): sys.exit()
      shutil.move('build', oldbuild)

   # create build dir
   if verbose:
      print("Continue (will create build dir)? [Y/n]")
      answer = input()
      if 'n' in answer.lower(): sys.exit()
   os.makedirs('build')

##########################################################################

def run_cmake(compiler, parallelization, flags, verbose, code):
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

   # add code compile flags
   comp_host = 'OFF'
   comp_imp = 'OFF'
   comp_susc = 'OFF'
   comp_pkkr = 'OFF'
   comp_voro = 'OFF'
   if code=='kkrhost':
      comp_host = 'ON'
   elif code=='kkrimp':
      comp_imp = 'ON'
   elif code=='kkrsusc':
      comp_susc = 'ON'
   elif code=='pkkprime':
      comp_pkkr = 'ON'
   elif code=='voronoi':
      comp_voro = 'ON'
   else:
      print("Error setting code option!")
      sys.exit()
   task += " -DCOMPILE_KKRHOST={} ".format(comp_host)
   task += " -DCOMPILE_KKRIMP={} ".format(comp_imp)
   task += " -DCOMPILE_KKRSUSC={} ".format(comp_susc)
   task += " -DCOMPILE_KKRPKKPRIME={} ".format(comp_pkkr)
   task += " -DCOMPILE_VORONOI={} ".format(comp_voro)

   # add additional flags if given
   for flag in flags:
      task += "-D{} ".format(flag)

   # finalize task and run
   task += '..'

   print("\nNow run cmake command: '{}'".format(task))
   print("\n***********************************************************\n")

   if verbose:
      print("Continue (will run cmake command)? [Y/n]")
      answer = input()
      if 'n' in answer.lower(): sys.exit()
   call(task, shell=True)

##########################################################################

def main(argv):
   # print greeting lines
   greeting()
   # process script options
   try:
      if len(argv)==0: usage()
      opts, args = getopt.getopt(argv, "ivhdm", ["interactive", "verbose", "help", "debug", "machine=", "compiler=", "parallelization=", "flags=","program="])
      #print(argv, len(argv), opts,args)
   except getopt.GetoptError:
      usage()

   # define defaults
   compiler = None
   parallelization = None
   flags = []
   verbose = False
   code = None

   # first check for vebosity level (determines additional printing)
   for opt, arg in opts:
      if opt in ("-v", "--verbose"):
         verbose = True

   # now deal with different input options:
   for opt, arg in opts:
      if opt in ("-h", "--help"):
         usage()

      elif opt in ("-i", "--interactive"):
         compiler, parallelization, flags, code = read_interactive(flags, codeslist)

      elif opt=="--machine":
         compiler, parallelization, flags = read_machine(flags, args)

      else: # read in options from compiler and parallelization info
         if opt=="--compiler":
            compiler = arg
            if compiler=='': compiler = None
         
         if opt=="--parallelization":
            parallelization = arg
            if parallelization=='': parallelization = None
         
      # finally allow to add flags additionally
      if "-i" not in opts and "--interactive" not in opts:
          if opt=="--flags" :
             for iarg in arg.split(','):
                 flags.append(iarg)

          if opt=="--program":
             code = arg
             if code not in codeslist:
                print("your input is not a valid code name. Please chose a valid scheme.")
                print("Supported codes: {}".format(codelist))

      if opt in ("-d", "--debug"):
          flags.append("CMAKE_BUILD_TYPE=Debug")

   # check dependencies
   create_build_dir(verbose)

   # create build dir
   check_dependencies(verbose)

   # run cmake with previously read-in options
   run_cmake(compiler, parallelization, flags, verbose, code)

   # print final messages
   goodbye()

##########################################################################

# run script
if __name__ == "__main__":
   main(sys.argv[1:])
