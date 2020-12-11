cmake_minimum_required(VERSION 2.6.4)

project(JuKKR)
enable_language(Fortran)

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# here the executables are written to
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ./ )

# maybe useful for future library dependecies
set(complib "")

# we define some compiler flags
include( ./cmake/compiler_config.txt )

#define build type and compile options
include( ./cmake/compile_options.txt )

# default values for compile options (see cmake/compile_options.txt for details)
option(TIMING_DETAIL "Write detailed timing information (e.g. in rllsll)" OFF)
option(ENABLE_DEBUG_SYMBOLS "Add debug symbols (-g) everywhere" OFF)
option(ENABLE_OPTRPT "Generate optimization reports (*.optrpt files) at compilation" OFF)
option(OPTFLAGS_xHOST "Turn on xHost optimization flag for release build" ON)
option(OPTFLAGS_AVX512 "Use AVX512 instruction set for release build" OFF)
option(OPTFLAGS_AMD "FLAGS for AMD64 node on iffslurm (see compiler_config.txt)" OFF)
option(ENABLE_MPI "Enable MPI parallelization" ON)
option(ENABLE_OMP "Enable OpenMP parallelization" OFF)
option(ENABLE_OMP_EVERYWHERE "Add openmp compile flags everywhere" OFF)


# default values of codes to compile:
option(COMPILE_KKRHOST "Compile KKRhost" ON)
option(COMPILE_KKRIMP "Compile KKRimp" OFF)
option(COMPILE_KKRSUSC "Compile KKRsusc" OFF)
option(COMPILE_PKKPRIME "Compile PKKprime" OFF)
option(COMPILE_VORONOI "Compile voronoi" ON)
option(COMPILE_RHOQ "Compile rhoq" OFF)
option(WRITE_RLLSLL_STANDALONE_INPUT "Set option to write out rllsll input files" OFF)
option(ENABLE_RLLSLL_STANDALONE "Compile rllsll standalone test exec" OFF)

if(ENABLE_RLLSLL_STANDALONE)
   include( ./cmake/source_list_rllsll_test.txt )
endif()

if(COMPILE_KKRHOST)
   include( ./cmake/source_list_KKRhost.txt )
endif()

if(COMPILE_KKRIMP)
   include( ./cmake/source_list_KKRimp.txt )
endif()

if(COMPILE_KKRSUSC)
   include( ./cmake/source_list_KKRsusc.txt )
endif()

if(COMPILE_PKKPRIME)
   include( ./cmake/source_list_PKKprime.txt )
endif()

if(COMPILE_VORONOI)
   include( ./cmake/source_list_voronoi.txt )
endif()

if(COMPILE_RHOQ)
   include( ./cmake/source_list_rhoq.txt )
endif()

