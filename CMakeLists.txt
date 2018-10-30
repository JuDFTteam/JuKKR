cmake_minimum_required(VERSION 2.6.4)

project(JuKKR)
enable_language(Fortran)

#call cmake like this:
#
# mkdir build
# cd build
# FC=<compiler-you-want-to-use> cmake -D<options> ..
#
# where with FC you coose the compiler, e.g. 
#  FC=gfortran, FC=mpif90 (gfortran with MPI), FC=ifort, FC=mpiifort
#
# List of -D<options>=default:
#  -DENABLE_MPI=ON
#  -DENABLE_OMP=OFF
#  -DENABLE_COV=OFF
#  -DCMAKE_BUILD_TYPE=Release
#  -DENABLE_BdG=OFF

# this is where we will place the Fortran module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

# here the executables are written to
set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ./ )


# we define some compiler flags
include( ./cmake/compiler_config.txt )

#define build type and compile options
include( ./cmake/compile_options.txt )


# maybe useful for future library dependecies
set(complib "")


# default values of codes to compile:
option(COMPILE_KKRHOST "Compile KKRhost" ON)
option(COMPILE_KKRIMP "Compile KKRimp" OFF)
option(COMPILE_KKRSUSC "Compile KKRsusc" OFF)
option(COMPILE_PKKPRIME "Compile PKKprime" OFF)
option(COMPILE_VORONOI "Compile voronoi" OFF)


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

