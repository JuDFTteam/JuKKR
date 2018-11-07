cmake_minimum_required(VERSION 2.6.4)

project(JuKKR)
enable_language(Fortran)

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
option(COMPILE_KKRIMP "Compile KKRimp" ON)
option(COMPILE_KKRSUSC "Compile KKRsusc" OFF)
option(COMPILE_PKKPRIME "Compile PKKprime" ON)
option(COMPILE_VORONOI "Compile voronoi" ON)
option(COMPILE_RHOQ "Compile rhoq" ON)


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

