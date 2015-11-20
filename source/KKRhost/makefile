####################################################
####################################################
#         Julich KKR code makefile
####################################################
#
# to compile use the following options:
#
# make           -    standart compilation using
#                     ifort on a single core
# make mpi       -    for using MPI parallelization 
#                     with ifort
# make clean     -    removes all *.o files
# make openmp    -    uses openmp
#
default: serial
####################################################
# specify the compiler which is going to be used   #
####################################################
# default value
FC =  ifort
# change default value if a flag is matched
mpi:       FC =  mpiifort
# mac:
mpi:       FC =  /usr/local/impi/bin/mpif90  
mpidebug:   FC = mpiifort


# information for version-tracking in executable:
# default
COMPVER = 'serial'
debug:     COMPVER = 'debug'
mpi:       COMPVER = 'mpi'
mpidebug:  COMPVER = 'mpidebug'
openmp:    COMPVER = 'openmp'


####################################################
# Objects to compile                               #
####################################################


BIN = .
OBJ = ./OBJ
SRA = ./ALL_SOURCE_FILES
COM = ./ALL_SOURCE_FILES
DRV = ./ALL_SOURCE_FILES
REL = ./ALL_SOURCE_FILES
VXC = ./ALL_SOURCE_FILES
VIR = ./ALL_SOURCE_FILES

INCP = ./inc.p
INCR = ./sprkkr_rmesh.dim


$(OBJ)/mod_mympi.o  : $(COM)/mod_mympi.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/mod_mympi.F90 -c -o $@
$(OBJ)/mod_types.o  : $(COM)/mod_types.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/mod_types.F90 -c -o $@
$(OBJ)/mod_timing.o  : $(COM)/mod_timing.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/mod_timing.F90 -c -o $@

FILENAME0 = kkr.x
include makefile0
OBJECTS0 = $(PROGOBJS0)

include makefile1a_omp
OBJECTS1a= $(PROGOBJS1a)

include makefile1b_omp
OBJECTS1b= $(PROGOBJS1b)

include makefile1c_omp
OBJECTS1c= $(PROGOBJS1c)

include makefile2
OBJECTS2 = $(PROGOBJS2)

$(OBJ)/main_all.o  : $(COM)/main_all.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/main_all.F90 -c -o $@

####################################################
# Libraries                                        #
####################################################

# standart libs
# serial:
#LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
# openmp:
#LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# mac libs:
#libfolder = /usr/local/intel-12.1.3/compiler/lib/intel64
#LLIBS=-L${libfolder} -lifcore -limf -Wl,-rpath,${libfolder} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
LLIBS = -mkl -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000

# openmp libs
#LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# mac:
#LLIBSOPENMP = -L${libfolder} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
LLIBSOPENMP = -mkl -Wl,-stack_size,0x40000000,-stack_addr,0xf0000000

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP=
# change to if flag is matched
openmp:   CPP = -D CPP_OMP
mpidebug: CPP = -D CPP_MPI
mpi:      CPP = -D CPP_MPI

####################################################
# Compiler flags                                   #
####################################################
# Default:
FFLAGS       =  -O2 -r8 -traceback -module $(OBJ) 

# Change to if flag is matched:
openmp:   FFLAGS = -O2 -r8 -traceback   -module $(OBJ)
debug:    FFLAGS = -CB -check all -traceback -g -fpe0 -warn -debug extended -mcmodel=medium -shared-intel
mpidebug: FFLAGS = -CB -check all -traceback -mcmodel=medium

OPENMPFLAGS = 
openmp: OPENMPFLAGS = -openmp

####################################################
# Linking                                          #
####################################################

SOURCEDIR = .

.PHONY: mkdirobj
.PHONY: copyinc
.PHONY: clean
.PHONY: clear
.PHONY: serial
.PHONY: debug
.PHONY: openmp
.PHONY: mpi
.PHONY: mpidebug
.PHONY: install
.PHONY: getversion

getversion:
	@echo $(COMPVER)>compver; \
	$(COM)/getversion.sh

mkdirobj:
	@if test ! -d $(OBJ); \
	then mkdir $(OBJ); fi

copyinc:
	cp -p inc.p   $(SRA)/inc.p; \
	cp -p sprkkr_rmesh.dim   $(REL)/sprkkr_rmesh.dim
clean:
	rm -rf *.mod *.o OBJ/ *__genmod.*
clear:
	rm -rf *.mod *.o OBJ/ *__genmod.*

serial:   mkdirobj getversion copyinc $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o 
	  $(FC)     -o $(FILENAME0) $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)
debug:    mkdirobj getversion copyinc $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o 
	  $(FC)     -o $(FILENAME0) $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)
openmp:   mkdirobj getversion copyinc $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBSOPENMP)
mpi:      mkdirobj getversion copyinc $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBSOPENMP)
mpidebug: mkdirobj getversion copyinc $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) $(OBJ)/mod_mympi.o $(OBJ)/mod_types.o $(OBJ)/mod_timing.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)

install: clean serial

%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: $(SOURCEDIR)/%.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

