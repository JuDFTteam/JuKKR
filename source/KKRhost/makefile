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
mpi:       FC = mpiifort
mpidebug:   FC = mpiifort

####################################################
# Objects to compile                               #
####################################################


BIN = .
OBJ = ./OBJ
SRA = ./SRC_SRA
COM = ./SRC_COM
DRV = ./SRC_DRV
REL = ./SRC_REL
VXC = ./SRC_VXC
VIR = ./SRC_VIRTUAL

INCP = ./inc.p
INCR = ./sprkkr_rmesh.dim


mod_mympi.o  : $(COM)/mod_mympi.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/mod_mympi.F90 -c -o $@
mod_types.o  : $(COM)/mod_types.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/mod_types.F90 -c -o $@

FILENAME0 = kkr.x  #kkr0.exe
include makefile0
OBJECTS0 = $(PROGOBJS0)

FILENAME1a = kkr1a.exe
include makefile1a_omp
OBJECTS1a = $(PROGOBJS1a)

FILENAME1b = kkr1b.exe
include makefile1b_omp
OBJECTS1b = $(PROGOBJS1b)

FILENAME1c = kkr1c.exe
include makefile1c_omp
OBJECTS1c = $(PROGOBJS1c)

FILENAME2 = kkr2.exe
include makefile2
OBJECTS2 = $(PROGOBJS2)

$(OBJ)/main_all.o  : $(COM)/main_all.F90                 ; \
                        $(FC) $(FFLAGS) $(CPP) $(COM)/main_all.F90 -c -o $@

####################################################
# Libraries                                        #
####################################################

# standart libs
LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

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
FFLAGS       =  -O2 -r8 -traceback
# Change to if flag is matched:
openmp:   FFLAGS = -O2 -r8 
debug:    FFLAGS = -CB -check all -traceback -g
mpidebug: FFLAGS = -CB -check all -traceback -g

OPENMPFLAGS = 
openmp: OPENMPFLAGS = -openmp
####################################################
# Linking                                          #
####################################################

SOURCEDIR = .

mkdirobj:
	@if test ! -d $(OBJ); \
	then mkdir $(OBJ); fi

copyinc:
	cp -p inc.p   $(SRA)/inc.p; \
	cp -p inc.p   $(COM)/inc.p; \
	cp -p inc.p   $(VXC)/inc.p; \
	cp -p sprkkr_rmesh.dim   $(REL)/sprkkr_rmesh.dim
clean:
	rm -rf *.mod *.o OBJ/
clear:
	rm -rf *.mod *.o OBJ/

serial:   mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o 
	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)
debug:    mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o 
	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)
openmp:   mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBSOPENMP)
#mpi:      mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
#	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)
mpi:      mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBSOPENMP)
mpidebug: mkdirobj copyinc mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	  $(FC)     -o $(FILENAME0) mod_mympi.o mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)

install: clean serial

%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: $(SOURCEDIR)/%.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

