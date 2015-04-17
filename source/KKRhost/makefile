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




FILENAME0 = kkr.x  #kkr0.exe
include makefile0
OBJECTS0 = $(PROGOBJS0)

mod_types.o  : $(COM)/mod_types.f90                 ; \
                        $(FC) $(FFLAGS) $(COM)/mod_types.f90 -c -o $@

$(OBJ)/main_all.o  : $(COM)/main_all.F90                 ; \
                        $(FC) $(FFLAGS) $(COM)/main_all.F90 -c -o $@

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

####################################################
# Libraries                                        #
####################################################

# standart libs
#LLIBS = -llapack_ifort -lblas_ifort  -L/usr/local/intel/lib/intel64 -lifcore -limf  -Wl,-rpath,/usr/local/intel/lib/intel64
LLIBS  = -L/usr/local/intel/lib/intel64 -lifcore -limf -Wl,-rpath,/usr/local/intel/lib/intel64 -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

# Special libraries for gfortran:
LLIBSGFORT = -L/usr/local/nag/lib -lnag -llapack_ifort -lblas_ifort  -L/usr/local/intel/current/lib/intel64 -lifcore -liomp5 -lintlc -limf -lirc -lpthread -Wl,-rpath,/usr/local/intel/current/lib/intel64

# MacOSX libraries:
LLIBSMAC  =  -L/usr/local/Library/Frameworks/Intel_MKL.framework/Libraries/32 -lmkl_intel \
              -lmkl_intel_thread -lmkl_core -lguide -lpthread -i_dynamic

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP=
# change to if flag is matched
mpidebug:      CPP = -D MPI
mpi:      CPP = -D MPI
gfort:    CPP = -D GFORT

####################################################
# Compiler flags                                   #
####################################################
# Default:
FFLAGS       =  -O2 -c -r8  -traceback
# Change to if flag is matched:
debug:    FFLAGS = -CB -check all -traceback -g
mpidebug: FFLAGS = -CB -check all -traceback -g
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

#serialall: mkdirobj copyinc $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS1c) $(OBJ)/main_all.o
#	$(FC)     -o $(FILENAME0) $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS1c) $(OBJ)/main_all.o $(LLIBS)
serialall: mkdirobj copyinc mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o
	$(FC)     -o $(FILENAME0) mod_types.o $(OBJECTS0) $(OBJECTS1a) $(OBJECTS1b) $(OBJECTS1c) $(OBJECTS2) $(OBJ)/main_all.o $(LLIBS)

serial0:   mkdirobj copyinc $(OBJECTS0)
	$(FC)     -o $(FILENAME0) $(OBJECTS0)  $(LLIBS)
openmp0:   mkdirobj copyinc $(OBJECTS0)
	$(FC)     -o $(FILENAME0) $(OBJECTS0)  $(LLIBSOPENMP)
debug0:    mkdirobj copyinc $(OBJECTS0)
	$(FC)     -o $(FILENAME0) $(OBJECTS0)  $(LLIBS)
mpi0:      mkdirobj copyinc $(OBJECTS0)
	$(FC)     -o $(FILENAME0) $(OBJECTS0)  $(LLIBS)
mpidebug0: mkdirobj copyinc $(OBJECTS0)
	$(FC)     -o $(FILENAME0) $(OBJECTS0)  $(LLIBS)
serial1a:   mkdirobj copyinc $(OBJECTS1a)
	$(FC)     -o $(FILENAME1a) $(OBJECTS1a)  $(LLIBS)
openmp1a:   mkdirobj copyinc $(OBJECTS1a)
	$(FC)     -o $(FILENAME1a) $(OBJECTS1a)  $(LLIBSOPENMP)
debug1a:    mkdirobj copyinc $(OBJECTS1a)
	$(FC)     -o $(FILENAME1a) $(OBJECTS1a)  $(LLIBS)
mpi1a:      mkdirobj copyinc $(OBJECTS1a)
	$(FC)     -o $(FILENAME1a) $(OBJECTS1a)  $(LLIBS)
mpidebug1a: mkdirobj copyinc $(OBJECTS1a)
	$(FC)     -o $(FILENAME1a) $(OBJECTS1a)  $(LLIBS)
serial1b:   mkdirobj copyinc $(OBJECTS1b)
	$(FC)     -o $(FILENAME1b) $(OBJECTS1b)  $(LLIBS)
openmp1b:   mkdirobj copyinc $(OBJECTS1b)
	$(FC)     -o $(FILENAME1b) $(OBJECTS1b)  $(LLIBSOPENMP)
debug1b:    mkdirobj copyinc $(OBJECTS1b)
	$(FC)     -o $(FILENAME1b) $(OBJECTS1b)  $(LLIBS)
mpi1b:      mkdirobj copyinc $(OBJECTS1b)
	$(FC)     -o $(FILENAME1b) $(OBJECTS1b)  $(LLIBS)
mpidebug1b: mkdirobj copyinc $(OBJECTS1b)
	$(FC)     -o $(FILENAME1b) $(OBJECTS1b)  $(LLIBS)
serial1c:   mkdirobj copyinc $(OBJECTS1c)
	$(FC)     -o $(FILENAME1c) $(OBJECTS1c)  $(LLIBS)
openmp1c:   mkdirobj copyinc $(OBJECTS1c)
	$(FC)     -o $(FILENAME1c) $(OBJECTS1c)  $(LLIBSOPENMP)
debug1c:    mkdirobj copyinc $(OBJECTS1c)
	$(FC)     -o $(FILENAME1c) $(OBJECTS1c)  $(LLIBS)
mpi1c:      mkdirobj copyinc $(OBJECTS1c)
	$(FC)     -o $(FILENAME1c) $(OBJECTS1c)  $(LLIBS)
mpidebug1c: mkdirobj copyinc $(OBJECTS1c)
	$(FC)     -o $(FILENAME1c) $(OBJECTS1c)  $(LLIBS)
serial2:   mkdirobj copyinc $(OBJECTS2)
	$(FC)     -o $(FILENAME2) $(OBJECTS2)  $(LLIBS)
openmp2:   mkdirobj copyinc $(OBJECTS2)
	$(FC)     -o $(FILENAME2) $(OBJECTS2)  $(LLIBSOPENMP)
debug2:    mkdirobj copyinc $(OBJECTS2)
	$(FC)     -o $(FILENAME2) $(OBJECTS2)  $(LLIBS)
mpi2:      mkdirobj copyinc $(OBJECTS2)
	$(FC)     -o $(FILENAME2) $(OBJECTS2)  $(LLIBS)
mpidebug2: mkdirobj copyinc $(OBJECTS2)
	$(FC)     -o $(FILENAME2) $(OBJECTS2)  $(LLIBS)


serial:   serial0 serial1a serial1b serial1c serial2
openmp:   openmp0 openmp1a openmp1b openmp1c openmp2
debug:    debug0 debug1a debug1b debug1c debug2
mpi:      mpi0 mpi1a mpi1b mpi1c mpi2
mpidebug: mpidebug0 mpidebug1a mpidebug1b mpidebug1c mpidebug2 

clean:
	rm -rf *.mod *.o OBJ/
clear:
	rm -rf *.mod *.o OBJ/

install: clean serial

%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: $(SOURCEDIR)/%.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

