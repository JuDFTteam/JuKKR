####################################################
####################################################
#         Julich KKR IMPURITY makefile
####################################################
#
# to compile use the following options:
#
# make           -    standart compilation using
#                     ifort on a single core
# make mpi       -    for using MPI parallelization 
#                     with ifort
# make gfort     -    uses gfort as the compiler
# make pg        -    uses portland group compiler
# make mac       -    uses compiler options for MacOSX
# make clean     -    removes all *.o files
# make openmp    -    uses openmp
#
####################################################
# specify the compiler which is going to be used   #
####################################################
# default value
FC =  ifort
# on Aachen cluster:
# $(FC) points to standard compiler in environment
#

# change default value if a flag is matched
mpi:        FC = mpiifort
# Aachen cluster:
#mpi:        FC = $(MPIFC)
mpitrace:  FC = vtf90 -vt:f90 $(MPIFC)
mpidebug:   FC = mpiifort
gfort:      FC = /usr/users/iff_th1/bauer/bin/gfortran
pg:         FC = pgf90 
hybrid:     FC= mpiifort -traceback 


# default
COMPVER = 'serial'
mpi:       COMPVER = 'mpi'
mpidebug:  COMPVER = 'mpidebug'
openmp:    COMPVER = 'openmp'
gfort:     COMPVER = 'gfort'
pg:        COMPVER = 'pg'
hybrid:    COMPVER = 'hybrid'


# default
COMPVER = 'serial'
mpi:       COMPVER = 'mpi'
mpidebug:  COMPVER = 'mpidebug'
openmp:    COMPVER = 'openmp'
gfort:     COMPVER = 'gfort'
pg:        COMPVER = 'pg'
hybrid:    COMPVER = 'hybrid'

####################################################
# Objects to compile                               #
####################################################


OBJECTS = mod_mympi.o mod_types.o mod_timing.o chebint.o rllsll.o 
#OBJECTS = rllsllsourceterms.o nrtype.o chebint.o type_config.o type_cellnew.o checknan.o mathtools.o config.o chebint.o physic_params.o beshank.o timing.o chebyshev.o rllslltools.o rllsllsourceterms.o rllsll.o 

####################################################
# Libraries                                        #
####################################################

# standart libs
LLIBS = -llapack_ifort -lblas_ifort  -L/usr/local/intel/lib/intel64 -lifcore -limf  -Wl,-rpath,/usr/local/intel/lib/intel64
# Aachen cluster
### for sequential libraries (no hybrid parallelisation) use this
#LLIBS = $(FLAGS_MKL_LINKER) #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# Aachen cluster
#LLIBSOPENMP = $(FLAGS_MKL_LINKER) -lmkl_intel_lp64 

### link parallel omp libraries for hybrid parallelisation
LLIBSOPENMP+=  -lmkl_intel_thread -liomp5


# Special libraries for gfortran:
LLIBSGFORT = -L/usr/local/nag/lib -lnag -llapack_ifort -lblas_ifort  -L/usr/local/intel/current/lib/intel64 -lifcore -liomp5 -lintlc -limf -lirc -lpthread -Wl,-rpath,/usr/local/intel/current/lib/intel64

# MacOSX libraries:
LLIBSMAC  =  -mkl #-L/usr/local/Library/Frameworks/Intel_MKL.framework/Libraries/32 -lmkl_intel \
              -lmkl_intel_thread -lmkl_core -lguide -lpthread -i_dynamic

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP=
# change to if flag is matched
mpidebug:      CPP = -D MPI
mpitrace:      CPP = -D MPI
mpi:      CPP = -D MPI
gfort:    CPP = -D GFORT
hybrid:   CPP = -D MPI

####################################################
# Compiler flags                                   #
####################################################
# Default:
FFLAGS       =  -O2 
# Change to if flag is matched:
gfort:    FFLAGS = -O0 -fcheck=all -frange-check  -fbounds-check -g -ff2c
pg:       FFLAGS = -C
debug:    FFLAGS = -CB -check all -traceback -g
#debug:    FFLAGS = -debug all -fno-omit-frame-pointer -ftrapuv -check all -traceback  -O0
mpidebug: FFLAGS = -CB -check all -traceback -g
hybrid:   FFLAGS = -O2 -openmp
####################################################
# Linking                                          #
####################################################
FILENAME = test_rllsll.exe
SOURCEDIR = .

.PHONY: cbg
.PHONY: hybrid
.PHONY: openmp
.PHONY: debug
.PHONY: pg
.PHONY: gfort
.PHONY: mpi
.PHONY: mpitrace
.PHONY: mpidebug
.PHONY: mac

cbg:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
hybrid:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBSOPENMP)
openmp:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBSOPENMP)
debug:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
pg:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
gfort:  $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBSGFORT)
mpi:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
mpitrace:  $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
mpidebug:  $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBS)
mac:    $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS)  $(LLIBSMAC)

%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: $(SOURCEDIR)/%.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

####################################################
# Clean up flags                                   #
####################################################

clean:
	rm -f *.mod *.o $(FILENAME)
cleanxc:
	rm -f exchangecorrel.o

