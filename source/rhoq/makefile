

default: install

FC =  ifort
mpi: FC = mpiifort
mpidebug: FC = mpiifort


FILENAME0 = rhoq.x
# mpi: FILENAME0 = rhoq.x_mpi
# mpidebug: FILENAME0 = rhoq.x_mpi

####################################################
# Libraries                                        #
####################################################

# standart libs
#LLIBS  = -mkl -liomp5  #-mkl #-lmkl_intel_thread -lmkl_core -lpthread -liomp5
LLIBS  = -mkl=sequential -liomp5  #-mkl #-lmkl_intel_thread -lmkl_core -lpthread -liomp5
#LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5


# LLIBS = $(LLIBSOPENMP)

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP= -D noCPP 
#-D noCPP 
# change to if flag is matched
openmp:   CPP = -D CPP_OMP
mpidebug: CPP = -D CPP_MPI
mpi:      CPP = -D CPP_MPI

####################################################
# Compiler flags                                   #
####################################################
# Default:
# use -mcmodel=medium -shared-intel if you get the error message 
# '...relocation truncated to fit...' when trying to compile, 
# this happens for example for natyp>=75, iemxd>=101, KORBIT=1
FFLAGS       =  -O3 -qopenmp -xhost -g -traceback
#-vec-report1 -xhost
#-vec-report5
# -qopt-report=5 -qopt-report-phase=loop,vec -qopt-report-file=stdout
#-CB -check all -check noarg_temp_created -g
# Change to if flag is matched:
openmp:   FFLAGS = -O3 -qopenmp -xhost -g -qopt-report=5 -qopt-report-phase=loop,vec -qopt-report-file=optrep.txt
#debug:    FFLAGS = -O0 -openmp -r8 -CB -check all,noarg_temp_created -traceback
debug:    FFLAGS = -O0 -openmp -r8 -CB -check all,noarg_temp_created,uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -g
#FFLAGS = -O0 -CB -check all,noarg_temp_created  -traceback -g -fpe0 -warn -debug extended -qopenmp #-mcmodel=medium -shared-intel -qopenmp
mpidebug: FFLAGS = -g -CB -check all,noarg_temp_created  -traceback  -O0 -fpe0 -warn -debug extended -qopenmp

OPENMPFLAGS = 
openmp: OPENMPFLAGS = -qopenmp


####################################################
# Linking                                          #
####################################################

SOURCEDIR = .

clean:
	rm -rf *.mod *.o *__genmod.f90 *.i90
clear:
	rm -rf *.mod *.o *__genmod.f90 *.i90

SRCS = version.o version_info.o  mod_mympi.o mod_types.o mod_timing.o ioinput.o splint_real.o spline_real.o interpolspline.o gaunt_new.o gaunt2.o grule.o rcstop.o cspher.o calc_jlk.o opt.o test.o cinit.o intcheb_cell.o interpol_shape.o read_shape.o mod_rhoq.o

COMPVER='serial'
debug:    COMPVER='debug'
openmp:   COMPVER='openmp'
mpi:      COMPVER='mpi'
mpidebug: COMPVER='mpidebug'

getversion:
	@echo "Wrinting Version file: "; git describe;\
	echo $(COMPVER)>compver; \
	echo $(FFLAGS) $(OPENMPFLAGS)>compflag; \
	echo $(LLIBS) >complib; \
	./getversion.sh

compile:  getversion $(SRCS)
	$(FC) $(LLIBS) -o $(FILENAME0) $(SRCS) $

serial: compile
debug: compile
mpi: compile
openmp: compile
mpidebug: compile
install: clean mpi


mod_%.o: %.f90
	$(FC) -c $(FFLAGS)  $<
mod_%.o: %.f
	$(FC) -c $(FFLAGS)  $<
mod_%.o: %.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
mod_%.o: %.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

%.o: %.f90
	$(FC) -c $(FFLAGS)  $<
%.o: %.f
	$(FC) -c $(FFLAGS)  $<
%.o: %.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: %.F
	$(FC) -c $(FFLAGS) $(CPP)  $<
