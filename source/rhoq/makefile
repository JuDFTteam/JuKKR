

default: serial

FC =  ifort

FILENAME0 = rhoq.x

####################################################
# Libraries                                        #
####################################################

# standart libs
LLIBS  = -mkl #-lmkl_intel_thread -lmkl_core -lpthread -liomp5
#LLIBS  = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

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
# use -mcmodel=medium -shared-intel if you get the error message 
# '...relocation truncated to fit...' when trying to compile, 
# this happens for example for natyp>=75, iemxd>=101, KORBIT=1
FFLAGS       =  -O3 -openmp -traceback -g #-CB -check all -check noarg_temp_created -g
# Change to if flag is matched:
openmp:   FFLAGS = -O2 #–Wl,-stack_size,0x1000000000,-stack_addr,0xc000000000 -heap-arrays # -r8 
debug:    FFLAGS = -O0 -CB -check all,noarg_temp_created  -traceback -g -fpe0 -warn -debug extended -openmp #-mcmodel=medium -shared-intel -openmp
mpidebug: FFLAGS = -g -CB -check all,noarg_temp_created  -traceback #-mcmodel=medium

OPENMPFLAGS = 
openmp: OPENMPFLAGS = -openmp
####################################################
# Linking                                          #
####################################################

SOURCEDIR = .

clean:
	rm -rf *.mod *.o *__genmod.f90 *.i90
clear:
	rm -rf *.mod *.o *__genmod.f90 *.i90

SRCS = ioinput.o splint_real.o spline_real.o interpolspline.o gaunt_new.o gaunt2.o grule.o rcstop.o cspher.o calc_jlk.o opt.o test.o cinit.o mod_mympi.o mod_types.o intcheb_cell.o interpol_shape.o read_shape.o mod_rhoq.o

serial:  $(SRCS) 
	  $(FC) $(LLIBS)    -o $(FILENAME0) $(SRCS) $
debug:   $(SRCS) 
	  $(FC) $(LLIBS)    -o $(FILENAME0) $(SRCS) $
install: clean serial

%.o: %.f90
	$(FC) -c $(FFLAGS)  $<
%.o: %.f
	$(FC) -c $(FFLAGS)  $<
%.o: %.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: %.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

