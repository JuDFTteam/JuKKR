#========================================================
# JUQUEEN: 
#FC = mpixlf90_r
#FFLAGS = -O2 -w -qarch=qp -qtune=qp -qsmp=omp -qnosave -traceback
#LDFLAGS = -L$(LAPACK_LIB) -L/bgsys/local/lib -lesslbg -llapack -lesslbg -qsmp=omp
#CPPFLAGS = -WF,-DCPP_MPI
#========================================================

#========================================================
# JUROPA
#FC = mpif90
#CPPFLAGS =  -D CPP_MPI
#FFLAGS = -O3 -r8
#LDFLAGS = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -lmkl -lguide -lmkl_lapack \
             -I/usr/local/fftw/3.2.1/include -L/usr/local/fftw/3.2.1/lib -lfftw3
#========================================================

#========================================================
# IFF cluster
FC = mpiifort
FFLAGS = -r8 -O3 -traceback
###FFLAGS += -openmp
CPPFLAGS =  -D CPP_MPI
###CPPFLAGS += -D CPP_DEBUG
###LDFLAGS= -L/usr/local/intel/Compiler/11.1/059/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 #OMP-Version of LApack
###LDFLAGS= -openmp -L/usr/local/intel/Compiler/11.1/059/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread #Standard LApack + separate OMP
###LDFLAGS= -L/usr/local/intel/Compiler/11.1/059/mkl/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread #Standard LApack
LDFLAGS= -llapack_ifort -lblas_ifort \
-L/usr/local/intel/lib/intel64 -lifcore -limf -Wl,-rpath,/usr/local/intel/lib/intel64 #Standard LApack with new Compiler V12
#========================================================

#========================================================
# RWTH cluster
#FC = $(MPIFC)
#FFLAGS = $(FLAGS_FAST)
#CPPFLAGS =  -DCPP_MPI
#LDFLAGS= $(FLAGS_MKL_LINKER) -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#========================================================


#######################################
###          PROGRAM NAMES          ###
#######################################
PKKR  = Pkkr.x
BAND  = band.x
MERG  = mergerefined.x
REFI  = refineBZparts.x
AMAT  = Amatprecalc.x
SPMX  = calculate_spinmixing.x
VISD  = visdata.x
VINT  = vis2int.x
TEST  = test.x
#TTLO = test_line_order.x


#######################################
###          SOURCES                ###
#######################################

PKKRobj =	type_inc.o mod_mympi.o mod_ioformat.o mod_ioinput.o mod_mathtools.o mod_parutils.o mod_vtkxml.o timing.o \
                mod_iohelp.o type_data.o mod_eigvects.o mod_dlke.o mod_spintools.o \
                mod_kkrmat.o mod_symmetries.o \
                mod_read.o \
                mod_fermisurf_basic.o mod_fermisurf_3D.o mod_fermisurf_2D.o mod_fermisurf.o mod_calconfs.o \
                mod_scattering.o \
                mod_routines.o

BANDobj =	type_inc.o mod_mympi.o mod_ioformat.o mod_ioinput.o mod_mathtools.o mod_parutils.o \
		type_data.o mod_eigvects.o mod_dlke.o \
		mod_kkrmat.o \
		mod_read.o \
		mod_bandstr.o \
                bands.o

MERGobj =	type_inc.o mod_mympi.o mod_ioformat.o mod_ioinput.o mod_mathtools.o mod_parutils.o mod_vtkxml.o timing.o \
                mod_iohelp.o type_data.o mod_eigvects.o mod_dlke.o mod_spintools.o \
                mod_kkrmat.o mod_symmetries.o \
                mod_read.o \
		mod_fermisurf_basic.o mod_fermisurf_3D.o mod_fermisurf_2D.o mod_fermisurf.o \
		mergerefined.o

REFIobj = 	mod_mathtools.o mod_vtkxml.o mod_ioinput.o mod_ioformat.o \
		refineBZparts.o

AMATobj =	type_inc.o mod_mympi.o mod_ioformat.o mod_ioinput.o mod_mathtools.o mod_parutils.o mod_vtkxml.o timing.o \
                mod_iohelp.o type_data.o mod_eigvects.o mod_dlke.o mod_spintools.o \
                mod_kkrmat.o mod_symmetries.o \
                mod_read.o \
                mod_fermisurf_basic.o mod_fermisurf_3D.o mod_fermisurf_2D.o mod_fermisurf.o mod_calconfs.o \
                mod_scattering.o \
                Amatprecalc.o

TESTobj =	test.o

#TTLOobj =       test_line_order.o

#######################################
###         COMPILING               ###
#######################################
.SUFFIXES: .f90 .F90
%.o: %.F90
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

#######################################
###           LINKING               ###
#######################################
$(PKKR):     $(PKKRobj) main.o
	$(FC) $^ -o $@ $(LDFLAGS)

$(BAND):     $(BANDobj)
	$(FC) $^ -o $@ $(LDFLAGS)

$(MERG):     $(MERGobj)
	$(FC) $^ -o $@ $(LDFLAGS)

$(REFI):     $(REFIobj)
	$(FC) $^ -o $@ $(LDFLAGS)

$(AMAT):     $(AMATobj)
	$(FC) $^ -o $@ $(LDFLAGS)

$(SPMX):     $(PKKRobj) calculate_spinmixing.o
	$(FC) $^ -o $@ $(LDFLAGS)

$(VISD):     $(PKKRobj) visdata.o
	$(FC) $^ -o $@ $(LDFLAGS)

$(VINT):     $(PKKRobj) vis2int.o
	$(FC) $^ -o $@ $(LDFLAGS)

$(TEST):     $(PKKRobj) $(TESTobj)
	$(FC) $^ -o $@ $(LDFLAGS)

#$(TTLO):     $(PKKRobj) test_line_order.o
#	$(FC) $^ -o $@ $(LDFLAGS)

##########################################################################
###                     COMMON DECLARATIONS                            ###
##########################################################################

default: $(PKKR)

clear:
	rm -f *.o *.mod *~ *.x 

clean:
	rm -f *.o *.mod *~ *.x 

all: $(PKKR) $(BAND) $(REFI) $(MERG) $(AMAT) $(SPMX) $(VISD) $(VINT)

install: clear all

##########################################################################
