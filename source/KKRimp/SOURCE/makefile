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
#FC =  ifort
# on Aachen cluster:
# $(FC) points to standard compiler in environment
#

# change default value if a flag is matched
#mpi:        FC = mpiifort
mpi:        FC = $(MPIFC)
mpitrace:  FC = vtf90 -vt:f90 $(MPIFC)
mpidebug:   FC = mpiifort
gfort:      FC = /usr/users/iff_th1/bauer/bin/gfortran
pg:         FC = pgf90 
hybrid:     FC= $(MPIFC)

####################################################
# Objects to compile                               #
####################################################
OBJECTS_TYPE  =   nrtype.o  type_gmatbulk.o type_tmat.o type_gmat.o type_gmatonsite.o type_gauntshape.o type_gauntcoeff.o type_corestate.o type_gref.o type_tbcluster.o type_density.o type_shapefun.o type_cell.o type_cellnew.o type_config.o type_energyparts.o type_wavefunction.o type_ldau.o

OBJECTS_MATH =  physic_params.o checknan.o log.o mathtools.o config.o beshank.o wrapper.o simp3.o simpk.o regsol.o soutk.o

OBJECTS_DIRAC  = DIRAC/DiracConfig.o DIRAC/Constants.o DIRAC/SpinSphericals.o DIRAC/Lebedev.o DIRAC/RelativisticGauntCoefficients.o DIRAC/Potential.o    DIRAC/SourceTerms.o  
# DIRAC/BasisTransform.o
OBJECTS_DIRAC2  = Constants.o  Lebedev.o SpinSphericals.o RelativisticGauntCoefficients.o Potential.o    SourceTerms.o  

OBJECTS_LDAU = rwldaupot.o averagewldau.o phicalc.o initldau.o rclm.o calcwldau.o

OBJECTS_MAIN  =  wavefunctodisc.o change_nrmin.o gauntharmonics.o  sinwk.o  ymy.o amngaunt.o amn2010.o convol.o shftvout.o vintras.o vinters2010.o timing.o arrayparams.o cinit.o dsort.o clustcomp.o  gauntshape.o read_spinorbit.o read_angle.o read_potential.o read_atominfo.o irwsol.o beshan.o  cradwf.o wfmesh.o zgeinv1.o wfint0.o wfint.o csout.o wftsca.o vllns.o csinwd.o regns.o pnstmat.o calctmat.o calctmatfull.o irwns.o csimpk.o pnsqns.o  BasisTransform.o chebyshev.o cheb2oldgridc.o  chebint.o mixbroydenspin.o rotaterealspace.o rotatespinframe.o rotatespin.o orbitalmoment.o spinorbit.o vllmatsra.o vllmat_new.o calccouplingconstants.o rllslltools.o rllsllsourceterms.o rllsll.o complexdos3.o rholm.o rhoin.o rhoout.o rhooutnew.o rhons.o rhoval.o rhoval_new.o rhovalfull.o gdyson.o mpienergy.o interpolpot-splint.o interpolpot.o interpolatecell.o checkinterpolation.o calcsph.o wronskian.o  calctmat_bauernew_testtools.o calctmat_bauernew.o energyloop.o wrmoms.o rhocoreint.o rhocore.o rhototb.o dysonvirtatom.o preconditioning.o ecoub.o exchangecorrel.o espcb.o epotinb.o etotb1.o vxcdrv.o mixstr.o mixbroyden.o rites.o force.o forceh.o forcxc.o calcforce.o cnwf011.o bessel1.o jmtrx.o umatrix.o utrafo.o cellproperties.o  kkrflex.o 

OBJECTS =  $(OBJECTS_TYPE) $(OBJECTS_MATH) $(OBJECTS_DIRAC)  $(OBJECTS_LDAU) $(OBJECTS_MAIN) 
OBJECTS2 = $(OBJECTS_TYPE) $(OBJECTS_MATH) $(OBJECTS_DIRAC2) $(OBJECTS_LDAU) $(OBJECTS_MAIN) 

####################################################
# Libraries                                        #
####################################################

# standart libs
#LLIBS = -llapack_ifort -lblas_ifort  -L/usr/local/intel/lib/intel64 -lifcore -limf  -Wl,-rpath,/usr/local/intel/lib/intel64
# Aachen cluster
#LLIBS =-L/usr/local/intel/lib/intel64 -lifcore -limf -Wl,-rpath,/usr/local/intel/lib/intel64 -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
### for sequential libraries (no hybrid parallelisation) use this
LLIBS = $(FLAGS_MKL_LINKER) #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# openmp libs
#LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
### link parallel omp libraries for hybrid parallelisation
LLIBSOPENMP = $(FLAGS_MKL_LINKER) -lmkl_intel_lp64 
LLIBSOPENMP+=  -lmkl_intel_thread -liomp5

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
####################################################
# Linking                                          #
####################################################
FILENAME = kkrflex.exe
SOURCEDIR = .

cbg:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
hybrid:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBSOPENMP)
openmp:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBSOPENMP)
debug:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
pg:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
gfort: $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBSGFORT)
mpi:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
mpitrace: $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
mpidebug: $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)
mac:   $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBSMAC)

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
	rm -f *.mod *.o ../$(FILENAME)
cleanxc:
	rm -f exchangecorrel.o

