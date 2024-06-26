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
default: install
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
rwthmpi:    FC = $(MPIFC)
mpitrace:   FC = vtf90 -vt:f90 $(MPIFC)
mpidebug:   FC = mpiifort
gfort:      FC = gfortran
pg:         FC = pgf90 
hybrid:     FC= mpiifort
rwthhybrid: FC= $(MPIFC) 
jureca:     FC = mpif90
hazelhen:   FC = ftn # run 'module swap PrgEnv-cray PrgEnv-intel' and 'module unload cray-libsci' in terminal before compiling


# default
COMPVER = 'serial'
mpi:       COMPVER = 'mpi'
rwthmpi:       COMPVER = 'rwthmpi'
rwht:       COMPVER = 'rwth'
mpidebug:  COMPVER = 'mpidebug'
openmp:    COMPVER = 'openmp'
gfort:     COMPVER = 'gfort'
pg:        COMPVER = 'pg'
hybrid:    COMPVER = 'hybrid'
rwthhybrid:    COMPVER = 'rwthhybrid'
jureca   : COMPVER='jureca'
hazelhen   : COMPVER='hazelhen'

####################################################
# Objects to compile                               #
####################################################
OBJECTS_TYPE  =   nrtype.o  type_gmatbulk.o type_tmat.o type_gmat.o type_gmatonsite.o type_gauntshape.o type_gauntcoeff.o type_corestate.o type_gref.o type_tbcluster.o type_density.o type_shapefun.o type_cell.o type_cellnew.o type_config.o type_energyparts.o type_wavefunction.o type_ldau.o type_cellorbit.o mod_version.o mod_version_info.o

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
LLIBS = -mkl
# Aachen cluster
### for sequential libraries (no hybrid parallelisation) use this
LLIBSRWTH = $(FLAGS_MKL_LINKER) #-lmkl_intel_lp64 -lmkl_sequential -lmkl_core

# openmp libs
LLIBSOPENMP = -L/usr/local/intel/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
# Aachen cluster
LLIBSOPENMPRWTH = $(FLAGS_MKL_LINKER) -lmkl_intel_lp64 
### link parallel omp libraries for hybrid parallelisation
LLIBSOPENMPRWTH+=  -lmkl_intel_thread -liomp5


# Special libraries for gfortran:
LLIBSGFORT = -L/usr/local/nag/lib -lnag -llapack_ifort -lblas_ifort  -L/usr/local/intel/current/lib/intel64 -lifcore -liomp5 -lintlc -limf -lirc -lpthread -Wl,-rpath,/usr/local/intel/current/lib/intel64

# MacOSX libraries:
LLIBSMAC  =  -mkl #-L/usr/local/Library/Frameworks/Intel_MKL.framework/Libraries/32 -lmkl_intel \
              -lmkl_intel_thread -lmkl_core -lguide -lpthread -i_dynamic

# JURECA:
JURECALIBS = -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread

# HazelHen:
HAZELHENLIBS = -mkl

# set library paths according to compile option
hybrid: LLIBS=$(LLIBSOPENMP)
rwthhybrid: LLIBS=$(LLIBSOPENMPRWTH)
rwth: LLIBS=$(LLIBSRWTH)
openmp: LLIBS=$(LLIBSOPENMP)
gfort: LLIBS=$(LLIBSGFORT)
rwthmpi: LLIBS=$(LLIBSRWTH)
mac: LLIBS=$(LLIBSMAC)
jureca: LLIBS=$(JURECALIBS)
hazelhen: LLIBS=$(HAZELHENLIBS)

####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP=
# change to if flag is matched
mpidebug:      CPP = -D MPI
mpitrace:      CPP = -D MPI
mpi:           CPP = -D MPI
rwthmpi:       CPP = -D MPI
gfort:         CPP = -D GFORT
hybrid:        CPP = -D MPI -D CPP_HYBRID
rwthhybrid:    CPP = -D MPI -D CPP_HYBRID
jureca:        CPP = -D MPI
hazelhen:      CPP = -D MPI

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
mpidebug: FFLAGS = -r8 -CB -check all -check uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -g
hybrid:   FFLAGS = -O2 -qopenmp
rwthhybrid:   FFLAGS = -O2 -qopenmp
####################################################
# Linking                                          #
####################################################
FILENAME = kkrflex.exe
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
.PHONY: getversion
.PHONY: print_rwthwarn
.PHONY: print_hazelwarn

getversion:
	@echo "Wrinting Version file: "; git describe;\
	echo $(COMPVER)>compver; \
	echo $(FFLAGS) $(OPENMPFLAGS)>compflag; \
	echo $(LLIBS) >complib; \
	$(SOURCEDIR)/getversion.sh

print_rwthwarn:
	@echo "****************************************\n"; \
	echo "!!!Warning!!!"; \
	echo "run 'module switch openmpi intelmpi'"; \
	echo "in the terminal before compiling and running the code!\n"; \
	echo "****************************************\n"

print_hazelwarn:
	@echo "****************************************\n"; \
	echo "!!!Warning!!!"; \
	echo "run 'module swap PrgEnv-cray PrgEnv-intel'"; \
	echo "and 'module unload cray-libsci'"; \
	echo "in the terminal before compiling and running the code!\n"; \
	echo "****************************************\n"


compile: getversion $(OBJECTS)
	$(FC)     -o $(FILENAME) $(OBJECTS2)  $(LLIBS)

# compile code
cbg:        compile
hybrid:     compile
openmp:     compile
debug:      compile
mpi:        compile
mpitrace:   compile
mpidebug:   compile
rwthhybrid: print_rwthwarn compile
rwth:       print_rwthwarn compile
rwthmpi:    print_rwthwarn compile
pg:         compile
gfort:      compile
mac:        compile
jureca:     compile
hazelhen:   print_hazelwarn compile

install: clean mpi

%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.F90
	$(FC) -c $(FFLAGS) $(CPP)  $<
%.o: $(SOURCEDIR)/%.F
	$(FC) -c $(FFLAGS) $(CPP)  $<

$(SOURCEDIR)/mod_version.o  : $(SOURCEDIR)/version.f90 ; \
                        $(FC) $(FFLAGS) $(SOURCEDIR)/version.f90 -c -o $@

$(SOURCEDIR)/mod_version_info.o  : $(SOURCEDIR)/version_info.F90 ; \
                        $(FC) $(FFLAGS) $(CPP) $(SOURCEDIR)/version_info.F90 -c -o $@


####################################################
# Clean up flags                                   #
####################################################

clean:
	rm -f *.mod *.o ../$(FILENAME)
cleanxc:
	rm -f exchangecorrel.o

