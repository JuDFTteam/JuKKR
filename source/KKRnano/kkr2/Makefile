# Makefile for kkr2

# default values
#PLATFORM = ifort_IFF_NOOMP
PLATFORM = ifort_IFF
TYPE = debug
#TYPE = nodebug

PROGRAM = kkr2.exe

# Path to put object files
# and module files
BUILDDIR = $(HOME)/build/kkr2


ifeq ($(PLATFORM),ifort_IFF)
FC = mpif77
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -O2 -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/em64t'  #The paths have changed
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/intel64'
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib'
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib/intel64'
LDFLAGS += -lfftw3
endif

ifeq ($(PLATFORM),JUROPA)
FC = mpif77
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces -I/usr/local/fftw/3.2.1/include
FC90FLAGS =
FCFLAGS += -I$(BUILDDIR)
LDFLAGS = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -L/usr/local/fftw/3.2.1/lib
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS += -lfftw3 -lmkl -lguide
endif


ifeq ($(PLATFORM),ifort_IFF_NOOMP)
FC = mpif90 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS +=  -llapack_ifort -lblas_ifort
LDFLAGS += -lfftw3
endif

ifeq ($(TYPE),debug)
FCFLAGS += -g -traceback -check all
endif

# ============================= JUGENE ========================================
ifeq ($(PLATFORM),JUGENE)

ifeq ($(USETOOL),scalasca)
  FC = scalasca -instrument -user #end of line
  TYPE = debug
else
  FC =
endif


FC += mpixlf77_r
FC90 = mpixlf90_r
FCFLAGS = -qsmp=omp -I/bgsys/local/fftw3/include
#-qsmallstack=dynlenonheap
FCFLAGS += -qmoddir=$(BUILDDIR) -I$(BUILDDIR)
LDFLAGS = -L/bgsys/local/lib -lesslbg -I/bgsys/local/fftw3/include -L/bgsys/local/fftw3/lib -lfftw3

ifeq ($(TYPE),debug)
  FCFLAGS += -g -qnosave -qinitauto=7FF7FFFF
#-qcheck -qflttrap=inv:zero:nanq:ov:und:enable -qarch=450 
else
  FCFLAGS += -O3 -qipa
endif

endif




# Paths to look for files
#DIRS = . dir1 dir2
DIRS = . DebugHelpers IterativeSolver Mixing XC quadrature madelung EnergyMesh SingleSiteRef

#also add BUILDDIR to VPATH to look if object files have changed
VPATH = $(DIRS) $(BUILDDIR)

SRCS = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f))
SRCS90 = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f90))

# notdir extracts only filename
OBJS =  $(notdir ${SRCS:.f=.o})
OBJS += $(notdir ${SRCS90:.f90=.o})

.PHONY: all
all:	$(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -o $(PROGRAM) $(addprefix $(BUILDDIR)/,$(OBJS)) $(LDFLAGS) 

%.o: %.f
	$(FC) $(FCFLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.f90
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F
	$(FC) $(FCFLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F90
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@

.PHONY: clean
clean:
	rm -f $(BUILDDIR)/*.o 
	rm -f $(BUILDDIR)/*.mod
	rm -f $(BUILDDIR)/*__genmod.f90
	rm -f *.exe

.PHONY: test
test:
	@echo $(SRCS)
	@echo $(SRCS90)
	@echo $(OBJS)
	@echo $(VPATH)
	
#Module dependencies

main2.o: main2_arrays_mod.o main2_aux_mod.o muffin_tin_zero_mod.o EnergyMesh_mod.o common_mpi.o common_testc.o common_optc.o lloyds_formula_mod.o KKRSelfConsistency_mod.o
main2_arrays_mod.o: error_messages_mod.o
kkrmat01.o: lloyds_formula_mod.o
lloyds_formula_mod.o: kkr_helpers_mod.o
lattice3d.o: common_mpi.o
test.o: common_testc.o
opt.o: common_optc.o

tbref.o: SingleSiteRef_mod.o
SingleSiteRef_mod.o: kkr_helpers_mod.o



