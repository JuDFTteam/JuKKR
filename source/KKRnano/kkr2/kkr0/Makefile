# Makefile for kkr0

# default values
PLATFORM = ifort
TYPE = debug

PROGRAM = kkr0.exe

# Path to put object files
# and module files
BUILDDIR = $(HOME)/build/kkr0

ifeq ($(PLATFORM),gfortran)
FC = gfortran
FC90 = gfortran
FCFLAGS = -g -M $(BUILDDIR)
FC90FLAGS = -Wall -std=f95 -fimplicit-none
FCFLAGS += -I $(BUILDDIR)
LDFLAGS = -llapack_gfortran -lblas_gfortran
endif

ifeq ($(PLATFORM),ifort)
FC = ifort
FC90 = ifort
# -module <path> specifies where to put .mod files
FCFLAGS = -warn -debug -traceback -gen-interfaces -warn interfaces -module $(BUILDDIR) -DTASKLOCAL_FILES
FC90FLAGS = -warn
FCFLAGS += -I $(BUILDDIR) 
LDFLAGS =-L/ -llapack_ifort -lblas_ifort
endif

ifeq ($(PLATFORM),rwth)
FC = ifort
FC90 = ifort
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -warn -debug -traceback -gen-interfaces -warn interfaces -module $(BUILDDIR)
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR) -DTASKLOCAL_FILES
LDFLAGS = -O2 -L/opt/intel/Compiler/12.1/3.293/rwthlnk/mkl/lib/intel64 -L/opt/intel/Compiler/12.1/3.293/rwthlnk/mkl/lib/ia32 -lrwthmkl -lpthread
endif

ifeq ($(PLATFORM),JUGENE)
FC = xlf_r
FC90 = xlf_r
FCLAGS = -O2 -qmoddir=$(BUILDDIR) -I $(BUILDDIR)
FC90FLAGS =
LDFLAGS = -lessl
endif

ifeq ($(PLATFORM),JUQUEEN)
FC = xlf_r
FC90 = xlf_r
FCFLAGS = -q64 -O2 -qmoddir=$(BUILDDIR) -I$(BUILDDIR)
FC90FLAGS =
LDFLAGS = -L/opt/ibmmath/essl/5.1/lib64 -lesslbg
endif

# Paths to look for files
#DIRS = . dir1 dir2
DIRS = . ../config_reader ../lattice ../kkr0common ../datastructures
#also add BUILDDIR to VPATH to look if object files have changed
VPATH = $(DIRS) $(BUILDDIR)

SRCS = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f))
SRCS90 = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f90))
SRCSFPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F))
SRCS90FPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F90))

# notdir extracts only filename
OBJS =  $(notdir ${SRCS:.f=.o})
OBJS += $(notdir ${SRCS90:.f90=.o})
OBJS += $(notdir ${SRCSFPP:.F=.o})
OBJS += $(notdir ${SRCS90FPP:.F90=.o})

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
	@echo $(PLATFORM)
	@echo $(TYPE)
	
#Module dependencies

config_reader.o: config_reader_dictionary.o
main0.o: config_reader.o InputParams_mod.o Main2Arrays_mod.o DimParams_mod.o
Main2Arrays_mod.o: DimParams_mod.o
DimParams_mod.o: config_reader.o

BasisAtom_mod.o: CellData_mod.o PotentialData_mod.o AtomicCoreData_mod.o RadialMeshData_mod.o
CellData_mod.o: ShapefunData_mod.o

startb1_wrapper.o: CellData_mod.o RadialMeshData_mod.o BasisAtom_mod.o
startb1.o: BasisAtom_mod.o
