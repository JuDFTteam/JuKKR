# Makefile for kkr2

# default values
PLATFORM = ifort_IFF
TYPE = debug

PROGRAM = kkr2.exe

# Path to put object files
# and module files
BUILDDIR = $(HOME)/build/kkr2


ifeq ($(PLATFORM),ifort_IFF)
FC = mpif77
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -O2 -module $(BUILDDIR)
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/em64t'
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib'
LDFLAGS += -lfftw3
endif

# Paths to look for files
#DIRS = . dir1 dir2
DIRS = .
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

.PHONY: clean
clean:
	rm -f $(BUILDDIR)/*.o 
	rm -f $(BUILDDIR)/*.mod
	rm -f *.exe

.PHONY: test
test:
	@echo $(SRCS)
	@echo $(SRCS90)
	@echo $(OBJS)
	@echo $(VPATH)
	
#Module dependencies

# only main2 and lattice3d use common_mpi
main2.o: inc_p_wrapper_module.o common_mpi.o
lattice3d.o: common_mpi.o

inc_p_wrapper_module.o: inc.p inc.cls
kloopz1.o: inc_p_wrapper_module.o
kkrmat01.o: inc_p_wrapper_module.o
eprdist.o: inc_p_wrapper_module.o
prinvit.o: inc_p_wrapper_module.o

