# Makefile for kkr0

# default values
PLATFORM = gfortran
TYPE = debug

PROGRAM = kkr0.exe

# Path to put object files
# and module files
BUILDDIR = ~/build/kkr0

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
FCFLAGS = -debug -module $(BUILDDIR)
FC90FLAGS = -warn
FCFLAGS += -I $(BUILDDIR) 
LDFLAGS =-L/ -llapack_ifort -lblas_ifort
endif

# Paths to look for files
#DIRS = . dir1 dir2
DIRS = .
VPATH = $(DIRS)

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
	rm $(BUILDDIR)/*.o 
	rm $(BUILDDIR)/*.mod
	rm *.exe

.PHONY: test
test:
	@echo $(SRCS)
	@echo $(SRCS90)
	@echo $(OBJS)
	@echo $(VPATH)
	
#Module dependencies

inc_p_wrapper_module.o: inc.p inc.cls
inc_p_replace_vars.o: inc.p inc.cls
inc_p_replace.o: inc_p_replace_vars.o
