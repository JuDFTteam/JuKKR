# Makefile for kkr0

BUILDDIR = ~/build/kkr0/
PROGRAM = kkr0.exe

FC = gfortran
#FC = ifort
FC90 = gfortran
#FC90 = ifort -warn 
FCFLAGS = -g
FC90FLAGS = -g -Wall -std=f95 -fimplicit-none
#LDFLAGS=-L/ -llapack_ifort -lblas_ifort
LDFLAGS= -llapack_gfortran -lblas_gfortran

# add -std=f95

SRCS = $(wildcard *.f)
SRCS90 = $(wildcard *.f90)
OBJS = ${SRCS:.f=.o}
OBJS90 = ${SRCS90:.f90=.o}

#MODULESRC = inc_p_replace_vars.f inc_p_replace.f
#MODULESRC90 = inputcard_reader.f90

.SUFFIXES:	.f .f90

inc_p_wrapper_module.o:
	$(FC) $(FCFLAGS) -c inc_p_wrapper_module.f

inc_replace_vars.o:
	$(FC) $(FCFLAGS) -c inc_p_replace_vars.f

inc_p_replace.o:	inc_replace_vars.o
	$(FC) $(FCFLAGS) -c inc_p_replace.f

inputcard_reader.o:	inc_p_replace.o
	$(FC90) $(FC90FLAGS) -c inputcard_reader.f90

.f90.o:
	$(FC90) $(FC90FLAGS) -c $<
	
.f.o:
	$(FC) $(FCFLAGS) -c $<

all:	$(OBJS) $(OBJS90) $(SRCS) $(SRCS90)
		$(FC) $(FCFLAGS) -o $(PROGRAM) $(OBJS) $(OBJS90) $(LDFLAGS) 

clean:
	rm *.o 
	rm *.mod
	rm *.exe
	
test:
	cd $(BUILDDIR)
	pwd
