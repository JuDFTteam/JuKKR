# Makefile for kkr0

BUILDDIR = ~/build/kkr0/
PROGRAM = kkr0.exe

FC = gfortran
#FC = ifort
FC90 = gfortran
#FC90 = ifort -warn 
FCFLAGS = -g
FC90FLAGS = -g -Wall -fimplicit-none
#LDFLAGS=-L/ -llapack_ifort -lblas_ifort
LDFLAGS= -llapack_gfortran -lblas_gfortran

# add -std=f95

SRCS = $(wildcard *.f)
SRCS90 = $(wildcard *.f90)
OBJS = ${SRCS:.f=.o}
OBJS90 = ${SRCS90:.f90=.o}

.SUFFIXES:	.f .f90

.f90.o:
	$(FC90) $(FC90FLAGS) -c $<
	
.f.o:
	$(FC) $(FCFLAGS) -c $<

all:	$(OBJS90) $(OBJS)
		$(FC) $(FCFLAGS) -o $(PROGRAM) $(OBJS) $(OBJS90) $(LDFLAGS) 

clean:
	rm *.o 
	rm *.mod
	rm *.exe
	
test:
	cd $(BUILDDIR)
	pwd
