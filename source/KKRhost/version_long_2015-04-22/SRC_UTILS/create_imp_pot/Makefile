PROG =	create_imp_pot

SRCS =	create_imp_pot.f

OBJS =	create_imp_pot.o

LIBS =	

CC = cc
CFLAGS = -O
FC = ifort
FFLAGS = -O -traceback -check all
F90 = ifort
F90FLAGS = -O
LDFLAGS = 

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

