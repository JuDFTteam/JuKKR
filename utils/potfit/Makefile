PROG =	potfit

SRCS =	potfit4tbkkr.f coreshell.f  rmtcalc.f seekreal.f spline.f \
	tbkkr_pformat.f

OBJS =	potfit4tbkkr.o coreshell.o  rmtcalc.o seekreal.o spline.o \
	tbkkr_pformat.o

LIBS =	

CC = cc
CFLAGS = -O
FC = g77
#FFLAGS = -O
F90 = f90
F90FLAGS = -O
LDFLAGS = -s

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f  $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

