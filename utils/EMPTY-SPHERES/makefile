FC = ifort
LD = $(FC)
FFLAGS = -O2  
#FFLAGS = -O0 -CB -check all -traceback -g
#LIBS = -llapack
#LIBS = -llapack_ifort -lblas_ifort  -L/usr/local/intel/lib/intel64 -lifcore -limf  -Wl,-rpath,/usr/local/intel/lib/intel64
LIBS = -mkl

SOURCEDIR = .

NAME = empty_spheres.exe

files = stdtypes.o mtprng.o ioinput.o empty_spheres.o 

#.f.o:
#	${FC} ${FFLAGS}  -c $<
#.suffixes: .f .o
#.suffixes: .f90 .mod


es:	${SOURCEDIR}/${files}
	${LD} ${files} ${LIBS}  -o empty_spheres.exe

#cbg: $(OBJECTS)
#        $(FC) -o $(NAME) $(OBJECTS) $(LIBS)
#mac: $(OBJECTS)
#        $(FC) -o $(NAME) $(OBJECTS) $(LIBSMAC)


%.o: $(SOURCEDIR)/%.f90
	$(FC) -c $(FFLAGS)  $<
%.o: $(SOURCEDIR)/%.f
	$(FC) -c $(FFLAGS)  $<


