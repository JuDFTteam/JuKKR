####################################################
####################################################
#         voronoi@JuKKR makefile
####################################################
#
# to compile use the following options:
#
# make           -    standart compilation using
#                     ifort on a single core
# make serial    -    like make without deleting all preciously compiled files (clean)
# make debug     -    compile with debug flags
# make rwth      -    compile on RWTH cluste
# make clean     -    removes all *.o files
#
default: install
####################################################
# specify the compiler which is going to be used   #
####################################################
# default value
#FC =g77
FC =ifort


# information for version-tracking in executable:
# default
COMPVER ='serial'
# change default if compile option is hit
debug: COMPVER ='debug'
rwth: COMPVER ='rwth'


####################################################
# C preprocessor options                           #
####################################################
# Default option
CPP= -D 

####################################################
# Objects to compile                               #
####################################################
dir      = .

host = `hostname`

voronoi = crospr.o vertex3d.o analyzevert3d.o halfspace.o polyhedron08.o voronoi12.o \
distplane.o normalplane0.o sortvertices.o ioinput.o \
shape12.o

newsubs = maindriver12.o clsgen_voronoi.o clsgen_tb.o elementdatabase.o lattix12.o \
scalevec2000.o readinput12.o rrgen2000.o spatpr.o dsort.o test.o \
clustcomp_voronoi.o clustcomp_tb.o veq.o vsub.o vadd.o ddet33.o rcstop.o vmul.o scalpr.o \
opt.o cross.o mtmesh.o sname.o genpotstart12.o length.o spline.o \
splint.o ritesone.o jellstart12.o pdbplot.o \
writeshape.o readimpatoms12.o clsgenimp12.o rationalbasis.o scalevecimp.o testpanel.o \
fakeshape.o findpanels.o suggestpts.o divpanels.o

extras = equivws.o wsclasses.o  

${voronoi}:	inc.geometry

${newsubs}:	inc.geometry

%.o: %.f
	$(FC) -c $(FFLAGS)  $<
%.o: %.f90
	$(FC) -c $(FFLAGS)  $<
%.o: %.F
	$(FC) $(CPP) -c $(FFLAGS)  $<
%.o: %.F90
	$(FC) $(CPP) -c $(FFLAGS)  $<

####################################################
# Libraries and compile flags                      #
####################################################

# standard libs
# serial:
LIBS=-L../ -lblas
# RWTH cluster:
rwth: LIBS= -mkl

# standard compile flags
FFLAGS=-O2 -r8 -traceback -i8 #-heap-arrays -mcmodel=large
# For debugging, checking array bounds, etc (makes code slower):
debug: FFLAGS=-v -O0 -r8 -CB -check all -check uninit -ftrapuv -gen-interfaces -warn all -warn notruncated_source -fpe0 -debug extended -traceback -g


####################################################
# Do stuff                                         #
####################################################

.PHONY: getversion
.PHONY: clear
.PHONY: clean
.PHONY: compile
.PHONY: debug
.PHONY: rwth
.PHONY: serial
.PHONY: install

getversion:
	@echo "Wrinting Version file: "; git describe;\
	echo $(COMPVER)>compver; \
	echo $(FFLAGS) >compflag; \
	echo $(LIBS) >complib; \
	${dir}/getversion.sh

clear:
	rm -f *.o ../voronoi.exe *.mod *__genmod.*

clean:
	rm -f *.o ../voronoi.exe *.mod *__genmod.*

compile: ${dir}/getversion ${dir}/version.o ${dir}/version_info.o ${dir}/${voronoi} ${dir}/${newsubs}
	${FC} ${dir}/version.o ${dir}/version_info.o ${voronoi} ${newsubs} ${LIBS}  -o ../voronoi.exe



serial: compile
debug: compile
rwth: compile
install: clean serial


