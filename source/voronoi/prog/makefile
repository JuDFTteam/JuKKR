dir      = .
#FC =g77
FC =ifort
#FFLAGS= -O0 
FFLAGS=-v -O -r8 -heap-arrays -mcmodel=large
# For debugging, checking array bounds, etc (makes code slower):
# FFLAGS=-v -O0 -r8 -CB -check all -traceback -g -heap-arrays -mcmodel=large
#LD=g77
LD =ifort
LIBS=-L../ -lblas
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


%.o: %.f
	$(FC) -c $(FFLAGS)  $<
%.o: %.f90
	$(FC) -c $(FFLAGS)  $<
voron:	${dir}/${voronoi} ${dir}/${newsubs} ${dir}/${newme}
	${LD} ${voronoi} ${newsubs} ${T} ${LIBS}  -o ../voronoi.exe

${voronoi}:	inc.geometry

${newsubs}:	inc.geometry

${newme}:	inc.geometry



