# FC = ifort
# FFLAGS = -debug -traceback -gen-interfaces -warn all -check all
FC = gfortran 
FFLAGS = -Wall -fno-second-underscore -fPIC -O3 -funroll-loops -g
FFLAGS += -DDEBUGSHAPEFUNCTIONS
LDFLAGS =


SRCS= ShapeGeometryHelpers_mod.f  ShapeIntegrationHelpers_mod.f

SRCS90= ShapeFunctions_mod.f90 \
        ShapeCriticalPoints_mod.f90  ShapeStandardMesh_mod.f90 \
        angles_common.f90 tetrahedra_common.f90 ShapeIntegration_mod.f90

SRCS90FPP = shape_constants_mod.F90

OBJS=${SRCS:.f=.o}
OBJS+=${SRCS90:.f90=.o}
OBJS+=${SRCS90FPP:.F90=.o}

LIBRARY=libshapefunctions.a

.PHONY: all
all: $(LIBRARY)

$(LIBRARY): $(OBJS)
	ar rc $(LIBRARY) $(OBJS)

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@
%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@
	
.PHONY: clean
clean:
	rm -f *.o 
	rm -f *.mod
	rm -f *__genmod.f90
	rm -f *.a
	
.PHONY: test
test:
	@echo $(OBJS)
	
#Module dependencies
angles_common.o: shape_constants_mod.o
tetrahedra_common.o: shape_constants_mod.o
ShapeFunctions_mod.o: shape_constants_mod.o ShapeCriticalPoints_mod.o ShapeStandardMesh_mod.o ShapeIntegration_mod.o
ShapeCriticalPoints_mod.o: shape_constants_mod.o angles_common.o tetrahedra_common.o ShapeGeometryHelpers_mod.o
ShapeGeometryHelpers_mod.o: shape_constants_mod.o
ShapeIntegrationHelpers_mod.o: shape_constants_mod.o
ShapeIntegration_mod.o: shape_constants_mod.o angles_common.o tetrahedra_common.o ShapeIntegrationHelpers_mod.o
ShapeStandardMesh_mod.o: shape_constants_mod.o
