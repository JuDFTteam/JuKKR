# Makefile for kkr2

# default values
PLATFORM = ifort
TYPE = nodebug

PROGRAM = kkr.exe

# Path to put object files
# and module files
BUILDDIR = $(HOME)/build

EXTRA_FLAGS =

# ========= IFF Workstations ============

ifeq ($(PLATFORM),ifort)
FC = mpif77 -warn all
FC90 = mpif90 -warn all
# -module <path> specifies where to put .mod files
FCFLAGS  = -mkl
FCFLAGS += -module $(BUILDDIR)
FCFLAGS += -I $(BUILDDIR)
# FCFLAGS += -openmp
FCFLAGS += -gen-interfaces -warn interfaces
FCFLAGS += -D TASKLOCAL_FILES 
FCFLAGS += -D USE_VOROWEIGHTS
# FCFLAGS += -D DEBUG1
FC90FLAGS =
# FCFLAGS += -C -traceback -CB
# FCFLAGS += -O0 -g
# LDFLAGS  = -O0 -g
FCFLAGS += -O3  -g 
LDFLAGS  = -O3  -g 
# LDFLAGS += -C -traceback
# FCFLAGS += -ipo
# LDFLAGS += -ipo
ifeq ($(TYPE),debug)
FCFLAGS += -g
else
FCFLAGS += -D NOLOGGING
endif
# set run paths to look for shared libraries (-Wl means pass to linker)
# LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/intel64'
# LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib/intel64'
endif

ifeq ($(PLATFORM),ifort_IFF)
FC = mpif77 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces -DTASKLOCAL_FILES -DUSE_VOROWEIGHTS
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS = -O3
ifeq ($(TYPE),debug)
FCFLAGS += -g
endif
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/intel64'
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib/intel64'
endif


ifeq ($(PLATFORM),ifort_IFF_NOOMP)
FC = mpif90 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS +=  -llapack_ifort -lblas_ifort
LDFLAGS += -lfftw3
endif

# ================== JURECA ======================


# ========= HPC Cluster - RWTH Aachen ============
ifeq ($(PLATFORM),rwth)
FC = $(MPIF77) -warn
FC90 = $(MPIFC) -warn
# -module <path> specifies where to put .mod files
FCFLAGS = $(FLAGS_OPENMP) -module $(BUILDDIR) -DUSE_VOROWEIGHTS
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I$(BUILDDIR) -DTASKLOCAL_FILES
LDFLAGS = -O3 $(FLAGS_MATH_LINKER)
# -lpthread
endif


# ========= MAC ============ no openmp
ifeq ($(PLATFORM),MAC)
FC = mpif77 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces -DUSE_VOROWEIGHTS
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/'
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/compiler/lib/intel64'
LDFLAGS = -L/usr/local/lib -llapack -lblas
LDFLAGS += -lfftw3
endif


# ============================= JUQUEEN========================================
ifeq ($(PLATFORM),JUQUEEN)

ifeq ($(USETOOL),scalasca)
FC = scalasca -instrument -user #end of line
FC90 = scalasca -instrument -user #end of line
TYPE = nodebug
else
FC =
FC90 = 
endif

FC += mpixlf77_r
FC90 += mpixlf90_r
#FCFLAGS = -qsmp=omp -I/bgsys/local/fftw3/3.3.2/fftw_g/include/ 
FCFLAGS = -qsmp=omp:noauto -qnosave -qdirective -qxlf2003=polymorphic -I/bgsys/local/fftw3/3.3.2/fftw_g/include/ 
FCFLAGS += -qmoddir=$(BUILDDIR) -I$(BUILDDIR) -WF,-DUSE_VOROWEIGHTS
#LDFLAGS = -L/opt/ibmmath/essl/5.1/lib64 -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lesslsmpbg -lfftw3 
LDFLAGS = -L/opt/ibmmath/essl/5.1/lib64 -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lesslbg -lfftw3 
#-lfftw3

ifeq ($(TYPE),debug)
FCFLAGS += -q64 -O3 -qstrict -g -qnosave -C -qinitauto=7FF7FFFF -WF,-DDEBUG1
#-qcheck -qflttrap=inv:zero:nanq:ov:und:enable -qarch=450 
else
FCFLAGS += -q64 -O3 -qstrict -WF,-DNOLOGGING
# -qipa
endif

ifeq ($(USETOOL),scalasca)
FCFLAGS += -g
endif

endif

#==============================================================================

# Paths to look for files
DIRS = . DebugHelpers IterativeSolver XC madelung \
         LDAU atomiccore \
         xccouplings energy force lattice KKRnanoParallel \
         MixingNew datastructures \
         shapefun IterativeSolverBCP debug_morgan  

#also add BUILDDIR to VPATH to look if object files have changed
VPATH = $(DIRS) $(BUILDDIR)

SRCS = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f))
SRCS90 = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f90))
SRCSFPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F))
SRCS90FPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F90))

# notdir extracts only filename
OBJS =  $(notdir ${SRCS:.f=.o})
OBJS += $(notdir ${SRCS90:.f90=.o})
OBJS += $(notdir ${SRCSFPP:.F=.o})
OBJS += $(notdir ${SRCS90FPP:.F90=.o})

.PHONY: all
all:	$(PROGRAM)

$(PROGRAM): $(OBJS)
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -o $(PROGRAM) $(addprefix $(BUILDDIR)/,$(OBJS)) $(LDFLAGS) 

ifeq ($(PLATFORM),ifort_IFF)
EXTRA_FLAGS = -heap-arrays
# these files need stack allocated automatic arrays
vbrmv_mat_mod.o: EXTRA_FLAGS =
sprszmm.o: EXTRA_FLAGS =
appblckcirc.o: EXTRA_FLAGS =
bcpwupper.o: EXTRA_FLAGS =
# TODO: unclear why this module does not work with -heap-arrays
ScatteringCalculation_mod.o: EXTRA_FLAGS =
# disable pointer runtime check only for the following file:
TFQMRSolver_mod.o: EXTRA_FLAGS = -check nopointer
endif

%.o: %.f
	$(FC) $(FCFLAGS) $(EXTRA_FLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.f90
	$(FC90) $(FCFLAGS) $(EXTRA_FLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F
	$(FC) $(FCFLAGS) $(EXTRA_FLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F90
	$(FC90) $(FCFLAGS) $(FC90FLAGS) $(EXTRA_FLAGS) -c $< -o $(BUILDDIR)/$@

.PHONY: clean
clean:
	rm -f $(BUILDDIR)/*.o 
	rm -f $(BUILDDIR)/*.mod
	rm -f $(BUILDDIR)/*__genmod.f90
	rm -f *.exe

.PHONY: test
test:
	@echo $(SRCS)
	@echo $(SRCS90)
	@echo $(OBJS)
	@echo $(VPATH)

.PHONY: depend 
depend:
		$(SRCS90)
		makedepend
	
#======================== Module dependencies ========================================

KKRnano.o: Logging_mod.o KKRzero_mod.o PotentialConverter_mod.o KKRnanoParallel_mod.o BasisAtom_mod.o wrappers_mod.o RadialMeshData_mod.o main2_aux_mod.o ScatteringCalculation_mod.o Main2Arrays_mod.o KKRnano_Comm_mod.o ProcessKKRresults_mod.o InputParams_mod.o EBalanceHandler_mod.o LDAUData_mod.o TimerMpi_mod.o EnergyMesh_mod.o DimParams_mod.o CalculationData_mod.o
TruncationZone_mod.o: Main2Arrays_mod.o
kkrmat01_new.o: arraytest2_mod.o Logging_mod.o dlke0_smat_mod.o KKROperator_mod.o one_sided_commZ_mod.o ClusterInfo_mod.o TFQMRSolver_mod.o MultScatData_mod.o jij_calc_mod.o TEST_lcutoff_mod.o SparseMatrixDescription_mod.o BCPOperator_mod.o SolverOptions_mod.o SolverStats_mod.o InitialGuess_mod.o fillKKRMatrix_mod.o
tbref.o: SingleSiteRef_mod.o
ProcessKKRresults_mod.o: arraytest2_mod.o Logging_mod.o lloyds_formula_mod.o wrappers_mod.o brydbm_new_com_mod.o Main2Arrays_mod.o CellData_mod.o EnergyResults_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o CalculationData_mod.o KKRnanoParallel_mod.o muffin_tin_zero_mod.o RadialMeshData_mod.o MadelungPotential_mod.o broyden_kkr_mod.o EnergyMesh_mod.o BasisAtom_mod.o DensityResults_mod.o lloyd0_new_mod.o KKRresults_mod.o InputParams_mod.o TimerMpi_mod.o DimParams_mod.o BroydenData_mod.o NearField_calc_mod.o debug_morgan_mod.o total_energy_mod.o
NearField_calc_mod.o: RadialMeshData_mod.o NearField_com_mod.o BasisAtom_mod.o KKRnanoParallel_mod.o CalculationData_mod.o
NearField_mod.o: Harmonics_mod.o MadelungCalculator_mod.o
KKRnano_Comm_mod.o: jij_calc_mod.o KKRnanoParallel_mod.o comm_patternsD_mod.o comm_patternsZ_mod.o comm_patternsC_mod.o
InputParams_mod.o: ConfigReader_mod.o
PolygonFaces_mod.o: Constants_mod.o
ShapeGeometryHelpers_mod.o: Exceptions_mod.o Constants_mod.o
ShapeStandardMesh_mod.o: Exceptions_mod.o Constants_mod.o
wrappers_mod.o: ValenceDensity_mod.o SingleSite_mod.o BasisAtom_mod.o RadialMeshData_mod.o CellData_mod.o EnergyMesh_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o
DensityResults_mod.o: DimParams_mod.o
BCPOperator_mod.o: OperatorT_mod.o ClusterInfo_mod.o SolverOptions_mod.o
CalculationData_mod.o: TruncationZone_mod.o Main2Arrays_mod.o CellData_mod.o EnergyResults_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o InitialGuess_mod.o KKRnanoParallel_mod.o MadelungCalculator_mod.o RadialMeshData_mod.o ConstructShapes_mod.o InterpolateBasisAtom_mod.o BasisAtom_mod.o ShapefunData_mod.o DensityResults_mod.o KKRresults_mod.o InputParams_mod.o RefCluster_mod.o DimParams_mod.o ClusterInfo_mod.o BroydenData_mod.o JijData_mod.o TEST_lcutoff_mod.o
NearField_com_mod.o: Logging_mod.o one_sided_commD_mod.o NearField_mod.o MadelungCalculator_mod.o NearField_kkr_mod.o
ShapeIntegrationHelpers_mod.o: Constants_mod.o
ConstructShapes_mod.o: RefCluster_mod.o ShapefunData_mod.o Voronoi_mod.o ShapeFunctions_mod.o
EnergyMesh_mod.o: KKRnanoParallel_mod.o EnergyMeshHelpers_mod.o
InterpolateBasisAtom_mod.o: RadialMeshData_mod.o BasisAtom_mod.o PotentialData_mod.o AtomicCoreData_mod.o
SingleSiteRef_mod.o: Harmonics_mod.o kkr_helpers_mod.o SingleSiteHelpers_mod.o Exceptions_mod.o
TFQMRSolver_mod.o: Solver_mod.o OperatorT_mod.o mminvmod_oop_mod.o SolverStats_mod.o
kloopz1_mod.o: KKROperator_mod.o ClusterInfo_mod.o BCPOperator_mod.o MultScatData_mod.o kkrmat01_new.o jij_calc_mod.o TEST_lcutoff_mod.o TFQMRSolver_mod.o SolverOptions_mod.o InitialGuess_mod.o
ShapeCriticalPoints_mod.o: Constants_mod.o ShapeGeometryHelpers_mod.o  PolygonFaces_mod.o
MadelungPotential_mod.o: BasisAtom_mod.o MadelungCalculator_mod.o RadialMeshData_mod.o DensityResults_mod.o EnergyResults_mod.o CalculationData_mod.o
EBalanceHandler_mod.o: KKRnanoParallel_mod.o Exceptions_mod.o
ShapeIntegration_mod.o: ShapeIntegrationHelpers_mod.o Constants_mod.o  PolygonFaces_mod.o
BasisAtom_mod.o: Exceptions_mod.o RadialMeshData_mod.o CellData_mod.o PotentialData_mod.o AtomicCoreData_mod.o
mminvmod_oop_mod.o: Logging_mod.o OperatorT_mod.o SolverStats_mod.o
KKRresults_mod.o: DimParams_mod.o
DimParams_mod.o: ConfigReader_mod.o Exceptions_mod.o
fillKKRMatrix_mod.o: SparseMatrixDescription_mod.o
Main2Arrays_mod.o: DimParams_mod.o
CellData_mod.o: ShapefunData_mod.o
main2_aux_mod.o: TimerMpi_mod.o
NearField_kkr_mod.o: NearField_mod.o
jij_calc_mod.o: JijData_mod.o
broyden_kkr_mod.o: brydbm_new_com_mod.o broyden_second_mod.o BroydenData_mod.o CalculationData_mod.o
ShapeFunctions_mod.o: PolygonFaces_mod.o  ShapeCriticalPoints_mod.o Constants_mod.o ShapeIntegration_mod.o ShapeStandardMesh_mod.o
ConfigReader_mod.o: ConfigReaderDictionary_mod.o
KKROperator_mod.o: MultScatData_mod.o OperatorT_mod.o ClusterInfo_mod.o SparseMatrixDescription_mod.o vbrmv_mat_mod.o
lloyd0_new_mod.o: KKRnanoParallel_mod.o BasisAtom_mod.o RadialMeshData_mod.o CellData_mod.o EnergyMesh_mod.o LDAUData_mod.o GauntCoefficients_mod.o
vbrmv_mat_mod.o: SparseMatrixDescription_mod.o
ClusterInfo_mod.o: TruncationZone_mod.o RefCluster_mod.o one_sided_commI_mod.o
ScatteringCalculation_mod.o: SingleSiteRef_mod.o arraytest2_mod.o Logging_mod.o InputParams_mod.o TruncationZone_mod.o wrappers_mod.o Main2Arrays_mod.o KKRnano_Comm_mod.o jij_calc_mod.o BCPOperator_mod.o LDAUData_mod.o one_sided_commZ_mod.o GauntCoefficients_mod.o InitialGuess_mod.o CalculationData_mod.o KKRnanoParallel_mod.o EBalanceHandler_mod.o EnergyMesh_mod.o SolverOptions_mod.o KKROperator_mod.o BasisAtom_mod.o KKRresults_mod.o TFQMRSolver_mod.o RefCluster_mod.o TimerMpi_mod.o DimParams_mod.o kloopz1_mod.o ClusterInfo_mod.o JijData_mod.o TEST_lcutoff_mod.o MultScatData_mod.o
TEST_lcutoff_mod.o: TruncationZone_mod.o Main2Arrays_mod.o lcutoff_mod.o
MultScatData_mod.o: TEST_lcutoff_mod.o ClusterInfo_mod.o fillKKRMatrix_mod.o SparseMatrixDescription_mod.o
total_energy_mod.o: Quadrature_mod.o BasisAtom_mod.o RadialMeshData_mod.o CellData_mod.o ShapeGauntCoefficients_mod.o
vxcgga.o: XCFunctionals_mod.o
KKRzero_mod.o: Exceptions_mod.o Startb1_mod.o EnergyMeshHelpers_mod.o Lattice_mod.o BrillouinZone_mod.o ConfigReader_mod.o ConfigReaderDictionary_mod.o InputParams_mod.o Main2Arrays_mod.o DimParams_mod.o BasisAtom_mod.o CellData_mod.o PotentialData_mod.o AtomicCoreData_mod.o RadialMeshData_mod.o ShapefunData_mod.o read_formatted_mod.o read_formatted_shapefun_mod.o ldauinfo_read.o ldaustart.o
SingleSite_mod.o: SingleSiteHelpers_mod.o
ldauphi.o: SingleSite_mod.o
ValenceDensity_mod.o: Quadrature_mod.o SingleSite_mod.o SingleSiteHelpers_mod.o
lloyd0_new_mod.o: ValenceDensity_mod.o
Voronoi_mod.o: Exceptions_mod.o Sorting_mod.o
Exceptions_mod.o: Errors_mod.o Warnings_mod.o StringHelpers_mod.o
Errors_mod.o: Warnings_mod.o
ldauovrlp.o: Quadrature_mod.o
ldauinit.o: Quadrature_mod.o
ldauphi.o: Quadrature_mod.o
force.o: Quadrature_mod.o
forceh.o: Quadrature_mod.o
forcxc.o: Quadrature_mod.o
corel.o: Quadrature_mod.o
mtzero_new.o: Quadrature_mod.o
epotinb_new.o: Quadrature_mod.o
vxclm.o: Quadrature_mod.o
rhototb_new.o: Quadrature_mod.o
vxcgga.o: Quadrature_mod.o
GauntCoefficients_mod.o: Harmonics_mod.o
ShapeGauntCoefficients_mod.o: Harmonics_mod.o
MadelungCalculator_mod.o: Harmonics_mod.o Lattice_mod.o
sphere_gga.o: Harmonics_mod.o
sphere_nogga.o: Harmonics_mod.o
strmat_new.o: Harmonics_mod.o
debug_morgan_mod.o: Harmonics_mod.o
BrillouinZone_mod.o: Symmetry_mod.o VectorMath_mod.o
Lattice_mod.o: VectorMath_mod.o Constants_mod.o
Symmetry_mod.o: VectorMath_mod.o Exceptions_mod.o
GaussWeights_mod.o: Exceptions_mod.o
EnergyMeshHelpers_mod.o: GaussWeights_mod.o Constants_mod.o
jij_calc_mod.o: Sorting_mod.o
PotentialConverter_mod.o: DimParams_mod.o BasisAtom_mod.o RadialMeshData_mod.o PotentialData_mod.o
Startb1_mod.o: RadialMeshData_mod.o BasisAtom_mod.o read_formatted_mod.o read_formatted_shapefun_mod.o
strmat_new.o: Exceptions_mod.o
one_sided_commZ_mod.o: ChunkIndex_mod.o 
one_sided_commI_mod.o: ChunkIndex_mod.o 
one_sided_commD_mod.o: ChunkIndex_mod.o 
ClusterInfo_mod.o: one_sided_commI_mod.o
RadialMeshData_mod.o: Exceptions_mod.o


# DO NOT DELETE
