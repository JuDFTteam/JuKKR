# Makefile for kkr2

# default values
#PLATFORM = ifort_IFF_NOOMP
PLATFORM = ifort_IFF
TYPE = debug
#TYPE = nodebug

PROGRAM = kkr2.exe

# Path to put object files
# and module files
BUILDDIR = $(HOME)/build/kkr2


ifeq ($(PLATFORM),ifort_IFF)
FC = mpif77 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/em64t'  #The paths have changed
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/intel64'
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib'
LDFLAGS +=  -Wl,-rpath '/usr/local/intel/lib/intel64'
LDFLAGS += -lfftw3
endif

ifeq ($(PLATFORM),rwth)
FC = mpif77 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I$(BUILDDIR) -I/usr/local/fftw/3.2.1/include -DTASKLOCAL_FILES
LDFLAGS = -O2 -L/opt/intel/Compiler/12.1/3.293/rwthlnk/mkl/lib/intel64 -L/opt/intel/Compiler/12.1/3.293/rwthlnk/mkl/lib/ia32 -lrwthmkl -lpthread
endif

# ========= MAC ============ no openmp
ifeq ($(PLATFORM),MAC)
FC = mpif77 -warn
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces
FC90FLAGS =
FCFLAGS += -I $(BUILDDIR)
LDFLAGS =
# set run paths to look for shared libraries (-Wl means pass to linker)
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/mkl/lib/'
#LDFLAGS +=  -Wl,-rpath '/usr/local/intel/compiler/lib/intel64'
LDFLAGS = -L/usr/local/lib -llapack -lblas
LDFLAGS += -lfftw3
endif

ifeq ($(PLATFORM),JUROPA)
FC = mpif77
FC90 = mpif90 -warn
# -module <path> specifies where to put .mod files
FCFLAGS = -openmp -mkl -module $(BUILDDIR)
FCFLAGS += -gen-interfaces -warn interfaces -I/usr/local/fftw/3.2.1/include
FC90FLAGS =
FCFLAGS += -I$(BUILDDIR)
LDFLAGS = -L/opt/intel/Compiler/11.0/074/mkl/lib/em64t -L/usr/local/fftw/3.2.1/lib
# set run paths to look for shared libraries (-Wl means pass to linker)
LDFLAGS += -lfftw3 -lmkl -lguide
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

ifeq ($(TYPE),debug)
FCFLAGS += -debug -g -O0 -traceback -check all -DDEBUG1
#-auto -ftrapuv -check uninit
else
FCFLAGS += -O2 -g -DNOLOGGING
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
FCFLAGS = -qsmp=omp:noauto -qnosave -qdirective -I/bgsys/local/fftw3/3.3.2/fftw_g/include/ 
FCFLAGS += -qmoddir=$(BUILDDIR) -I$(BUILDDIR)
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
#DIRS = . dir1 dir2
DIRS = . DebugHelpers IterativeSolver XC quadrature madelung \
         EnergyMesh SingleSiteRef SingleSite LDAU atomiccore rhovalence \
         xccouplings energy force harmonics lattice KKRnanoParallel \
         timing kkr0common MixingNew datastructures config_reader \
         cluster shapefun

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

$(PROGRAM): globals $(OBJS)
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -o $(PROGRAM) $(addprefix $(BUILDDIR)/,$(OBJS)) $(LDFLAGS) 

%.o: %.f
	$(FC) $(FCFLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.f90
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F
	$(FC) $(FCFLAGS) -c $< -o $(BUILDDIR)/$@
	
%.o: %.F90
	$(FC90) $(FCFLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@

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
	
# Compile these files first
globals: arraytest2_mod.o Logging_mod.o
	
#Module dependencies

main2.o: main2_aux_mod.o  EnergyMesh_mod.o \
         MadelungCalculator_mod.o  \
         KKRnanoParallel_mod.o KKRnano_Comm_mod.o \
         GauntCoefficients_mod.o ShapeGauntCoefficients_mod.o TimerMpi_mod.o \
         EBalanceHandler_mod.o \
         RadialMeshData_mod.o CellData_mod.o wrappers_mod.o BasisAtom_mod.o LDAUData_mod.o \
         JijData_mod.o BroydenData_mod.o Main2Arrays_mod.o DimParams_mod.o InputParams_mod.o \
         ScatteringCalculation_mod.o KKRresults_mod.o DensityResults_mod.o ProcessKKRresults_mod.o \
         CalculationData_mod.o
         
ProcessKKRresults_mod.o: KKRnanoParallel_mod.o lloyds_formula_mod.o main2_aux_mod.o \
            muffin_tin_zero_mod.o EnergyMesh_mod.o MadelungCalculator_mod.o \
            lloyd0_new_mod.o GauntCoefficients_mod.o ShapeGauntCoefficients_mod.o \
            RadialMeshData_mod.o CellData_mod.o BasisAtom_mod.o LDAUData_mod.o \
            TimerMpi_mod.o BroydenData_mod.o brydbm_new_com_mod.o \
            wrappers_mod.o DimParams_mod.o InputParams_mod.o Main2Arrays_mod.o \
            KKRresults_mod.o DensityResults_mod.o EnergyResults_mod.o \
            MadelungPotential_mod.o
         
ScatteringCalculation_mod.o: EnergyMesh_mod.o KKRnanoParallel_mod.o KKRnano_Comm_mod.o \
                             GauntCoefficients_mod.o TimerMpi_mod.o \
                             EBalanceHandler_mod.o TEST_lcutoff_mod.o \
                             wrappers_mod.o BasisAtom_mod.o LDAUData_mod.o \
                             JijData_mod.o Main2Arrays_mod.o DimParams_mod.o InputParams_mod.o \
                             KKRresults_mod.o CalculationData_mod.o kloopz1_mod.o \
                             one_sided_commZ_mod.o
         
kkrmat01.o: lloyds_formula_mod.o initialGuess_store_mod.o TEST_lcutoff_mod.o
lloyds_formula_mod.o: kkr_helpers_mod.o
test.o: common_testc.o
opt.o: common_optc.o
kloopz1.o: kkrmat01.o kkrmat01_new.o TEST_lcutoff_mod.o

tbref.o: SingleSiteRef_mod.o
SingleSiteRef_mod.o: kkr_helpers_mod.o

KKRnano_Comm_mod.o: KKRnanoParallel_mod.o comm_patternsZ_mod.o comm_patternsC_mod.o comm_patternsD_mod.o

mminvmod_mod.o: vbrmv_mat_mod.o SparseMatrixDescription_mod.o
vbrmv_mat_mod.o: SparseMatrixDescription_mod.o
kkrmat01_new.o: lloyds_formula_mod.o initialGuess_store_mod.o mminvmod_mod.o fillKKRMatrix_mod.o \
                dlke0_smat_mod.o TEST_lcutoff_mod.o SparseMatrixDescription_mod.o \
                one_sided_commZ_mod.o
                
fillKKRMatrix_mod.o: SparseMatrixDescription_mod.o
EBalanceHandler_mod.o: KKRnanoParallel_mod.o

#l-cutoff TEST
TEST_lcutoff_mod.o: lcutoff_mod.o TruncationZone_mod.o Main2Arrays_mod.o

BasisAtom_mod.o: CellData_mod.o PotentialData_mod.o AtomicCoreData_mod.o RadialMeshData_mod.o
CellData_mod.o: ShapefunData_mod.o

wrappers_mod.o: BasisAtom_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o LDAUData_mod.o EnergyMesh_mod.o

main2_aux_mod.o: BasisAtom_mod.o RadialMeshData_mod.o

EnergyMesh_mod.o: EnergyMeshHelpers_mod.o KKRnanoParallel_mod.o
lloyd0_new_mod.o: BasisAtom_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o \
                  LDAUData_mod.o EnergyMesh_mod.o KKRnanoParallel_mod.o

Main2Arrays_mod.o: DimParams_mod.o
KKRresults_mod.o: DimParams_mod.o
DensityResults_mod.o: DimParams_mod.o

# TODO
CalculationData_mod.o: main2_aux_mod.o  EnergyMesh_mod.o \
         MadelungCalculator_mod.o  \
         KKRnanoParallel_mod.o KKRnano_Comm_mod.o \
         GauntCoefficients_mod.o ShapeGauntCoefficients_mod.o TimerMpi_mod.o \
         EBalanceHandler_mod.o \
         RadialMeshData_mod.o CellData_mod.o wrappers_mod.o BasisAtom_mod.o LDAUData_mod.o \
         JijData_mod.o BroydenData_mod.o Main2Arrays_mod.o DimParams_mod.o InputParams_mod.o \
         KKRresults_mod.o DensityResults_mod.o \
         EnergyResults_mod.o TruncationZone_mod.o RefCluster_mod.o TEST_lcutoff_mod.o \
         ClusterInfo_mod.o InterpolateBasisAtom_mod.o ConstructShapes_mod.o
         
ClusterInfo_mod.o: RefCluster_mod.o TruncationZone_mod.o one_sided_commI_mod.o

MadelungPotential_mod.o: CalculationData_mod.o
kloopz1_mod.o: kkrmat01_new.o kkrmat01.o TEST_lcutoff_mod.o
TruncationZone_mod.o: Main2Arrays_mod.o
config_reader.o: config_reader_dictionary.o
DimParams_mod.o: config_reader.o

#Module dependencies shape-functions
ConstructShapes_mod.o: RefCluster_mod.o ShapefunData_mod.o
voronoi08.o: voronoi08_mod.o
shapewrapper.o: ShapeFunctions_mod.o

angles_common.o: shape_constants_mod.o
tetrahedra_common.o: shape_constants_mod.o
ShapeFunctions_mod.o: shape_constants_mod.o ShapeCriticalPoints_mod.o ShapeStandardMesh_mod.o ShapeIntegration_mod.o
ShapeCriticalPoints_mod.o: shape_constants_mod.o angles_common.o tetrahedra_common.o ShapeGeometryHelpers_mod.o
ShapeGeometryHelpers_mod.o: shape_constants_mod.o
ShapeIntegrationHelpers_mod.o: shape_constants_mod.o
ShapeIntegration_mod.o: shape_constants_mod.o angles_common.o tetrahedra_common.o ShapeIntegrationHelpers_mod.o
ShapeStandardMesh_mod.o: shape_constants_mod.o

broyden_kkr_mod.o: CalculationData_mod.o broyden_second_mod.o
