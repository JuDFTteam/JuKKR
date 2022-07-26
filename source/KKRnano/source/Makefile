# Makefile for kkr2

### default values
PLATFORM ?= ifort
### the user may specify TYPE=debug as a command line argument to make
TYPE ?= nodebug
### the user may specify SMP=openmp as a command line argument to make
SMP ?= none

PROGRAM = kkr.exe

# Path to put object files and module files
BUILDDIR = ./build

EXTRA_FLAGS ?=

### can we make use of the tfQMRgpu library?
tfQMRgpu ?= no
ifeq ($(tfQMRgpu),yes)
	EXTRA_FLAGS += -D has_tfQMRgpu
	TFQMRGPU_PATH = $(HOME)/tfQMRgpu
endif


FC90FLAGS = 
FCFLAGS =
PPFLAGS =
PPFLAGS += -D USE_VOROWEIGHTS
ifeq ($(TYPE),debug)
	PPFLAGS += -D DEBUG
endif
# PPFLAGS += -D NDEBUG

# PPFLAGS += -D BENCHMARK_tfQMR
# PPFLAGS += -D EXPORT_tfQMR_PROBLEM

# PPFLAGS += -D TRANSPOSE_TO_ROW_MAJOR
### TRANSPOSE_TO_ROW_MAJOR seems to work, however, there was a bug in col_norms about the non-linearity of DZNRM2, so we have to see which references went wrong
#PPFLAGS += -D SUPERCELL_ELECTROSTATICS=8
LDFLAGS =

# ========= IFF Workstations ============
ifeq ($(PLATFORM),ifort)
###> module load intel impi mkl
	FC = mpiifort -warn all
	FC90 = mpiifort -warn all
	# -module <path> specifies where to put .mod files
	ifeq ($(SMP),openmp)
		FCFLAGS += -qopenmp
	endif
	FCFLAGS += -module $(BUILDDIR)
	FCFLAGS += -I $(BUILDDIR)
	# PPFLAGS += -D TASKLOCAL_FILES
	FCFLAGS += -g
	LDFLAGS += -g -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl 
	ifeq ($(TYPE),debug_light)
		FCFLAGS += -O0 -g
	else ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
	else ifeq ($(TYPE),mtradii)
		FCFLAGS += -O3
		PPFLAGS += -D USE_MTRADII 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),debug_mtradii)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all
		PPFLAGS += -D USE_MTRADII -D USE_OLD_MESH -D NOLOGGING
 	else ifeq ($(TYPE),voronoi_mesh)
		FCFLAGS += -O2
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug_light)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),voronoi_mesh_print_mtradii)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH -D PRINT_MTRADII
	else
		FCFLAGS += -O2 
		PPFLAGS += -D NOLOGGING
	endif
endif

ifeq ($(PLATFORM),gfortran)
	FC = mpifort
	FC90 = mpifort 
	FCFLAGS += -B $(BUILDDIR)
	# -J <path> specifies where to put .mod files
	FCFLAGS += -J $(BUILDDIR)  
	FCFLAGS += -I $(BUILDDIR) -I ./ 
	ifeq ($(SMP),openmp)
		FCFLAGS += -fopenmp
	endif		
	# FCFLAGS += -I /usr/include ## for fftw3
	# PPFLAGS += -D HAS_FFTW3    ## for fftw3
	ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -Wall -Warray-bounds
	else
		FCFLAGS += -O3
	endif
	FC90FLAGS = -ffree-line-length-0
	PPFLAGS += -D NOLOGGING
	# FCFLAGS += -D TASKLOCAL_FILES ### breaks when TASKLOCAL_FILES is used
# 	LDFLAGS += -L /usr/local/atlas/lib -llapack -latlas
	LDFLAGS += -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 ###-lpthread
# 	# 	LDFLAGS += -lgfortran
	### -lf77blas -lcblas -latlas
endif


# ================== JURECA ======================
# The following software stages were tested and work:
# module --force purge
# module use /usr/local/software/jureca/OtherStages
# module load Stages/2016a
# module load intel-para
ifeq ($(PLATFORM),jureca)
	FC = mpif77
	FC90 = mpif90
	# -module <path> specifies where to put .mod files
  	ifeq ($(SMP),openmp)
 		FCFLAGS += -qopenmp -qopt-report
        endif
	FCFLAGS += -module $(BUILDDIR)
	FCFLAGS += -I $(BUILDDIR)
	# PPFLAGS += -D TASKLOCAL_FILES
	FCFLAGS += -g -march=core-avx2
	LDFLAGS += -g -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl 
	ifeq ($(TYPE),debug_light)
		FCFLAGS += -O0 -g
	else ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
	else ifeq ($(TYPE),mtradii)
		FCFLAGS += -O2
		PPFLAGS += -D USE_MTRADII 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),debug_mtradii)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all
		PPFLAGS += -D USE_MTRADII
 	else ifeq ($(TYPE),voronoi_mesh)
		FCFLAGS += -O2
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug_light)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else
		FCFLAGS += -O2
		PPFLAGS += -D NOLOGGING
	endif
endif

# ================== JURECA KNL ======================
# optimal setup on one node seems to be 8MPI/8OMP for 64 atoms
# code does not run on standard JURECA nodes because of AVX512
# module purge
# module use /usr/local/software/mod_environment
# module load Architecture/KNL
# module load Intel IntelMPI
ifeq ($(PLATFORM),jureca_knl)
	FC = mpiifort
	FC90 = mpiifort
	# -module <path> specifies where to put .mod files
  	ifeq ($(SMP),openmp)
 		FCFLAGS += -qopenmp -qopt-report
        endif
	FCFLAGS += -module $(BUILDDIR)
	FCFLAGS += -I $(BUILDDIR) -I${MKLROOT}/include
	# PPFLAGS += -D TASKLOCAL_FILES
	FCFLAGS += -g
	LDFLAGS += -g -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl 
	ifeq ($(TYPE),debug_light)
		FCFLAGS += -O0 -g
	else ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
	else ifeq ($(TYPE),mtradii)
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xhost 
		PPFLAGS += -D USE_MTRADII 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),debug_mtradii)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all
		PPFLAGS += -D USE_MTRADII
 	else ifeq ($(TYPE),voronoi_mesh)
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xMIC-AVX512 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug_light)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xMIC-AVX512
		PPFLAGS += -D NOLOGGING
	endif
endif
# ================== HAWK ===========================
# Usage of Intel compiler and Intel MKL seems to give best results:
# module swap PrgEnv-cray PrgEnv-intel
# module unload cray-libsci
# Remarks: - Bad GGA performance -> Use LDA!
#  	   - Home directory slow -> Use Workspace (see HLRS wiki)
#  	   - If 'seg 0' error occurs, add 'module unload craype-hugepages16M' to job file
ifeq ($(PLATFORM),hawk)
	FC = ftn
	FC90 = ftn
	# -module <path> specifies where to put .mod files
  	ifeq ($(SMP),openmp)
 		FCFLAGS += -qopenmp 
        endif
	FCFLAGS += -I $(BUILDDIR) 
	FCFLAGS += -I . 
	# PPFLAGS += -D TASKLOCAL_FILES
	FCFLAGS += -g
	LDFLAGS += -g 
	LDFLAGS += -mkl # use Intel MKL library
	ifeq ($(TYPE),debug_light)
		FCFLAGS += -O0 -g
	else ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
	else ifeq ($(TYPE),mtradii)
		FCFLAGS += -O2
		PPFLAGS += -D USE_MTRADII 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),debug_mtradii)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all
		PPFLAGS += -D USE_MTRADII
 	else ifeq ($(TYPE),voronoi_mesh)
		FCFLAGS += -O2
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
 	else ifeq ($(TYPE),voronoi_mesh_debug_light)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else
		FCFLAGS += -O2  
		PPFLAGS += -D NOLOGGING
	endif
endif



# ================== QPACE3 ======================
ifeq ($(PLATFORM),qpace3)
	FC = mpiifort
	FC90 = mpiifort
	# -module <path> specifies where to put .mod files
  	ifeq ($(SMP),openmp)
 		FCFLAGS += -qopenmp -qopt-report
        endif
	FCFLAGS += -module $(BUILDDIR)
	FCFLAGS += -I $(BUILDDIR)
	# PPFLAGS += -D TASKLOCAL_FILES
	FCFLAGS += -g
	LDFLAGS += -g -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl 
	ifeq ($(TYPE),debug_light)
		FCFLAGS += -O0 -g
	else ifeq ($(TYPE),debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
	else ifeq ($(TYPE),mtradii)
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xhost 
		PPFLAGS += -D USE_MTRADII 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else ifeq ($(TYPE),debug_mtradii)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all
		PPFLAGS += -D USE_MTRADII
 	else ifeq ($(TYPE),voronoi_mesh)
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xMIC-AVX512 
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH #-D CPP_OMP
 	else ifeq ($(TYPE),voronoi_mesh_debug)
		FCFLAGS += -O0 -g -check bounds -traceback -debug all -check all -fpe1
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH #-D CPP_OMP
 	else ifeq ($(TYPE),voronoi_mesh_debug_light)
		FCFLAGS += -O0 -g
		PPFLAGS += -D NOLOGGING -D USE_OLD_MESH
	else
		FCFLAGS += -ipo -O2 -no-prec-div -fp-model fast=2 -xMIC-AVX512
		PPFLAGS += -D NOLOGGING
	endif
endif


# ========= HPC Cluster - RWTH Aachen ============
# MIGHT BE OUTDATED!!!
ifeq ($(PLATFORM),rwth)
	FC = $(MPIF77) -warn
	FC90 = $(MPIFC) -warn
	# -module <path> specifies where to put .mod files
	FCFLAGS += $(FLAGS_OPENMP) -module $(BUILDDIR) 
	FCFLAGS += -gen-interfaces -warn interfaces
	FC90FLAGS +=
	FCFLAGS += -I$(BUILDDIR) 
	PPFLAGS += -D TASKLOCAL_FILES
	LDFLAGS += -O3 $(FLAGS_MATH_LINKER)
	# -lpthread
endif

# ============================= JUQUEEN========================================
ifeq ($(PLATFORM),juqueen)

#ifeq ($(USETOOL),scalasca)
#FC = scalasca -instrument -user #end of line
#FC90 = scalasca -instrument -user #end of line
#TYPE = nodebug
#else
FC =
FC90 = 
#endif

FC += mpixlf77_r
FC90 += mpixlf90_r
FCFLAGS = -qnosave -qdirective -qxlf2003=polymorphic  -I$(PWD) 
FCFLAGS += -qmoddir=$(BUILDDIR) -I$(BUILDDIR) -I/bgsys/local/fftw3/3.3.2/fftw_g/include/
PPFLAGS = -WF,-DUSE_VOROWEIGHTS -WF,-DCOMPUTE_tref_LOCALLY -WF,-DNOLOGGING 

ifeq ($(SMP),openmp)
	FCFLAGS += -qsmp=omp
	PPFLAGS += -WF,-DCPP_hybrid
endif

#ESSL
#LDFLAGS = -L/bgsys/local/lib -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lesslbg -lfftw3

#ESSLSMP
LDFLAGS = -L/bgsys/local/lib -lesslsmpbg -L/opt/ibmcmp/xlsmp/3.1/bglib64 -lxlsmp -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lfftw3 

ifeq ($(TYPE),debug)
	FCFLAGS += -q64 -O0 -qstrict -g -qnosave -C -qinitauto=7FF7FFFF -WF,-DDEBUG1
else ifeq ($(TYPE),scorep)
	FC = scorep mpixlf77_r	
	FC90 = scorep mpixlf90_r
	FCFLAGS += -q64 -O3 -qstrict
else ifeq ($(TYPE),voronoi_mesh)
	FCFLAGS += -q64 -O3 -qstrict
	PPFLAGS += -WF,-DUSE_OLD_MESH 
else
	FCFLAGS += -q64 -O3 -qstrict
endif

#ifeq ($(USETOOL),scalasca)
#FCFLAGS += -g
#endif

endif
# ============================= JUQUEENgfort ========================================
ifeq ($(PLATFORM),JUQUEENgfort)
	ifeq ($(USETOOL),scalasca)
		FC = scalasca -instrument -user #end of line
		FC90 = scalasca -instrument -user #end of line
		TYPE = nodebug
	else
		FC =
		FC90 = 
	endif

	FC += mpigfortran
	FC90 += mpigfortran
	# FCFLAGS += -qsmp=omp -I/bgsys/local/fftw3/3.3.2/fftw_g/include/ 
	FCFLAGS += -qsmp=omp:noauto -qnosave -qdirective -qxlf2003=polymorphic -I/bgsys/local/fftw3/3.3.2/fftw_g/include/ 
	FCFLAGS += -qmoddir=$(BUILDDIR) -I$(BUILDDIR) 
	PPFLAGS += -WF,-DUSE_VOROWEIGHTS
	# LDFLAGS += -L/opt/ibmmath/essl/5.1/lib64 -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lesslsmpbg -lfftw3 
	LDFLAGS += -L/opt/ibmmath/essl/5.1/lib64 -L/bgsys/local/fftw3/3.3.2/fftw_g/lib/ -lesslbg -lfftw3 
endif

#==============================================================================

# Paths to look for files
DIRS = . DebugHelpers IterativeSolver XC madelung LDAU atomiccore xccouplings energy lattice \
       parallel MixingNew datastructures shapefun IterativeSolverBCP debug_morgan bfield 

#also add BUILDDIR to VPATH to look if object files have changed
VPATH = $(DIRS) $(BUILDDIR)

SRCS = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f))
SRCS90 = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.f90))
SRCSFPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F))
SRCS90FPP = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F90))
SRCS95TMPL = $(foreach DIR,$(DIRS),$(wildcard $(DIR)/*.F95))


OBJS=

ifeq ($(tfQMRgpu),yes)
	LDFLAGS += -L$(TFQMRGPU_PATH)/lib64 -ltfQMRgpu_Fortran -ltfQMRgpu
	OBJS += tfQMRgpu_mod.o
endif


LINKER= $(FC90)


# notdir extracts only filename
OBJS += $(notdir ${SRCS:.f=.o})
OBJS += $(notdir ${SRCS90:.f90=.o})
OBJS += $(notdir ${SRCSFPP:.F=.o})
OBJS += $(notdir ${SRCS90FPP:.F90=.o})
OBJS += $(notdir ${SRCS95TMPL:.F95=.o})


.PHONY: all
all:	$(PROGRAM)

$(PROGRAM): $(OBJS)
	$(LINKER) -o $(PROGRAM) $(addprefix $(BUILDDIR)/,$(OBJS)) $(LDFLAGS) 

%.o: %.f
	$(FC) $(FCFLAGS) $(EXTRA_FLAGS) -c $< -o $(BUILDDIR)/$@

%.o: %.f90
	$(FC90) $(FCFLAGS) $(EXTRA_FLAGS) $(FC90FLAGS) -c $< -o $(BUILDDIR)/$@

%.o: %.F
	$(FC) $(FCFLAGS) $(EXTRA_FLAGS) $(PPFLAGS) $(PPFLAGS) -c $< -o $(BUILDDIR)/$@

%.o: %.F90
	$(FC90) $(FCFLAGS) $(EXTRA_FLAGS) $(FC90FLAGS) $(PPFLAGS) -c $< -o $(BUILDDIR)/$@

ifeq ($(tfQMRgpu),yes)
# add a target tfQMRgpu
tfQMRgpu: $(TFQMRGPU_PATH)/example/tfqmrgpu_Fortran_example.F90
	$(FC90) -I $(TFQMRGPU_PATH)/tfQMRgpu/include $(FCFLAGS) $(EXTRA_FLAGS) \
	$(FC90FLAGS) $(PPFLAGS) -c $< -o $(BUILDDIR)/tfQMRgpu_mod.o
endif


### Fortran-templates:
%.o: %.F95
	sed -e 's/_TYPE/D/' $< >  $(BUILDDIR)/$*.F90
	sed -e 's/_TYPE/C/' $< >> $(BUILDDIR)/$*.F90
	sed -e 's/_TYPE/Z/' $< >> $(BUILDDIR)/$*.F90
	sed -e 's/_TYPE/I/' $< >> $(BUILDDIR)/$*.F90
	$(FC90) $(FCFLAGS) $(EXTRA_FLAGS) $(FC90FLAGS) $(PPFLAGS) -c $(BUILDDIR)/$*.F90 -o $(BUILDDIR)/$@

.PHONY: clean
clean:
	rm -f $(BUILDDIR)/*.o 
	rm -f $(BUILDDIR)/*.mod
	rm -f $(BUILDDIR)/*__genmod.f90
	rm -f *.exe

.PHONY: test
test:
	@echo $(OBJS)
	@echo $(VPATH)
	@echo $(SRCS90FPP)

.PHONY: depend 
depend:
	$(SRCS90)
	makedepend

ifeq ($(tfQMRgpu),yes)
	# manually add a dependency
	KKRmat_mod.o: tfQMRgpu_mod.o
endif


#======================== Module dependencies ========================================

KKRnano.o: Logging_mod.o KKRzero_mod.o PotentialConverter_mod.o KKRnanoParallel_mod.o BasisAtom_mod.o AtomicCore_mod.o RadialMeshData_mod.o main2_aux_mod.o ScatteringCalculation_mod.o Main2Arrays_mod.o KKRnano_Comm_mod.o ProcessKKRresults_mod.o InputParams_mod.o EBalanceHandler_mod.o LDAUData_mod.o TimerMpi_mod.o EnergyMesh_mod.o DimParams_mod.o CalculationData_mod.o
TruncationZone_mod.o: Main2Arrays_mod.o
KKRmat_mod.o: ExchangeTable_mod.o two_sided_comm_mod.o arraytest2_mod.o Logging_mod.o TimerMpi_mod.o KKROperator_mod.o one_sided_comm_mod.o ClusterInfo_mod.o IterativeSolver_mod.o jij_calc_mod.o Truncation_mod.o SparseMatrixDescription_mod.o BCPOperator_mod.o SolverStats_mod.o InitialGuess_mod.o fillKKRMatrix_mod.o DirectSolver_mod.o
ProcessKKRresults_mod.o: arraytest2_mod.o KKRnanoParallel_mod.o AtomicForce_mod.o Logging_mod.o Lloyds_formula_mod.o wrappers_mod.o brydbm_new_com_mod.o Main2Arrays_mod.o ShapefunData_mod.o EnergyResults_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o CalculationData_mod.o muffin_tin_zero_mod.o RadialMeshData_mod.o MadelungPotential_mod.o broyden_kkr_mod.o EnergyMesh_mod.o BasisAtom_mod.o DensityResults_mod.o Lloyd0_new_mod.o KKRresults_mod.o InputParams_mod.o TimerMpi_mod.o DimParams_mod.o BroydenData_mod.o NearField_calc_mod.o debug_morgan_mod.o total_energy_mod.o ChebMeshData_mod.o
NearField_kkr_mod.o: Constants_mod.o
NearField_calc_mod.o: RadialMeshData_mod.o NearField_com_mod.o BasisAtom_mod.o CalculationData_mod.o
NearField_mod.o: Harmonics_mod.o MadelungCalculator_mod.o Constants_mod.o NearField_kkr_mod.o
KKRnano_Comm_mod.o: KKRnanoParallel_mod.o jij_calc_mod.o comm_patterns_mod.o Exceptions_mod.o bfield.o
InputParams_mod.o: ConfigReader_mod.o
PolygonFaces_mod.o: Constants_mod.o
ShapeGeometryHelpers_mod.o: Exceptions_mod.o Constants_mod.o
ShapeStandardMesh_mod.o: Exceptions_mod.o Constants_mod.o
wrappers_mod.o: Warnings_mod.o ValenceDensity_mod.o SingleSite_mod.o BasisAtom_mod.o RadialMeshData_mod.o ShapefunData_mod.o EnergyMesh_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o NonCollinearMagnetism_mod.o ##AtomicCore_mod.o
BCPOperator_mod.o: ClusterInfo_mod.o Exceptions_mod.o
CalculationData_mod.o: Exceptions_mod.o ExchangeTable_mod.o KKRnanoParallel_mod.o JelliumPotentials_mod.o LatticeVectors_mod.o TruncationZone_mod.o Main2Arrays_mod.o EnergyResults_mod.o LDAUData_mod.o ShapeGauntCoefficients_mod.o GauntCoefficients_mod.o InitialGuess_mod.o MadelungCalculator_mod.o RadialMeshData_mod.o ConstructShapes_mod.o InterpolateBasisAtom_mod.o BasisAtom_mod.o ShapefunData_mod.o DensityResults_mod.o KKRresults_mod.o InputParams_mod.o RefCluster_mod.o DimParams_mod.o ClusterInfo_mod.o BroydenData_mod.o JijData_mod.o Truncation_mod.o ChebMeshData_mod.o NonCollinearMagnetismData_mod.o bfield.o
NearField_com_mod.o: Logging_mod.o one_sided_comm_mod.o NearField_mod.o MadelungCalculator_mod.o NearField_kkr_mod.o
ShapeIntegrationHelpers_mod.o: Constants_mod.o
ConstructShapes_mod.o: RefCluster_mod.o ShapefunData_mod.o Voronoi_mod.o ShapeFunctions_mod.o LatticeVectors_mod.o
EnergyMesh_mod.o: EnergyMeshHelpers_mod.o
InterpolateBasisAtom_mod.o: RadialMeshData_mod.o BasisAtom_mod.o PotentialData_mod.o AtomicCoreData_mod.o
SingleSiteRef_mod.o: Harmonics_mod.o kkr_helpers_mod.o SingleSiteHelpers_mod.o Exceptions_mod.o
IterativeSolver_mod.o: tfQMR_mod.o SolverStats_mod.o TimerMpi_mod.o
kloopz1_mod.o: KKROperator_mod.o ClusterInfo_mod.o BCPOperator_mod.o KKRmat_mod.o jij_calc_mod.o IterativeSolver_mod.o InitialGuess_mod.o Constants_mod.o TimerMpi_mod.o
ShapeCriticalPoints_mod.o: Constants_mod.o ShapeGeometryHelpers_mod.o PolygonFaces_mod.o
MadelungPotential_mod.o: BasisAtom_mod.o MadelungCalculator_mod.o RadialMeshData_mod.o DensityResults_mod.o EnergyResults_mod.o CalculationData_mod.o
EBalanceHandler_mod.o: KKRnanoParallel_mod.o Exceptions_mod.o
ShapeIntegration_mod.o: ShapeIntegrationHelpers_mod.o Constants_mod.o PolygonFaces_mod.o
BasisAtom_mod.o: Exceptions_mod.o RadialMeshData_mod.o ShapefunData_mod.o PotentialData_mod.o AtomicCoreData_mod.o ChebMeshData_mod.o
tfQMR_mod.o: Logging_mod.o SolverStats_mod.o KKROperator_mod.o BCPOperator_mod.o TimerMpi_mod.o
KKRresults_mod.o: DimParams_mod.o
DimParams_mod.o: ConfigReader_mod.o Exceptions_mod.o
fillKKRMatrix_mod.o: SparseMatrixDescription_mod.o Exceptions_mod.o
Main2Arrays_mod.o: DimParams_mod.o
main2_aux_mod.o: TimerMpi_mod.o
jij_calc_mod.o: JijData_mod.o
broyden_kkr_mod.o: brydbm_new_com_mod.o broyden_second_mod.o BroydenData_mod.o CalculationData_mod.o
ShapeFunctions_mod.o: PolygonFaces_mod.o ShapeCriticalPoints_mod.o Constants_mod.o ShapeIntegration_mod.o ShapeStandardMesh_mod.o
ConfigReader_mod.o: ConfigReaderDictionary_mod.o
KKROperator_mod.o: Truncation_mod.o ClusterInfo_mod.o fillKKRMatrix_mod.o SparseMatrixDescription_mod.o bsrmm_mod.o
Lloyd0_new_mod.o: BasisAtom_mod.o RadialMeshData_mod.o ShapefunData_mod.o EnergyMesh_mod.o LDAUData_mod.o GauntCoefficients_mod.o ValenceDensity_mod.o 
ClusterInfo_mod.o: Statistics_mod.o TruncationZone_mod.o RefCluster_mod.o two_sided_comm_mod.o ExchangeTable_mod.o
ScatteringCalculation_mod.o: KKRnanoParallel_mod.o Constants_mod.o SingleSiteRef_mod.o arraytest2_mod.o Logging_mod.o InputParams_mod.o TruncationZone_mod.o wrappers_mod.o Main2Arrays_mod.o KKRnano_Comm_mod.o jij_calc_mod.o BCPOperator_mod.o LDAUData_mod.o one_sided_comm_mod.o GauntCoefficients_mod.o InitialGuess_mod.o CalculationData_mod.o EBalanceHandler_mod.o EnergyMesh_mod.o KKROperator_mod.o BasisAtom_mod.o KKRresults_mod.o IterativeSolver_mod.o RefCluster_mod.o TimerMpi_mod.o DimParams_mod.o kloopz1_mod.o ClusterInfo_mod.o JijData_mod.o Truncation_mod.o two_sided_comm_mod.o ChebMeshData_mod.o NonCollinearMagnetism_mod.o NonCollinearMagnetism_Helpers_mod.o
total_energy_mod.o: Quadrature_mod.o BasisAtom_mod.o RadialMeshData_mod.o ShapefunData_mod.o ShapeGauntCoefficients_mod.o
vxcgga.o: XCFunctionals_mod.o
KKRzero_mod.o: Exceptions_mod.o BrillouinZoneMesh_mod.o PositionReader_mod.o MadelungCalculator_mod.o Startb1_mod.o EnergyMesh_mod.o EnergyMeshHelpers_mod.o Lattice_mod.o BrillouinZone_mod.o ConfigReader_mod.o ConfigReaderDictionary_mod.o InputParams_mod.o Main2Arrays_mod.o DimParams_mod.o BasisAtom_mod.o ShapefunData_mod.o PotentialData_mod.o AtomicCoreData_mod.o RadialMeshData_mod.o ShapefunData_mod.o read_formatted_mod.o read_formatted_shapefun_mod.o ldauinfo_read.o ldaustart.o
SingleSite_mod.o: SingleSiteHelpers_mod.o
ldauphi.o: SingleSite_mod.o
SingleSiteHelpers_mod.o: Constants_mod.o
ValenceDensity_mod.o: Quadrature_mod.o SingleSite_mod.o SingleSiteHelpers_mod.o
Voronoi_mod.o: Exceptions_mod.o Sorting_mod.o
Exceptions_mod.o: Errors_mod.o Warnings_mod.o StringHelpers_mod.o
Errors_mod.o: Warnings_mod.o
ldauovrlp.o: Quadrature_mod.o
ldauinit.o: Quadrature_mod.o
ldauphi.o: Quadrature_mod.o
AtomicForce_mod.o: Quadrature_mod.o Constants_mod.o
mtzero_new.o: Quadrature_mod.o
epotinb_new.o: Quadrature_mod.o
vxclm.o: Quadrature_mod.o
rhototb_new.o: Quadrature_mod.o
vxcgga.o: Quadrature_mod.o
GauntCoefficients_mod.o: Harmonics_mod.o
ShapeGauntCoefficients_mod.o: Harmonics_mod.o
MadelungCalculator_mod.o: Harmonics_mod.o Lattice_mod.o Sorting_mod.o
sphere_gga.o: Harmonics_mod.o
sphere_nogga.o: Harmonics_mod.o
debug_morgan_mod.o: Harmonics_mod.o
BrillouinZone_mod.o: Symmetry_mod.o VectorMath_mod.o
Lattice_mod.o: VectorMath_mod.o Constants_mod.o
Symmetry_mod.o: VectorMath_mod.o Exceptions_mod.o
GaussWeights_mod.o: Exceptions_mod.o
EnergyMeshHelpers_mod.o: GaussWeights_mod.o Constants_mod.o
jij_calc_mod.o: Sorting_mod.o
PotentialConverter_mod.o: DimParams_mod.o BasisAtom_mod.o RadialMeshData_mod.o PotentialData_mod.o
Startb1_mod.o: RadialMeshData_mod.o BasisAtom_mod.o read_formatted_mod.o read_formatted_shapefun_mod.o InputParams_mod.o ChebMeshData_mod.o
RadialMeshData_mod.o: Exceptions_mod.o
read_formatted_mod.o: Exceptions_mod.o
BrillouinZone_mod.o: Exceptions_mod.o BrillouinZoneMesh_mod.o 
KKRnanoParallel_mod.o: Exceptions_mod.o
MadelungCalculator_mod.o: Exceptions_mod.o
PositionReader_mod.o: Exceptions_mod.o
RefCluster_mod.o: Sorting_mod.o
LatticeVectors_mod.o: Sorting_mod.o
jij_calc_mod.o: Symmetry_mod.o
Harmonics_mod.o: Constants_mod.o
KKRmat_mod.o: Exceptions_mod.o InitialGuess_mod.o jij_calc_mod.o SolverStats_mod.o IterativeSolver_mod.o BCPOperator_mod.o KKROperator_mod.o TimerMpi_mod.o fillKKRMatrix_mod.o ClusterInfo_mod.o
BrillouinZoneMesh_mod.o: Exceptions_mod.o
AtomicCore_mod.o: Exceptions_mod.o Quadrature_mod.o Constants_mod.o
JelliumPotentials_mod.o: Exceptions_mod.o
DirectSolver_mod.o: Exceptions_mod.o KKROperator_mod.o fillKKRMatrix_mod.o
ExchangeTable_mod.o: ChunkIndex_mod.o
one_sided_comm_mod.o: ChunkIndex_mod.o 
two_sided_comm_mod.o: ExchangeTable_mod.o
LocalAtomData_mod.o: Exceptions_mod.o RadialMeshData_mod.o ShapefunData_mod.o BasisAtom_mod.o LDAUData_mod.o JijData_mod.o RefCluster_mod.o MadelungCalculator_mod.o DensityResults_mod.o EnergyResults_mod.o KKRresults_mod.o
bsrmm_mod.o: CacheOverlap_mod.o
CacheOverlap_mod.o:
ChebMeshData_mod.o: InputParams_mod.o RadialMeshData_mod.o
Truncation_mod.o: Logging_mod.o Exceptions_mod.o TruncationZone_mod.o
NonCollinearMagnetism_mod.o: RadialMeshData_mod.o ChebMeshData_mod.o read_formatted_shapefun_mod.o NonCollinearMagnetism_Helpers_mod.o bfield.o torque.o
NonCollinearMagnetismData_mod.o: Exceptions_mod.o
vintras_new.o: SingleSiteHelpers_mod.o
bfield.o: NonCollinearMagnetism_Helpers_mod.o
torque.o: NonCollinearMagnetism_Helpers_mod.o bfield.o

# DO NOT DELETE
