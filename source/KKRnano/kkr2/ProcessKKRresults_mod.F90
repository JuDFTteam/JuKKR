#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

module ProcessKKRresults_mod

  implicit none

  public :: processKKRresults
  private :: calculateDensities
  private :: calculatePotentials

CONTAINS

!------------------------------------------------------------------------------
!> Returns 1 when target rms error has been reached,
!> master rank adds 2 if STOP-file present.
integer function processKKRresults(iter, calc_data, my_mpi, emesh, dims, params, arrays, &
                                   program_timer)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use EnergyMesh_mod

  use CalculationData_mod

  use RadialMeshData_mod
  use BasisAtom_mod

  use TimerMpi_mod

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use DensityResults_mod
  use EnergyResults_mod

  implicit none

  include 'mpif.h'

  integer, intent(in)                                  :: iter
  type (CalculationData), intent(inout)                :: calc_data
  type (KKRnanoParallel), intent(in)                   :: my_mpi
  type (EnergyMesh), intent(inout)                     :: emesh
  type (Main2Arrays), intent(in)                       :: arrays
  type (DimParams)  , intent(in)                       :: dims
  type (InputParams), intent(in)                       :: params
  type (TimerMpi), intent(in)                          :: program_timer

  ! locals
  type (BasisAtom) , pointer                           :: atomdata
  type (DensityResults), pointer                       :: densities
  type (EnergyResults), pointer                        :: energies

  type (RadialMeshData), pointer :: mesh
  integer :: I1
  integer :: ierr
  integer :: ilocal
  integer :: num_local_atoms

  processKKRresults = 0
  num_local_atoms = getNumLocalAtoms(calc_data)

  densities    => getDensities(calc_data, 1)
  energies     => getEnergies(calc_data, 1)

  ! kkr
  !  |
  !  v
  call calculateDensities(iter, calc_data, my_mpi, dims, params, &
                          program_timer, arrays, emesh)
  ! |
  ! v
  ! modified: densities, emesh, energies (ESPV only)
  ! |
  ! v
  call calculatePotentials(iter, calc_data, my_mpi, dims, params, &
                           program_timer, arrays)
  ! |
  ! v
  ! modified: atomdata, energies
  !
  ! |
  ! v

  ! mix_potential returns 1 when target_rms reached, otherwise 0
  processKKRresults = mix_potential(calc_data, iter, params, dims, my_mpi)

  ! |
  ! v
  ! modified: atomdata


  ! use any atomdata to open file - use reclen stored in calc_data
  atomdata => getAtomData(calc_data, 1)
  call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew", &
                                    getMaxReclenPotential(calc_data))

  do ilocal = 1, num_local_atoms ! no OpenMP
    atomdata => getAtomData(calc_data, ilocal)
    energies => getEnergies(calc_data, ilocal)
    mesh => atomdata%mesh_ptr
    I1 = getAtomIndexOfLocal(calc_data, ilocal)

    TESTARRAYLOG(3, atomdata%potential%VINS)
    TESTARRAYLOG(3, atomdata%potential%VISP)
    TESTARRAYLOG(3, atomdata%potential%VONS)

  !----------------------------------------------------------------------
  ! -->    reset to start new iteration
  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    call resetPotentials(atomdata)

  ! ----------------------------------------------------- output_potential
    call writeBasisAtomPotentialDA(atomdata, 37, I1)
  ! ----------------------------------------------------- output_potential

  end do ! ilocal

  call closeBasisAtomPotentialDAFile(37)

  call OUTTIME(isMasterRank(my_mpi),'potential written .. ', &
               getElapsedTime(program_timer), iter)

! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS
  call OUTTIME(isMasterRank(my_mpi),'barrier begin.. ', &
               getElapsedTime(program_timer), iter)

  call MPI_BARRIER(getMySEcommunicator(my_mpi),IERR)

  call OUTTIME(isMasterRank(my_mpi),'barrier end ... ', &
               getElapsedTime(program_timer), iter)

! -----------------------------------------------------------------
! BEGIN: only MASTERRANK is working here
! -----------------------------------------------------------------
  if(isMasterRank(my_mpi)) then

    ! DOS was written to file 'results1' and read out here just
    ! to be written in routine wrldos
    ! also other stuff is read from results1 (and results2)

    ! TODO: note: title written to DOS files is not correct
    ! - taken from 1st local atom only
    ! TODO: Fermi energy written to DOS files is not correct
    call RESULTS(dims%LRECRES2,densities%IEMXD,ITER,dims%LMAXD, &
    arrays%NAEZ,emesh%NPOL, &
    dims%NSPIND,params%KPRE,params%KTE,atomdata%potential%LPOT, &
    emesh%E1,emesh%E2,emesh%TK,emesh%EFERMI, &
    params%ALAT,atomdata%core%ITITLE(:,1:dims%NSPIND), &
    densities%total_charge_neutrality, &
    arrays%ZAT,emesh%EZ,emesh%WEZ,params%LDAU, &
    arrays%iemxd)

    call OUTTIME(isMasterRank(my_mpi),'results......',getElapsedTime(program_timer), iter)

    ! manual exit possible by creation of file 'STOP' in home directory
    processKKRresults = processKKRresults + 2*stopfile_flag()

  endif
! -----------------------------------------------------------------
! END: only MASTERRANK is working here
! -----------------------------------------------------------------

end function


!------------------------------------------------------------------------------
!> Performs potential mixing.
!>
!> Returns 1 when params%target_rms is reached, otherwise 0
!>
!> Selection of algorithm according to value of params%imix
!> imix = 0     straight mixing
!> imix = 1, 2  not used - defaults to straight mixing
!> imix = 3     Broyden's 1st method
!> imix = 4     Broyden's 2nd method
!> imix = 5     generalised Anderson
!> imix = 6     Broyden's 2nd method with support for num_local_atoms > 1
integer function mix_potential(calc_data, iter, params, dims, my_mpi)

  use KKRnanoParallel_mod
  use CalculationData_mod
  use RadialMeshData_mod
  use BasisAtom_mod
  use BroydenData_mod
  use BRYDBM_new_com_mod
  use wrappers_mod
  use DimParams_mod
  use InputParams_mod
  use broyden_kkr_mod

  implicit none

  integer, intent(in)                                  :: iter
  type (CalculationData), intent(inout)                :: calc_data
  type (KKRnanoParallel), intent(in)                   :: my_mpi
  type (DimParams)  , intent(in)                       :: dims
  type (InputParams), intent(in)                       :: params

  type (BasisAtom) , pointer                           :: atomdata
  type (BroydenData), pointer                          :: broyden

  type (RadialMeshData), pointer :: mesh
  double precision :: RMSAVQ ! rms error charge dens. (contribution of all local sites)
  double precision :: RMSAVM ! rms error mag. density (contribution of all local sites)
  double precision :: RMSAVQ_single
  double precision :: RMSAVM_single
  integer :: ilocal
  integer :: num_local_atoms

  mix_potential = 0

  num_local_atoms = getNumLocalAtoms(calc_data)

  ! -->   calculation of RMS and final construction of the potentials (straight mixing)
  RMSAVQ = 0.0d0
  RMSAVM = 0.0d0

  ! straight/simple mixing
  !$omp parallel do reduction(+: RMSAVQ, RMSAVM) private(ilocal, atomdata, RMSAVQ_single, RMSAVM_single)
  do ilocal = 1, num_local_atoms
    atomdata => getAtomData(calc_data, ilocal)

    call MIXSTR_wrapper(atomdata, RMSAVQ_single, RMSAVM_single, params%MIXING, params%FCM)

    RMSAVQ = RMSAVQ + RMSAVQ_single
    RMSAVM = RMSAVM + RMSAVM_single
  end do
  !$omp end parallel do

  ! summation and output of RMS error
  call RMSOUT_com(RMSAVQ,RMSAVM,ITER,dims%NSPIND,dims%NAEZ, &
                 getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi))

  ! check if target rms error has been reached and set abort flag
  if (rmsavq <= params%target_rms) then
    mix_potential = 1
    if (isMasterRank(my_mpi)) write(*,*) "TARGET RMS ERROR REACHED..."
  end if

  ! straight mixing is called in any case before
  ! Broyden mixing - it is undone in Broyden routines

! -->   potential mixing procedures: Broyden or Andersen updating schemes
  if (params%IMIX>=3 .and. params%IMIX<6) then

    if (num_local_atoms > 1) then
      if (isMasterRank(my_mpi)) write(*,*) "Broyden mixing and num_local_atoms > 1 not supported."
      STOP
    end if

    ! Take data from 1st local atom, since only one local atom is supported
    atomdata     => getAtomData(calc_data, 1)
    broyden      => getBroyden(calc_data, 1)
    mesh => atomdata%mesh_ptr

    call BRYDBM_new_com(atomdata%potential%VISP,atomdata%potential%VONS, &
    atomdata%potential%VINS, &
    atomdata%potential%LMPOT,mesh%R,mesh%DRDI,broyden%MIXING, &
    mesh%IRC,mesh%IRMIN,atomdata%potential%NSPIN, &
    broyden%IMIX,ITER, &
    broyden%UI2,broyden%VI2,broyden%WIT,broyden%SM1S,broyden%FM1S, &
    getMyAtomRank(my_mpi), &
    getMySEcommunicator(my_mpi), &
    broyden%itdbryd, mesh%irmd, atomdata%potential%irnsd, &
    atomdata%potential%nspin)

  ! this method supports num_local_atoms > 1
  else if (params%imix == 6) then
    ! use Broyden mixing that supports num_local_atoms > 1
    call mix_broyden2_com(calc_data, iter, getMySEcommunicator(my_mpi))
  end if

end function

!------------------------------------------------------------------------------
!> Write forces to file 'forces'.
!>
!> Gather all forces at rank 'master', only this rank writes the file.
!> Since the amount of data for forces is low this is a reasonable approach.
subroutine output_forces(calc_data, master, rank, comm)
  use CalculationData_mod
  use DensityResults_mod
  implicit none
  include 'mpif.h'
  type(CalculationData), intent(in) :: calc_data
  integer, intent(in) :: master
  integer, intent(in) :: rank
  integer, intent(in) :: comm

  integer :: ilocal, I1, num_local_atoms
  type(DensityResults), pointer :: densities
  double precision, allocatable :: force_buffer(:,:)
  double precision, allocatable :: local_buffer(:,:)
  integer :: nranks
  integer :: ierr
  integer, parameter :: NUM = 3

  num_local_atoms = getNumLocalAtoms(calc_data) ! must be the same for all ranks

  if (rank == master) then
    call MPI_Comm_size(comm, nranks, ierr)
    allocate(force_buffer(NUM, num_local_atoms * nranks))
  else
    allocate(force_buffer(1, 1))
  end if

  allocate(local_buffer(NUM, num_local_atoms))

  do ilocal = 1, num_local_atoms
    densities    => getDensities(calc_data, ilocal)

    local_buffer(:, ilocal) = densities%force_flm(-1:1)
  end do

  call MPI_Gather(local_buffer, NUM*num_local_atoms, MPI_DOUBLE_PRECISION, &
                  force_buffer, NUM*num_local_atoms, MPI_DOUBLE_PRECISION, &
                  master, comm, ierr)

  if (rank == master) then
    call openForceFile()
    do I1 = 1, num_local_atoms * nranks
      call writeForces(force_buffer(:,I1), I1)
    end do
    call closeForceFile()
  end if

  deallocate(force_buffer)
  deallocate(local_buffer)

end subroutine

!==============================================================================

!------------------------------------------------------------------------------
!> Calculate densities.
!>
!> output: emesh (Fermi-energy updated, renormalized weights), densities, ldau_data?
!> files written: 'results1'
subroutine calculateDensities(iter, calc_data, my_mpi, dims, params, &
                              program_timer, arrays, emesh)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use lloyds_formula_mod, only: renormalizeDOS

  use EnergyMesh_mod

  use lloyd0_new_mod

  use CalculationData_mod
  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod

  use RadialMeshData_mod
  use CellData_mod
  use BasisAtom_mod

  use LDAUData_mod

  use TimerMpi_mod

  use wrappers_mod

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use KKRresults_mod
  use DensityResults_mod
  use EnergyResults_mod

  implicit none

  integer, intent(in)                        :: iter
  type (CalculationData), intent(inout)      :: calc_data
  type (KKRnanoParallel), intent(in)         :: my_mpi
  type (EnergyMesh), intent(inout)           :: emesh
  type (Main2Arrays), intent(in)             :: arrays
  type (DimParams), intent(in)               :: dims
  type (InputParams), intent(in)             :: params
  type (TimerMpi), intent(in)                :: program_timer

  ! locals
  type (ShapeGauntCoefficients), pointer     :: shgaunts  ! const ref.
  type (GauntCoefficients), pointer          :: gaunts    ! const ref
  type (BasisAtom), pointer                  :: atomdata  ! not const
  type (LDAUData), pointer                   :: ldau_data ! not const
  type (KKRresults), pointer                 :: kkr       ! const ref
  type (DensityResults), pointer             :: densities ! not const
  type (EnergyResults), pointer              :: energies  ! not const
  type (RadialMeshData), pointer             :: mesh
  type (CellData), pointer                   :: cell

  double complex, parameter      :: CZERO = (0.0d0, 0.0d0)
  logical :: LdoRhoEF
  integer :: I1
  double precision :: denEf !< charge density at Fermi level
  double precision :: chrgNt !< charge neutrality
  double precision :: denEf_local
  double precision :: chrgNt_local
  double precision :: new_fermi
  integer :: ilocal
  integer :: num_local_atoms

  num_local_atoms = getNumLocalAtoms(calc_data)

  shgaunts     => getShapeGaunts(calc_data)
  gaunts       => getGaunts(calc_data)
  atomdata     => getAtomData(calc_data, 1)
  ldau_data    => getLDAUData(calc_data, 1)
  kkr          => getKKR(calc_data, 1)
  densities    => getDensities(calc_data, 1)
  energies     => null()
  mesh         => null()
  cell         => null()

  I1 = 0

  if (dims%LLY /= 0 .and. num_local_atoms > 1) then
    if (isMasterRank(my_mpi)) write(*,*) "Lloyd's formula and num_local_atoms > 1 not supported."
    STOP
  endif

  ! out: emesh, RNORM
  call lloyd0_wrapper_com(atomdata, my_mpi, kkr%LLY_GRDT, &
                          emesh, densities%RNORM, &
                          dims%LLY, params%ICST, params%NSRA, &
                          kkr%GMATN, gaunts, ldau_data)

  if (dims%LLY == 1) then
    TESTARRAYLOG(3, emesh%WEZRN)
    TESTARRAYLOG(3, densities%RNORM)
    call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......', &
                 getElapsedTime(program_timer),ITER)
  endif

  ! now WEZRN stores the weights for E-integration

  DENEF = 0.0D0
  CHRGNT = 0.0D0

  LDORHOEF = emesh%NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

!------------------------------------------------------------------------------
  !$omp parallel do reduction(+: chrgnt, denef) private(ilocal, atomdata, &
  !$omp densities, energies, kkr, ldau_data, denef_local, chrgnt_local)
  do ilocal = 1, num_local_atoms
    atomdata  => getAtomData(calc_data, ilocal)
    densities => getDensities(calc_data, ilocal)
    energies  => getEnergies(calc_data, ilocal)
    kkr       => getKKR(calc_data, ilocal)
    ldau_data => getLDAUData(calc_data, ilocal)
!------------------------------------------------------------------------------

    if (params%LDAU) then
      ldau_data%DMATLDAU = CZERO
    endif

    ! has to be done after Lloyd
    ! output: RHO2NS, R2NEF, DEN, ESPV
    densities%DEN = CZERO

    call RHOVAL_wrapper(atomdata, LdoRhoEF, params%ICST, params%NSRA, &
                        densities%RHO2NS, densities%R2NEF, &
                        densities%DEN, energies%ESPV, kkr%GMATN, &
                        gaunts, emesh, ldau_data)

  ! ----------------------------------------------------------------------
  ! -->   determine total charge expanded in spherical harmonics
  ! -------------------------------------------------------------- density
    ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
    !         core charge density added to densities%RHO2NS
    call RHOTOTB_wrapper(densities%CATOM, densities%RHO2NS, atomdata)

    CHRGNT_local = densities%CATOM(1) - atomdata%Z_nuclear

    if (dims%LLY == 1) then
      call renormalizeDOS(densities%DEN,densities%RNORM, &
                          densities%LMAXD+1,densities%IEMXD, &
                          densities%NSPIND,densities%IEMXD)
    end if

    ! calculate DOS at Fermi level
    DENEF_local = calcDOSatFermi(densities%DEN, densities%IEMXD, &
                                   densities%IEMXD, densities%LMAXD+1, &
                                   densities%NSPIND)

    ! ---> l/m_s/atom-resolved charges, output -> CHARGE
    ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
    ! CHARGE -> written to result file
    call calcChargesLres(densities%CHARGE, densities%DEN, emesh%ielast, &
                         densities%LMAXD+1, densities%NSPIND, emesh%WEZ, &
                         densities%IEMXD)

    DENEF = DENEF + DENEF_local
    CHRGNT = CHRGNT + CHRGNT_local

!------------------------------------------------------------------------------
  end do ! ilocal
  !$omp end parallel do
!------------------------------------------------------------------------------


  call sumNeutralityDOSFermi_com(CHRGNT, DENEF, getMySEcommunicator(my_mpi))

  ! write to 'results1' - only to be read in in results.f
  ! necessary for density of states calculation, otherwise
  ! only for informative reasons
  if (params%KTE >= 0) then
    call openResults1File(dims%IEMXD, dims%LMAXD, emesh%NPOL)

    do ilocal = 1, num_local_atoms
      atomdata  => getAtomData(calc_data, ilocal)
      densities => getDensities(calc_data, ilocal)
      I1 = getAtomIndexOfLocal(calc_data, ilocal)

      call writeResults1File(densities%CATOM, densities%CHARGE, densities%DEN, &
                             atomdata%core%ECORE, I1, emesh%NPOL, &
                             atomdata%core%QC_corecharge)
    end do

    call closeResults1File()
  endif

  call OUTTIME(isMasterRank(my_mpi),'density calculated ..', &
               getElapsedTime(program_timer),ITER)

!------------------------------------------------------------------------------
  ! OMP had problems here
  !!!$omp parallel do private(ilocal, atomdata, densities, energies, mesh, cell) lastprivate(new_fermi)
  do ilocal = 1, num_local_atoms
    atomdata  => getAtomData(calc_data, ilocal)
    densities => getDensities(calc_data, ilocal)
    energies  => getEnergies(calc_data, ilocal)
    mesh         => atomdata%mesh_ptr
    cell         => atomdata%cell_ptr
!------------------------------------------------------------------------------
    densities%total_charge_neutrality = CHRGNT

    new_fermi = emesh%E2

    ! allow only a maximal Fermi Energy shift of 0.03 Ry
    call doFermiEnergyCorrection(atomdata, isMasterRank(my_mpi), arrays%naez, &
                                 0.03d0, CHRGNT, DENEF, densities%R2NEF, &
                                 energies%ESPV, densities%RHO2NS, new_fermi)

    !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)

    ! calculate multipole moments
    densities%CMOM   = 0.0D0
    densities%CMINST = 0.0D0
    call RHOMOM_NEW_wrapper(densities%CMOM,densities%CMINST, &
                            densities%RHO2NS(:,:,1), cell, mesh, shgaunts)

!------------------------------------------------------------------------------
  end do ! ilocal
  !!!$omp end parallel do
!------------------------------------------------------------------------------

  emesh%E2 = new_fermi  ! Assumes that for every atom the same Fermi correction
                        ! was calculated !!!

  call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......', &
               getElapsedTime(program_timer),ITER)


end subroutine


!------------------------------------------------------------------------------
!> Calculate potentials.
!>
!> Output: atomdata, ldau_data
!> Files written: 'results2'
subroutine calculatePotentials(iter, calc_data, my_mpi, dims, params, &
                               program_timer, &
                               arrays)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use CalculationData_mod
  use MadelungPotential_mod
  use ShapeGauntCoefficients_mod
  use muffin_tin_zero_mod

  use BasisAtom_mod

  use LDAUData_mod

  use TimerMpi_mod

  use wrappers_mod

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use DensityResults_mod
  use EnergyResults_mod
  use RadialMeshData_mod

  use NearField_calc_mod

  implicit none

  integer, intent(in)                        :: iter
  type (CalculationData), intent(inout)      :: calc_data
  type (KKRnanoParallel), intent(in)         :: my_mpi
  type (Main2Arrays), intent(in)             :: arrays
  type (DimParams), intent(in)               :: dims
  type (InputParams), intent(in)             :: params
  type (TimerMpi), intent(in)                :: program_timer

  ! locals
  type (ShapeGauntCoefficients), pointer     :: shgaunts     ! const ref
  type (BasisAtom), pointer                  :: atomdata     ! not const
  type (LDAUData), pointer                   :: ldau_data    ! not const
  type (EnergyResults), pointer              :: energies     ! not const
  type (DensityResults), pointer             :: densities    ! not const
  type (RadialMeshData), pointer             :: mesh

  integer :: I1
  integer :: lcoremax
  double precision :: VAV0, VOL0
  double precision :: VAV0_local, VOL0_local
  double precision :: VBC_new(2)
  integer :: ilocal
  integer :: num_local_atoms
  logical :: calc_force
  double precision :: force_flmc(-1:1)
    character(len=:), allocatable :: str


  num_local_atoms = getNumLocalAtoms(calc_data)

  shgaunts     => getShapeGaunts(calc_data)
  atomdata     => getAtomData(calc_data, 1)
  ldau_data    => getLDAUData(calc_data, 1)
  densities    => null()
  energies     => getEnergies(calc_data, 1)


  I1 = 0

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

  calc_force = (params%KFORCE == 1) ! calculate force at each iteration

!------------------------------------------------------------------------------
  !!!$omp parallel do private(ilocal, atomdata, densities)
  ! OpenMP problems?
  do ilocal = 1, num_local_atoms
    atomdata     => getAtomData(calc_data, ilocal)
    densities    => getDensities(calc_data, ilocal)
!------------------------------------------------------------------------------

    !output: VONS
    call VINTRAS_wrapper(densities%RHO2NS(:,:,1), shgaunts, atomdata)

!------------------------------------------------------------------------------
  end do ! ilocal
  !!!$omp end parallel do
!------------------------------------------------------------------------------

  call OUTTIME(isMasterRank(my_mpi),'VINTRAS ......',&
               getElapsedTime(program_timer),ITER)

  ! perform near field correction for ALL local atoms
  if (params%near_field > 0) then
    call add_near_field_corr(calc_data, arrays, params%alat, my_mpi)

    call OUTTIME(isMasterRank(my_mpi),'near field....',&
                 getElapsedTime(program_timer),ITER)
  end if

  ! output: VONS (changed), VMAD
  ! operation on all atoms! O(N**2)
  call addMadelungPotentialnew_com(calc_data, arrays%ZAT, getMyAtomRank(my_mpi), &
                                dims%atoms_per_proc, &
                                getMySEcommunicator(my_mpi))

  call OUTTIME(isMasterRank(my_mpi),'VMADELBLK ......', &
               getElapsedTime(program_timer),ITER)

  VAV0 = 0.0d0
  VOL0 = 0.0d0

  if (params%KTE >= 0) call openResults2File(dims%LRECRES2)

!------------------------------------------------------------------------------
  !!!$omp parallel do reduction(+: VAV0, VOL0) &
  !!!$omp private(ilocal, atomdata, densities, energies, ldau_data, I1, &
  !!!$omp         lcoremax, VAV0_local, VOL0_local, mesh, force_flmc)
  do ilocal = 1, num_local_atoms
    atomdata     => getAtomData(calc_data, ilocal)
    densities    => getDensities(calc_data, ilocal)
    energies     => getEnergies(calc_data, ilocal)
    ldau_data    => getLDAUData(calc_data, ilocal)
    mesh => atomdata%mesh_ptr
    I1 = getAtomIndexOfLocal(calc_data, ilocal)
!------------------------------------------------------------------------------

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

    if (calc_force) then

      call FORCEH(densities%force_FLM,atomdata%potential%LPOT, &
                  densities%RHO2NS,atomdata%potential%VONS, &
                  mesh%R,mesh%DRDI,mesh%IMT,atomdata%Z_nuclear,mesh%irmd)

      force_FLMC = 0.0d0 ! temporary needed later in forcxc

      call FORCE(densities%force_FLM,force_FLMC,atomdata%potential%LPOT, &
                 atomdata%potential%NSPIN, atomdata%core%RHOCAT, &
                 atomdata%potential%VONS, mesh%R, mesh%DRDI, &
                 mesh%IMT, mesh%irmd)

    end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

    lcoremax = 0
    if (params%KTE==1) then
      ! These energies have to be calculated BEFORE the XC-potential is added!
      ! calculate total energy and individual contributions if requested
      ! core electron contribution
      call ESPCB_wrapper(energies%ESPC, LCOREMAX, atomdata)

      ! output: EPOTIN
      call EPOTINB_wrapper(energies%EPOTIN,densities%RHO2NS,atomdata)

      ! output: ECOU - l resolved Coulomb energy
      call ECOUB_wrapper(densities%CMOM, energies%ECOU, densities%RHO2NS, &
                         shgaunts, atomdata)

    end if

  ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

  ! =====================================================================
    ! TODO: OpenMP critical !!! VXCDRV is most likely not threadsafe!
    ! output: VONS (changed), EXC (exchange energy) (l-resolved)

    !!!$omp critical
    call VXCDRV_wrapper(energies%EXC, params%KXC, densities%RHO2NS, &
                        shgaunts, atomdata)
    !!!$omp end critical

  ! =====================================================================

  ! Force calculation continues here

    if (calc_force) then
      call FORCXC(densities%force_FLM,force_FLMC,atomdata%potential%LPOT, &
                  atomdata%potential%NSPIN, atomdata%core%RHOCAT, &
                  atomdata%potential%VONS, mesh%R, &
                  mesh%DRDI, mesh%IMT, mesh%irmd)
    end if

    ! writes some results (mostly energies) to direct access file 'results2',
    ! those are read in again in routine 'RESULTS'
    if (params%KTE >= 0) then
      call writeResults2File(densities%CATOM, energies%ECOU, ldau_data%EDCLDAU, &
                             energies%EPOTIN, energies%ESPC, energies%ESPV, ldau_data%EULDAU, &
                             energies%EXC, I1, LCOREMAX, energies%VMAD)
    end if

    ! calculate new muffin-tin zero. output: VAV0, VOL0
    VAV0_local = 0.0d0
    VOL0_local = 0.0d0
    call MTZERO_wrapper(VAV0_local, VOL0_local, atomdata)
    VAV0 = VAV0 + VAV0_local
    VOL0 = VOL0 + VOL0_local

!------------------------------------------------------------------------------
  end do ! ilocal
  !!!$omp end parallel do
!------------------------------------------------------------------------------

  if (params%KTE >= 0) call closeResults2File()

  call OUTTIME(isMasterRank(my_mpi),'before CONVOL.....', &
               getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

  call allreduceMuffinTinShift_com(getMySEcommunicator(my_mpi), VAV0, VBC_new, VOL0)

  if(isMasterRank(my_mpi)) then
    call printMuffinTinShift(VAV0, VBC_new, VOL0)
  end if

!------------------------------------------------------------------------------
  !$omp parallel do private(ilocal, atomdata, energies)
  do ilocal = 1, num_local_atoms
    atomdata     => getAtomData(calc_data, ilocal)
    energies    => getEnergies(calc_data, ilocal)
!------------------------------------------------------------------------------

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

! -->   shift potential by VBC and multiply with shape functions - output: VONS
!       add also an optional muffin-tin-zero shift 'mt_zero_shift'
    energies%VBC = VBC_new + params%mt_zero_shift
    call CONVOL_wrapper(energies%VBC, shgaunts, atomdata)

!------------------------------------------------------------------------------
  end do
  !$omp end parallel do
!------------------------------------------------------------------------------

! LDAU
  ldau_data%EULDAU = 0.0D0
  ldau_data%EDCLDAU = 0.0D0

  if (ldau_data%LDAU.and.ldau_data%NLDAU>=1) then
    call LDAUWMAT(I1,ldau_data%NSPIND,ITER,params%MIXING,ldau_data%DMATLDAU, &
                  ldau_data%NLDAU,ldau_data%LLDAU, &
                  ldau_data%ULDAU,ldau_data%JLDAU, &
                  ldau_data%UMLDAU,ldau_data%WMLDAU, &
                  ldau_data%EULDAU,ldau_data%EDCLDAU, &
                  ldau_data%lmaxd)
  endif
! LDAU

  call OUTTIME(isMasterRank(my_mpi),'calculated pot ......',getElapsedTime(program_timer),ITER)

end subroutine

!------------------------------------------------------------------------------
subroutine openForceFile()
  implicit none
  integer :: reclen
  double precision :: dummy(-1:1)

  inquire (iolength = reclen) dummy ! get reclen for 3 doubles
  open(91, access='direct', file='forces', recl=reclen, form='unformatted', &
           status='replace')
end subroutine

!------------------------------------------------------------------------------
subroutine writeForces(forces_flm, recnr)
  implicit none
  double precision, intent(in) :: forces_flm(-1:1)
  integer, intent(in) :: recnr

  write(91, rec=recnr) forces_flm

end subroutine

!------------------------------------------------------------------------------
subroutine closeForceFile()
  implicit none
  close(91)
end subroutine

  !----------------------------------------------------------------------------
  !> Open Results2 File
  subroutine openResults2File(LRECRES2)
    implicit none
    integer :: LRECRES2
    open (72,access='direct',recl=LRECRES2,file='results2', &
    form='unformatted')
  end subroutine

  !----------------------------------------------------------------------------
  !> Write calculated stuff into historical 'results2' file
  subroutine writeResults2File(CATOM, ECOU, EDCLDAU, EPOTIN, ESPC, ESPV, EULDAU, EXC, I1, LCOREMAX, VMAD)
    implicit none
    double precision :: CATOM(:)
    double precision :: ECOU(:)
    double precision :: EDCLDAU
    double precision :: EPOTIN
    double precision :: ESPC(:,:)
    double precision :: ESPV(:,:)
    double precision :: EULDAU
    double precision :: EXC(:)
    integer :: I1
    integer :: LCOREMAX
    double precision :: VMAD

    write(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
                     EULDAU,EDCLDAU
  end subroutine

  !----------------------------------------------------------------------------
  !> Close the file 'results2'
  subroutine closeResults2File()
    implicit none
    close(72)
  end subroutine

  !----------------------------------------------------------------------------
  !> Open the file 'results1'
  subroutine openResults1File(IEMXD, LMAXD, NPOL)
    implicit none
    integer :: IEMXD
    integer :: LMAXD
    integer :: NPOL

    integer :: LRECRES1

    LRECRES1 = 8*43 + 16*(LMAXD+2)
    if (NPOL==0) then
      LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
    end if

    open (71,access='direct',recl=LRECRES1,file='results1', &
    form='unformatted')
  end subroutine

  !----------------------------------------------------------------------------
  !> Write some stuff to the 'results1' file
  subroutine writeResults1File(CATOM, CHARGE, DEN, ECORE, I1, NPOL, QC)
    implicit none
    double precision :: CATOM(:)
    double precision :: CHARGE(:,:)
    double complex :: DEN(:,:,:)
    double precision :: ECORE(20,2)
    integer :: I1
    integer :: NPOL
    double precision :: QC

    if (NPOL==0) then
      write(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN  ! write density of states (DEN) only when certain options set
    else
      write(71,rec=I1) QC,CATOM,CHARGE,ECORE
    end if
  end subroutine

  !---------------------------------------------------------------------------
  !> Closes the file 'results1'
  subroutine closeResults1File()
    implicit none
    close(71)
  end subroutine

  !----------------------------------------------------------------------------
  !> Communicate and sum up contributions for charge neutrality and
  !> density of states at Fermi level.
  subroutine sumNeutralityDOSFermi_com(CHRGNT, DENEF, communicator)
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: CHRGNT
    double precision, intent(inout) :: DENEF
    integer, intent(in) :: communicator
    !----------------------------------

    double precision :: WORK1(2)
    double precision :: WORK2(2)
    integer :: ierr

    WORK1(1) = CHRGNT
    WORK1(2) = DENEF

    call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
    communicator,IERR)

    CHRGNT = WORK2(1)
    DENEF  = WORK2(2)

  end subroutine

  !----------------------------------------------------------------------------
  !> Extracts the DOS at the Fermi level and stores result in DENEF
  !> @param IELAST energy index of Fermi level
  !> @param LMAXD1 lmax+1
  double precision function calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIN)
    implicit none

    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    integer, intent(in) :: IELAST
    integer, intent(in) :: IEMXD
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: NSPIN

    integer :: ISPIN
    integer :: L
    double precision :: PI
    double precision :: DENEF

    PI = 4.0D0*ATAN(1.0D0)

    ! get density of states at Fermi-level
    DENEF = 0.0d0
    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        DENEF = DENEF - 2.0D0 * &
        DIMAG(DEN(L,IELAST,ISPIN))/PI/DBLE(NSPIN)
      end do
    end do

    calcDOSatFermi = DENEF

  end function calcDOSatFermi

  !----------------------------------------------------------------------------
  !> Calculates the l-resolved charges.
  !> This is done by energy integration in the complex plane over the imaginary
  !> part of the diagonal of the structural Green's function (DEN)

  subroutine calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIN, WEZ, IEMXD)
    implicit none

    double precision, intent(out) :: CHARGE(0:LMAXD1,2)
    double complex, intent(in) :: DEN(0:LMAXD1,IEMXD,NSPIN)
    doublecomplex, intent(in) :: WEZ(IEMXD)

    integer, intent(in) :: IEMXD
    integer, intent(in) :: NSPIN
    integer, intent(in) :: LMAXD1
    integer, intent(in) :: IELAST

    integer :: ISPIN
    integer :: L
    integer :: IE

    ! ---> l/m_s/atom-resolved charges

    do ISPIN = 1,NSPIN
      do L = 0,LMAXD1
        CHARGE(L,ISPIN) = 0.0D0

        do IE = 1,IELAST
          CHARGE(L,ISPIN) = CHARGE(L,ISPIN) + &
          DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))/ &
          DBLE(NSPIN)
        end do

      end do
    end do
  end subroutine calcChargesLres

  !----------------------------------------------------------------------------
  !> correct Fermi-energy (charge neutrality).
  subroutine doFermiEnergyCorrection(atomdata, output, naez, max_shift, CHRGNT, DENEF, R2NEF, &
                                     ESPV, RHO2NS, E2)
    use BasisAtom_mod
    use RadialMeshData_mod
    implicit none

    type (BasisAtom), intent(in) :: atomdata
    logical, intent(in) :: output !< output to stdout - yes/no
    integer, intent(in) :: naez
    double precision, intent(in) :: max_shift !< maximally allowed Fermi-Energy shift (good: 0.03d0)
    double precision, intent(in) :: CHRGNT
    double precision, intent(in) :: DENEF
    double precision, dimension(:,:,:), intent(in) ::  R2NEF

    ! in,out - arguments
    double precision, dimension(0:,:), intent(inout) :: ESPV
    double precision, dimension(:,:,:), intent(inout) ::  RHO2NS
    double precision, intent(inout) :: E2

    !-------- locals
    integer :: nspind
    type (RadialMeshData), pointer :: mesh

    double precision :: E2SHIFT
    double precision :: EFold
    double precision :: DF
    double precision :: PI
    integer :: ispin, lm, lmpotd

    PI = 4.0D0*ATAN(1.0D0)

    nspind = atomdata%nspin
    lmpotd = atomdata%potential%lmpot

    mesh => atomdata%mesh_ptr
! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),max_shift)*DSIGN(1.0D0,E2SHIFT)
        EFOLD = E2

        E2 = E2 - E2SHIFT

        if( output ) then
          call printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
        end if

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIND)
! ----------------------------------------------------------------------

        do ISPIN = 1,NSPIND

! -->     get correct density and valence band energies

          ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
          EFOLD*CHRGNT/DBLE(NSPIND*NAEZ)

          do LM = 1,LMPOTD
            call DAXPY(mesh%IRC,DF,R2NEF(1,LM,ISPIN),1, &
            RHO2NS(1,LM,ISPIN),1)
          end do
        end do
! ----------------------------------------------------------------------
  end subroutine

  !---------------------------------------------------------------------------
  !> Checks for file 'STOP' in current working directory, returns 1 if it
  !> exists, otherwise 0.
  !>
  !> Should be called only by master rank!!!
  integer function stopfile_flag()
    implicit none
    logical :: stopfile_exists
    stopfile_exists = .false.
    inquire(file='STOP',exist=stopfile_exists)
    stopfile_flag = 0
    if (stopfile_exists) then
      stopfile_flag = 1
    end if
  end function

  !----------------------------------------------------------------------------
  !> Print Fermi-Energy information to screen.
  subroutine printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
    implicit none
    double precision :: DENEF
    double precision :: E2
    double precision :: E2SHIFT
    double precision :: EFOLD
    integer :: NAEZ

    write (6,fmt=9020) EFOLD,E2SHIFT

    ! --> divided by NAEZ because the weight of each atom has been already
    !     taken into account in 1c

    write (6,fmt=9030) E2,DENEF/DBLE(NAEZ)
    write(6,'(79(1H+),/)')
9020 format ('                old', &
    ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
9030 format ('                new', &
    ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)
  end subroutine

end module ProcessKKRresults_mod
