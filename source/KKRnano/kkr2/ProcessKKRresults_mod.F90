#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

module ProcessKKRresults_mod

  implicit none

  public :: processKKRresults
  private :: calculateDensities
  private :: calculatePotentials

CONTAINS

subroutine processKKRresults(iter, calc_data, my_mpi, emesh, dims, params, arrays, &
                             program_timer)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use main2_aux_mod
  use EnergyMesh_mod

  use MadelungCalculator_mod, only: MadelungLatticeSum

  use CalculationData_mod

  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod

  use RadialMeshData_mod
  use BasisAtom_mod

  use LDAUData_mod

  use TimerMpi_mod
  use BroydenData_mod
  use BRYDBM_new_com_mod

  use wrappers_mod

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use KKRresults_mod
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
  type (MadelungLatticeSum), pointer                   :: madelung_sum
  type (ShapeGauntCoefficients), pointer               :: shgaunts
  type (GauntCoefficients), pointer                    :: gaunts
  type (BasisAtom) , pointer                           :: atomdata
  type (LDAUData)  , pointer                           :: ldau_data
  type (BroydenData), pointer                          :: broyden
  type (KKRresults) , pointer                          :: kkr
  type (DensityResults), pointer                       :: densities
  type (EnergyResults), pointer                        :: energies

  type (RadialMeshData), pointer :: mesh
  integer :: I1
  integer :: ierr
  double precision :: RMSAVQ ! rms error magnetisation dens. (contribution of single site)
  double precision :: RMSAVM ! rms error charge density (contribution of single site)
  logical, external :: testVFORM

  madelung_sum => getMadelungSum(calc_data, 1)
  shgaunts     => getShapeGaunts(calc_data)
  gaunts       => getGaunts(calc_data)
  atomdata     => getAtomData(calc_data, 1)
  ldau_data    => getLDAUData(calc_data, 1)
  broyden      => getBroyden(calc_data, 1)
  kkr          => getKKR(calc_data, 1)
  densities    => getDensities(calc_data, 1)
  energies     => getEnergies(calc_data, 1)

  mesh => atomdata%mesh_ptr

  I1 = atomdata%atom_index

  ! kkr
  !  |
  !  v
  call calculateDensities(iter, calc_data, my_mpi, dims, params, &
                          program_timer, arrays, emesh)
  ! |
  ! v
  ! densities, emesh, energies (ESPV only)
  ! |
  ! v
  call calculatePotentials(iter, my_mpi, dims, params, madelung_sum, shgaunts, &
                           program_timer, densities, energies, arrays, &
                           ldau_data, atomdata)
  ! |
  ! v
  ! atomdata, energies

! -->   calculation of RMS and final construction of the potentials (straight mixing)
  call MIXSTR_wrapper(atomdata, RMSAVQ, RMSAVM, params%MIXING, params%FCM)

  ! output of RMS error
  call RMSOUT_com(RMSAVQ,RMSAVM,ITER,dims%NSPIND,dims%NAEZ, &
                 getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi))

  ! it is weird that straight mixing is called in any case before
! -->   potential mixing procedures: Broyden or Andersen updating schemes
  if (params%IMIX>=3) then
    call BRYDBM_new_com(atomdata%potential%VISP,atomdata%potential%VONS,atomdata%potential%VINS, &
    atomdata%potential%LMPOT,mesh%R,mesh%DRDI,broyden%MIXING, &
    mesh%IRC,mesh%IRMIN,atomdata%potential%NSPIN, &
    broyden%IMIX,ITER, &
    broyden%UI2,broyden%VI2,broyden%WIT,broyden%SM1S,broyden%FM1S, &
    getMyAtomRank(my_mpi), &
    getMySEcommunicator(my_mpi), &
    broyden%itdbryd, mesh%irmd, atomdata%potential%irnsd, atomdata%potential%nspin)
  endif

  TESTARRAYLOG(3, atomdata%potential%VINS)
  TESTARRAYLOG(3, atomdata%potential%VISP)
  TESTARRAYLOG(3, atomdata%potential%VONS)

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  call resetPotentials(mesh%IRC, mesh%IRMD, mesh%IRMIN, &
                       atomdata%potential%IRMIND, atomdata%potential%LMPOT, &
                       atomdata%potential%NSPIN, atomdata%potential%VINS, &
                       atomdata%potential%VISP, atomdata%potential%VONS)

! ----------------------------------------------------- output_potential
  call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
  call writeBasisAtomPotentialDA(atomdata, 37, I1)
  call closeBasisAtomPotentialDAFile(37)
! ----------------------------------------------------- output_potential

! write formatted potential if file VFORM exists - contains bad inquire
! - bad check deactivated when KTE<0
  if (ITER == params%SCFSTEPS .and. params%KTE >= 0) then
    if (testVFORM()) then
      call writeFormattedPotential(emesh%E2, params%ALAT, energies%VBC, params%KXC, atomdata)
    endif
  endif

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
    call RESULTS(dims%LRECRES2,params%IELAST,ITER,arrays%LMAXD,arrays%NAEZ,emesh%NPOL, &
    dims%NSPIND,params%KPRE,params%KTE,arrays%LPOT,emesh%E1,emesh%E2,emesh%TK,emesh%EFERMI, &
    params%ALAT,atomdata%core%ITITLE(:,1:arrays%NSPIND),densities%total_charge_neutrality, &
    arrays%ZAT,emesh%EZ,emesh%WEZ,params%LDAU, &
    arrays%iemxd)

    call OUTTIME(isMasterRank(my_mpi),'results......',getElapsedTime(program_timer), iter)

    ! only MASTERRANK updates, other ranks get it broadcasted later
    ! (although other processes could update themselves)
    call updateEnergyMesh(emesh)

    ! write file 'energy_mesh'
    if (emesh%NPOL /= 0) emesh%EFERMI = emesh%E2  ! if not a DOS-calculation E2 coincides with Fermi-Energy

    call writeEnergyMesh(emesh)

    call printDoubleLineSep()

  endif
! -----------------------------------------------------------------
! END: only MASTERRANK is working here
! -----------------------------------------------------------------

end subroutine

!==============================================================================

!------------------------------------------------------------------------------
!> Calculate densities.
!>
!> output: emesh (Fermi-energy updated, renormalized weights), densities, ldau_data?, arrays
!> files written: 'results1'
subroutine calculateDensities(iter, calc_data, my_mpi, dims, params, &
                              program_timer, arrays, emesh)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use lloyds_formula_mod, only: renormalizeDOS

  use main2_aux_mod
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

  double complex, parameter      :: CZERO = (0.0d0, 0.0d0)
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  logical :: LdoRhoEF
  integer :: I1
  double precision :: denEf !< charge density at Fermi level
  double precision :: chrgNt !< charge neutrality

  shgaunts     => getShapeGaunts(calc_data)
  gaunts       => getGaunts(calc_data)
  atomdata     => getAtomData(calc_data, 1)
  ldau_data    => getLDAUData(calc_data, 1)
  kkr          => getKKR(calc_data, 1)
  densities    => getDensities(calc_data, 1)
  energies     => getEnergies(calc_data, 1)

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  I1 = atomdata%atom_index

  densities%CMOM   = 0.0D0
  densities%CMINST = 0.0D0

  ! out: emesh, RNORM
  call lloyd0_wrapper_com(atomdata, my_mpi, kkr%LLY_GRDT, emesh, densities%RNORM, &
                          dims%LLY, params%ICST, params%NSRA, kkr%GMATN, gaunts, ldau_data)

  if (dims%LLY == 1) then
    TESTARRAYLOG(3, emesh%WEZRN)
    TESTARRAYLOG(3, densities%RNORM)
    call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......',getElapsedTime(program_timer),ITER)
  endif

  ! now WEZRN stores the weights for E-integration

  densities%DEN = CZERO
  DENEF = 0.0D0
  CHRGNT = 0.0D0

  if (params%LDAU) then
    ldau_data%DMATLDAU = CZERO
  endif

  LDORHOEF = emesh%NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

  ! has to be done after Lloyd
  ! output: RHO2NS, R2NEF, DEN, ESPV
  call RHOVAL_wrapper(atomdata, LdoRhoEF, params%ICST, params%NSRA, densities%RHO2NS, densities%R2NEF, &
                      densities%DEN, energies%ESPV, kkr%GMATN, gaunts, emesh, ldau_data)

! ----------------------------------------------------------------------
! -->   determine total charge expanded in spherical harmonics
! -------------------------------------------------------------- density
  ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
  call RHOTOTB_wrapper(densities%CATOM, densities%RHO2NS, atomdata)

  CHRGNT = CHRGNT + densities%CATOM(1) - atomdata%Z_nuclear

  if (dims%LLY == 1) then
    call renormalizeDOS(densities%DEN,densities%RNORM,densities%LMAXD+1,densities%IEMXD,arrays%NSPIND,densities%IEMXD)
  end if

  ! calculate DOS at Fermi level
  DENEF = DENEF + calcDOSatFermi(densities%DEN, params%IELAST, densities%IEMXD, densities%LMAXD+1, densities%NSPIND)

  ! ---> l/m_s/atom-resolved charges, output -> CHARGE
  ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
  ! CHARGE -> written to result file
  call calcChargesLres(densities%CHARGE, densities%DEN, params%IELAST, densities%LMAXD+1, densities%NSPIND, emesh%WEZ, densities%IEMXD)

  call sumNeutralityDOSFermi_com(CHRGNT, DENEF, getMySEcommunicator(my_mpi))
  densities%total_charge_neutrality = CHRGNT

  ! write to 'results1' - only to be read in in results.f
  ! necessary for density of states calculation, otherwise
  ! only for informative reasons
  if (params%KTE >= 0) then
    call openResults1File(arrays%IEMXD, arrays%LMAXD, emesh%NPOL)
    call writeResults1File(densities%CATOM, densities%CHARGE, densities%DEN, &
                           atomdata%core%ECORE, I1, emesh%NPOL, atomdata%core%QC_corecharge)
    call closeResults1File()
  endif

  call OUTTIME(isMasterRank(my_mpi),'density calculated ..',getElapsedTime(program_timer),ITER)

  call doFermiEnergyCorrection(atomdata, isMasterRank(my_mpi), arrays%naez, &
                               0.03d0, CHRGNT, DENEF, densities%R2NEF, &
                               energies%ESPV, densities%RHO2NS, emesh%E2)

  !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)
  call RHOMOM_NEW_wrapper(densities%CMOM,densities%CMINST,densities%RHO2NS(:,:,1), cell, mesh, shgaunts)

  call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......',getElapsedTime(program_timer),ITER)

end subroutine


!------------------------------------------------------------------------------
!> Calculate potentials.
!>
!> Output: atomdata, ldau_data, arrays
!> Files written: 'results2'
subroutine calculatePotentials(iter, my_mpi, dims, params, madelung_sum, &
                               shgaunts, program_timer, densities, energies, &
                               arrays, ldau_data, atomdata)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use main2_aux_mod

  use MadelungCalculator_mod
  use ShapeGauntCoefficients_mod
  use muffin_tin_zero_mod

  use RadialMeshData_mod
  use BasisAtom_mod

  use LDAUData_mod

  use TimerMpi_mod

  use wrappers_mod

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use DensityResults_mod
  use EnergyResults_mod

  implicit none

  integer, intent(in)                        :: iter
  type (ShapeGauntCoefficients), intent(in)  :: shgaunts
  type (MadelungLatticeSum), intent(in)      :: madelung_sum
  type (KKRnanoParallel), intent(in)         :: my_mpi
  type (BasisAtom), intent(inout)            :: atomdata
  type (LDAUData), intent(inout)             :: ldau_data
  type (Main2Arrays), intent(in)             :: arrays
  type (DimParams), intent(in)               :: dims
  type (InputParams), intent(in)             :: params
  type (EnergyResults), intent(inout)        :: energies
  type (DensityResults), intent(inout)       :: densities
  type (TimerMpi), intent(in)                :: program_timer

  ! locals
  type (RadialMeshData), pointer :: mesh
  integer :: I1
  double precision :: VMAD
  integer :: lcoremax
  double precision :: VAV0, VOL0


  mesh => atomdata%mesh_ptr

  I1 = atomdata%atom_index
  VMAD = 0.0d0
  VAV0 = 0.0d0
  VOL0 = 0.0d0
  lcoremax = 0

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
  !output: VONS
  call VINTRAS_wrapper(densities%RHO2NS(:,:,1), shgaunts, atomdata)

  TESTARRAYLOG(3, atomdata%potential%VONS)
  TESTARRAYLOG(3, densities%RHO2NS)

  call OUTTIME(isMasterRank(my_mpi),'VINTRAS ......',getElapsedTime(program_timer),ITER)

  ! output: VONS (changed), VMAD
  call addMadelungPotential_com(madelung_sum, densities%CMOM, densities%CMINST, arrays%NSPIND, &
       atomdata%potential%VONS, arrays%ZAT, mesh%R, mesh%IRCUT, mesh%IPAN, VMAD, &
       getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi), getNumAtomRanks(my_mpi), mesh%irmd, mesh%ipand)

  call OUTTIME(isMasterRank(my_mpi),'VMADELBLK ......',getElapsedTime(program_timer),ITER)

! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

!            if (KFORCE==1 .and. ITER==SCFSTEPS) then
! !---------------------------------------------------------------------
!              call FORCEH(CMOM,FLM,LPOT,I1,RHO2NS,VONS, &
!              R,DRDI,IMT,ZAT,irmd)  ! TODO: get rid of atom parameter I1
!              call FORCE(FLM,FLMC,LPOT,NSPIND,I1,RHOCAT,VONS,R, &
!              DRDI,IMT,naez,irmd)
! !---------------------------------------------------------------------
!            end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

  if (params%KTE==1) then
    ! calculate total energy and individual contributions if requested
    ! core electron contribution
    call ESPCB_wrapper(energies%ESPC, LCOREMAX, atomdata)
    ! output: EPOTIN
    call EPOTINB_wrapper(energies%EPOTIN,densities%RHO2NS,atomdata)
    ! output: ECOU - l resolved Coulomb energy
    call ECOUB_wrapper(densities%CMOM, energies%ECOU, densities%RHO2NS, shgaunts, atomdata)

    call OUTTIME(isMasterRank(my_mpi),'KTE ......',getElapsedTime(program_timer),ITER)
  end if

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
  ! output: VONS (changed), EXC (exchange energy) (l-resolved)
  call VXCDRV_wrapper(energies%EXC, params%KXC, densities%RHO2NS, shgaunts, atomdata)

  call OUTTIME(isMasterRank(my_mpi),'VXCDRV ......',getElapsedTime(program_timer),ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

!            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
!              call FORCXC_com(FLM,FLMC,LPOT,NSPIND,I1,RHOCAT,VONS,R, &
!              ALAT,DRDI,IMT,ZAT, &
!              getMyAtomRank(my_mpi), &
!              getMySEcommunicator(my_mpi), &
!              naez, irmd)
! ---------------------------------------------------------------------
!            end if

  ! unnecessary I/O? see results.f
  if (params%KTE >= 0) then
    call openResults2File(dims%LRECRES2)
    call writeResults2File(densities%CATOM, energies%ECOU, ldau_data%EDCLDAU, &
                           energies%EPOTIN, energies%ESPC, energies%ESPV, ldau_data%EULDAU, &
                           energies%EXC, I1, LCOREMAX, VMAD)
    call closeResults2File()
  end if

  ! calculate new muffin-tin zero. output: VAV0, VOL0
  call MTZERO_wrapper(VAV0, VOL0, atomdata)

  call OUTTIME(isMasterRank(my_mpi),'MTZERO ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

  call allreduceMuffinTinShift_com(getMySEcommunicator(my_mpi), VAV0, energies%VBC, VOL0)

  if(isMasterRank(my_mpi)) then
    call printMuffinTinShift(VAV0, energies%VBC, VOL0)
  end if

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

! -->   shift potential by VBC and multiply with shape functions - output: VONS
  call CONVOL_wrapper(energies%VBC, shgaunts, atomdata)

! LDAU
  ldau_data%EULDAU = 0.0D0
  ldau_data%EDCLDAU = 0.0D0

  if (ldau_data%LDAU.and.ldau_data%NLDAU>=1) then
    call LDAUWMAT(I1,ldau_data%NSPIND,ITER,params%MIXING,ldau_data%DMATLDAU,ldau_data%NLDAU,ldau_data%LLDAU, &
                  ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%UMLDAU,ldau_data%WMLDAU,ldau_data%EULDAU,ldau_data%EDCLDAU, &
                  ldau_data%lmaxd)
  endif
! LDAU

  call OUTTIME(isMasterRank(my_mpi),'calculated pot ......',getElapsedTime(program_timer),ITER)

end subroutine

end module ProcessKKRresults_mod
