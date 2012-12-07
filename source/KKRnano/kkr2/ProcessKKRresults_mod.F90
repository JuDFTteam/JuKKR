#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

module ProcessKKRresults_mod

  implicit none

  public :: processKKRresults
  private :: calculateDensities

CONTAINS

subroutine processKKRresults(iter, kkr, my_mpi, atomdata, emesh, dims, params, arrays, gaunts, shgaunts, madelung_calc, program_timer, &
                             densities, broyden, ldau_data)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use lloyds_formula_mod, only: renormalizeDOS

  use main2_aux_mod
  use muffin_tin_zero_mod
  use EnergyMesh_mod

  use MadelungCalculator_mod
  use lloyd0_new_mod

  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod

  use RadialMeshData_mod
  use CellData_mod
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

  implicit none

  include 'mpif.h'

  integer, intent(in)           :: iter
  type (MadelungCalculator)     :: madelung_calc
  type (ShapeGauntCoefficients) :: shgaunts
  type (GauntCoefficients)      :: gaunts
  type (KKRnanoParallel)        :: my_mpi
  type (BasisAtom)              :: atomdata
  type (EnergyMesh)             :: emesh
  type (LDAUData)               :: ldau_data
  type (BroydenData)            :: broyden
  type (Main2Arrays)            :: arrays
  type (DimParams)              :: dims
  type (InputParams)            :: params
  type (KKRresults)             :: kkr
  type (DensityResults)         :: densities
  type (TimerMpi)               :: program_timer

  ! locals
  double complex, parameter      :: CZERO = (0.0d0, 0.0d0)
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  integer :: I1
  double precision :: VMAD
  integer :: ierr
  integer :: lcoremax
  double precision :: EPOTIN, VAV0, VOL0
  double precision :: RMSAVQ ! rms error magnetisation dens. (contribution of single site)
  double precision :: RMSAVM ! rms error charge density (contribution of single site)
  logical, external :: testVFORM

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  VMAD = 0.0d0
  I1 = atomdata%atom_index

  call calculateDensities(iter, my_mpi, atomdata, dims, params, gaunts, &
                          shgaunts, kkr, program_timer, &
                          ldau_data, arrays, emesh, densities)

  call calculatePotentials(iter, my_mpi, dims, params, madelung_calc, shgaunts, &
                           program_timer, densities, arrays, &
                           ldau_data, atomdata)

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

  call resetPotentials(mesh%IRC, mesh%IRMD, mesh%IRMIN, atomdata%potential%IRMIND, atomdata%potential%LMPOT, &
                       atomdata%potential%NSPIN, atomdata%potential%VINS, atomdata%potential%VISP, atomdata%potential%VONS) ! Note: only LMPIC=1 processes

! ----------------------------------------------------- output_potential
  call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
  call writeBasisAtomPotentialDA(atomdata, 37, I1)
  call closeBasisAtomPotentialDAFile(37)
! ----------------------------------------------------- output_potential

! write formatted potential if file VFORM exists - contains bad inquire
! - bad check deactivated when KTE<0
  if (ITER == params%SCFSTEPS .and. params%KTE >= 0) then
    if (testVFORM()) then
      call writeFormattedPotential(emesh%E2, params%ALAT, arrays%VBC, params%KXC, atomdata)
    endif
  endif

! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS

  call MPI_BARRIER(getMySEcommunicator(my_mpi),IERR)

! -----------------------------------------------------------------
! BEGIN: only MASTERRANK is working here
! -----------------------------------------------------------------
  if(isMasterRank(my_mpi)) then

    ! DOS was written to file 'results1' and read out here just
    ! to be written in routine wrldos
    ! also other stuff is read from results1 (and results2)
    call RESULTS(dims%LRECRES2,params%IELAST,ITER,arrays%LMAXD,arrays%NAEZ,emesh%NPOL, &
    dims%NSPIND,params%KPRE,params%KTE,arrays%LPOT,emesh%E1,emesh%E2,emesh%TK,emesh%EFERMI, &
    params%ALAT,atomdata%core%ITITLE(:,1:arrays%NSPIND),densities%CHRGNT,arrays%ZAT,emesh%EZ,emesh%WEZ,params%LDAU, &
    arrays%iemxd)

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
subroutine calculateDensities(iter, my_mpi, atomdata, dims, params, gaunts, shgaunts, kkr, program_timer, &
                              ldau_data, arrays, emesh, densities)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use KKRnanoParallel_mod

  use lloyds_formula_mod, only: renormalizeDOS

  use main2_aux_mod
  use EnergyMesh_mod

  use lloyd0_new_mod

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

  implicit none

  integer, intent(in)                        :: iter
  type (ShapeGauntCoefficients), intent(in)  :: shgaunts
  type (GauntCoefficients), intent(in)       :: gaunts
  type (KKRnanoParallel), intent(in)         :: my_mpi
  type (BasisAtom), intent(inout)            :: atomdata
  type (EnergyMesh), intent(inout)           :: emesh
  type (LDAUData), intent(inout)             :: ldau_data
  type (Main2Arrays), intent(inout)          :: arrays
  type (DimParams), intent(in)               :: dims
  type (InputParams), intent(in)             :: params
  type (KKRresults), intent(inout)           :: kkr  ! should be in only
  type (DensityResults), intent(inout)       :: densities
  type (TimerMpi), intent(in)                :: program_timer

  ! locals
  double complex, parameter      :: CZERO = (0.0d0, 0.0d0)
  type (RadialMeshData), pointer :: mesh
  type (CellData), pointer       :: cell
  logical :: LdoRhoEF
  integer :: I1

  mesh => atomdata%mesh_ptr
  cell => atomdata%cell_ptr

  I1 = atomdata%atom_index

  ! out: emesh, RNORM
  call lloyd0_wrapper_com(atomdata, my_mpi, kkr%LLY_GRDT, emesh, arrays%RNORM, &
                          dims%LLY, params%ICST, params%NSRA, kkr%GMATN, gaunts, ldau_data)

  if (dims%LLY == 1) then
    TESTARRAYLOG(3, emesh%WEZRN)
    TESTARRAYLOG(3, arrays%RNORM)
    call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......',getElapsedTime(program_timer),ITER)
  endif

  ! now WEZRN stores the weights for E-integration

  densities%DEN = CZERO
  densities%DENEF = 0.0D0

  if (params%LDAU) then
    ldau_data%DMATLDAU = CZERO
  endif

  LDORHOEF = emesh%NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

  ! has to be done after Lloyd
  ! output: RHO2NS, R2NEF, DEN, ESPV
  call RHOVAL_wrapper(atomdata, LdoRhoEF, params%ICST, params%NSRA, densities%RHO2NS, densities%R2NEF, &
                      densities%DEN, arrays%ESPV, kkr%GMATN, gaunts, emesh, ldau_data)

! ----------------------------------------------------------------------
! -->   determine total charge expanded in spherical harmonics
! -------------------------------------------------------------- density
  ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
  call RHOTOTB_wrapper(densities%CATOM, densities%RHO2NS, atomdata)

  densities%CHRGNT = densities%CHRGNT + densities%CATOM(1) - atomdata%Z_nuclear

  if (dims%LLY == 1) then
    call renormalizeDOS(densities%DEN,arrays%RNORM,densities%LMAXD+1,densities%IEMXD,arrays%NSPIND,densities%IEMXD)
  end if

  ! calculate DOS at Fermi level
  densities%DENEF = calcDOSatFermi(densities%DEN, params%IELAST, densities%IEMXD, densities%LMAXD+1, densities%NSPIND)

  ! ---> l/m_s/atom-resolved charges, output -> CHARGE
  ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
  ! CHARGE -> written to result file
  call calcChargesLres(densities%CHARGE, densities%DEN, params%IELAST, densities%LMAXD+1, densities%NSPIND, emesh%WEZ, densities%IEMXD)

  call sumNeutralityDOSFermi_com(densities%CHRGNT, densities%DENEF, getMySEcommunicator(my_mpi))

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
                               0.03d0, densities%CHRGNT, densities%DENEF, densities%R2NEF, &
                               arrays%ESPV, densities%RHO2NS, emesh%E2)

  !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)
  call RHOMOM_NEW_wrapper(densities%CMOM,densities%CMINST,densities%RHO2NS(:,:,1), cell, mesh, shgaunts)

  call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......',getElapsedTime(program_timer),ITER)

end subroutine


!------------------------------------------------------------------------------
subroutine calculatePotentials(iter, my_mpi, dims, params, madelung_calc, shgaunts, program_timer, densities, arrays, ldau_data, atomdata)

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

  implicit none

  integer, intent(in)                        :: iter
  type (ShapeGauntCoefficients), intent(in)  :: shgaunts
  type (MadelungCalculator), intent(inout)   :: madelung_calc  ! should be 'in' only
  type (KKRnanoParallel), intent(in)         :: my_mpi
  type (BasisAtom), intent(inout)            :: atomdata
  type (LDAUData), intent(inout)             :: ldau_data
  type (Main2Arrays), intent(inout)          :: arrays
  type (DimParams), intent(in)               :: dims
  type (InputParams), intent(in)             :: params

  type (DensityResults), intent(inout)       :: densities
  type (TimerMpi), intent(in)                :: program_timer

  ! locals
  type (RadialMeshData), pointer :: mesh
  integer :: I1
  double precision :: VMAD
  integer :: lcoremax
  double precision :: EPOTIN, VAV0, VOL0

  mesh => atomdata%mesh_ptr

  I1 = atomdata%atom_index
  VMAD = 0.0d0
  EPOTIN = 0.0d0
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
  call addMadelungPotential_com(madelung_calc, densities%CMOM, densities%CMINST, arrays%NSPIND, &
       arrays%NAEZ, atomdata%potential%VONS, arrays%ZAT, mesh%R, mesh%IRCUT, mesh%IPAN, VMAD, &
       arrays%SMAT, getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi), getNumAtomRanks(my_mpi), mesh%irmd, mesh%ipand)

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
    call ESPCB_wrapper(arrays%ESPC, LCOREMAX, atomdata)
    ! output: EPOTIN
    call EPOTINB_wrapper(EPOTIN,densities%RHO2NS,atomdata)
    ! output: ECOU - l resolved Coulomb energy
    call ECOUB_wrapper(densities%CMOM, arrays%ECOU, densities%RHO2NS, shgaunts, atomdata)
  end if

  call OUTTIME(isMasterRank(my_mpi),'KTE ......',getElapsedTime(program_timer),ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
  ! output: VONS (changed), EXC (exchange energy) (l-resolved)
  call VXCDRV_wrapper(arrays%EXC,params%KXC,densities%RHO2NS, shgaunts, atomdata)

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
    call writeResults2File(densities%CATOM, arrays%ECOU, ldau_data%EDCLDAU, &
                           EPOTIN, arrays%ESPC, arrays%ESPV, ldau_data%EULDAU, &
                           arrays%EXC, I1, LCOREMAX, VMAD)
    call closeResults2File()
  end if

  ! calculate new muffin-tin zero. output: VAV0, VOL0
  call MTZERO_wrapper(VAV0, VOL0, atomdata)

  call OUTTIME(isMasterRank(my_mpi),'MTZERO ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

  call OUTTIME(isMasterRank(my_mpi),'calculated pot ......',getElapsedTime(program_timer),ITER)

  call allreduceMuffinTinShift_com(getMySEcommunicator(my_mpi), VAV0, arrays%VBC, VOL0)

  if(isMasterRank(my_mpi)) then
    call printMuffinTinShift(VAV0, arrays%VBC, VOL0)
  end if

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

! -->   shift potential by VBC and multiply with shape functions - output: VONS
  call CONVOL_wrapper(arrays%VBC, shgaunts, atomdata)

! LDAU
  ldau_data%EULDAU = 0.0D0
  ldau_data%EDCLDAU = 0.0D0

  if (ldau_data%LDAU.and.ldau_data%NLDAU>=1) then
    call LDAUWMAT(I1,ldau_data%NSPIND,ITER,params%MIXING,ldau_data%DMATLDAU,ldau_data%NLDAU,ldau_data%LLDAU, &
                  ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%UMLDAU,ldau_data%WMLDAU,ldau_data%EULDAU,ldau_data%EDCLDAU, &
                  ldau_data%lmaxd)
  endif
! LDAU

end subroutine

end module ProcessKKRresults_mod
