! KKRnano
! massive parallel KKR for nanoscaled systems

#include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"

program MAIN2

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use common_testc
  use common_optc

  use KKRnanoParallel_mod
  use KKRnano_Comm_mod

  use lloyds_formula_mod

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

  use JijData_mod
  use LDAUData_mod

  use TimerMpi_mod
  use EBalanceHandler_mod
  use BroydenData_mod
  use BRYDBM_new_com_mod

  use wrappers_mod

  use TEST_lcutoff_mod !TODO: remove

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use KKRresults_mod

  use ScatteringCalculation_mod, only: energyloop

  implicit none
  include 'mpif.h'

  type (MadelungCalculator) :: madelung_calc
  type (ShapeGauntCoefficients) :: shgaunts
  type (GauntCoefficients) :: gaunts

  !     .. Parameters ..
  double complex, parameter :: CZERO = (0.0D0,0.0D0)

  !     ..
  !     .. Local Scalars ..

  double precision::DENEF
  double precision::CHRGNT

  double precision::RMSAVM      ! rms error magnetisation dens. (contribution of single site)
  double precision::RMSAVQ      ! rms error charge density (contribution of single site)

  type (TimerMpi) :: program_timer
  type (TimerMpi) :: iteration_timer

  type (EBalanceHandler) :: ebalance_handler

  integer::ITER
  integer::NOITER
  integer::NOITER_ALL
  integer::I1
  logical::XCCPL
  logical::LDORHOEF

  integer :: BCP ! TODO: remove - is dummy
  integer :: IGUESS ! TODO: remove - is dummy

  !     .. Local Arrays ..

  double precision::VAV0
  double precision::VOL0

  double precision::EPOTIN
  double precision::VMAD

  integer::LCOREMAX

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

  type(KKRnanoParallel) :: my_mpi

  integer :: flag
  logical, external :: testVFORM

  type (RadialMeshData), target :: mesh
  type (CellData), target       :: cell
  type (BasisAtom), target      :: atomdata
  type (EnergyMesh), target     :: emesh
  type (LDAUData), target       :: ldau_data
  type (JijData), target        :: jij_data
  type (BroydenData), target    :: broyden
  type (Main2Arrays), target    :: arrays
  type (DimParams), target      :: dims
  type (InputParams)            :: params
  type (KKRresults)             :: kkr

  call createDimParams(dims) ! read dim. parameters from 'inp0.unf'

! -----------------------------------------------------------------------------
  call createKKRnanoParallel(my_mpi, dims%NAEZ, dims%SMPID, dims%EMPID)
  call setKKRnanoNumThreads(dims%nthrds)
  call printKKRnanoInfo(my_mpi, dims%nthrds)
!------------------------------------------------------------------------------

  if (getMyWorldRank(my_mpi) < 128) then ! max. 128 logfiles
    OPENLOG(getMyWorldRank(my_mpi), 3)
  else
    OPENLOG(getMyWorldRank(my_mpi), 0)
  endif

! ========= TIMING =========================================================
    call resetTimer(program_timer)
    if (isMasterRank(my_mpi)) then
      open (2,file='time-info',form='formatted')
    endif
!========= TIMING END ======================================================

!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call createMain2Arrays(arrays, dims)
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  !every process does this!
  call readKKR0InputNew(dims%NSYMAXD, params%ALAT, arrays%ATOM, BCP, arrays%BRAVAIS, &
                        arrays%CLS, arrays%DSYMLL, arrays%EZOA, params%FCM, params%GMAX, params%ICST, &
                        IGUESS, params%IMIX, arrays%INDN0, &
                        arrays%ISYMINDEX, &
                        params%JIJ, params%KFORCE, arrays%KMESH, params%KPRE, params%KTE, params%KXC, &
                        params%LDAU, params%MAXMESH, &
                        params%MIXING, arrays%NACLS, params%NCLS, params%NR, params%NREF, &
                        params%NSRA, params%NSYMAT, arrays%NUMN0, OPTC, params%QMRBOUND, &
                        arrays%RBASIS, arrays%RCLS, params%RCUTJIJ, arrays%REFPOT, params%RMAX, arrays%RMTREF, &
                        arrays%RR, params%SCFSTEPS, TESTC, arrays%VREF, arrays%ZAT)


  !if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  call consistencyCheck03(arrays%ATOM, arrays%CLS, arrays%EZOA, arrays%INDN0, &
                          arrays%NACLS, arrays%NACLSD, arrays%NAEZ, arrays%NCLSD, params%NR, arrays%NUMN0)

  if ((params%JIJ .eqv. .true.) .and. (arrays%nspind /= 2)) then
    write(*,*) "ERROR: Jij calculation not possible for spin-unpolarized calc."
    stop
  end if

!=====================================================================
!     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
! ... and wait after SC-ITER loop
!=====================================================================

  ! This if closes several hundreds of lines later!
  if (isActiveRank(my_mpi)) then

!+++++++++++ pre self-consistency preparation

    I1 = getMyAtomId(my_mpi) !assign atom number for the rest of the program

    call createKKRresults(kkr, dims)

    ! ---------------------------------------------------------- k_mesh
    call readKpointsFile(arrays%BZKP, params%MAXMESH, arrays%NOFKS, arrays%VOLBZ, arrays%VOLCUB)  !every process does this!

    call OUTTIME(isMasterRank(my_mpi),'input files read.....', &
                                       getElapsedTime(program_timer), 0)

    call createBasisAtom(atomdata, I1, dims%lpot, dims%nspind, dims%irmind, dims%irmd)
    call openBasisAtomDAFile(atomdata, 37, "atoms")
    call readBasisAtomDA(atomdata, 37, I1)
    call closeBasisAtomDAFile(37)

    if (isInMasterGroup(my_mpi)) then
      call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
      call readBasisAtomPotentialDA(atomdata, 37, I1)
      call closeBasisAtomPotentialDAFile(37)
    end if

    call createCellData(cell, dims%irid, (2*dims%LPOT+1)**2, dims%nfund)
    call openCellDataDAFile(cell, 37 , "cells")
    call readCellDataDA(cell, 37, getCellIndex(atomdata))
    call closeCellDataDAFile(37)

    call associateBasisAtomCell(atomdata, cell)

    call createRadialMeshData(mesh, dims%irmd, dims%ipand)
    call openRadialMeshDataDAFile(mesh, 37 , "meshes")
    call readRadialMeshDataDA(mesh, 37, I1)
    call closeRadialMeshDataDAFile(37)

    call associateBasisAtomMesh(atomdata, mesh)

    call createLDAUData(ldau_data, params%ldau, dims%irmd, dims%lmaxd, dims%nspind)
    call createJijData(jij_data, params%jij, params%rcutjij, dims%nxijd,dims%lmmaxd,dims%nspind)
    call createBroydenData(broyden, dims%ntird, dims%itdbryd, params%imix, params%mixing)

    call createEnergyMesh(emesh, dims%iemxd)
    params%ielast = dims%iemxd

    call readEnergyMesh(emesh)  !every process does this!

    call createMadelungCalculator(madelung_calc, dims%lmaxd, params%ALAT, params%RMAX, params%GMAX, &
                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)

    call calculateMadelungLatticeSum(madelung_calc, dims%naez, I1, arrays%rbasis, arrays%smat)

    call OUTTIME(isMasterRank(my_mpi),'Madelung sums calc...',getElapsedTime(program_timer), 0)

    call createGauntCoefficients(gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(shgaunts, dims%lmaxd)

    call createEBalanceHandler(ebalance_handler, params%ielast)
    call initEBalanceHandler(ebalance_handler, my_mpi)
    call setEqualDistribution(ebalance_handler, (emesh%NPNT1 == 0))

    call initLcutoff(arrays%rbasis, arrays%bravais, arrays%lmmaxd, I1) !TODO: remove
    WRITELOG(3, *) "lm-array: ", lmarray

!+++++++++++
    ASSERT( arrays%ZAT(I1) == atomdata%Z_nuclear )

   !flag = 0
   !99 continue
   !if (flag == 0) goto 99

! ######################################################################
! ######################################################################
    do ITER = 1, params%SCFSTEPS
! ######################################################################
! ######################################################################

      call resetTimer(iteration_timer)

      NOITER = 0

      if (isMasterRank(my_mpi)) then
        call printDoubleLineSep(unit_number = 2)
        call OUTTIME(isMasterRank(my_mpi),'started at ..........', &
                     getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      arrays%CMOM   = 0.0D0
      arrays%CMINST = 0.0D0
             CHRGNT = 0.0D0

      WRITELOG(2, *) "Iteration Atom ", ITER, I1

!=======================================================================
! xccpl

      XCCPL = .false.

      ! calculate exchange couplings only at last self-consistency step and when Jij=true
      if ((ITER==params%SCFSTEPS).and.params%JIJ) XCCPL = .true.

      if (XCCPL) then

        call CLSJIJ(I1,dims%NAEZ,arrays%RR,params%NR,arrays%RBASIS,jij_data%RCUTJIJ,params%NSYMAT,arrays%ISYMINDEX, &
                    jij_data%IXCP,jij_data%NXCP,jij_data%NXIJ,jij_data%RXIJ,jij_data%RXCCLS,jij_data%ZKRXIJ, &
                    arrays%nrd, jij_data%nxijd)

        jij_data%JXCIJINT = CZERO
        jij_data%GMATXIJ = CZERO

      endif

      ! New: instead of reading potential every time, communicate it
      call communicatePotential(my_mpi, atomdata%potential%VISP, atomdata%potential%VINS, atomdata%core%ECORE)

      ! Core relaxation - only mastergroup needs results
      if (isInMasterGroup(my_mpi)) then
        call RHOCORE_wrapper(emesh%E1, params%NSRA, atomdata)
      endif

! LDA+U
      if (params%LDAU) then

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

        call LDAUINIT(I1,ITER,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU,ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)

      endif
! LDA+U

      call OUTTIME(isMasterRank(my_mpi),'initialized .........', &
                   getElapsedTime(program_timer),ITER)

      ! Scattering calculations - that is what KKR is all about
      call energyLoop(iter, atomdata, emesh, params, dims, gaunts, &
                      ebalance_handler, my_mpi, arrays, kkr, jij_data, ldau_data)

      call OUTTIME(isMasterRank(my_mpi),'G obtained ..........',getElapsedTime(program_timer),ITER)

!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
      if (isInMasterGroup(my_mpi)) then

        ! out: emesh, RNORM
        call lloyd0_wrapper_com(atomdata, my_mpi, kkr%LLY_GRDT, emesh, arrays%RNORM, &
                                dims%LLY, params%ICST, params%NSRA, arrays%GMATN, gaunts, ldau_data)

        if (dims%LLY == 1) then
          TESTARRAYLOG(3, emesh%WEZRN)
          TESTARRAYLOG(3, arrays%RNORM)
          call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......',getElapsedTime(program_timer),ITER)
        endif

        ! now WEZRN stores the weights for E-integration

        arrays%DEN = CZERO
        DENEF = 0.0D0

        if (params%LDAU) then
          ldau_data%DMATLDAU = CZERO
        endif

        LDORHOEF = emesh%NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

        ! has to be done after Lloyd
        ! output: RHO2NS, R2NEF, DEN, ESPV
        call RHOVAL_wrapper(atomdata, LdoRhoEF, params%ICST, params%NSRA, arrays%RHO2NS, arrays%R2NEF, &
                            arrays%DEN, arrays%ESPV, arrays%GMATN, gaunts, emesh, ldau_data)

! ----------------------------------------------------------------------
! -->   determine total charge expanded in spherical harmonics
! -------------------------------------------------------------- density
        ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
        call RHOTOTB_wrapper(arrays%CATOM, arrays%RHO2NS, atomdata)

        CHRGNT = CHRGNT + arrays%CATOM(1) - atomdata%Z_nuclear

        if (dims%LLY == 1) then
          call renormalizeDOS(arrays%DEN,arrays%RNORM,arrays%LMAXD+1,params%IELAST,arrays%NSPIND,arrays%IEMXD)
        end if

        ! calculate DOS at Fermi level
        DENEF = calcDOSatFermi(arrays%DEN, params%IELAST, arrays%IEMXD, arrays%LMAXD+1, arrays%NSPIND)

        ! ---> l/m_s/atom-resolved charges, output -> CHARGE
        ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
        ! CHARGE -> written to result file
        call calcChargesLres(arrays%CHARGE, arrays%DEN, params%IELAST, arrays%LMAXD+1, arrays%NSPIND, emesh%WEZ, arrays%IEMXD)

        call sumNeutralityDOSFermi_com(CHRGNT, DENEF, getMySEcommunicator(my_mpi))

        ! write to 'results1' - only to be read in in results.f
        ! necessary for density of states calculation, otherwise
        ! only for informative reasons
        if (params%KTE >= 0) then
          call openResults1File(arrays%IEMXD, arrays%LMAXD, emesh%NPOL)
          call writeResults1File(arrays%CATOM, arrays%CHARGE, arrays%DEN, atomdata%core%ECORE, I1, emesh%NPOL, atomdata%core%QC_corecharge)
          call closeResults1File()
        endif

        call OUTTIME(isMasterRank(my_mpi),'density calculated ..',getElapsedTime(program_timer),ITER)

        call doFermiEnergyCorrection(atomdata, isMasterRank(my_mpi), arrays%naez, 0.03d0, CHRGNT, DENEF, arrays%R2NEF, &
                                     arrays%ESPV, arrays%RHO2NS, emesh%E2)

        !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)
        call RHOMOM_NEW_wrapper(arrays%CMOM,arrays%CMINST,arrays%RHO2NS(:,:,1), cell, mesh, shgaunts)

        call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
        !output: VONS
        call VINTRAS_wrapper(arrays%RHO2NS(:,:,1), shgaunts, atomdata)

        TESTARRAYLOG(3, atomdata%potential%VONS)
        TESTARRAYLOG(3, arrays%RHO2NS)

        call OUTTIME(isMasterRank(my_mpi),'VINTRAS ......',getElapsedTime(program_timer),ITER)

        ! output: VONS (changed), VMAD
        call addMadelungPotential_com(madelung_calc, arrays%CMOM, arrays%CMINST, arrays%NSPIND, &
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
          call EPOTINB_wrapper(EPOTIN,arrays%RHO2NS,atomdata)
          ! output: ECOU - l resolved Coulomb energy
          call ECOUB_wrapper(arrays%CMOM, arrays%ECOU, arrays%RHO2NS, shgaunts, atomdata)
        end if

        call OUTTIME(isMasterRank(my_mpi),'KTE ......',getElapsedTime(program_timer),ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
        ! output: VONS (changed), EXC (exchange energy) (l-resolved)
        call VXCDRV_wrapper(arrays%EXC,params%KXC,arrays%RHO2NS, shgaunts, atomdata)

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
          call writeResults2File(arrays%CATOM, arrays%ECOU, ldau_data%EDCLDAU, EPOTIN, arrays%ESPC, arrays%ESPV, ldau_data%EULDAU, arrays%EXC, I1, LCOREMAX, VMAD)
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

      endif
! -----------------------------------------------------------------
! END: only processes in MASTERGROUP are working here
!-----------------------------------------------------------------

! -----------------------------------------------------------------
! BEGIN: only MASTERRANK is working here
! -----------------------------------------------------------------
      if(isMasterRank(my_mpi)) then

        ! DOS was written to file 'results1' and read out here just
        ! to be written in routine wrldos
        ! also other stuff is read from results1 (and results2)
        call RESULTS(dims%LRECRES2,params%IELAST,ITER,arrays%LMAXD,arrays%NAEZ,emesh%NPOL, &
        dims%NSPIND,params%KPRE,params%KTE,arrays%LPOT,emesh%E1,emesh%E2,emesh%TK,emesh%EFERMI, &
        params%ALAT,atomdata%core%ITITLE(:,1:arrays%NSPIND),CHRGNT,arrays%ZAT,emesh%EZ,emesh%WEZ,params%LDAU, &
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
      call broadcastEnergyMesh_com(my_mpi, emesh)

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      getMyActiveCommunicator(my_mpi),IERR) ! TODO: allreduce not necessary, only master rank needs NOITER_ALL, use reduce instead

      if(isMasterRank(my_mpi)) then
        call printSolverIterationNumber(ITER, NOITER_ALL)
        call writeIterationTimings(ITER, getElapsedTime(program_timer), getElapsedTime(iteration_timer))
      endif

! manual exit possible by creation of file 'STOP' in home directory
      if (isManualAbort_com(getMyWorldRank(my_mpi), &
          getMyActiveCommunicator(my_mpi)) .eqv. .true.) exit

! ######################################################################
! ######################################################################
    enddo          ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################

    if (isMasterRank(my_mpi)) close(2)    ! TIME

    !if (KFORCE==1) close(54)

    call destroyMadelungCalculator(madelung_calc)
    call destroyGauntCoefficients(gaunts)
    call destroyShapeGauntCoefficients(shgaunts)

    call destroyEBalanceHandler(ebalance_handler)

    call destroyBasisAtom(atomdata)
    call destroyCellData(cell)
    call destroyRadialMeshData(mesh)

    call destroyBroydenData(broyden)
    call destroyLDAUData(ldau_data)
    call destroyJijData(jij_data)
    call destroyKKRresults(kkr)

! ======================================================================

  endif ! active Ranks

  CLOSELOG

!------------------------------------------------------------------------------
  call destroyMain2Arrays(arrays)
  call destroyDimParams(dims)
!------------------------------------------------------------------------------

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================
! Free KKRnano mpi resources

  call destroyKKRnanoParallel(my_mpi)

end program MAIN2
