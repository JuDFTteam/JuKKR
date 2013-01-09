! KKRnano
! massive parallel KKR for nanoscaled systems

#include "DebugHelpers/test_macros.h"
#include "DebugHelpers/logging_macros.h"

program MAIN2

  USE_LOGGING_MOD

  use common_testc
  use common_optc

  use KKRnanoParallel_mod
  use KKRnano_Comm_mod

  use main2_aux_mod
  use EnergyMesh_mod

  use MadelungCalculator_mod

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

  use wrappers_mod, only: rhocore_wrapper

  use TEST_lcutoff_mod !TODO: remove

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod
  use KKRresults_mod
  use DensityResults_mod

  use ScatteringCalculation_mod, only: energyloop
  use ProcessKKRresults_mod

  use CalculationData_mod

  implicit none

  type (CalculationData) :: calc_data

  type (MadelungCalculator), target :: madelung_calc
  type (MadelungLatticeSum) :: madelung_sum
  type (ShapeGauntCoefficients) :: shgaunts
  type (GauntCoefficients) :: gaunts

  !     .. Local Scalars ..

  type (TimerMpi) :: program_timer
  type (TimerMpi) :: iteration_timer

  type (EBalanceHandler) :: ebalance_handler

  integer::ITER
  integer::I1

  type(KKRnanoParallel) :: my_mpi

  integer :: flag

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
  type (DensityResults)         :: densities

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
  call readKKR0InputNew(dims%NSYMAXD, params%ALAT, arrays%ATOM, arrays%BRAVAIS, &
                        arrays%CLS, arrays%DSYMLL, arrays%EZOA, params%FCM, params%GMAX, params%ICST, &
                        params%IMIX, arrays%INDN0, &
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

  call consistencyCheck03(arrays%ATOM, arrays%CLS, arrays%EZOA, &
                          arrays%INDN0, arrays%NACLS, arrays%NACLSD, &
                          arrays%NAEZ, arrays%NCLSD, params%NR, arrays%NUMN0)

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

    call createCalculationData(calc_data, dims, params, arrays, my_mpi)
    call prepareMadelung(calc_data, arrays)

    call createKKRresults(kkr, dims)
    call createDensityResults(densities, dims)

    ! ---------------------------------------------------------- k_mesh
    call readKpointsFile(arrays%BZKP, params%MAXMESH, arrays%NOFKS, &
                         arrays%VOLBZ, arrays%VOLCUB)  !every process does this!

    I1 = getMyAtomId(my_mpi) !assign atom number for the rest of the program

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

    call createEnergyMesh(emesh, dims%iemxd) !!!!
    params%ielast = dims%iemxd !!!!
    call readEnergyMesh(emesh)  !every process does this!  !!!!

    call OUTTIME(isMasterRank(my_mpi),'input files read.....', &
                                       getElapsedTime(program_timer), 0)

    call createMadelungCalculator(madelung_calc, dims%lmaxd, params%ALAT, params%RMAX, params%GMAX, &
                                  arrays%BRAVAIS, dims%NMAXD, dims%ISHLD)

    call createMadelungLatticeSum(madelung_sum, madelung_calc, dims%naez)

    call calculateMadelungLatticeSum(madelung_sum, I1, arrays%rbasis)

    call OUTTIME(isMasterRank(my_mpi),'Madelung sums calc...',getElapsedTime(program_timer), 0)

    call createGauntCoefficients(gaunts, dims%lmaxd)
    call createShapeGauntCoefficients(shgaunts, dims%lmaxd)

    call createEBalanceHandler(ebalance_handler, params%ielast)
    call initEBalanceHandler(ebalance_handler, my_mpi)
    call setEqualDistribution(ebalance_handler, (emesh%NPNT1 == 0))

    call initLcutoff(arrays%rbasis, arrays%bravais, arrays%lmmaxd, I1)
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

      if (isMasterRank(my_mpi)) then
        call printDoubleLineSep(unit_number = 2)
        call OUTTIME(isMasterRank(my_mpi),'started at ..........', &
                     getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      WRITELOG(2, *) "Iteration Atom ", ITER, I1

      ! New: instead of reading potential every time, communicate it
      ! between energy and spin processes of same atom
      call communicatePotential(my_mpi, atomdata%potential%VISP, &
                                atomdata%potential%VINS, atomdata%core%ECORE)

      ! Core relaxation - only mastergroup needs results
      if (isInMasterGroup(my_mpi)) then
        call RHOCORE_wrapper(emesh%E1, params%NSRA, atomdata)
      endif

! LDA+U
      if (params%LDAU) then

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

        call LDAUINIT(I1,ITER,params%NSRA,ldau_data%NLDAU,ldau_data%LLDAU, &
                      ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,ldau_data%NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      ldau_data%lmaxd, mesh%irmd, mesh%ipand)

      endif
! LDA+U

      call OUTTIME(isMasterRank(my_mpi),'initialized .........', &
                   getElapsedTime(program_timer),ITER)

      ! Scattering calculations - that is what KKR is all about
      ! output: ebalance_handler, kkr (!), jij_data, ldau_data
      call energyLoop(iter, atomdata, emesh, params, dims, gaunts, &
                      ebalance_handler, my_mpi, arrays, kkr, jij_data, ldau_data)

      call OUTTIME(isMasterRank(my_mpi),'G obtained ..........', &
                   getElapsedTime(program_timer),ITER)

!----------------------------------------------------------------------
! BEGIN only processes in master-group are working
!----------------------------------------------------------------------
      if (isInMasterGroup(my_mpi)) then
        ! output: atomdata, densities, broyden, ldau_data, emesh (only correct for master)
        call processKKRresults(iter, kkr, my_mpi, atomdata, emesh, dims, &
                               params, arrays, gaunts, shgaunts, &
                               madelung_sum, program_timer, &
                               densities, broyden, ldau_data)
      endif
!----------------------------------------------------------------------
! END only processes in master-group are working
!----------------------------------------------------------------------

      call broadcastEnergyMesh_com(my_mpi, emesh)

      !call MPI_ALLREDUCE(kkr%NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      !getMyActiveCommunicator(my_mpi),IERR)
      ! TODO: allreduce not necessary, only master rank needs NOITER_ALL, use reduce instead

      ! TODO
      if(isMasterRank(my_mpi)) then
        !call printSolverIterationNumber(ITER, NOITER_ALL)
        call writeIterationTimings(ITER, getElapsedTime(program_timer), &
                                         getElapsedTime(iteration_timer))
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

    call destroyMadelungLatticeSum(madelung_sum)
    call destroyMadelungCalculator(madelung_calc)
    call destroyGauntCoefficients(gaunts)
    call destroyShapeGauntCoefficients(shgaunts)

    call destroyBasisAtom(atomdata)
    call destroyCellData(cell)
    call destroyRadialMeshData(mesh)

    call destroyBroydenData(broyden)
    call destroyLDAUData(ldau_data)
    call destroyJijData(jij_data)
    call destroyDensityResults(densities)
    call destroyKKRresults(kkr)

    !--------------
    call destroyEBalanceHandler(ebalance_handler)
    call destroyEnergyMesh(emesh)

    call destroyCalculationData(calc_data)

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
