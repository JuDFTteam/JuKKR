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

  use main2_arrays_mod ! If you can't find one of the unbelievable amount of arrays, look here
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

  implicit none
  include 'mpif.h'

  type (MadelungCalculator) :: madelung_calc
  type (ShapeGauntCoefficients) :: shgaunts
  type (GauntCoefficients) :: gaunts

  !     .. Parameters ..

  integer, parameter ::   NSYMAXD = 48
  double complex, parameter :: CZERO = (0.0D0,0.0D0)

  !     ..
  !     .. Local Scalars ..
  double complex:: JSCAL        ! scaling factor for Jij calculation

  double precision::DENEF
  double precision::CHRGNT
  double precision::RCUTJIJ
  double precision::ALAT
  double precision::FCM
  double precision::MIXING
  double precision::RMAX
  double precision::GMAX

  double precision::RMSAVM      ! rms error magnetisation dens. (contribution of single site)
  double precision::RMSAVQ      ! rms error charge density (contribution of single site)

  type (TimerMpi) :: program_timer
  type (TimerMpi) :: iteration_timer
  type (TimerMpi) :: mult_scattering_timer
  type (TimerMpi) :: single_site_timer

  type (EBalanceHandler) :: ebalance_handler

  integer::ITER
  integer::SCFSTEPS
  integer::IMIX
  integer::NOITER
  integer::NOITER_ALL
  integer::KPRE
  integer::KTE
  integer::KXC
  integer::KFORCE
  integer::ISPIN
  integer::I1
  integer::NR
  integer::EKM
  logical::XCCPL
  logical::JIJ
  logical::LDORHOEF
  logical::LDAU

  ! static arrays
  double precision::BRAVAIS(3,3)
  double precision::VBC(2)
  integer::ISYMINDEX(NSYMAXD)

  !     .. Local Arrays ..

  double precision::VAV0
  double precision::VOL0

  integer::NMESH
  integer::NSYMAT
  integer::MAXMESH

  double precision:: QMRBOUND
  integer::IGUESS
  integer::BCP
  integer::PRSPIN

  double precision::EPOTIN
  double precision::VMAD

  integer::LCOREMAX

  integer::ICST
  integer::NSRA
  integer::NCLS
  integer::NREF
  integer::RF
  !integer::NXIJ  !DEL

  integer::IE
  integer::IELAST

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

  type(KKRnanoParallel) :: my_mpi

  integer :: flag
  logical, external :: testVFORM

  type (RadialMeshData), target :: mesh
  type (CellData), target :: cell
  type (BasisAtom), target :: atomdata
  type (EnergyMesh) :: emesh
  type (LDAUData) :: ldau_data
  type (JijData) :: jij_data
  type (BroydenData) :: broyden

  call read_dimension_parameters()

  ! consistency check of some dimension parameters
  call consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)

! -----------------------------------------------------------------------------
  call createKKRnanoParallel(my_mpi, NAEZ, SMPID, EMPID)
  call setKKRnanoNumThreads(nthrds)
  call printKKRnanoInfo(my_mpi, nthrds)
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

  ! from dimension parameters - calculate some derived parameters
  call getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, LMAXD, &
                            LMAXD1, LMMAXD, LMPOTD, LMXSPD, &
                            LRECRES2, MMAXD, NAEZ, NGUESSD, NSPIND, NTIRD)


!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call allocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  !every process does this!
  call readKKR0InputNew(NSYMAXD, ALAT, ATOM, BCP, BRAVAIS, &
                        CLS, DSYMLL, EZOA, FCM, GMAX, ICST, &
                        IGUESS, IMIX, INDN0, &
                        ISYMINDEX, &
                        JIJ, KFORCE, KMESH, KPRE, KTE, KXC, &
                        LDAU, MAXMESH, &
                        MIXING, NACLS, NCLS, NR, NREF, &
                        NSRA, NSYMAT, NUMN0, OPTC, QMRBOUND, &
                        RBASIS, RCLS, RCUTJIJ, REFPOT, RMAX, RMTREF, &
                        RR, SCFSTEPS, TESTC, VREF, ZAT)


  !if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  call consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)

  if ((JIJ .eqv. .true.) .and. (nspind /= 2)) then
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

    ! ---------------------------------------------------------- k_mesh
    call readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)  !every process does this!

    call OUTTIME(isMasterRank(my_mpi),'input files read.....', &
                                       getElapsedTime(program_timer), 0)

    call createBasisAtom(atomdata, I1, lpot, nspind, irmind, irmd)
    call openBasisAtomDAFile(atomdata, 37, "atoms")
    call readBasisAtomDA(atomdata, 37, I1)
    call closeBasisAtomDAFile(37)

    if (isInMasterGroup(my_mpi)) then
      call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
      call readBasisAtomPotentialDA(atomdata, 37, I1)
      call closeBasisAtomPotentialDAFile(37)
    end if

    call createCellData(cell, irid, (2*LPOT+1)**2, nfund)
    call openCellDataDAFile(cell, 37 , "cells")
    call readCellDataDA(cell, 37, getCellIndex(atomdata))
    call closeCellDataDAFile(37)

    call associateBasisAtomCell(atomdata, cell)

    call createRadialMeshData(mesh, irmd, ipand)
    call openRadialMeshDataDAFile(mesh, 37 , "meshes")
    call readRadialMeshDataDA(mesh, 37, I1)
    call closeRadialMeshDataDAFile(37)

    call associateBasisAtomMesh(atomdata, mesh)

    call createLDAUData(ldau_data, ldau, irmd, lmaxd, nspind)
    call createJijData(jij_data, jij, rcutjij, nxijd,lmmaxd,nspind)
    call createBroydenData(broyden, ntird, itdbryd, imix, mixing)

    call createEnergyMesh(emesh, iemxd)
    ielast = iemxd

    call readEnergyMesh(emesh)  !every process does this!

    call createMadelungCalculator(madelung_calc, lmaxd, ALAT, RMAX, GMAX, &
                                  BRAVAIS, NMAXD, ISHLD)

    call calculateMadelungLatticeSum(madelung_calc, naez, I1, rbasis, smat)

    call OUTTIME(isMasterRank(my_mpi),'Madelung sums calc...',getElapsedTime(program_timer), 0)

    call createGauntCoefficients(gaunts, lmaxd)
    call createShapeGauntCoefficients(shgaunts, lmaxd)

    call createEBalanceHandler(ebalance_handler, ielast)
    call initEBalanceHandler(ebalance_handler, my_mpi)
    call setEqualDistribution(ebalance_handler, (emesh%NPNT1 == 0))

    call initLcutoff(rbasis, bravais, lmmaxd, I1) !TODO: remove
    WRITELOG(3, *) "lm-array: ", lmarray

!+++++++++++
    ASSERT( ZAT(I1) == atomdata%Z_nuclear )
    ASSERT( ielast == iemxd )

   !flag = 0
   !99 continue
   !if (flag == 0) goto 99

! ######################################################################
! ######################################################################
    do ITER = 1, SCFSTEPS
! ######################################################################
! ######################################################################

      call resetTimer(iteration_timer)
      call resetTimer(mult_scattering_timer)
      call stopTimer(mult_scattering_timer)

      EKM    = 0
      NOITER = 0

      if (isMasterRank(my_mpi)) then
        call printDoubleLineSep(unit_number = 2)
        call OUTTIME(isMasterRank(my_mpi),'started at ..........', &
                     getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      CMOM   = 0.0D0
      CMINST = 0.0D0
      CHRGNT = 0.0D0

      WRITELOG(2, *) "Iteration Atom ", ITER, I1

!=======================================================================
! xccpl

      XCCPL = .false.

      ! calculate exchange couplings only at last self-consistency step and when Jij=true
      if ((ITER==SCFSTEPS).and.JIJ) XCCPL = .true.

      if (XCCPL) then

        call CLSJIJ(I1,NAEZ,RR,NR,RBASIS,jij_data%RCUTJIJ,NSYMAT,ISYMINDEX, &
                    jij_data%IXCP,jij_data%NXCP,jij_data%NXIJ,jij_data%RXIJ,jij_data%RXCCLS,jij_data%ZKRXIJ, &
                    nrd, jij_data%nxijd)

        jij_data%JXCIJINT = CZERO
        jij_data%GMATXIJ = CZERO

      endif

      ! New: instead of reading potential every time, communicate it
      call communicatePotential(my_mpi, atomdata%potential%VISP, atomdata%potential%VINS, atomdata%core%ECORE)

      ! Core relaxation - only mastergroup needs results
      if (isInMasterGroup(my_mpi)) then
        call RHOCORE_wrapper(emesh%E1, NSRA, atomdata)
      endif

! LDA+U
      if (LDAU) then

        ldau_data%EREFLDAU = emesh%EFERMI
        ldau_data%EREFLDAU = 0.48    ! ???

        call LDAUINIT(I1,ITER,NSRA,ldau_data%NLDAU,ldau_data%LLDAU,ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%EREFLDAU, &
                      atomdata%potential%VISP,NSPIND,mesh%R,mesh%DRDI, &
                      atomdata%Z_nuclear,mesh%IPAN,mesh%IRCUT, &
                      ldau_data%PHILDAU,ldau_data%UMLDAU,ldau_data%WMLDAU, &
                      lmaxd, irmd, ipand)

      endif
! LDA+U

! TIME
      call OUTTIME(isMasterRank(my_mpi),'initialized .........', &
                   getElapsedTime(program_timer),ITER)
      call resetTimer(single_site_timer)
! TIME

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================

      do IE = 1,IELAST
! IE ====================================================================
        if (getMyEnergyId(my_mpi)==ebalance_handler%EPROC(IE)) then
! IE ====================================================================
          call startEBalanceTiming(ebalance_handler, IE)

          WRITELOG(2, *) "Working on energy point ", IE

          do RF = 1,NREF
            call TREF(emesh%EZ(IE),VREF(RF),LMAXD,RMTREF(RF), &
                      TREFLL(1,1,RF),DTREFLL(1,1,RF), LLY)
          end do

          TESTARRAYLOG(3, TREFLL)
          TESTARRAYLOG(3, DTREFLL)

          call GREF_com(emesh%EZ(IE),ALAT,gaunts%IEND,NCLS,NAEZ, &
                        gaunts%CLEB,RCLS,ATOM,CLS,gaunts%ICLEB, &
                        gaunts%LOFLM,NACLS, &
                        REFPOT, &
                        TREFLL,DTREFLL,GREFN,DGREFN, &
                        LLY_G0TR(:,IE), &
                        getMyAtomRank(my_mpi),getMySEcommunicator(my_mpi),&
                        getNumAtomRanks(my_mpi), &
                        lmaxd, naclsd, gaunts%ncleb, nrefd, nclsd, &
                        LLY)

          TESTARRAYLOG(3, GREFN)
          TESTARRAYLOG(3, DGREFN)

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================
!------------------------------------------------------------------------------
!     beginning of SMPID-parallel section
!------------------------------------------------------------------------------
spinloop: do ISPIN = 1,NSPIND
            if(isWorkingSpinRank(my_mpi, ispin)) then

              if (SMPID==1) then
                PRSPIN   = ISPIN
              else
                PRSPIN   = 1
              endif

              call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, ICST, NSRA, gaunts, TMATN, TR_ALPH, ldau_data)

              jij_data%DTIXIJ(:,:,ISPIN) = TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

              if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
                call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, ICST, NSRA, gaunts, DTDE, TR_ALPH, ldau_data)
              end if

              RF = REFPOT(I1)
              call substractReferenceTmatrix(TMATN(:,:,ISPIN), TREFLL(:,:,RF), LMMAXD)
              ! do the same for derivative of T-matrix
              call substractReferenceTmatrix(DTDE(:,:,ISPIN), DTREFLL(:,:,RF), LMMAXD)

              ! TMATN now contains Delta t = t - t_ref !!!
              ! DTDE now contains Delta dt !!!

              ! renormalize TR_ALPH
              TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(CLS(I1), IE)

              NMESH = KMESH(IE)

              if( getMyAtomRank(my_mpi)==0 ) then
                if (KTE >= 0) call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, NMESH)
              end if

              call stopTimer(single_site_timer)
              call resumeTimer(mult_scattering_timer)

! <<>> Multiple scattering part

!                  if (I1 == 1 .and. IE == IELAST .and. ISPIN == 1) then
!                     DEBUG_dump_matrix = .true.
!                  else
!                     DEBUG_dump_matrix = .false.
!                  endif

              TESTARRAYLOG(3, TMATN(:,:,ISPIN))

              call KLOOPZ1( &
              GMATN(1,1,1,ISPIN), &
              ALAT,IE,ITER,NAEZ, &
              NOFKS(NMESH),VOLBZ(NMESH), &
              BZKP(1,1,NMESH),VOLCUB(1,NMESH), &
              CLS,NACLS,RR, &
              EZOA,ATOM,GREFN,DGREFN, &
              NSYMAT,DSYMLL, &
              TMATN(:,:,ISPIN),DTDE(:,:,ISPIN), &
              NUMN0,INDN0,I1, &
              PRSC(1,1,PRSPIN), &
              EKM,NOITER, &
              QMRBOUND,IGUESS,BCP, &
              jij_data%NXIJ,XCCPL,jij_data%IXCP,jij_data%ZKRXIJ, &
              LLY_GRDT(IE,ISPIN),TR_ALPH(ISPIN), &
              jij_data%GMATXIJ(1,1,1,ISPIN), &
              getMySEcommunicator(my_mpi),getNumAtomRanks(my_mpi), &
              iemxd, &
              lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
              jij_data%nxijd, nguessd, kpoibz, nrd, ekmd)

              TESTARRAYLOG(3, GMATN(:,:,IE,ISPIN))

              call stopTimer(mult_scattering_timer)
              call resumeTimer(single_site_timer)

            endif
          end do spinloop                          ! ISPIN = 1,NSPIN
!------------------------------------------------------------------------------
!        End of SMPID-parallel section
!------------------------------------------------------------------------------
! SPIN ==================================================================
!     END do loop over spins
! SPIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! xccpl

          if (XCCPL) then
            call jijSpinCommunication_com(my_mpi, jij_data%GMATXIJ, jij_data%DTIXIJ)

            ! calculate DTIXIJ = T_down - T_up
            call calcDeltaTupTdown(jij_data%DTIXIJ)

            JSCAL = emesh%WEZ(IE)/DBLE(NSPIND)

            call jijLocalEnergyIntegration(my_mpi, JSCAL, jij_data%GMATXIJ, &
                                            jij_data%DTIXIJ(:,:,1), jij_data%RXIJ, jij_data%NXIJ, jij_data%IXCP, &
                                            jij_data%RXCCLS, jij_data%JXCIJINT)
          end if

! xccpl
! End of Jij calculation
! =====================================================================

          call stopEBalanceTiming(ebalance_handler, ie)

! IE ====================================================================
        endif
! IE ====================================================================

! for initial guess calculate sparse indices combining IE.KPT
        EKM = EKM + NOFKS(KMESH(IE))

      end do                   ! IE = 1,IELAST

! IE ====================================================================
!     END do loop over energies (EMPID-parallel)
! IE ====================================================================

      call stopTimer(single_site_timer)

!=======================================================================
!communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
      call collectMSResults_com(my_mpi, GMATN, LLY_GRDT, ebalance_handler%EPROC)
!=======================================================================

! TIME
      call OUTTIME(isMasterRank(my_mpi),'Single Site took.....',getElapsedTime(single_site_timer),ITER)
      call OUTTIME(isMasterRank(my_mpi),'Mult. Scat. took.....',getElapsedTime(mult_scattering_timer),ITER)
      call OUTTIME(isMasterRank(my_mpi),'G obtained ..........',getElapsedTime(program_timer),ITER)
! TIME

!=======================================================================
!     output of Jij's
!=======================================================================
      if (XCCPL) then

        call jijReduceIntResults_com(my_mpi, jij_data%JXCIJINT)

        if (isInMasterGroup(my_mpi)) then
          call writeJiJs(I1,jij_data%RXIJ,jij_data%NXIJ,jij_data%IXCP,jij_data%RXCCLS,jij_data%JXCIJINT, jij_data%nxijd)
        end if
      endif

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors - all processes SYNCED
!=======================================================================
      call updateEBalance_com(ebalance_handler, my_mpi)

!=======================================================================
!     in case of IGUESS and EMPID > 1 initial guess arrays might
!     have to be adjusted to new distributions
!=======================================================================
      if ((IGUESS==1).and.(EMPID>1)) then

        do ISPIN = 1,NSPIND
          if(isWorkingSpinRank(my_mpi, ispin)) then

            if (SMPID==1) then
              PRSPIN   = ISPIN
            else
              PRSPIN   = 1
            endif

            WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
            WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old
            call redistributeInitialGuess_com(my_mpi, PRSC(:,:,PRSPIN), &
                 ebalance_handler%EPROC, ebalance_handler%EPROC_old, KMESH, NofKs)

          endif
        enddo

      endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================

      TESTARRAYLOG(3, GMATN)
      TESTARRAYLOG(3, LLY_GRDT)

!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
      if (isInMasterGroup(my_mpi)) then

        ! out: emesh, RNORM
        call lloyd0_wrapper_com(atomdata, my_mpi, LLY_GRDT, emesh, RNORM, LLY, ICST, NSRA, GMATN, gaunts, ldau_data)

        if (LLY == 1) then
          TESTARRAYLOG(3, emesh%WEZRN)
          TESTARRAYLOG(3, RNORM)
          call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......',getElapsedTime(program_timer),ITER)
        endif

        ! now WEZRN stores the weights for E-integration

        DEN = CZERO
        DENEF = 0.0D0

        if (LDAU) then
          ldau_data%DMATLDAU = CZERO
        endif

        LDORHOEF = emesh%NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

        ! has to be done after Lloyd
        ! output: RHO2NS, R2NEF, DEN, ESPV
        call RHOVAL_wrapper(atomdata, LdoRhoEF, ICST, NSRA, RHO2NS, R2NEF, &
                            DEN, ESPV, GMATN, gaunts, emesh, ldau_data)

! ----------------------------------------------------------------------
! -->   determine total charge expanded in spherical harmonics
! -------------------------------------------------------------- density
        ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
        call RHOTOTB_wrapper(CATOM, RHO2NS, atomdata)

        CHRGNT = CHRGNT + CATOM(1) - atomdata%Z_nuclear

        if (LLY == 1) then
          call renormalizeDOS(DEN,RNORM,LMAXD1,IELAST,NSPIND,IEMXD)
        end if

        ! calculate DOS at Fermi level
        DENEF = calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIND)

        ! ---> l/m_s/atom-resolved charges, output -> CHARGE
        ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
        ! CHARGE -> written to result file
        call calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIND, emesh%WEZ, IEMXD)

        call sumNeutralityDOSFermi_com(CHRGNT, DENEF, getMySEcommunicator(my_mpi))

        ! write to 'results1' - only to be read in in results.f
        ! necessary for density of states calculation, otherwise
        ! only for informative reasons
        if (KTE >= 0) then
          call openResults1File(IEMXD, LMAXD, emesh%NPOL)
          call writeResults1File(CATOM, CHARGE, DEN, atomdata%core%ECORE, I1, emesh%NPOL, atomdata%core%QC_corecharge)
          call closeResults1File()
        endif

        call OUTTIME(isMasterRank(my_mpi),'density calculated ..',getElapsedTime(program_timer),ITER)

        call doFermiEnergyCorrection(atomdata, isMasterRank(my_mpi), naez, 0.03d0, CHRGNT, DENEF, R2NEF, &
                                     ESPV, RHO2NS, emesh%E2)

        !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)
        call RHOMOM_NEW_wrapper(CMOM,CMINST,RHO2NS(:,:,1), cell, mesh, shgaunts)

        call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
        !output: VONS
        call VINTRAS_wrapper(RHO2NS(:,:,1), shgaunts, atomdata)

        TESTARRAYLOG(3, atomdata%potential%VONS)
        TESTARRAYLOG(3, RHO2NS)

        call OUTTIME(isMasterRank(my_mpi),'VINTRAS ......',getElapsedTime(program_timer),ITER)

        ! output: VONS (changed), VMAD
        call addMadelungPotential_com(madelung_calc, CMOM, CMINST, NSPIND, &
             NAEZ, atomdata%potential%VONS, ZAT, mesh%R, mesh%IRCUT, mesh%IPAN, VMAD, &
             SMAT, getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi), getNumAtomRanks(my_mpi), irmd, ipand)

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

        if (KTE==1) then
          ! calculate total energy and individual contributions if requested
          ! core electron contribution
          call ESPCB_wrapper(ESPC, LCOREMAX, atomdata)
          ! output: EPOTIN
          call EPOTINB_wrapper(EPOTIN,RHO2NS,atomdata)
          ! output: ECOU - l resolved Coulomb energy
          call ECOUB_wrapper(CMOM, ECOU, RHO2NS, shgaunts, atomdata)
        end if

        call OUTTIME(isMasterRank(my_mpi),'KTE ......',getElapsedTime(program_timer),ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
        ! output: VONS (changed), EXC (exchange energy) (l-resolved)
        call VXCDRV_wrapper(EXC,KXC,RHO2NS, shgaunts, atomdata)

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
        if (KTE >= 0) then
          call openResults2File(LRECRES2)
          call writeResults2File(CATOM, ECOU, ldau_data%EDCLDAU, EPOTIN, ESPC, ESPV, ldau_data%EULDAU, EXC, I1, LCOREMAX, VMAD)
          call closeResults2File()
        end if

        ! calculate new muffin-tin zero. output: VAV0, VOL0
        call MTZERO_wrapper(VAV0, VOL0, atomdata)

        call OUTTIME(isMasterRank(my_mpi),'MTZERO ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

        call OUTTIME(isMasterRank(my_mpi),'calculated pot ......',getElapsedTime(program_timer),ITER)

        call allreduceMuffinTinShift_com(getMySEcommunicator(my_mpi), VAV0, VBC, VOL0)

        if(isMasterRank(my_mpi)) then
          call printMuffinTinShift(VAV0, VBC, VOL0)
        end if

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

! -->   shift potential by VBC and multiply with shape functions - output: VONS
        call CONVOL_wrapper(VBC, shgaunts, atomdata)

! LDAU
        ldau_data%EULDAU = 0.0D0
        ldau_data%EDCLDAU = 0.0D0

        if (ldau_data%LDAU.and.ldau_data%NLDAU>=1) then
          call LDAUWMAT(I1,NSPIND,ITER,MIXING,ldau_data%DMATLDAU,ldau_data%NLDAU,ldau_data%LLDAU, &
                        ldau_data%ULDAU,ldau_data%JLDAU,ldau_data%UMLDAU,ldau_data%WMLDAU,ldau_data%EULDAU,ldau_data%EDCLDAU, &
                        lmaxd)
        endif
! LDAU

! -->   calculation of RMS and final construction of the potentials (straight mixing)
        call MIXSTR_wrapper(atomdata, RMSAVQ, RMSAVM, MIXING, FCM)

        ! output of RMS error
        call RMSOUT_com(RMSAVQ,RMSAVM,ITER,NSPIND,NAEZ, &
                       getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi))

        ! it is weird that straight mixing is called in any case before
! -->   potential mixing procedures: Broyden or Andersen updating schemes
        if (IMIX>=3) then
          call BRYDBM_new_com(atomdata%potential%VISP,atomdata%potential%VONS,atomdata%potential%VINS, &
          LMPOTD,mesh%R,mesh%DRDI,broyden%MIXING, &
          mesh%IRC,mesh%IRMIN,NSPIND, &
          broyden%IMIX,ITER, &
          broyden%UI2,broyden%VI2,broyden%WIT,broyden%SM1S,broyden%FM1S, &
          getMyAtomRank(my_mpi), &
          getMySEcommunicator(my_mpi), &
          broyden%itdbryd, irmd, irnsd, nspind)
        endif

        TESTARRAYLOG(3, atomdata%potential%VINS)
        TESTARRAYLOG(3, atomdata%potential%VISP)
        TESTARRAYLOG(3, atomdata%potential%VONS)

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        call resetPotentials(mesh%IRC, IRMD, mesh%IRMIN, IRMIND, LMPOTD, &
                             NSPIND, atomdata%potential%VINS, atomdata%potential%VISP, atomdata%potential%VONS) ! Note: only LMPIC=1 processes

! ----------------------------------------------------- output_potential
        call openBasisAtomPotentialDAFile(atomdata, 37, "vpotnew")
        call writeBasisAtomPotentialDA(atomdata, 37, I1)
        call closeBasisAtomPotentialDAFile(37)
! ----------------------------------------------------- output_potential

! write formatted potential if file VFORM exists - contains bad inquire
! - bad check deactivated when KTE<0
        if (ITER == SCFSTEPS .and. KTE >= 0) then
          if (testVFORM()) then
            call writeFormattedPotential(emesh%E2, ALAT, VBC, KXC, atomdata)
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
        call RESULTS(LRECRES2,IELAST,ITER,LMAXD,NAEZ,emesh%NPOL, &
        NSPIND,KPRE,KTE,LPOT,emesh%E1,emesh%E2,emesh%TK,emesh%EFERMI, &
        ALAT,atomdata%core%ITITLE(:,1:NSPIND),CHRGNT,ZAT,emesh%EZ,emesh%WEZ,LDAU, &
        iemxd)

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

! ======================================================================

  endif ! active Ranks

  CLOSELOG

!------------------------------------------------------------------------------
  call deallocate_main2_arrays()
!------------------------------------------------------------------------------

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================
! Free KKRnano mpi resources

  call destroyKKRnanoParallel(my_mpi)

end program MAIN2
