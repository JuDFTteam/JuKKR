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
  use KKRSelfConsistency_mod

  use main2_aux_mod
  use muffin_tin_zero_mod
  use EnergyMesh_mod

  use MadelungCalculator_mod
  use lloyd0_new_mod

  use GauntCoefficients_mod
  use ShapeGauntCoefficients_mod

  use TimerMpi_mod

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
  double complex:: Delta_E_z

  double precision::DENEF
  double precision::E1
  double precision::E2
  double precision::TK
  double precision::EFERMI
  double precision::EFOLD
  double precision::CHRGNT
  double precision::E2SHIFT
  double precision::RCUTJIJ
  double precision::DF
  double precision::ALAT
  double precision::FCM
  double precision::MIXING
  double precision::RMAX
  double precision::GMAX
  double precision::PI
  double precision::RMSAVM      ! rms error magnetisation dens. (contribution of single site)
  double precision::RMSAVQ      ! rms error charge density (contribution of single site)
  double precision::EREFLDAU    ! LDA+U

  type (TimerMpi) :: program_timer
  type (TimerMpi) :: iteration_timer
  type (TimerMpi) :: mult_scattering_timer
  type (TimerMpi) :: single_site_timer

  real::TIME_E
  real::TIME_EX

  integer::ICELL
  integer::NPNT1
  integer::NPNT2
  integer::NPNT3
  integer::NPOL
  integer::ITER
  integer::SCFSTEPS
  integer::IMIX
  integer::NOITER
  integer::NOITER_ALL
  integer::ISHIFT
  integer::KPRE
  integer::KTE
  integer::KVMAD
  integer::KXC
  integer::KFORCE
  integer::IPOT
  integer::ISPIN
  integer::I1
  integer::I1BRYD
  integer::LM
  integer::NR
  integer::EKM
  logical::XCCPL
  logical::JIJ
  logical::LDORHOEF
  logical::LDAU
  logical::ERESJIJ


  ! static arrays
  double precision::BRAVAIS(3,3)
  double precision::VBC(2)
  integer::ISYMINDEX(NSYMAXD)
  double precision::ECORE(20,2)

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
  double precision::EULDAU
  double precision::EDCLDAU

  double precision::QC                ! core charge
  double precision::VMAD

  integer::LCOREMAX

  integer::ICST
  integer::NSRA
  integer::NCLS
  integer::NREF
  integer::RF
  integer::NLDAU
  integer::NXIJ

  integer::IE
  integer::IELAST

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

  type(KKRnanoParallel) :: my_mpi

  integer :: flag

 !============================================================= CONSTANTS
  PI = 4.0D0*ATAN(1.0D0)
!=============================================================

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
                            LRECRES2, MMAXD, NAEZ, NGUESSD, NPOTD, NSPIND, NTIRD)


!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call allocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  !every process does this!
  call readKKR0Input      (NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, ICST, &
                           IFUNM, IGUESS, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMSP, MAXMESH, &
                           MIXING, NACLS, NCLS, NCORE, NFU, NR, NREF, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VREF, ZAT)

  ! ---------------------------------------------------------- k_mesh
  call readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)  !every process does this!

  call readEnergyMesh(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ) !every process does this!

  if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

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

    call OUTTIME(isMasterRank(my_mpi),'input files read.....', &
                                       getElapsedTime(program_timer), 0)

! initialise the arrays for (gen. Anderson/Broyden) potential mixing
    UI2 = 0.00
    VI2 = 0.00
    SM1S = 0.00
    FM1S = 0.00

    call createMadelungCalculator(madelung_calc, lmaxd, ALAT, RMAX, GMAX, &
                                  BRAVAIS, NMAXD, ISHLD)

!+++++++++++ atom - parallel  TODO: replace with better construct
    do I1 = 1,NAEZ
      if(getMyAtomRank(my_mpi)==MAPBLOCK(I1,1,NAEZ,1,0,getNumAtomRanks(my_mpi)-1)) then

        call calculateMadelungLatticeSum(madelung_calc, naez, I1, rbasis, smat)

        call OUTTIME(isMasterRank(my_mpi),'Madelung sums calc...',getElapsedTime(program_timer), 0)

        if (isInMasterGroup(my_mpi)) then
          call openPotentialFile(LMPOTD, IRNSD, IRMD)
          call readPotential(I1, VISP, VINS, ECORE)
          call closePotentialFile()
        end if

        call initLcutoff(rbasis, bravais, lmmaxd, I1) !TODO: remove
        WRITELOG(3, *) "lm-array: ", lmarray

      end if
    end do
!+++++++++++ end atom - parallel ++++++++++++++++++++++++++++++++

    call createGauntCoefficients(gaunts, lmaxd)
    call createShapeGauntCoefficients(shgaunts, lmaxd)

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
        call OUTTIME(isMasterRank(my_mpi),'started at ..........',getElapsedTime(program_timer),ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      CMOM   = 0.0D0
      CMINST = 0.0D0
      CHRGNT = 0.0D0

      ! needed for results.f - find better solution - unnecessary I/O
      ! actually only LMPIC==1 process needed to open these files!!!
      if (KTE >= 0) call openResults1File(IEMXD, LMAXD, NPOL)

!N ====================================================================
!     BEGIN do loop over atoms (NMPID-parallel)
!N ====================================================================

      do I1 = 1,NAEZ
        if(getMyAtomRank(my_mpi)==MAPBLOCK(I1,1,NAEZ,1,0,getNumAtomRanks(my_mpi)-1)) then

        WRITELOG(2, *) "Iteration Atom ", ITER, I1
!=======================================================================
! xccpl

          XCCPL = .false.
          ERESJIJ = .false.

          ! calculate exchange couplings only at last self-consistency step and when Jij=true
          if ((ITER==SCFSTEPS).and.JIJ) XCCPL = .true.

          if (XCCPL) then

            !inquire(file='ERESJIJ',exist=ERESJIJ)  ! deactivated, doesn't work anyway

            call CLSJIJ(I1,NAEZ,RR,NR,RBASIS,RCUTJIJ,NSYMAT,ISYMINDEX, &
                        IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ, &
                        nrd, nxijd)

            JXCIJINT = CZERO
            GMATXIJ = CZERO

          endif

          ! New: instead of reading potential every time, communicate it
          !call readPotential(I1, VISP, VINS, ECORE)
          call communicatePotential(my_mpi, VISP, VINS, ECORE)

! LDA+U
          if (LDAU) then

            EREFLDAU = EFERMI
            EREFLDAU = 0.48    ! ???

            call LDAUINIT(I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU, &
                          VISP,NSPIND,R(1,I1),DRDI(1,I1), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          PHILDAU,UMLDAU,WMLDAU, &
                          lmaxd, irmd, ipand)

          endif
! LDA+U

! TIME
          call OUTTIME(isMasterRank(my_mpi),'initialized .........',getElapsedTime(program_timer),ITER)
          call resetTimer(single_site_timer)
! TIME

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================

          do IE = 1,IELAST

            call CPU_TIME(TIME_E)

            ETIME(IE) = 0.0D0

            if (ITER==1.and.IE==1) then

              call EBALANCE1(IELAST, EPROC, EPROCO, empid, iemxd)

            endif

! IE ====================================================================
            if (getMyEnergyId(my_mpi)==EPROC(IE)) then
! IE ====================================================================

              WRITELOG(2, *) "Working on energy point ", IE

              do RF = 1,NREF
                call TREF(EZ(IE),VREF(RF),LMAXD,RMTREF(RF), &
                          TREFLL(1,1,RF),DTREFLL(1,1,RF), LLY)
              end do

              TESTARRAYLOG(4, TREFLL)
              TESTARRAYLOG(4, DTREFLL)

              call GREF_com(EZ(IE),ALAT,gaunts%IEND,NCLS,NAEZ, &
                            gaunts%CLEB,RCLS,ATOM,CLS,gaunts%ICLEB,gaunts%LOFLM,NACLS, &
                            REFPOT, &
                            TREFLL,DTREFLL,GREFN,DGREFN, &
                            LLY_G0TR(:,IE), &
                            getMyAtomRank(my_mpi),getMySEcommunicator(my_mpi),getNumAtomRanks(my_mpi), &
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
spinloop:     do ISPIN = 1,NSPIND
                if(isWorkingSpinRank(my_mpi, ispin)) then

                  if (SMPID==1) then
                    PRSPIN   = ISPIN
                  else
                    PRSPIN   = 1
                  endif

                  call CALCTMAT(LDAU,NLDAU,ICST, &
                                NSRA,EZ(IE), &
                                DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                IRCUT(0,I1),gaunts%CLEB,gaunts%LOFLM,gaunts%ICLEB,gaunts%IEND, &
                                TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAXD, &
                                LLDAU,WMLDAU(1,1,1,ISPIN), &
                                gaunts%ncleb, ipand, irmd, irnsd)

                  DTIXIJ(:,:,ISPIN) = TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

                  if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula

                    call calcdtmat_DeltaEz(delta_E_z, IE, NPNT1, NPNT2, NPNT3, TK)

                    call CALCDTMAT(LDAU,NLDAU,ICST, &
                                  NSRA,EZ(IE),delta_E_z, &
                                  DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                  VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                  IRCUT(0,I1),gaunts%CLEB,gaunts%LOFLM,gaunts%ICLEB,gaunts%IEND, &
                                  DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAXD, &
                                  LLDAU,WMLDAU(1,1,1,ISPIN), &
                                  gaunts%ncleb, ipand, irmd, irnsd)
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
                    if (KTE >= 0) call printEnergyPoint(EZ(IE), IE, ISPIN, NMESH)
                  end if

                  call stopTimer(single_site_timer)
                  call resumeTimer(mult_scattering_timer)

! <<>> Multiple scattering part

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
                  NXIJ,XCCPL,IXCP,ZKRXIJ, &
                  LLY_GRDT(IE,ISPIN),TR_ALPH(ISPIN), &
                  GMATXIJ(1,1,1,ISPIN), &
                  getMySEcommunicator(my_mpi),getNumAtomRanks(my_mpi), &
                  iemxd, &
                  lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
                  nxijd, nguessd, kpoibz, nrd, ekmd)

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

                 call jijSpinCommunication_com(my_mpi, GMATXIJ, DTIXIJ)

                 ! calculate DTIXIJ = T_down - T_up
                 call calcDeltaTupTdown(DTIXIJ)

                 JSCAL = WEZ(IE)/DBLE(NSPIND)

                 call jijLocalEnergyIntegration(my_mpi, JSCAL, GMATXIJ, &
                                                DTIXIJ(:,:,1), RXIJ, NXIJ, IXCP, &
                                                RXCCLS, JXCIJINT)

              end if

! xccpl
! End of Jij calculation
! =====================================================================

              call CPU_TIME(TIME_EX)
              ETIME(IE) = TIME_EX-TIME_E

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
          TESTARRAYLOG(4, TMATN)
          TESTARRAYLOG(4, TR_ALPH)
          TESTARRAYLOG(4, DTDE)
          TESTARRAYLOG(4, GMATN)
          TESTARRAYLOG(4, LLY_GRDT)

!=======================================================================
!communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
          call collectMSResults_com(my_mpi, GMATN, LLY_GRDT, EPROC)
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

            call jijReduceIntResults_com(my_mpi, JXCIJINT)

            if (isInMasterGroup(my_mpi)) then
              call writeJiJs(I1,RXIJ,NXIJ,IXCP,RXCCLS,JXCIJINT, nxijd)
            end if
          endif

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors
!=======================================================================
          if (empid > 1) call EBALANCE2(IELAST,NPNT1, &
                         getMyWorldRank(my_mpi),getMyActiveCommunicator(my_mpi), &
                         ETIME,EPROC,EPROCO, &
                         empid, iemxd) ! should be communicated - dependence on floating point ops!!

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

                !call EPRDIST
                call redistributeInitialGuess_com(my_mpi, PRSC(:,:,PRSPIN), EPROC, EPROCO, KMESH, NofKs)

              endif
            enddo

          endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================

          TESTARRAYLOG(3, GMATN)
          TESTARRAYLOG(3, LLY_G0TR)
          TESTARRAYLOG(3, LLY_GRDT)

!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
          if (isInMasterGroup(my_mpi)) then

            if (LLY==1) then
              ! get WEZRN and RNORM, the important input from previous
              ! calculations is LLY_GRDT_ALL

              ICELL = NTCELL(I1)
              call LLOYD0_NEW(EZ,WEZ,gaunts%CLEB,DRDI(:,I1),R(:,I1),IRMIN(I1), &
                              VINS,VISP,THETAS(:,:,ICELL),ZAT(I1),gaunts%ICLEB, &
                              IFUNM(:,ICELL),IPAN(I1),IRCUT(:,I1),LMSP(:,ICELL), &
                              gaunts%JEND,gaunts%LOFLM,ICST,IELAST,gaunts%IEND,NSPIND,NSRA, &
                              WEZRN,RNORM, &
                              GMATN, &
                              LLY_GRDT, &
                              LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,DMATLDAU, &
                              getMySEcommunicator(my_mpi), &
                              lmaxd, irmd, irnsd, iemxd, &
                              irid, nfund, ipand, gaunts%ncleb)

              TESTARRAYLOG(3, WEZRN)
              TESTARRAYLOG(3, RNORM)

! IME
              call OUTTIME(isMasterRank(my_mpi),'Lloyd processed......',getElapsedTime(program_timer),ITER)
! IME
            else ! no Lloyd

              do IE=1,IELAST
                WEZRN(IE,1) = WEZ(IE)
                WEZRN(IE,2) = WEZ(IE)
              enddo
            endif

            ! now WEZRN stores the weights for E-integration

            DEN = CZERO
            DENEF = 0.0D0

            if (LDAU) then
              DMATLDAU = CZERO
            endif

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN ==================================================================

            do ISPIN = 1,NSPIND
              ICELL = NTCELL(I1)
              IPOT = (I1-1) * NSPIND + ISPIN

              LDORHOEF = NPOL/=0  ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

              ! contains a loop over energies, TODO: remove spin dependence
              ! output: RHO2NS, R2NEF, DEN, ESPV
              call RHOVAL(LDORHOEF,ICST,IELAST, &
                          NSRA,ISPIN,NSPIND,EZ,WEZRN(1,ISPIN), &   ! unfortunately spin-dependent
                          DRDI(1,I1),R(1,I1),IRMIN(I1), &
                          VINS(IRMIND,1,ISPIN),VISP(1,ISPIN), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          THETAS(1,1,ICELL),IFUNM(1,ICELL),LMSP(1,ICELL), &
                          RHO2NS,R2NEF, &
                          DEN(0,1,ISPIN),ESPV(0,ISPIN), &
                          gaunts%CLEB,gaunts%LOFLM,gaunts%ICLEB,gaunts%IEND,gaunts%JEND, &
                          GMATN, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          iemxd, &
                          lmaxd, irmd, irnsd, irid, ipand, nfund, gaunts%ncleb)

              ! output: ECORE, NCORE, LCORE, RHOCAT?, QC
              call RHOCORE(E1,NSRA,ISPIN,NSPIND,I1, &  ! I1 is used only for debugging output
                           DRDI(1,I1),R(1,I1),VISP(1,ISPIN), &
                           A(I1),B(I1),ZAT(I1), &
                           IRCUT(0,I1),RHOCAT,QC, &
                           ECORE(1,ISPIN),NCORE(IPOT),LCORE(1,IPOT), &
                           irmd, ipand)

            end do

! SPIN ==================================================================
!      END do loop over spins
! SPIN ===================================================================

            if (LLY == 1) then
              call renormalizeDOS(DEN,RNORM,LMAXD1,IELAST,NSPIND,IEMXD)
            end if

            ! calculate DOS at Fermi level
            DENEF = calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIND)

            ! ---> l/m_s/atom-resolved charges, output -> CHARGE
            ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
            ! CHARGE -> written to result file
            call calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIND, WEZ, IEMXD)

! LDAU

            EULDAU = 0.0D0
            EDCLDAU = 0.0D0

            if (LDAU.and.NLDAU>=1) then
              call LDAUWMAT(I1,NSPIND,ITER,MIXING,DMATLDAU,NLDAU,LLDAU, &
                            ULDAU,JLDAU,UMLDAU,WMLDAU,EULDAU,EDCLDAU, &
                            lmaxd)
            endif
! LDAU

! ----------------------------------------------------------------------
! -->   determine total charge density expanded in spherical harmonics
! -------------------------------------------------------------- density
            ! output: CATOM
            call RHOTOTB_NEW(NSPIND,RHO2NS,RHOCAT, &
                         DRDI(:,I1),IRCUT(:,I1), &
                         LPOT,NFU(ICELL),LLMSP(1,ICELL),THETAS(:,:,ICELL),IPAN(I1), &
                         CATOM, &
                         irmd, irid, ipand, nfund)

            CHRGNT = CHRGNT + CATOM(1) - ZAT(I1)

            ! write to 'results1' - only to be read in in results.f
            ! necessary for density of states calculation, otherwise
            ! only for informative reasons
            if (KTE >= 0) call writeResults1File(CATOM, CHARGE, DEN, ECORE, I1, NPOL, QC)

          endif
!----------------------------------------------------------------------
! END L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
        end if
      end do
!N ====================================================================
!     END do loop over atoms (NMPID-parallel)
!N ====================================================================

      !call closePotentialFile()
      if (KTE >= 0) call closeResults1File()

      call OUTTIME(isMasterRank(my_mpi),'density calculated ..',getElapsedTime(program_timer),ITER)


! TODO: Only necessary for non-DOS calculation - otherwise proceed to RESULTS
!       Drawback: RESULTS has to be modified - extract DOS part
!----------------------------------------------------------------------
! BEGIN L-MPI: only processes in Master-Group are working
!----------------------------------------------------------------------
      if (isInMasterGroup(my_mpi)) then

        call sumNeutralityDOSFermi_com(CHRGNT, DENEF, getMySEcommunicator(my_mpi))

! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT) !FIXME: hardcoded maximal shift of 0.03
        EFOLD = E2

        if (ISHIFT < 2) E2 = E2 - E2SHIFT

        if( getMyAtomRank(my_mpi) == 0 ) then
          call printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
        end if

        if (KTE >= 0) call openResults2File(LRECRES2)

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIND)
! ----------------------------------------------------------------------

        VAV0 = 0.0D0
        VOL0 = 0.0D0

! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================
        do I1 = 1,NAEZ
          if(getMyAtomRank(my_mpi) == MAPBLOCK(I1,1,NAEZ,1,0,getNumAtomRanks(my_mpi)-1)) then

            ICELL = NTCELL(I1)

            do ISPIN = 1,NSPIND

! -->     get correct density and valence band energies

              ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
              EFOLD*CHRGNT/DBLE(NSPIND*NAEZ)

              do LM = 1,LMPOTD
                call DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1, &
                RHO2NS(1,LM,ISPIN),1)
              end do
            end do
! ----------------------------------------------------------------------
            !output: CMOM, CMINST
            call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),shgaunts%ILM,IFUNM(1,ICELL),shgaunts%IMAXSH,shgaunts%GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ipand, shgaunts%ngshd)

            call OUTTIME(isMasterRank(my_mpi),'RHOMOM ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
            !output: VONS
            call VINTRAS_NEW(LPOT,NSPIND,RHO2NS,VONS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),shgaunts%ILM,IFUNM(1,ICELL),shgaunts%IMAXSH,shgaunts%GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, shgaunts%ngshd, ipand)

            TESTARRAYLOG(3, VONS)
            TESTARRAYLOG(3, RHO2NS)

            call OUTTIME(isMasterRank(my_mpi),'VINTRAS ......',getElapsedTime(program_timer),ITER)

            ! output: VONS (changed), VMAD
            call addMadelungPotential_com(madelung_calc, CMOM, CMINST, NSPIND, &
                 NAEZ, VONS, ZAT, R(:,I1), IRCUT(:,I1), IPAN(I1), VMAD, &
                 SMAT, getMyAtomRank(my_mpi), getMySEcommunicator(my_mpi), getNumAtomRanks(my_mpi), irmd, ipand)

            call OUTTIME(isMasterRank(my_mpi),'VMADELBLK ......',getElapsedTime(program_timer),ITER)

! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

            if (KFORCE==1 .and. ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCEH(CMOM,FLM,LPOT,I1,RHO2NS,VONS, &
              R,DRDI,IMT,ZAT,irmd)  ! TODO: get rid of atom parameter I1
              call FORCE(FLM,FLMC,LPOT,NSPIND,I1,RHOCAT,VONS,R, &
              DRDI,IMT,naez,irmd)
! ---------------------------------------------------------------------
            end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

            if (KTE==1) then
              ! calculate total energy and individual contributions if requested
              IPOT = NSPIND* (I1-1) + 1
              call ESPCB_NEW(ESPC,NSPIND,ECORE,LCORE(:,IPOT),LCOREMAX,NCORE(IPOT))

              ! output: EPOTIN
              call EPOTINB_NEW(EPOTIN,NSPIND,RHO2NS,VISP,R(:,I1),DRDI(:,I1), &
              IRMIN(I1),IRWS(I1),LPOT,VINS,IRCUT(:,I1),IPAN(I1),ZAT(I1), &
              irmd, irnsd, ipand)

              ! output: ECOU - l resolved Coulomb energy
              call ECOUB_NEW(CMOM,ECOU,LPOT,NSPIND,RHO2NS, &
              VONS,ZAT(I1),R(:,I1), &
              DRDI(:,I1),KVMAD,IRCUT(:,I1),IPAN(I1),shgaunts%IMAXSH,IFUNM(1,ICELL), &
              shgaunts%ILM,shgaunts%GSH,THETAS(:,:,ICELL),LMSP(1,ICELL), &
              irmd, irid, nfund, ipand, shgaunts%ngshd)

            end if

            call OUTTIME(isMasterRank(my_mpi),'KTE ......',getElapsedTime(program_timer),ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
            ! output: VONS (changed), EXC (exchange energy)
            call VXCDRV_NEW(EXC,KTE,KXC,LPOT,NSPIND,RHO2NS, &
            VONS,R(:,I1),DRDI(:,I1),A(I1), &
            IRWS(I1),IRCUT(:,I1),IPAN(I1),shgaunts%GSH,shgaunts%ILM,shgaunts%IMAXSH,IFUNM(1,ICELL), &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, shgaunts%ngshd, ipand)

            call OUTTIME(isMasterRank(my_mpi),'VXCDRV ......',getElapsedTime(program_timer),ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCXC_com(FLM,FLMC,LPOT,NSPIND,I1,RHOCAT,VONS,R, &
              ALAT,DRDI,IMT,ZAT, &
              getMyAtomRank(my_mpi), &
              getMySEcommunicator(my_mpi), &
              naez, irmd)
! ---------------------------------------------------------------------
            end if

            ! unnecessary I/O? see results.f
            if (KTE >= 0) call writeResults2File(CATOM, ECOU, EDCLDAU, EPOTIN, ESPC, ESPV, EULDAU, EXC, I1, LCOREMAX, VMAD)

! Force calculation ends
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! =====================================================================

            !output: VAV0, VOL0
            call MTZERO_NEW(LMPOTD,NSPIND,VONS,ZAT(I1),R(:,I1),DRDI(:,I1),IMT(I1),IRCUT(:,I1), &
                            IPAN(I1),LMSP(1,ICELL),IFUNM(1,ICELL), &
                            THETAS(:,:,ICELL),IRWS(I1),VAV0,VOL0, &
                            irmd, irid, nfund, ipand)

            call OUTTIME(isMasterRank(my_mpi),'MTZERO ......',getElapsedTime(program_timer),ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
          end if
        end do
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================

        call OUTTIME(isMasterRank(my_mpi),'calculated pot ......',getElapsedTime(program_timer),ITER)

        call allreduceMuffinTinShift_com(getMySEcommunicator(my_mpi), VAV0, VBC, VOL0)

        call shiftMuffinTinZero(ISHIFT, VBC, E2SHIFT) ! purpose? ISHIFT usually=0

        if(isMasterRank(my_mpi)) then
          call printMuffinTinShift(VAV0, VBC, VOL0)
        end if

! =====================================================================
! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

!+++++++++++++++++ BEGIN ATOM PARALLEL +++++++++++++++++++++++++++++++
        do I1 = 1,NAEZ
          if(getMyAtomRank(my_mpi)== MAPBLOCK(I1,1,NAEZ,1,0,getNumAtomRanks(my_mpi)-1)) then

            ICELL = NTCELL(I1)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
            do ISPIN = 1,NSPIND
              IPOT = NSPIND* (I1-1) + ISPIN

              call shiftPotential(VONS(:,:,ISPIN), IRCUT(IPAN(I1),I1), VBC(ISPIN))
              !output: VONS (changed)
              call CONVOL_NEW(IRCUT(1,I1),IRC(I1), &
                          shgaunts%IMAXSH(shgaunts%LMPOTD),shgaunts%ILM,IFUNM(1,ICELL),LMPOTD,shgaunts%GSH, &
                          THETAS(:,:,ICELL),ZAT(I1), &
                          R(1,I1),VONS(1,1,ISPIN),LMSP(1,ICELL), &
                          irid, nfund, irmd, shgaunts%ngshd)

            end do

! -->   calculation of RMS and final construction of the potentials (straight mixing)

            call MIXSTR_NEW(RMSAVQ,RMSAVM,LMPOTD,NSPIND,MIXING,FCM, &
                            IRC(I1),IRMIN(I1),R(:,I1),DRDI(:,I1),VONS,VISP,VINS, &
                            irmd, irnsd)

            I1BRYD=I1
          end if
        end do
!+++++++++++++++++ END ATOM PARALLEL +++++++++++++++++++++++++++++++

       ! it is weird that straight mixing is called in any case before
! -->  potential mixing procedures: Broyden or Andersen updating schemes
        if (IMIX>=3) then
          call BRYDBM_com(VISP,VONS,VINS, &
          LMPOTD,R,DRDI,MIXING, &
          IRC,IRMIN,NSPIND,I1BRYD, &
          IMIX,ITER, &
          UI2,VI2,WIT,SM1S,FM1S, &
          getMyAtomRank(my_mpi), &
          getMySEcommunicator(my_mpi), &
          itdbryd, irmd, irnsd, nspind)
        endif

        TESTARRAYLOG(4, VINS)
        TESTARRAYLOG(4, VISP)
        TESTARRAYLOG(4, VONS)

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!--------- BEGIN Atom-parallel ----------------------------------------
        do I1 = 1,NAEZ
          if(getMyAtomRank(my_mpi) == MAPBLOCK(I1,1,NAEZ,1,0,getNumAtomRanks(my_mpi)-1)) then

            call resetPotentials(IRC(I1), IRMD, IRMIN(I1), IRMIND, LMPOTD, &
                                 NSPIND, VINS, VISP, VONS) ! Note: only LMPIC=1 processes

! ----------------------------------------------------- output_potential
            call openPotentialFile(LMPOTD, IRNSD, IRMD)
            call writePotential(I1, VISP, VINS, ECORE)
            call closePotentialFile()
! ----------------------------------------------------- output_potential

          end if
        end do
! -------- END Atom-parallel ------------------------------------------

        if (KTE >= 0) call closeResults2File()

! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
! ====== write RMS convergency data - written by
! MYRANK=0   (RMSOUT) =================================================
! write formatted potential if file VFORM exists
! Note: LMPIC is always 1 here!

        call RMSOUT_com(RMSAVQ,RMSAVM,ITER,E2, &
        SCFSTEPS,VBC,NSPIND,NAEZ, &
        KXC,LPOT,A,B,IRC, &
        VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
        ECORE,LCORE,NCORE,ZAT,ITITLE, &
        getMyAtomRank(my_mpi), &
        getMySEcommunicator(my_mpi),getNumAtomRanks(my_mpi), &
        irmd, irnsd)

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
        ! also other stuff is read from results1
        call RESULTS(LRECRES2,IELAST,ITER,LMAXD,NAEZ,NPOL, &
        NSPIND,KPRE,KTE,LPOT,E1,E2,TK,EFERMI, &
        ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,LDAU, &
        iemxd)

        ! only MASTERRANK updates, other ranks get it broadcasted later
        ! (although other processes could update themselves)
        call updateEnergyMesh(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3)

        ! write file 'energy_mesh'
        if (NPOL /= 0) EFERMI = E2  ! if not a DOS-calculation E2 coincides with Fermi-Energy

        call writeEnergyMesh(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)

        call printDoubleLineSep()

      endif
! -----------------------------------------------------------------
! END: only MASTERRANK is working here
! -----------------------------------------------------------------

      call broadcastEnergyMesh_com(getMyActiveCommunicator(my_mpi), 0, E1, E2, EZ, IEMXD, WEZ) ! BCRANK = 0

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      getMyActiveCommunicator(my_mpi),IERR) ! TODO: allreduce not necessary, only master rank needs NOITER_ALL, use reduce instead

      if(isMasterRank(my_mpi)) then

        call printSolverIterationNumber(ITER, NOITER_ALL)
        call writeIterationTimings(ITER, getElapsedTime(program_timer), getElapsedTime(iteration_timer))

      endif

! manual exit possible by creation of file 'STOP' in home directory
      if (isManualAbort_com(getMyWorldRank(my_mpi), getMyActiveCommunicator(my_mpi)) .eqv. .true.) exit

! ######################################################################
! ######################################################################
    enddo          ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################

    if (isMasterRank(my_mpi)) close(2)    ! TIME

    if (KFORCE==1) close(54)

    call destroyMadelungCalculator(madelung_calc)
    call destroyGauntCoefficients(gaunts)
    call destroyShapeGauntCoefficients(shgaunts)

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
