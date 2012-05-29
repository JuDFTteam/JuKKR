! KKRnano
! massive parallel KKR for nanoscaled systems

program MAIN2

  !use mpi
  use common_testc
  use common_optc
  use common_mpi

  use main2_arrays_mod ! If you can't find one of the unbelievable amount of arrays, look here
  use lloyds_formula_mod
  use KKRSelfConsistency_mod

  use main2_aux_mod
  use muffin_tin_zero_mod
  use EnergyMesh_mod

  implicit none
  include 'mpif.h'

  !     .. Parameters ..

  integer::   NSYMAXD
  parameter (NSYMAXD=48)
  double complex:: CZERO
  parameter      (CZERO=(0.0D0,0.0D0))

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
  double precision::MIX
  double precision::MIXING
  double precision::QBOUND
  double precision::VOLUME0
  double precision::RMAX
  double precision::GMAX
  double precision::FPI
  double precision::PI
  double precision::RFPI
  double precision::RMSAVM      ! rms error magnetisation dens. (contribution of single site)
  double precision::RMSAVQ      ! rms error charge density (contribution of single site)
  double precision::EREFLDAU    ! LDA+U

  real::TIME_I
  real::TIME_S
  real::TIME_E
  real::TIME_EX

  integer::ICELL
  integer::IR
  integer::LMAX
  integer::NPNT1
  integer::NPNT2
  integer::NPNT3
  integer::NPOL
  integer::NSPIN
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
  integer::LPOT
  integer::LMPOT
  integer::NAEZ
  integer::IEND
  integer::NCLEBD
  integer::LM1
  integer::LM2
  integer::IEND1
  integer::IPOT
  integer::ISPIN
  integer::I1
  integer::I1BRYD
  integer::LM
  integer::NR
  integer::EKM
  logical::XCCPL
  logical::JIJ
  logical::STOPIT
  logical::LDORHOEF
  logical::LDAU
  logical::ERESJIJ


  ! static arrays
  double precision::BRAVAIS(3,3)
  double precision::RECBV(3,3)
  double precision::VBC(2)
  integer::ISYMINDEX(NSYMAXD)
  double precision::ECORE(20,2)
  double precision:: WORK1(2)

  !     .. Local Arrays ..

  double precision::VAV0
  double precision::VOL0

  integer::NMESH
  integer::NSYMAT
  integer::MAXMESH
  integer::NGMAX
  integer::NRMAX
  integer::NSHLG
  integer::NSHLR

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

  !     .. L-MPI
  integer::LMPIC
  integer::LSMPIC
  integer::LSMPIB

  ! S-MPI
  integer::SMPIC
  integer::SMPIB
  integer::MAPSPIN

  ! E-MPI
  integer::EMPIC
  integer::EMPIB
  integer::IE
  integer::IELAST

  !     .. ACTV-MPI
  integer::MYACTVRANK
  integer::ACTVCOMM
  integer::ACTVGROUP
  integer::ACTVSIZE
  integer::MYBCRANK
  integer::BCRANK

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

 !============================================================= CONSTANTS
  PI = 4.0D0*ATAN(1.0D0)
  FPI = 4.0D0*PI
  RFPI = SQRT(FPI)
!=============================================================

  call read_dimension_parameters()

  ! consistency check of some dimension parameters
  call consistencyCheck01(IEMXD, LMAXD, NSPIND, SMPID)

  ! from dimension parameters - calculate some derived parameters
  call getDerivedParameters(IGUESSD, IRMD, IRMIND, IRNSD, LASSLD, LM2D, LMAXD, &
                            LMAXD1, LMMAXD, LMPOTD, LMXSPD, LPOTD, LRECPOT, &
                            LRECRES2, MMAXD, NAEZD, NCLEB, NGUESSD, NPOTD, NSPIND, NTIRD)


!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call allocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  call readKKR0Input      (NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, CLEB1C, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, GSH, ICLEB1C, ICST, &
                           IEND1, IFUNM, IGUESS, ILM, IMAXSH, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JEND, JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMAX, LMPOT, LMSP, LOFLM1C, LPOT, MAXMESH, &
                           MIXING, NACLS, NAEZ, NCLS, NCORE, NFU, NR, NREF, NSPIN, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QBOUND, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, RECBV, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VOLUME0, VREF, ZAT)

  ! ---------------------------------------------------------- k_mesh
  call readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)

  call readEnergyMesh(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)

  if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  call consistencyCheck02(IELAST, IEMXD, IGUESS, IGUESSD, LMAX, LMAXD, NAEZ, NAEZD, &
                          NPNT1, NPNT2, NPNT3, NPOL, NR, NRD, NSPIN, NSPIND)

  call consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)

! ------------------------------------------------------------------

  ! TODO: get rid of LSMYRANK, LSRANK, LSMPIB, LSMPIC
  call IMPI(NAEZ,MYRANK,NROFNODES, &
            LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &
            LSMPIB,LSMPIC,LSRANK,LSMYRANK, &
            SMPIB,SMPIC,SRANK,SMYRANK, &
            EMPIB,EMPIC,ERANK,EMYRANK, &
            MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE, &
            lmpid, smpid, empid, nthrds)

!====================================================================

!=====================================================================
!     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
! ... and wait after SC-ITER loop
!=====================================================================

  ! This if closes several hundreds of lines later!
  if (LMPIC/=0.or.LSMPIC/=0) then   !     ACTVGROUP could also test EMPIC

    MYBCRANK = 0

! ========= TIMING ======================================================
    if (MYLRANK(1) == 0) then
      call CPU_TIME(TIME_I)
      open (2,file='time-info',form='formatted')
    endif
!========= TIMING END ======================================================

    CNVFAC = 1000.0D0

! initialise the arrays for (gen. Anderson/Broyden) potential mixing
    UI2 = 0.00
    VI2 = 0.00
    SM1S = 0.00
    FM1S = 0.00

! ######################################################################
! ######################################################################
    do ITER = 1, SCFSTEPS
! ######################################################################
! ######################################################################

      call CPU_TIME(TIME_S)

      EKM    = 0
      NOITER = 0

      if (MYLRANK(1)==0) then
        write(2,'(79(1H=))')
        call OUTTIME(MYLRANK(1),'started at ..........',TIME_I,ITER)
        write(2,'(79(1H=))')
      endif

      call GAUNT2(WG,YRG,LMAX)

      call MADELUNG3D(LPOT,YRG,WG,ALAT, &
      RMAX,GMAX,BRAVAIS,RECBV, &
      LMXSPD,LASSLD,LPOTD,LMPOTD, &
      NMAXD,ISHLD, &
      LMPOT,CLEB,ICLEB,IEND, &
      NCLEBD,LOFLM,DFAC, &
      NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM) ! does it have to be in SCF-loop?

      do LM = 1,LMPOTD
        CMOM(LM) = 0.0D0
        CMINST(LM) = 0.0D0
      end do

      CHRGNT = 0.0D0

      ! needed for results.f - find better solution - unnecessary I/O
      call openResults1File(IEMXD, LMAXD, NPOL)

      open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
      form='unformatted')


!N ====================================================================
!     BEGIN do loop over atoms (NMPID-parallel)
!N ====================================================================

      do I1 = 1,NAEZ
        if(MYLRANK(LMPIC)==MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

!=======================================================================
! xccpl

          XCCPL = .false.
          ERESJIJ = .false.

          ! calculate exchange couplings only at last self-consistency step and when Jij=true
          if ((ITER==SCFSTEPS).and.JIJ) XCCPL = .true.

          if (XCCPL) then

            inquire(file='ERESJIJ',exist=ERESJIJ)

            call CLSJIJ(I1,NAEZ,RR,NR,RBASIS,RCUTJIJ,NSYMAT,ISYMINDEX, &
                        IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ, &
                        nrd, nxijd)

            GMATXIJ = CZERO

          endif

! xccpl
!=======================================================================

          ! This read is probably NOT NECESSARY (except for 1st iteration) !!!
          ! TURNS out it is necessary - otherwise wrong results - why?
          read(66,rec=I1) VINS,VISP,ECORE  ! Read potential from file!!!

! LDA+U
          if (LDAU) then

            EREFLDAU = EFERMI
            EREFLDAU = 0.48       ! FIXME: hardcoded

            call LDAUINIT(I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU, &
                          VISP,NSPIN,R(1,I1),DRDI(1,I1), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          PHILDAU,UMLDAU,WMLDAU, &
                          lmax, irmd, ipand)

          endif
! LDA+U

! TIME
          call OUTTIME(MYLRANK(1),'initialized .........',TIME_I,ITER)
! TIME

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================

          do IE = 1,IELAST

            call CPU_TIME(TIME_E)

            ETIME(IE) = 0.0D0

            if (ITER==1.and.IE==1) then

              call EBALANCE('I',ITER,SCFSTEPS, &
                            IELAST,NPNT1, &
                            MYACTVRANK,ACTVCOMM, &
                            ETIME,EPROC,EPROCO, &
                            empid, iemxd)

            endif

! IE ====================================================================
            if (EMPIB==EPROC(IE)) then
! IE ====================================================================


              do RF = 1,NREF

                call TREF(EZ(IE),VREF(RF),LMAX,RMTREF(RF), &
                          TREFLL(1,1,RF),DTREFLL(1,1,RF), LLY)

              end do

              call GREF(EZ(IE),ALAT,IEND1,NCLS,NAEZ, &
                        CLEB1C,RCLS,ATOM,CLS,ICLEB1C,LOFLM1C,NACLS, &
                        REFPOT, &
                        TREFLL(1,1,1),DTREFLL(1,1,1),GREFN,DGREFN, &
                        IE, &
                        LLY_G0TR, &
                        LMPIC,MYLRANK,LCOMM,LSIZE, &
                        naez, lmax, naclsd, ncleb, nrefd, iemxd, nclsd, &
                        LLY, LMPID*SMPID*EMPID)

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================

spinloop:     do ISPIN = 1,NSPIN
!------------------------------------------------------------------------------
!         beginning of SMPID-parallel section
!------------------------------------------------------------------------------
                if (SMPID==1) then
                  MAPSPIN = 0
                  PRSPIN   = ISPIN
                else
                  MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                  PRSPIN   = 1
                endif

                if(SRANK(SMPIB,SMPIC)==MAPSPIN) then

                  call CALCTMAT(LDAU,NLDAU,ICST, &
                                NSRA,EZ(IE), &
                                DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                                LLDAU,WMLDAU(1,1,1,ISPIN), &
                                nspind, ncleb, ipand, irmd, irnsd)

                  if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula

                    call calcdtmat_DeltaEz(delta_E_z, IE, NPNT1, NPNT2, NPNT3, TK)

                    call CALCDTMAT(LDAU,NLDAU,ICST, &
                                  NSRA,EZ(IE),delta_E_z, &
                                  DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                  VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                  IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                  DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                                  LLDAU,WMLDAU(1,1,1,ISPIN), &
                                  nspind, ncleb, ipand, irmd, irnsd)
                  end if

                  ! calculate DTIXIJ = T_down - T_up
                  if (XCCPL) then
                    if (ISPIN==1) then
                      do LM1 = 1,LMMAXD
                        do LM2 = 1,LMMAXD
                          DTIXIJ(LM1,LM2) = TMATN(LM1,LM2,ISPIN)
                        enddo
                      enddo
                    else
                      do LM1 = 1,LMMAXD
                        do LM2 = 1,LMMAXD
                          DTIXIJ(LM1,LM2) = DTIXIJ(LM1,LM2)-TMATN(LM1,LM2,ISPIN)
                        enddo
                      enddo
                    endif
                  endif


                  RF = REFPOT(I1)
                  do LM1 = 1,LMMAXD
                    TMATN(LM1,LM1,ISPIN) =  TMATN(LM1,LM1,ISPIN) &
                    - TREFLL(LM1,LM1,RF)
                    DTDE(LM1,LM1,ISPIN) =  DTDE(LM1,LM1,ISPIN) &
                    - DTREFLL(LM1,LM1,RF)
                  end do

                  ! TMATN now contains Delta t = t - t_ref !!!
                  ! DTDE now contains Delta dt !!!

                  ! renormalize TR_ALPH
                  TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(IE,CLS(I1))

                  NMESH = KMESH(IE)

                  if( MYLRANK(LMPIC)==0 ) then
                    call printEnergyPoint(EZ(IE), IE, ISPIN, NMESH)
                  end if

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
                  SPRS(1,1,PRSPIN),PRSC(1,1,PRSPIN), &
                  EKM,NOITER, &
                  QMRBOUND,IGUESS,BCP,CNVFAC(1,PRSPIN), &
                  NXIJ,XCCPL,IXCP,ZKRXIJ, &
                  LLY_GRDT(IE,ISPIN),TR_ALPH(ISPIN), &
                  GMATXIJ(1,1,1,ISPIN), &
                  LCOMM(LMPIC),LSIZE(LMPIC), &
                  iemxd, nthrds, &
                  lmmaxd, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
                  nxijd, nguessd, kpoibz, nrd, ekmd)

                endif
!------------------------------------------------------------------------------
!        End of SMPID-parallel section
!------------------------------------------------------------------------------

              end do spinloop                          ! ISPIN = 1,NSPIN

! SPIN ==================================================================
!     END do loop over spins
! SPIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! xccpl

              if (XCCPL) then

                call SREDGX( NSPIN, &
                             MYRANK, &
                             SMPIC,SMYRANK, &
                             GMATXIJ, &
                             GXIJ_ALL, &
                             naez, lmax, lmpid, empid, smpid, nxijd)

                JSCAL = WEZ(IE)/DBLE(NSPIN)

                call XCCPLJIJ_START(I1,IE,JSCAL, &
                               RXIJ,NXIJ,IXCP,RXCCLS, &
                               GXIJ_ALL,DTIXIJ, &
                               LMPIC,LCOMM, &
                               JXCIJINT,ERESJIJ, &
                               naez, lmmaxd, nxijd, nspind, &
                               lmpid*smpid*empid)

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


!=======================================================================
!     "allreduce" information of 1 .. EMPID and 1 .. SMPID processors
!=======================================================================
          call SREDGM(NSPIN,IELAST, &
                      MYRANK, &
                      SMPIC,SMYRANK, &
                      EMPIC,EMYRANK,EPROC, &
                      GMATN,LLY_GRDT, &
                      GMATN_ALL,LLY_GRDT_ALL, &
                      naez, lmax, lmpid, smpid, empid, iemxd)
!=======================================================================
!=======================================================================

! TIME
          call OUTTIME(MYLRANK(1),'G obtained ..........',TIME_I,ITER)
! TIME

!=======================================================================
!     output of Jij's
!=======================================================================
          if (XCCPL) then

            call XCCPLJIJ_OUT(I1, &  ! I1 needed for filenames
                          RXIJ,NXIJ,IXCP,RXCCLS, &
                          LMPIC, &
                          MYRANK,EMPIC,EMYRANK, &
                          JXCIJINT, &
                          naez, nxijd, &
                          lmpid, smpid, empid)
          endif
!=======================================================================
!=======================================================================

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors
!=======================================================================
          call EBALANCE('R',ITER,SCFSTEPS, &
                        IELAST,NPNT1, &
                        MYACTVRANK,ACTVCOMM, &
                        ETIME,EPROC,EPROCO, &
                        empid, iemxd)
!=======================================================================
!=======================================================================


!=======================================================================
!     in case of IGUESS and EMPID > 1 initial guess arrays might
!     have to be adjusted to new distributions
!=======================================================================
          if ((IGUESS==1).and.(EMPID>1)) then

            do ISPIN = 1,NSPIN

              if (SMPID==1) then
                MAPSPIN = 0
                PRSPIN   = ISPIN
              else
                MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                PRSPIN   = 1
              endif

!       true beginning of SMPID-parallel section

              if(SRANK(SMPIB,SMPIC) == MAPSPIN) then

                call EPRDIST(IELAST,KMESH,NOFKS, &
                             PRSC(1,1,PRSPIN), &
                             SPRS(1,1,PRSPIN), &
                             CNVFAC(1,PRSPIN), &
                             MYRANK,EMPIC,EMYRANK, &
                             EPROC,EPROCO, &
                             lmpid, smpid, empid, naez, lmax, nguessd, ekmd, iemxd)

              endif
            enddo

          endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================
!=======================================================================


!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
          if (LMPIC==1) then

            if (LLY==1) then
              ! get WEZRN and RNORM, the important input from previous
              ! calculations is LLY_GRDT_ALL
              ! TODO: all THETAS passed, but only 1 needed
              ! here atom processes communicate with each other
              call LLOYD0(EZ,WEZ,CLEB1C,DRDI,R,IRMIN,VINS,VISP, &
                          THETAS,ZAT,ICLEB1C, &
                          IFUNM,IPAN,IRCUT,LMSP,JEND,LOFLM1C, &
                          NTCELL,ICST, &
                          IELAST,IEND1,NAEZ,NSPIN,NSRA, &
                          WEZRN,RNORM, &
                          GMATN_ALL, &
                          LLY_GRDT_ALL, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          LMPIC,MYLRANK, &
                          LCOMM,LSIZE, &
                          lmpid*smpid*empid, lmax, irmd, irnsd, iemxd, &
                          irid, nfund, ncelld, ipand, ncleb)

! IME
              call OUTTIME(MYLRANK(1),'Lloyd processed......',TIME_I,ITER)
! IME
            else ! no Lloyd

              do IE=1,IELAST
                WEZRN(IE,1) = WEZ(IE)
                WEZRN(IE,2) = WEZ(IE)
              enddo
            endif

            ! now WEZRN stores the weights for E-integration

            call CINIT(IEMXD*(LMAXD+2)*NSPIND,DEN)
            DENEF = 0.0D0

            if (LDAU) then
              call CINIT(MMAXD*MMAXD*NSPIND*LMAXD1,DMATLDAU(1,1,1,1))
            endif

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN ==================================================================

            do ISPIN = 1,NSPIN
              ICELL = NTCELL(I1)
              IPOT = (I1-1) * NSPIN + ISPIN

              LDORHOEF = NPOL/=0  ! needed in RHOVAL
              call RHOVAL(LDORHOEF,ICST,IELAST, &
                          NSRA,ISPIN,NSPIN,EZ,WEZRN(1,ISPIN), &
                          DRDI(1,I1),R(1,I1),IRMIN(I1), &
                          VINS(IRMIND,1,ISPIN),VISP(1,ISPIN), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          THETAS(1,1,ICELL),IFUNM(1,ICELL),LMSP(1,ICELL), &
                          RHO2NS,R2NEF, &
                          DEN(0,1,ISPIN),ESPV(0,ISPIN), &
                          CLEB1C,LOFLM1C,ICLEB1C,IEND1,JEND, &
                          GMATN_ALL, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          iemxd, &
                          lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

              call RHOCORE(E1,NSRA,ISPIN,NSPIN,I1, &  ! I1 is used only for debugging output
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
              call renormalizeDOS(DEN,RNORM,LMAXD1,IELAST,NSPIN,IEMXD)
            end if

            ! calculate DOS at Fermi level
            DENEF = calcDOSatFermi(DEN, IELAST, IEMXD, LMAXD1, NSPIN)

            ! ---> l/m_s/atom-resolved charges, output -> CHARGE
            ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
            ! CHARGE -> written to result file
            call calcChargesLres(CHARGE, DEN, IELAST, LMAXD1, NSPIN, WEZ, IEMXD)

! LDAU

            EULDAU = 0.0D0
            EDCLDAU = 0.0D0

            if (LDAU.and.NLDAU>=1) then

              call LDAUWMAT(I1,NSPIN,ITER,MIXING,DMATLDAU,NLDAU,LLDAU, &
                            ULDAU,JLDAU,UMLDAU,WMLDAU,EULDAU,EDCLDAU, &
                            lmaxd)

            endif

! LDAU

! ----------------------------------------------------------------------
! -->   determine total charge density expanded in spherical harmonics
! -------------------------------------------------------------- density

            call RHOTOTB_NEW(NSPIN,RHO2NS,RHOCAT, &
                         DRDI(:,I1),IRCUT(:,I1), &
                         LPOT,NFU(ICELL),LLMSP(1,ICELL),THETAS(:,:,ICELL),IPAN(I1), &
                         CATOM, &
                         irmd, irid, ipand, nfund)

            CHRGNT = CHRGNT + CATOM(1) - ZAT(I1)

            ! write to 'results1' - only to be read in in results.f
            ! necessary for density of states calculation, otherwise
            ! only for informative reasons
            call writeResults1File(CATOM, CHARGE, DEN, ECORE, I1, NPOL, QC)

          endif
!----------------------------------------------------------------------
! END L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
        end if
      end do

!N ====================================================================
!     END do loop over atoms (NMPID-parallel)
!N ====================================================================

      close(66)
      call closeResults1File()

      call OUTTIME(MYLRANK(1),'density calculated ..',TIME_I,ITER)


!----------------------------------------------------------------------
! BEGIN L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
      if (LMPIC==1) then

!****************************************************** MPI COLLECT DATA

        call MPI_ALLREDUCE(CHRGNT,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        call DCOPY(1,WORK1,1,CHRGNT,1)

        call MPI_ALLREDUCE(DENEF,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        call DCOPY(1,WORK1,1,DENEF,1)

!****************************************************** MPI COLLECT DATA

! ----------------------------------------------------------------------

! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT) !FIXME: hardcoded maximal shift of 0.03
        EFOLD = E2

        if (ISHIFT < 2) E2 = E2 - E2SHIFT

        if( MYLRANK(LMPIC) == 0 ) then
          call printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
        end if

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIN)
! ----------------------------------------------------------------------

        open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
        form='unformatted')
        call openResults2File(LRECRES2)

! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================
        do I1 = 1,NAEZ
          if(MYLRANK(LMPIC) == MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            ICELL = NTCELL(I1)

            do ISPIN = 1,NSPIN

! -->     get correct density and valence band energies

              ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
              EFOLD*CHRGNT/DBLE(NSPIN*NAEZ)

              do LM = 1,LMPOT
                call DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1, &
                RHO2NS(1,LM,ISPIN),1)
              end do

! ----------------------------------------------------------------------
            end do

            call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ipand, ngshd)

            call OUTTIME(MYLRANK(1),'RHOMOM ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

            call VINTRAS_NEW(LPOT,NSPIN,RHO2NS,VONS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(MYLRANK(1),'VINTRAS ......',TIME_I,ITER)

            call STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM, &
            RBASIS,SMAT,VOLUME0,LASSLD,LMXSPD,NAEZD,I1)

            call OUTTIME(MYLRANK(1),'STRMAT ......',TIME_I,ITER)

            call VMADELBLK(CMOM,CMINST,LPOT,NSPIN, &
            NAEZ,VONS,ZAT,R,IRCUT,IPAN, &
            VMAD, &
            LMPOT,SMAT,CLEB,ICLEB,IEND, &
            LMXSPD,NCLEBD,LOFLM,DFAC,I1, &
            LMPIC,MYLRANK, &
            LCOMM,LSIZE, &
            irmd, ipand, lmpid*smpid*empid)

            call OUTTIME(MYLRANK(1),'VMADELBLK ......',TIME_I,ITER)

! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

            if (KFORCE==1 .and. ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCEH(CMOM,FLM,LPOT,I1,RHO2NS,VONS, &
              R,DRDI,IMT,ZAT,irmd)
              call FORCE(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              DRDI,IMT,naez,irmd)
! ---------------------------------------------------------------------
            end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

            if (KTE==1) then
              call ESPCB(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE) !TODO

              call EPOTINB_NEW(EPOTIN,NSPIN,RHO2NS,VISP,R(:,I1),DRDI(:,I1), &
              IRMIN(I1),IRWS(I1),LPOT,VINS,IRCUT(:,I1),IPAN(I1),ZAT(I1), &
              irmd, irnsd, ipand)


              call ECOUB_NEW(CMOM,ECOU,LPOT,NSPIN,RHO2NS, &
              VONS,ZAT(I1),R(:,I1), &
              DRDI(:,I1),KVMAD,IRCUT(:,I1),IPAN(I1),IMAXSH,IFUNM(1,ICELL), &
              ILM,GSH,THETAS(:,:,ICELL),LMSP(1,ICELL), &
              irmd, irid, nfund, ipand, ngshd)

            end if
            call OUTTIME(MYLRANK(1),'KTE ......',TIME_I,ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
            call VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,I1,RHO2NS, &
            VONS,R,DRDI,A, &
            IRWS,IRCUT,IPAN,ICELL,GSH,ILM,IMAXSH,IFUNM(1,ICELL), &
            THETAS,LMSP(1,ICELL), &
            naez, irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(MYLRANK(1),'VXCDRV ......',TIME_I,ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCXC(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              ALAT,DRDI,IMT,ZAT, &
              LMPIC,MYLRANK, &
              LCOMM, &
              naez, irmd, lmpid*smpid*empid)
! ---------------------------------------------------------------------
            end if

            ! unnecessary I/O? see results.f
            call writeResults2File(CATOM, ECOU, EDCLDAU, EPOTIN, ESPC, ESPV, EULDAU, EXC, I1, LCOREMAX, VMAD)

! Force calculation ends
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! =====================================================================
            VAV0 = 0.0D0
            VOL0 = 0.0D0

            !TODO: all THETAS passed, only one used
            call MTZERO(LMPOT,I1,NSPIN,VONS,ZAT,R,DRDI,IMT,IRCUT, &
                        IPAN,ICELL,LMSP(1,ICELL),IFUNM(1,ICELL), &
                        THETAS,IRWS,VAV0,VOL0, &
                        irmd, irid, nfund, ipand)

            call OUTTIME(MYLRANK(1),'MTZERO ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
          end if
        end do
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================

        call OUTTIME(MYLRANK(1),'calculated pot ......',TIME_I,ITER)

        call allreduceMuffinTinShift_com(LCOMM(LMPIC), VAV0, VBC, VOL0)

        call shiftMuffinTinZero(ISHIFT, VBC, E2SHIFT) ! purpose? ISHIFT usually=0

        if(MYRANK==0) then
          call printMuffinTinShift(VAV0, VBC, VOL0)
        end if

! =====================================================================

! ---------------------------------------------------------------------

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

!+++++++++++++++++ BEGIN ATOM PARALLEL +++++++++++++++++++++++++++++++
        do I1 = 1,NAEZ

          if(MYLRANK(LMPIC)== &
          MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            ICELL = NTCELL(I1)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
            do ISPIN = 1,NSPIN
              IPOT = NSPIN* (I1-1) + ISPIN

              do IR = 1,IRCUT(IPAN(I1),I1)
                VONS(IR,1,ISPIN) = VONS(IR,1,ISPIN) + RFPI*VBC(ISPIN)
              end do

              ! TODO: all THETAS passed - only one used
              call CONVOL(IRCUT(1,I1),IRC(I1),ICELL, &
                          IMAXSH(LMPOT),ILM,IFUNM(1,ICELL),LMPOT,GSH, &
                          THETAS,ZAT(I1),RFPI, &
                          R(1,I1),VONS(1,1,ISPIN),LMSP(1,ICELL), &
                          irid, nfund, irmd, ngshd)

            end do

! -->   final construction of the potentials (straight mixing)
            MIX = MIXING
            RMSAVQ = 0.0D0
            RMSAVM = 0.0D0

            call MIXSTR(RMSAVQ,RMSAVM,LPOT,LMPOT, &
            I1,NSPIN, &
            ITER,RFPI,FPI, &
            MIX, &
            FCM,IRC,IRMIN,R,DRDI,VONS, &
            VISP,VINS, &
            naez, irmd, irnsd)

            I1BRYD=I1
          end if
        end do
!+++++++++++++++++ END ATOM PARALLEL +++++++++++++++++++++++++++++++

! -->  potential mixing procedures: Broyden or Andersen updating schemes
        if (IMIX>=3) then
          call BRYDBM(VISP,VONS,VINS, &
          LMPOT,R,DRDI,MIX, &
          IRC,IRMIN,NSPIN,I1BRYD, &
          IMIX,ITER, &
          UI2,VI2,WIT,SM1S,FM1S, &
          LMPIC,MYLRANK, &
          LCOMM, &
          itdbryd, irmd, irnsd, nspind, &
          LMPID * SMPID * EMPID)
        endif

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!--------- BEGIN Atom-parallel ----------------------------------------
        do I1 = 1,NAEZ
          if(MYLRANK(LMPIC)== &
          MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            call resetPotentials(IRC(I1), IRMD, IRMIN(I1), IRMIND, LMPOTD, &
                                 NSPIN, VINS, VISP, VONS) ! not sure if correct?

! ----------------------------------------------------- output_potential
            write(66,rec=I1) VINS,VISP,ECORE
! ----------------------------------------------------- output_potential
          end if
        end do
! -------- END Atom-parallel ------------------------------------------

        call closeResults2File()

! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
! ====== write RMS convergency data - not parallelized, written by
! MYRANK=0   (RMSOUT) =================================================
! write formatted potential if file VFORM exists
        call RMSOUT_com(RMSAVQ,RMSAVM,ITER,E2,EFOLD, &
        SCFSTEPS,VBC,QBOUND,NSPIN,NAEZ, &
        KXC,LPOT,A,B,IRC, &
        VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
        ECORE,LCORE,NCORE,ZAT,ITITLE, &
        LMPIC,MYLRANK, &
        LCOMM,LSIZE, &
        irmd, irnsd, lmpid*smpid*empid)

! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS

        call MPI_BARRIER(LCOMM(LMPIC),IERR)

! -----------------------------------------------------------------
! BEGIN: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------
        if(MYLRANK(LMPIC)==0) then

          ! DOS was written to file 'results1' and read out here just
          ! to be written in routine wrldos
          ! also other stuff is read from results1
          call RESULTS(LRECRES2,IELAST,ITER,LMAX,NAEZ,NPOL, &
          NSPIN,KPRE,KTE,LPOT,E1,E2,TK,EFERMI, &
          ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,LDAU, &
          iemxd)

          ! only ranks with MYLRANK(LMPIC)==0 update, other ranks get it broadcasted later
          call updateEnergyMesh(EZ,WEZ,IELAST,E1,E2,TK,NPOL,NPNT1,NPNT2,NPNT3)

          write(6,'(79(1H=))')

! .. get info on MYACTVRANK of this processor: to be used in
!    subsequent reduce-commands
          MYBCRANK = MYACTVRANK
! ..
        endif
! -----------------------------------------------------------------
! END: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------

! ..
      endif
! -----------------------------------------------------------------
! END: L-MPI: only processes with LMPIC = 1 are working here
! -----------------------------------------------------------------
      close(66)  ! close 'vpotnew'  ! this may be misplaced - could belong into previous if clause
! -----------------------------------------------------------------

      ! why? all processes except 1 have MYBCRANK = 0, this allreduce
      ! tells all the other processes who is the root
      ! not really necessary
      call MPI_ALLREDUCE(MYBCRANK,BCRANK,1,MPI_INTEGER,MPI_MAX, &
      ACTVCOMM,IERR)

      call broadcastEnergyMesh_com(ACTVCOMM, BCRANK, E1, E2, EZ, IEMXD, WEZ)

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      ACTVCOMM,IERR) ! TODO: allreduce not necessary, only master rank needs NOITER_ALL, use reduce instead

      if(MYLRANK(1)==0) then

        ! write file 'energy_mesh'
        call writeEnergyMesh(E1, E2, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ)

        call printSolverIterationNumber(ITER, NOITER_ALL)
        call writeIterationTimings(ITER, TIME_I, TIME_S)
      endif

! manual exit possible by creation of file 'STOP' in home directory

      inquire(file='STOP',exist=STOPIT)
      if (STOPIT) goto 200


! ######################################################################
! ######################################################################
    enddo          ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################


200 continue

    if (MYLRANK(1)==0) close(2)    ! TIME

    if (KFORCE==1) close(54)
! ======================================================================
! ======================================================================

! Free communicators and groups ..
! ..
    if (MYLRANK(LMPIC)>=0) then
      call MPI_COMM_FREE(LCOMM(LMPIC),IERR)
    endif

    call MPI_GROUP_FREE(LGROUP(LMPIC),IERR)

    call MPI_COMM_FREE(ACTVCOMM,IERR)
    call MPI_GROUP_FREE(ACTVGROUP,IERR)
! .. .

  endif ! ACTVGROUP

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================


!      WRITE(6,*) 'BARRIER i:',MYRANK
  call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!      WRITE(6,*) 'BARRIER f:',MYRANK
  call MPI_FINALIZE(IERR)

!-----------------------------------------------------------------------------
! Array DEallocations BEGIN
!-----------------------------------------------------------------------------
  call deallocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array DEallocations END
!-----------------------------------------------------------------------------

end program MAIN2
