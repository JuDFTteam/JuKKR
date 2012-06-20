! KKRnano
! massive parallel KKR for nanoscaled systems

program MAIN2

  !use mpi
  use common_testc
  use common_optc

  use KKRnano_mpi_mod

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
  double precision::MIXING
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
  integer::IEND
  integer::NCLEBD
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

  integer::MAPSPIN

  integer::IE
  integer::IELAST

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
                            LMAXD1, LMMAXD, LMPOTD, LMXSPD, &
                            LRECRES2, MMAXD, NAEZ, NCLEB, NGUESSD, NPOTD, NSPIND, NTIRD)


!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
  call allocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

  !every process does this!
  call readKKR0Input      (NSYMAXD, A, ALAT, ATOM, B, BCP, BRAVAIS, CLEB1C, &
                           CLS, DRDI, DSYMLL, EZOA, FCM, GMAX, GSH, ICLEB1C, ICST, &
                           IEND1, IFUNM, IGUESS, ILM, IMAXSH, IMIX, IMT, INDN0, IPAN, &
                           IRC, IRCUT, IRMIN, IRNS, IRWS, ISHIFT, ISYMINDEX, ITITLE, &
                           JEND, JIJ, KFORCE, KMESH, KPRE, KTE, KVMAD, KXC, LCORE, &
                           LDAU, LLMSP, LMSP, LOFLM1C, MAXMESH, &
                           MIXING, NACLS, NCLS, NCORE, NFU, NR, NREF, &
                           NSRA, NSYMAT, NTCELL, NUMN0, OPTC, QMRBOUND, R, &
                           RBASIS, RCLS, RCUTJIJ, RECBV, REFPOT, RMAX, RMT, RMTREF, &
                           RR, RWS, SCFSTEPS, TESTC, THETAS, VOLUME0, VREF, ZAT)

  ! ---------------------------------------------------------- k_mesh
  call readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)  !every process does this!

  call readEnergyMesh(E1, E2, EFERMI, EZ, IELAST, NPNT1, NPNT2, NPNT3, NPOL, TK, WEZ) !every process does this!

  if (KFORCE==1) open (54,file='force',form='formatted')   ! every process opens file 'force' !!!

 ! ======================================================================
 ! =                     End read in variables                          =
 ! ======================================================================

  call consistencyCheck03(ATOM, CLS, EZOA, INDN0, NACLS, NACLSD, NAEZ, NCLSD, NR, NUMN0)

  if ((JIJ .eqv. .true.) .and. (SMPID /= 1)) then
    write(*,*) "ERROR: Jij calculation is broken for spin-parallel calc. Set SMPID=1"
    stop
  end if

! -----------------------------------------------------------------------------
   call initialiseKKRnano_mpi_com(SMPID, EMPID, NAEZ, nthrds)
!------------------------------------------------------------------------------

!=====================================================================
!     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
! ... and wait after SC-ITER loop
!=====================================================================

  ! This if closes several hundreds of lines later!
  if (LMPIC/=0) then   !     ACTVGROUP could also test EMPIC or MYACTVRANK (preferred)

! ========= TIMING ======================================================
    if (is_Masterrank) then
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

    call GAUNT2(WG,YRG,LMAXD)

    call MADELUNG3D(LPOT,YRG,WG,ALAT, &
    RMAX,GMAX,BRAVAIS,RECBV, &
    LMXSPD,LASSLD,LPOT,LMPOTD, &
    NMAXD,ISHLD, &
    LMPOTD,CLEB,ICLEB,IEND, &
    NCLEBD,LOFLM,DFAC, &
    NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM, &
    MYRANK)

! ######################################################################
! ######################################################################
    do ITER = 1, SCFSTEPS
! ######################################################################
! ######################################################################

      call CPU_TIME(TIME_S)

      EKM    = 0
      NOITER = 0

      if (is_Masterrank) then
        call printDoubleLineSep(unit_number = 2)
        call OUTTIME(is_Masterrank,'started at ..........',TIME_I,ITER)
        call printDoubleLineSep(unit_number = 2)
      endif

      do LM = 1,LMPOTD
        CMOM(LM) = 0.0D0
        CMINST(LM) = 0.0D0
      end do

      CHRGNT = 0.0D0

      ! needed for results.f - find better solution - unnecessary I/O
      ! actually only LMPIC==1 process needed to open these files!!!
      call openResults1File(IEMXD, LMAXD, NPOL)

      call openPotentialFile(LMPOTD, IRNSD, IRMD)


!N ====================================================================
!     BEGIN do loop over atoms (NMPID-parallel)
!N ====================================================================

      do I1 = 1,NAEZ
        if(my_SE_rank==MAPBLOCK(I1,1,NAEZ,1,0,my_SE_comm_size-1)) then

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

            GMATXIJ = CZERO

          endif

          ! TODO: Find a better solution than exchanging potential by file I/O
          call readPotential(I1, VISP, VINS, ECORE)

! LDA+U
          if (LDAU) then

            EREFLDAU = EFERMI
            EREFLDAU = 0.48       ! FIXME: hardcoded

            call LDAUINIT(I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU, &
                          VISP,NSPIND,R(1,I1),DRDI(1,I1), &
                          ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                          PHILDAU,UMLDAU,WMLDAU, &
                          lmaxd, irmd, ipand)

          endif
! LDA+U

! TIME
          call OUTTIME(is_Masterrank,'initialized .........',TIME_I,ITER)
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
            if (my_energy_rank==EPROC(IE)) then
! IE ====================================================================

              do RF = 1,NREF
                call TREF(EZ(IE),VREF(RF),LMAXD,RMTREF(RF), &
                          TREFLL(1,1,RF),DTREFLL(1,1,RF), LLY)
              end do

              call GREF_com(EZ(IE),ALAT,IEND1,NCLS,NAEZ, &
                            CLEB1C,RCLS,ATOM,CLS,ICLEB1C,LOFLM1C,NACLS, &
                            REFPOT, &
                            TREFLL,DTREFLL,GREFN,DGREFN, &
                            LLY_G0TR(:,IE), &
                            my_SE_rank,my_SE_communicator,my_SE_comm_size, &
                            lmaxd, naclsd, ncleb, nrefd, nclsd, &
                            LLY)

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================

spinloop:     do ISPIN = 1,NSPIND
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

                if(my_spin_rank==MAPSPIN) then

                  call CALCTMAT(LDAU,NLDAU,ICST, &
                                NSRA,EZ(IE), &
                                DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAXD, &
                                LLDAU,WMLDAU(1,1,1,ISPIN), &
                                ncleb, ipand, irmd, irnsd)

                  if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula

                    call calcdtmat_DeltaEz(delta_E_z, IE, NPNT1, NPNT2, NPNT3, TK)

                    call CALCDTMAT(LDAU,NLDAU,ICST, &
                                  NSRA,EZ(IE),delta_E_z, &
                                  DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                  VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                  IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                  DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAXD, &
                                  LLDAU,WMLDAU(1,1,1,ISPIN), &
                                  ncleb, ipand, irmd, irnsd)
                  end if

                  ! calculate DTIXIJ = T_down - T_up
                  if (XCCPL) then
                    call calcDeltaTupTdown(DTIXIJ, ISPIN, LMMAXD, TMATN)
                  endif


                  RF = REFPOT(I1)
                  call substractReferenceTmatrix(TMATN(:,:,ISPIN), TREFLL(:,:,RF), LMMAXD)
                  ! do the same for derivative of T-matrix
                  call substractReferenceTmatrix(DTDE(:,:,ISPIN), DTREFLL(:,:,RF), LMMAXD)

                  ! TMATN now contains Delta t = t - t_ref !!!
                  ! DTDE now contains Delta dt !!!

                  ! renormalize TR_ALPH
                  TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(CLS(I1), IE)

                  NMESH = KMESH(IE)

                  if( my_SE_rank==0 ) then
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
                  my_SE_communicator,my_SE_comm_size, &
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

                call SREDGX( NSPIND, &
                             MYRANK, &
                             SMPIC,SMYRANK, &
                             GMATXIJ, &
                             GXIJ_ALL, &
                             naez, lmaxd, lmpid, empid, smpid, nxijd)

                JSCAL = WEZ(IE)/DBLE(NSPIND)

                call XCCPLJIJ_START(I1,IE,JSCAL, &
                               RXIJ,NXIJ,IXCP,RXCCLS, &
                               GXIJ_ALL,DTIXIJ, &
                               my_SE_communicator, &
                               JXCIJINT,ERESJIJ, &
                               naez, lmmaxd, nxijd, nspind)

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
          call SREDGM(NSPIND,IELAST, &
                      MYRANK, &
                      SMPIC,SMYRANK, &
                      EMPIC,EMYRANK,EPROC, &
                      GMATN,LLY_GRDT, &
                      GMATN_ALL,LLY_GRDT_ALL, &
                      naez, lmmaxd, lmpid, smpid, empid, iemxd)
!=======================================================================
!=======================================================================

! TIME
          call OUTTIME(is_Masterrank,'G obtained ..........',TIME_I,ITER)
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

            do ISPIN = 1,NSPIND

              if (SMPID==1) then
                MAPSPIN = 0
                PRSPIN   = ISPIN
              else
                MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                PRSPIN   = 1
              endif

!       true beginning of SMPID-parallel section

              if(my_spin_rank == MAPSPIN) then

                call EPRDIST(IELAST,KMESH,NOFKS, &
                             PRSC(1,1,PRSPIN), &
                             SPRS(1,1,PRSPIN), &
                             CNVFAC(1,PRSPIN), &
                             MYRANK,EMPIC,EMYRANK, &
                             EPROC,EPROCO, &
                             lmpid, smpid, empid, naez, lmaxd, nguessd, ekmd, iemxd)

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
              ! TODO: all THETAS passed, but only 1 needed, also ZAT
              ! here atom processes communicate with each other - LLOYD0_NEW already written
              call LLOYD0(EZ,WEZ,CLEB1C,DRDI,R,IRMIN,VINS,VISP, &
                          THETAS,ZAT,ICLEB1C, &
                          IFUNM,IPAN,IRCUT,LMSP,JEND,LOFLM1C, &
                          NTCELL,ICST, &
                          IELAST,IEND1,NAEZ,NSPIND,NSRA, &
                          WEZRN,RNORM, &
                          GMATN_ALL, &
                          LLY_GRDT_ALL, &
                          LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                          DMATLDAU, &
                          my_SE_rank, &
                          my_SE_communicator,my_SE_comm_size, &
                          lmaxd, irmd, irnsd, iemxd, &
                          irid, nfund, ncelld, ipand, ncleb)

! IME
              call OUTTIME(is_Masterrank,'Lloyd processed......',TIME_I,ITER)
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
              call RHOVAL(LDORHOEF,ICST,IELAST, &
                          NSRA,ISPIN,NSPIND,EZ,WEZRN(1,ISPIN), &   ! unfortunately spin-dependent
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

            call RHOTOTB_NEW(NSPIND,RHO2NS,RHOCAT, &
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

      call closePotentialFile()
      call closeResults1File()

      call OUTTIME(is_Masterrank,'density calculated ..',TIME_I,ITER)


! TODO: Only necessary for non-DOS calculation - otherwise proceed to RESULTS
!       Drawback: RESULTS has to be modified - extract DOS part
!----------------------------------------------------------------------
! BEGIN L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
      if (LMPIC==1) then

!****************************************************** MPI COLLECT DATA

        call MPI_ALLREDUCE(CHRGNT,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        my_SE_communicator,IERR)
        call DCOPY(1,WORK1,1,CHRGNT,1)

        call MPI_ALLREDUCE(DENEF,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM, &
        my_SE_communicator,IERR)
        call DCOPY(1,WORK1,1,DENEF,1)

!****************************************************** MPI COLLECT DATA

! ----------------------------------------------------------------------

! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT) !FIXME: hardcoded maximal shift of 0.03
        EFOLD = E2

        if (ISHIFT < 2) E2 = E2 - E2SHIFT

        if( my_SE_rank == 0 ) then
          call printFermiEnergy(DENEF, E2, E2SHIFT, EFOLD, NAEZ)
        end if

        call openPotentialFile(LMPOTD, IRNSD, IRMD)
        call openResults2File(LRECRES2)

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIND)
! ----------------------------------------------------------------------
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================
        do I1 = 1,NAEZ
          if(my_SE_rank == MAPBLOCK(I1,1,NAEZ,1,0,my_SE_comm_size-1)) then

            ICELL = NTCELL(I1)

            do ISPIN = 1,NSPIND

! -->     get correct density and valence band energies

              ESPV(0,ISPIN) = ESPV(0,ISPIN) - &
              EFOLD*CHRGNT/DBLE(NSPIND*NAEZ)

              do LM = 1,LMPOTD
                call DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1, &
                RHO2NS(1,LM,ISPIN),1)
              end do

! ----------------------------------------------------------------------
            end do

            call RHOMOM_NEW(CMOM,CMINST,LPOT,RHO2NS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ipand, ngshd)

            call OUTTIME(is_Masterrank,'RHOMOM ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

            call VINTRAS_NEW(LPOT,NSPIND,RHO2NS,VONS, &
            R(:,I1),DRDI(:,I1),IRCUT(:,I1),IPAN(I1),ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(is_Masterrank,'VINTRAS ......',TIME_I,ITER)

            call STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM, &
            RBASIS,SMAT,VOLUME0,LASSLD,LMXSPD,naez,I1)

            call OUTTIME(is_Masterrank,'STRMAT ......',TIME_I,ITER)

            call VMADELBLK_com(CMOM,CMINST,LPOT,NSPIND, &
            NAEZ,VONS,ZAT,R,IRCUT,IPAN, &
            VMAD, &
            LMPOTD,SMAT,CLEB,ICLEB,IEND, &
            LMXSPD,NCLEBD,LOFLM,DFAC,I1, &
            my_SE_rank, &
            my_SE_communicator,my_SE_comm_size, &
            irmd, ipand)

            call OUTTIME(is_Masterrank,'VMADELBLK ......',TIME_I,ITER)

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
              call ESPCB(ESPC,NSPIND,I1,ECORE,LCORE,LCOREMAX,NCORE) !TODO

              call EPOTINB_NEW(EPOTIN,NSPIND,RHO2NS,VISP,R(:,I1),DRDI(:,I1), &
              IRMIN(I1),IRWS(I1),LPOT,VINS,IRCUT(:,I1),IPAN(I1),ZAT(I1), &
              irmd, irnsd, ipand)

              call ECOUB_NEW(CMOM,ECOU,LPOT,NSPIND,RHO2NS, &
              VONS,ZAT(I1),R(:,I1), &
              DRDI(:,I1),KVMAD,IRCUT(:,I1),IPAN(I1),IMAXSH,IFUNM(1,ICELL), &
              ILM,GSH,THETAS(:,:,ICELL),LMSP(1,ICELL), &
              irmd, irid, nfund, ipand, ngshd)

            end if

            call OUTTIME(is_Masterrank,'KTE ......',TIME_I,ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================

            !call VXCDRV(EXC,KTE,KXC,LPOT,NSPIND,I1,RHO2NS, &
            !VONS,R,DRDI,A, &
            !IRWS,IRCUT,IPAN,ICELL,GSH,ILM,IMAXSH,IFUNM(1,ICELL), &
            !THETAS,LMSP(1,ICELL), &
            !naez, irmd, irid, nfund, ngshd, ipand)
            call VXCDRV_NEW(EXC,KTE,KXC,LPOT,NSPIND,RHO2NS, &
            VONS,R(:,I1),DRDI(:,I1),A(I1), &
            IRWS(I1),IRCUT(:,I1),IPAN(I1),GSH,ILM,IMAXSH,IFUNM(1,ICELL), &
            THETAS(:,:,ICELL),LMSP(1,ICELL), &
            irmd, irid, nfund, ngshd, ipand)

            call OUTTIME(is_Masterrank,'VXCDRV ......',TIME_I,ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCXC_com(FLM,FLMC,LPOT,NSPIND,I1,RHOCAT,VONS,R, &
              ALAT,DRDI,IMT,ZAT, &
              my_SE_rank, &
              my_SE_communicator, &
              naez, irmd)
! ---------------------------------------------------------------------
            end if

            ! unnecessary I/O? see results.f
            call writeResults2File(CATOM, ECOU, EDCLDAU, EPOTIN, ESPC, ESPV, EULDAU, EXC, I1, LCOREMAX, VMAD)

! Force calculation ends
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! =====================================================================
            VAV0 = 0.0D0
            VOL0 = 0.0D0

            call MTZERO_NEW(LMPOTD,NSPIND,VONS,ZAT(I1),R(:,I1),DRDI(:,I1),IMT(I1),IRCUT(:,I1), &
                            IPAN(I1),LMSP(1,ICELL),IFUNM(1,ICELL), &
                            THETAS(:,:,ICELL),IRWS(I1),VAV0,VOL0, &
                            irmd, irid, nfund, ipand)

            call OUTTIME(is_Masterrank,'MTZERO ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
          end if
        end do
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================

        call OUTTIME(is_Masterrank,'calculated pot ......',TIME_I,ITER)

        call allreduceMuffinTinShift_com(my_SE_communicator, VAV0, VBC, VOL0)

        call shiftMuffinTinZero(ISHIFT, VBC, E2SHIFT) ! purpose? ISHIFT usually=0

        if(is_Masterrank) then
          call printMuffinTinShift(VAV0, VBC, VOL0)
        end if

! =====================================================================

! ---------------------------------------------------------------------

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

!+++++++++++++++++ BEGIN ATOM PARALLEL +++++++++++++++++++++++++++++++
        do I1 = 1,NAEZ
          if(my_SE_rank== MAPBLOCK(I1,1,NAEZ,1,0,my_SE_comm_size-1)) then

            ICELL = NTCELL(I1)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
            do ISPIN = 1,NSPIND
              IPOT = NSPIND* (I1-1) + ISPIN

              do IR = 1,IRCUT(IPAN(I1),I1)
                VONS(IR,1,ISPIN) = VONS(IR,1,ISPIN) + RFPI*VBC(ISPIN)
              end do

              call CONVOL_NEW(IRCUT(1,I1),IRC(I1), &
                          IMAXSH(LMPOTD),ILM,IFUNM(1,ICELL),LMPOTD,GSH, &
                          THETAS(:,:,ICELL),ZAT(I1), &
                          R(1,I1),VONS(1,1,ISPIN),LMSP(1,ICELL), &
                          irid, nfund, irmd, ngshd)

            end do

! -->   final construction of the potentials (straight mixing)
            RMSAVQ = 0.0D0
            RMSAVM = 0.0D0

            call MIXSTR(RMSAVQ,RMSAVM,LPOT,LMPOTD, &
            I1,NSPIND, &
            ITER,RFPI,FPI, &
            MIXING, &
            FCM,IRC,IRMIN,R,DRDI,VONS, &
            VISP,VINS, &
            naez, irmd, irnsd)

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
          my_SE_rank, &
          my_SE_communicator, &
          itdbryd, irmd, irnsd, nspind)
        endif

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!--------- BEGIN Atom-parallel ----------------------------------------
        do I1 = 1,NAEZ
          if(my_SE_rank == MAPBLOCK(I1,1,NAEZ,1,0,my_SE_comm_size-1)) then

            call resetPotentials(IRC(I1), IRMD, IRMIN(I1), IRMIND, LMPOTD, &
                                 NSPIND, VINS, VISP, VONS) ! Note: only LMPIC=1 processes

! ----------------------------------------------------- output_potential
            call writePotential(I1, VISP, VINS, ECORE)
! ----------------------------------------------------- output_potential

          end if
        end do
! -------- END Atom-parallel ------------------------------------------

        call closeResults2File()

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
        my_SE_rank, &
        my_SE_communicator,my_SE_comm_size, &
        irmd, irnsd)

! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS

        call MPI_BARRIER(my_SE_communicator,IERR)

! -----------------------------------------------------------------
        call closePotentialFile()  ! close 'vpotnew' this was originally misplaced
! -----------------------------------------------------------------
      endif
! -----------------------------------------------------------------
! END: L-MPI: only processes with LMPIC = 1 are working here

! -----------------------------------------------------------------
! BEGIN: only MASTERRANK is working here
! -----------------------------------------------------------------
      if(is_Masterrank) then

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

        if (MYACTVRANK /= 0) then   !DEBUG
          write(*,*) "main2: Assertion MYACTVRANK==0 for MASTERRANK failed."
          write(*,*) MYACTVRANK, my_SE_rank
          STOP
        end if

        call printDoubleLineSep()

      endif
! -----------------------------------------------------------------
! END: only MASTERRANK is working here
! -----------------------------------------------------------------

      call broadcastEnergyMesh_com(ACTVCOMM, 0, E1, E2, EZ, IEMXD, WEZ) ! BCRANK = 0

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      ACTVCOMM,IERR) ! TODO: allreduce not necessary, only master rank needs NOITER_ALL, use reduce instead

      if(is_Masterrank) then

        call printSolverIterationNumber(ITER, NOITER_ALL)
        call writeIterationTimings(ITER, TIME_I, TIME_S)

      endif

! manual exit possible by creation of file 'STOP' in home directory
      if (isManualAbort_com(MYACTVRANK, ACTVCOMM) .eqv. .true.) exit


! ######################################################################
! ######################################################################
    enddo          ! SC ITERATION LOOP ITER=1, SCFSTEPS
! ######################################################################
! ######################################################################

    if (is_Masterrank) close(2)    ! TIME

    if (KFORCE==1) close(54)
! ======================================================================
! ======================================================================

  endif ! ACTVGROUP

!=====================================================================
!     processors not fitting in NAEZ*LMPID do nothing ...
! ... and wait here
!=====================================================================
! Free communicators and groups ..

  call finaliseKKRnano_mpi_com()

!-----------------------------------------------------------------------------
! Array DEallocations BEGIN
!-----------------------------------------------------------------------------
  call deallocate_main2_arrays()
!-----------------------------------------------------------------------------
! Array DEallocations END
!-----------------------------------------------------------------------------

end program MAIN2
