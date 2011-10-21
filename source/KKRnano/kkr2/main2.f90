
! KKRnano
! massive parallel KKR for nanoscaled systems


program MAIN2

  use mpi
  use common_mpi
  use inc_p_wrapper_module
  implicit none

  !     .. Parameters ..
  integer::   LMMAXD
  parameter (LMMAXD= (LMAXD+1)**2)
  integer::   NPOTD
  parameter (NPOTD=NSPIND*NAEZD)
  integer::   LMAXD1
  parameter (LMAXD1=LMAXD+1)
  integer::    MMAXD
  parameter (MMAXD  = 2*LMAXD + 1)
  integer::   LM2D
  parameter (LM2D= (2*LMAXD+1)**2)
  integer::   LMXSPD
  parameter (LMXSPD= (2*LPOTD+1)**2)
  integer::   LASSLD
  parameter (LASSLD=4*LMAXD)
  integer::   LMPOTD
  parameter (LMPOTD= (LPOTD+1)**2)
  integer::   IRMIND
  parameter (IRMIND=IRMD-IRNSD)
  integer::   LRECPOT
  parameter (LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20))
  integer::   LRECRES1
  integer::LRECRES2
  parameter (LRECRES2=4+8*(NSPIND*(LMAXD+7)+2*LPOTD+4+2))
  integer::   NTIRD
  parameter (NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND)
  integer::   MAXMSHD
  parameter (MAXMSHD=8)
  integer::   LLYALM
  parameter (LLYALM=LLY*(NAEZD*LMMAXD-1)+1)
  integer::   LLYNGD
  parameter (LLYNGD=LLY*(NACLSD*LMMAXD-1)+1)
  integer::   NSYMAXD
  parameter (NSYMAXD=48)
  double complex :: CZERO
  parameter      (CZERO=(0.0D0,0.0D0))
  !     ..
  !     .. Local Scalars ..

  double complex :: JSCAL        ! scaling factor for Jij calculation
  double precision::DENEF
  double precision::E1
  double precision::E2
  double precision::TK
  double precision::EFERMI
  double precision::EBOT
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
  double precision::RMSAVM
  double precision::RMSAVQ
  double precision::EREFLDAU    ! LDA+U

  double precision::WALLCLOCK_I ! time management
  double precision::WALLCLOCK_F
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
  integer::IPF
  integer::ISHIFT
  integer::KPRE
  integer::KTE
  integer::KVMAD
  integer::KXC
  integer::KFORCE
  integer::L
  integer::LPOT
  integer::LMPOT
  integer::NAEZ
  integer::IEND
  integer::NCLEBD
  integer::LM1
  integer::LM2
  integer::IEND1
  integer::I
  integer::J
  integer::IPOT
  integer::ISPIN
  integer::I1
  integer::I1BRYD
  integer::IH
  integer::IRC1
  integer::IRMIN1
  integer::LM
  integer::NR
  integer::EKM
  integer::SYSTEM_I
  integer::SYSTEM_F
  integer::RATETIME
  integer::MAXTIME
  logical::TEST
  logical::LCORDENS
  logical::XCCPL
  logical::JIJ
  logical::STOPIT
  logical::LDORHOEF
  logical::LDAU
  logical::ERESJIJ


  !     .. Local Arrays ..
  double complex :: EZ(IEMXD)
  double complex :: WEZ(IEMXD)
  double complex :: DEZ(IEMXD)
  double complex :: WEZRN(IEMXD,2)
  double complex :: DEN(0:LMAXD1,IEMXD,NSPIND)
  double complex :: DSYMLL(LMMAXD,LMMAXD,48)
  double complex :: PHILDAU(IRMD,LMAXD1)
  double complex :: DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1) ! LDA+U
  double precision::WG(LASSLD)
  double precision::YRG(LASSLD,0:LASSLD,0:LASSLD)
  double precision::BRAVAIS(3,3)
  double precision::RBASIS(3,NAEZD)
  double precision::RECBV(3,3)
  double precision::VBC(2)
  double precision::VAV0
  double precision::VOL0
  double precision::SMAT(LMXSPD,NAEZD)
  double precision::CLEB(LMXSPD*LMPOTD)
  double precision::CLEB1C(NCLEB,2)
  double precision::RNORM(IEMXD,2)
  double precision::DFAC(0:LPOTD,0:LPOTD)
  double precision::RWS(NAEZD)
  double precision::RMT(NAEZD)
  double precision::GN(3,NMAXD)
  double precision::RM(3,NMAXD)
  double precision::BZKP(3,KPOIBZ,MAXMSHD)
  double precision::VOLCUB(KPOIBZ,MAXMSHD)
  double precision::VOLBZ(MAXMSHD)
  double precision::RR(3,0:NRD)
  double precision::VINS(IRMIND:IRMD,LMPOTD,2) ! .. input potential (spherical VISP, nonspherical VINS)
  double precision::VISP(IRMD,2)
  double precision::VONS(IRMD,LMPOTD,2)        !     .. output potential (nonspherical VONS)
  double precision::ULDAU(LMAXD1)              ! LDA+U
  double precision::JLDAU(LMAXD1)              ! LDA+U
  double precision::UMLDAU(MMAXD,MMAXD,MMAXD,MMAXD,LMAXD1) ! LDA+U
  double precision::WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)

  integer::NMESH
  integer::KMESH(IEMXD)
  integer::NOFKS(MAXMSHD)
  integer::EZOA(NACLSD,NAEZD)
  integer::NSYMAT
  integer::NUMN0(NAEZD)
  integer::INDN0(NAEZD,NACLSD)
  integer::ID
  integer::MAXMESH
  integer::NSG(ISHLD)
  integer::NSR(ISHLD)
  integer::LLDAU(LMAXD1)    ! LDA+U
  integer::NGMAX
  integer::NRMAX
  integer::NSHLG
  integer::NSHLR
  !     ..
  ! ----------------------------------------------------------------------
  !   ECOU(0:LPOTD,NAEZD)    ! Coulomb energy
  !   EPOTIN(NAEZD),         ! energy of input potential (EPOTINB
  !   ESPC(0:3,NPOTD),        ! energy single particle core
  !   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence
  !   EXC(0:LPOTD,NAEZD),    ! E_xc
  ! ----------------------------------------------------------------------
  double complex :: TMATN(LMMAXD,LMMAXD,NSPIND)
  double complex :: DTDE(LMMAXD,LMMAXD,NSPIND)
  double complex :: TREFLL(LMMAXD,LMMAXD,NREFD)
  double complex :: DTREFLL(LMMAXD,LMMAXD,NREFD)
  double complex :: DGREFN(LMMAXD,LMMAXD,NACLSD,NCLSD)
  double complex :: GREFN(LMMAXD,LMMAXD,NACLSD,NCLSD)
  double complex :: GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND)
  double complex :: GMATN_ALL(LMMAXD,LMMAXD,IEMXD,NSPIND)
  double complex :: DTIXIJ(LMMAXD,LMMAXD)
  double complex :: GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND)
  double complex :: GXIJ_ALL(LMMAXD,LMMAXD,NXIJD,NSPIND)

  !----- Initial guess ---------------------------------------------------
  complex::          PRSC(NGUESSD*LMMAXD,EKMD,NSPIND-SMPID+1)
  double precision:: QMRBOUND
  double precision:: CNVFAC(EKMD,NSPIND-SMPID+1)
  integer::          IGUESS
  integer::BCP
  integer::PRSPIN
  integer::          SPRS(NGUESSD*LMMAXD+1,EKMD+1,NSPIND-SMPID+1)

  !----- Lloyd -----------------------------------------------------------
  double complex :: LLY_G0TR(IEMXD,NCLSD)
  double complex :: LLY_GRDT(IEMXD,NSPIND)
  double complex :: LLY_GRDT_ALL(IEMXD,NSPIND)
  double complex :: TR_ALPH(NSPIND)
  double complex :: JXCIJINT(NXIJD)
  ! ----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  double precision::ECOU(0:LPOTD)
  double precision::EPOTIN
  double precision::ESPC(0:3,NSPIND)
  double precision::ESPV(0:LMAXD1,NSPIND)
  double precision::EXC(0:LPOTD)
  double precision::EULDAU
  double precision::EDCLDAU
  double precision::A(NAEZD)
  double precision::B(NAEZD)
  double precision::DRDI(IRMD,NAEZD)
  double precision::R(IRMD,NAEZD)
  double precision::ECORE(20,2)
  double precision::THETAS(IRID,NFUND,NCELLD)
  double precision::ZAT(NAEZD)

  ! ----------------------------------------------------------------------
  double precision:: RHOCAT(IRMD,2)
  double precision:: R2NEF(IRMD,LMPOTD,2)
  double precision:: RHO2NS(IRMD,LMPOTD,2)
  double precision:: CHARGE(0:LMAXD1,2)
  !     ..
  double precision::GSH(NGSHD)
  double precision::CMINST(LMPOTD)    ! charge moment of interstitial
  double precision::CMOM(LMPOTD)      ! LM moment of total charge
  double precision::CATOM(NSPIND)     ! total charge per atom
  double precision::QC                ! core charge
  double precision::VMAD
  !     ,,
  !     .. FORCES
  double precision::FLM(-1:1,NAEZD)
  double precision::FLMC(-1:1,NAEZD)
  !     .. MIXING
  double precision::SM1S(NTIRD)
  double precision::FM1S(NTIRD)
  double precision::UI2(NTIRD,2:ITDBRYD)
  double precision::VI2(NTIRD,2:ITDBRYD)
  double precision::WIT(2:ITDBRYD)

  ! ----------------------------------------------------------------------
  integer::IMT(NAEZD)
  integer::IPAN(NAEZD)
  integer::IRC(NAEZD)
  integer::IRNS(NAEZD)
  integer::IRCUT(0:IPAND,NAEZD)
  integer::IRMIN(NAEZD)
  integer::IRWS(NAEZD)
  integer::ITITLE(20,NPOTD)
  integer::LCORE(20,NPOTD)
  integer::LCOREMAX
  integer::LLMSP(NFUND,NAEZD)
  integer::NCORE(NPOTD)
  integer::NFU(NAEZD)
  integer::NTCELL(NAEZD)
  integer::ISYMINDEX(48)
  integer::ILM(NGSHD,3)
  integer::IMAXSH(0:LMPOTD)
  integer::ICLEB(LMXSPD*LMPOTD,3)
  integer::LOFLM(LMXSPD)
  integer::ICLEB1C(NCLEB,3)
  integer::LOFLM1C(LM2D)
  integer::IFUNM(LMXSPD,NAEZD)
  integer::ICST
  integer::NSRA
  integer::LMSP(LMXSPD,NAEZD)
  integer::JEND(LMPOTD,0:LMAXD,0:LMAXD)
  integer::NCLS
  integer::NREF
  integer::RF
  integer::NLDAU
  double precision::RCLS(3,NACLSD,NCLSD)
  double precision::RMTREF(NREFD)
  double precision::VREF(NAEZD)

  double precision::RXIJ(NXIJD)          ! interatomic distance Ri-Rj
  double precision::RXCCLS(3,NXIJD)      ! position relative of j rel. to i (sorted)
  double precision::ZKRXIJ(48,3,NXIJD)   ! set up in clsjij, used in kkrmat01
  integer::IXCP(NXIJD)                   ! index to atom in elem/cell at site in cluster
  integer::NXCP(NXIJD)                   ! index to bravais lattice at site in cluster
  integer::XIJ
  integer::NXIJ
  integer::ATOM(NACLSD,NAEZD)
  integer::CLS(NAEZD)
  integer::NACLS(NCLSD)
  integer::REFPOT(NAEZD)

  integer::NUTRC                         ! number of inequivalent atoms in the cluster
  integer::INTRC(NATRCD)                 ! pointer to atoms in the unit cell
  integer::NATRC                         ! number of atoms in cluster
  integer::ATTRC(NATRCD)                 ! index to atom in elem/cell at site in cluster
  integer::EZTRC(NATRCD)                 ! index to bravais lattice  at site in cluster

  !-----------------------------------------------------------------------
  !     ..
  !-----------------------------------------------------------------------
  !     .. MPI ..
  !     .. N-MPI

  integer::   IERR
  integer::   MAPBLOCK
  external     MAPBLOCK

  !    Now in common_mpi
  !    COMMON       /MPI/MYRANK,NROFNODES
  !    INTEGER ::      MYRANK,NROFNODES

  !     .. L-MPI
  integer::MYLRANK(LMPID*SMPID*EMPID)
  integer::LCOMM(LMPID*SMPID*EMPID)
  integer::LGROUP(LMPID*SMPID*EMPID)
  integer::LSIZE(LMPID*SMPID*EMPID)
  integer::LMPIC

  !     .. LS-MPI
  integer::LSRANK(LMPID,NAEZD*SMPID*EMPID)
  integer::LSMYRANK(LMPID,NAEZD*SMPID*EMPID)
  integer::LSMPIC
  integer::LSMPIB
  !     .. S-MPI
  integer::SRANK(SMPID,NAEZD*LMPID*EMPID)
  integer::SMYRANK(SMPID,NAEZD*LMPID*EMPID)
  integer::SMPIC
  integer::SMPIB
  integer::MAPSPIN

  !     .. E-MPI
  real::ETIME(IEMXD)
  integer::EMYRANK(EMPID,NAEZD*LMPID*SMPID)
  integer::ERANK(EMPID,NAEZD*LMPID*SMPID)
  integer::EMPIC
  integer::EMPIB
  integer::IE
  integer::IELAST
  integer::EPROC(IEMXD)
  integer::EPROCO(IEMXD)

  !     .. ACTV-MPI
  integer::MYACTVRANK
  integer::ACTVCOMM
  integer::ACTVGROUP
  integer::ACTVSIZE
  integer::MYBCRANK
  integer::BCRANK
  double precision :: WORK1(2)
  double precision :: WORK2(2)

  external MPI_ALLREDUCE,MPI_FINALIZE,MPI_SEND
  !     ..
  !     .. Arrays in Common ..
  character(len=8)::OPTC(8)
  character(len=8)::TESTC(16)
!     ..
!     .. Common blocks ..
  common /OPTC/OPTC
  common /TESTC/TESTC
!     ..
!     .. External Subroutines ..
  external CINIT,CONVOL,ECOUB,EPOTINB,ESPCB, &
  FORCE,FORCEH,FORCXC,MTZERO,KLOOPZ1, &
  TEST,VMADELBLK,VXCDRV,RMSOUT, &
  OUTTIME,TIME

!     .. Intrinsic Functions ..
  intrinsic DABS,ATAN,DMIN1,DSIGN,SQRT,MAX,DBLE
!     ..
!============================================================= CONSTANTS
  LCORDENS=.true.
  PI = 4.0D0*ATAN(1.0D0)
  FPI = 4.0D0*PI
  RFPI = SQRT(FPI)
  IPF = 74

! ======================================================================
! =             read in variables from unformatted files               =
! ======================================================================
  open (67,file='inp.unf',form='unformatted')
  read (67) KMESH,MAXMESH,RR,EZOA,NUMN0,INDN0,NSYMAT,DSYMLL
  read (67) NAEZ,NSPIN,IPAN,IRNS,IRCUT,LCORE,NCORE,NTCELL, &
  LMAX,LPOT,LMPOT
  read (67) IMIX,MIXING,QBOUND,FCM,KPRE,KTE, &
  KVMAD,KXC,ISHIFT,KFORCE,IGUESS,BCP,QMRBOUND
  read (67) A,B,DRDI,R,THETAS,ZAT,IMT,IRC, &
  IRMIN,IRWS,RWS,RMT,ITITLE,LLMSP,NFU
  read (67) ALAT,GSH,ILM,IMAXSH,TESTC,OPTC
  read (67) BRAVAIS,RBASIS,RECBV,VOLUME0,RMAX,GMAX
  read (67) IEND1,CLEB1C,ICLEB1C,LOFLM1C,JEND,IFUNM,LMSP,NSRA,ICST
  read (67) NCLS,NREF,RMTREF,VREF,RCLS, &
  ATOM,CLS,NACLS,REFPOT
  read (67) NR,RCUTJIJ,JIJ,LDAU
  read (67) ISYMINDEX,SCFSTEPS
  close (67)
! ---------------------------------------------------------- k_mesh

  open (52,file='kpoints',form='formatted')
  rewind (52)

  do L = 1,MAXMESH
    read (52,fmt='(I8,f15.10)') NOFKS(L),VOLBZ(L)
    read (52,fmt=*) (BZKP(ID,1,L),ID=1,3),VOLCUB(1,L)
    do I=2,NOFKS(L)
      read (52,fmt=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
    end do
  end do

  close (52)
  open (67,file='energy_mesh',form='unformatted')
  read (67) IELAST,EZ,WEZ,E1,E2
  read (67) NPOL,TK,NPNT1,NPNT2,NPNT3

  if ( NPOL==0 ) read(67) EFERMI
  close (67)
  if (KFORCE==1) open (54,file='force',form='formatted')

! ======================================================================
! =                     End read in variables                          =
! ======================================================================

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


  if (LMPIC/=0.or.LSMPIC/=0) then   !     ACTVGROUP

    MYBCRANK = 0

    open (IPF,access='direct',recl=NSPIN*8,file='RMS.unf', &
    form='unformatted')

! ======================================================================
! ========= TIMING ======================================================
    if (MYLRANK(1) == 0) then
      RATETIME = 100
      MAXTIME  = 100000

      call SYSTEM_CLOCK(SYSTEM_I,RATETIME,MAXTIME)

      WALLCLOCK_I = MPI_WTIME()

      call CPU_TIME(TIME_I)

      open (2,file='time-info',form='formatted')
    endif
!========= TIMING ======================================================

! -------------------------------------------------------- not.converged

    open (28,file='not.converged',form='formatted')
    read (28,'(1P,4D17.10)') EFOLD,VBC
    close (28)

    do IH=1,EKMD
      do ISPIN=1,NSPIND-SMPID+1
        CNVFAC(IH,ISPIN) = 1000.0D0
      enddo
    enddo

! initialise the arrays for (gen. Anderson/Broyden) potential mixing
    do IH=1,NTIRD
      do LM1=2,ITDBRYD
        UI2(IH,LM1)=0.00
        VI2(IH,LM1)=0.00
      enddo
      SM1S(IH)=0.00
      FM1S(IH)=0.00
    enddo


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

      LDORHOEF = NPOL/=0  ! needed in RHOVAL

      call GAUNT2(WG,YRG,LMAX)

      call MADELUNG3D(LPOT,YRG,WG,ALAT, &
      RMAX,GMAX,BRAVAIS,RECBV, &
      LMXSPD,LASSLD,LPOTD,LMPOTD, &
      NMAXD,ISHLD, &
      LMPOT,CLEB,ICLEB,IEND, &
      NCLEBD,LOFLM,DFAC, &
      NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM)

      do LM = 1,LMPOTD
        CMOM(LM) = 0.0D0
        CMINST(LM) = 0.0D0
      end do

      CHRGNT = 0.0D0

      LRECRES1 = 8*43 + 16*(LMAXD+2)
      if (NPOL==0 .or. TEST('DOS     ')) then
        LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
      end if

      open (71,access='direct',recl=LRECRES1,file='results1', &
      form='unformatted')
      open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
      form='unformatted')

      if (TRC==1) then
        open (37,access='direct',recl=LRECTRC, &
              file='trnc.unf',form='unformatted')
      end if


!N ====================================================================
!     BEGIN do loop over atoms (NMPID-parallel)
!N ====================================================================

      do I1 = 1,NAEZ
        if(MYLRANK(LMPIC)==MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

!=======================================================================
! ccpl

          XCCPL = .false.
          ERESJIJ = .false.

          ! calculate exchange couplings only at last self-consistency step and when Jij=true
          if ((ITER==SCFSTEPS).and.JIJ) XCCPL = .true.

          if (XCCPL) then

            inquire(file='ERESJIJ',exist=ERESJIJ)

            call CLSJIJ(I1,NAEZ,RR,NR,RBASIS,RCUTJIJ,LMPIC,NSYMAT,ISYMINDEX, &
                        IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ, &
                        nrd, nxijd)

            do ISPIN = 1, NSPIN
              do XIJ = 1, NXIJ
                do LM1 = 1,LMMAXD
                  do LM2 = 1,LMMAXD
                    GMATXIJ(LM1,LM2,XIJ,ISPIN) = CZERO
                  enddo
                enddo
              enddo
            enddo

          endif

! ccpl
!=======================================================================


          read(66,rec=I1) VINS,VISP,ECORE
          if (TRC==1) read(37,rec=I1) NATRC,ATTRC,EZTRC,NUTRC,INTRC


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

! IME
          call OUTTIME(MYLRANK(1),'initialized .........',TIME_I,ITER)
! IME

! E ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! E ====================================================================

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

! E ====================================================================
            if (EMPIB==EPROC(IE)) then
! E ====================================================================


              do RF = 1,NREF

                call TREF(EZ(IE),VREF(RF),LMAX,RMTREF(RF), &
                          TREFLL(1,1,RF),DTREFLL(1,1,RF))

              end do

              call GREF(EZ(IE),ALAT,IEND1,NCLS,NAEZ, &
                        CLEB1C,RCLS,ATOM,CLS,ICLEB1C,LOFLM1C,NACLS, &
                        REFPOT, &
                        TREFLL(1,1,1),DTREFLL(1,1,1),GREFN,DGREFN, &
                        IE, &
                        LLY_G0TR,I1, &
                        LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE)

spinloop:     do ISPIN = 1,NSPIN

                call CALCTMAT(LDAU,NLDAU,ICST, &
                              NSRA,EZ(IE), &
                              DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                              VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                              IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                              TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                              LLDAU,WMLDAU, &
                              nspind, ncleb, ipand, irmd, irnsd)

                if(LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
                  call CALCDTMAT(LDAU,NLDAU,ICST, &
                                NSRA,EZ(IE),IE,NPNT1,NPNT2,NPNT3,PI,TK, &
                                DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN), &
                                VISP(1,ISPIN),ZAT(I1),IPAN(I1), &
                                IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1, &
                                DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN, &
                                LLDAU,WMLDAU, &
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


! PIN ==================================================================
!     BEGIN do loop over spins (SMPID-parallel)
! PIN===================================================================

!         initialize LLY_GRDT

                TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(IE,CLS(I1))
                LLY_GRDT(IE,ISPIN) = CZERO

                if (SMPID==1) then
                  MAPSPIN = 0
                  PRSPIN   = ISPIN
                else
                  MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
                  PRSPIN   = 1
                endif

!         true beginning of SMPID-parallel section

                if(SRANK(SMPIB,SMPIC)==MAPSPIN) then

                  NMESH = KMESH(IE)

                  if( MYLRANK(LMPIC)==0 ) then
                    write (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)')  &
                    ' ** IE = ',IE,' ENERGY =',EZ(IE), &
                    ' KMESH = ', NMESH,' ISPIN = ',ISPIN
                  end if

! <<>>

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
                  LMPIC,LCOMM,LSIZE, &
                  iemxd, lmpid * smpid * empid, nthrds, &
                  lmax, naclsd, nclsd, xdim, ydim, zdim, natbld, LLY, &
                  nxijd, nguessd, kpoibz, nrd, ekmd)

                endif

              end do spinloop                          ! ISPIN = 1,NSPIN

! PIN ==================================================================
!     END do loop over spins (SMPID-parallel)
! PIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! ccpl

              if (XCCPL) then

                call SREDGX( &
                ISPIN,NSPIN, &
                MYRANK,NROFNODES, &
                SMPIB,SMPIC,SMYRANK,SRANK, &
                GMATXIJ, &
                GXIJ_ALL)

                JSCAL = WEZ(IE)/DBLE(NSPIN)

                call XCCPLJIJ( &
                'R',I1,IE,JSCAL, &
                RXIJ,NXIJ,IXCP,RXCCLS, &
                GXIJ_ALL,DTIXIJ, &
                LMPIC,LCOMM, &
                MYRANK,EMPIC,EMYRANK, &
                JXCIJINT,ERESJIJ)

              end if

! ccpl
! End of Jij calculation
! =====================================================================

              call CPU_TIME(TIME_EX)
              ETIME(IE) = TIME_EX-TIME_E

! E ====================================================================
            endif
! E ====================================================================

! for preconditioning purposes calculate sparse indices combining IE.KPT
            EKM = EKM + NOFKS(KMESH(IE))

          end do                   ! IE = 1,IELAST

          if (ERESJIJ) close(75)   ! FIXME: is opened in XCCPLJIJ, but not closed?!

! E ====================================================================
!     END do loop over energies (EMPID-parallel) to be implemented
! E ====================================================================


!=======================================================================
!     "allreduce" information of 1 .. EMPID and 1 .. SMPID processors
!=======================================================================
          call SREDGM(NSPIN,IELAST, &
                      MYRANK, &
                      SMPIC,SMYRANK,SRANK, &
                      EMPIC,EMYRANK,ERANK,EPROC, &
                      GMATN,LLY_GRDT, &
                      GMATN_ALL,LLY_GRDT_ALL)
!=======================================================================
!=======================================================================

! IME
          call OUTTIME(MYLRANK(1),'G obtained ..........',TIME_I,ITER)
! IME

!=======================================================================
!     output of Jij's .. calling xccpljij with flag 'F'
!=======================================================================
          if (XCCPL) then

            call XCCPLJIJ('F',I1,IE,JSCAL, &
                          RXIJ,NXIJ,IXCP,RXCCLS, &
                          GXIJ_ALL,DTIXIJ, &
                          LMPIC,LCOMM, &
                          MYRANK,EMPIC,EMYRANK, &
                          JXCIJINT,ERESJIJ)
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
!     in case of IGUESS and EMPID > 1 preconditioning arrays might
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

          endif
!=======================================================================
!=======================================================================


!----------------------------------------------------------------------
! BEGIN only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
          if (LMPIC==1) then

            if (LLY==1) then
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
                          DMATLDAU)

              if (LLY==1) then
                do IE=1,IELAST
                  do L=0,LMAXD1
                    DEN(L,IE,ISPIN)=DEN(L,IE,ISPIN)*RNORM(IE,ISPIN)
                  end do
                end do
              end if

              EBOT = E1

              call RHOCORE(EBOT,NSRA,ISPIN,NSPIN,I1, &
                           DRDI(1,I1),R(1,I1),VISP(1,ISPIN), &
                           A(I1),B(I1),ZAT(I1), &
                           IRCUT(0,I1),RHOCAT,QC, &
                           ECORE(1,ISPIN),NCORE(IPOT),LCORE(1,IPOT))

            end do

! SPIN ==================================================================
!      END do loop over spins
! SPIN ===================================================================

            do ISPIN = 1,NSPIN
              do L = 0,LMAXD1
                DENEF = DENEF - 2.0D0 * &
                DIMAG(DEN(L,IELAST,ISPIN))/PI/DBLE(NSPIN)
              end do
            end do

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


! DAU

            EULDAU = 0.0D0
            EDCLDAU = 0.0D0

            if (LDAU.and.NLDAU>=1) then

              call LDAUWMAT(I1,NSPIN,ITER,MIXING,DMATLDAU,NLDAU,LLDAU, &
                            ULDAU,JLDAU,UMLDAU,WMLDAU,EULDAU,EDCLDAU, &
                            lmaxd)

            endif

! DAU


! ----------------------------------------------------------------------

! -->   determine total charge density expanded in spherical harmonics

! -------------------------------------------------------------- density

            call RHOTOTB(NAEZ,I1,NSPIN,RHO2NS,RHOCAT, &
                         ZAT,DRDI,IRCUT, &
                         LPOT,NFU,LLMSP(1,ICELL),THETAS,ICELL,IPAN, &
                         CATOM)

            CHRGNT = CHRGNT + CATOM(1) - ZAT(I1)

            ! write to 'results1'
            if (NPOL==0 .or. TEST('DOS     ')) then
              write(71,rec=I1) QC,CATOM,CHARGE,ECORE,DEN  ! write density of states (DEN) only when certain options set
            else
              write(71,rec=I1) QC,CATOM,CHARGE,ECORE
            end if

          endif
!----------------------------------------------------------------------
! END L-MPI: only processes with LMPIC = 1 are working
!----------------------------------------------------------------------
        end if
      end do

!N ====================================================================
!     END do loop over atoms (NMPID-parallel)
!N ====================================================================

      if (TRC==1) close(37)
      close(66)
      close(71)

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

! **************************  ITERATION NUMBER  ************************
        open (28,file='not.converged',form='formatted')
        read (28,'(1P,4D17.10)') EFOLD,VBC
        close (28)
! **********************************************************************

! ----------------------------------------------------------------------

! --> determine new Fermi level due to valence charge up to
!     old Fermi level E2 and density of states DENEF

        E2SHIFT = CHRGNT/DENEF
        E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT) !FIXME: hardcoded
        EFOLD = E2

        if (ISHIFT < 2) E2 = E2 - E2SHIFT

        if( MYLRANK(LMPIC) == 0 ) then
          write (6,fmt=9020) EFOLD,E2SHIFT

! --> divided by NAEZ because the weight of each atom has been already
!     taken into account in 1c

          write (6,fmt=9030) E2,DENEF/DBLE(NAEZ)
          write(6,'(79(1H+),/)')
        end if

! ----------------------------------------------------------------------
        DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIN)
! ----------------------------------------------------------------------

        open (66,access='direct',recl=LRECPOT*2,file='vpotnew', &
        form='unformatted')
        open (72,access='direct',recl=LRECRES2,file='results2', &
        form='unformatted')

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
              if (LCORDENS) then

                do LM = 1,LMPOT
                  call DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1, &
                  RHO2NS(1,LM,ISPIN),1)
                end do

              end if
! ----------------------------------------------------------------------
            end do

            call RHOMOM(CMOM,CMINST,LPOT,I1,RHO2NS, &
            R,DRDI,IRWS,IRCUT,IPAN,ICELL,ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS,LMSP(1,ICELL))

            call OUTTIME(MYLRANK(1),'RHOMOM ......',TIME_I,ITER)

! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================

            call VINTRAS(LPOT,NSPIN,I1,RHO2NS,VONS, &
            R,DRDI,IRWS,IRCUT,IPAN,ICELL,ILM,IFUNM(1,ICELL),IMAXSH,GSH, &
            THETAS,LMSP(1,ICELL))

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
            LGROUP,LCOMM,LSIZE)

            call OUTTIME(MYLRANK(1),'VMADELBLK ......',TIME_I,ITER)

! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

            if (KFORCE==1 .and. ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCEH(CMOM,FLM,LPOT,NSPIN,I1,RHO2NS,VONS, &
              R,DRDI,IMT,ZAT,irmd)
              call FORCE(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              DRDI,IMT,naez,irmd)
! ---------------------------------------------------------------------
            end if

! Force Calculation stops here look after VXCDRV

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES

            if (KTE==1) then
              call ESPCB(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE)

              call EPOTINB(EPOTIN,NSPIN,I1,RHO2NS,VISP,R,DRDI, &
              IRMIN,IRWS,LPOT,VINS,IRCUT,IPAN,ZAT, &
              irmd, irnsd, ipand)

              call ECOUB(CMOM,ECOU,LPOT,NSPIN,I1,RHO2NS, &
              VONS,ZAT,R, &
              DRDI,IRWS,KVMAD,IRCUT,IPAN,IMAXSH,IFUNM(1,ICELL), &
              ILM,ICELL,GSH,THETAS,LMSP(1,ICELL), &
              irmd, irid, nfund, ipand, ngshd)

            end if
            call OUTTIME(MYLRANK(1),'KTE ......',TIME_I,ITER)
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

! =====================================================================
            call VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,I1,RHO2NS, &
            VONS,R,DRDI,A, &
            IRWS,IRCUT,IPAN,ICELL,GSH,ILM,IMAXSH,IFUNM(1,ICELL), &
            THETAS,LMSP(1,ICELL))

            call OUTTIME(MYLRANK(1),'VXCDRV ......',TIME_I,ITER)
! =====================================================================

! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

! Force calculation continues here

            if (KFORCE==1.and.ITER==SCFSTEPS) then
! ---------------------------------------------------------------------
              call FORCXC(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R, &
              ALAT,DRDI,IMT,ZAT, &
              LMPIC,MYLRANK, &
              LGROUP,LCOMM,LSIZE, &
              naez, irmd, lmpid*smpid*empid)
! ---------------------------------------------------------------------
            end if

            write(72,rec=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX, &
            EULDAU,EDCLDAU

! Force calculation ends
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! =====================================================================
            VAV0 = 0.0D0
            VOL0 = 0.0D0
            call MTZERO(LMPOT,I1,NSPIN,VONS,ZAT,R,DRDI,IMT,IRCUT, &
            IPAN,ICELL,LMSP(1,ICELL),IFUNM(1,ICELL), &
            THETAS,IRWS,VAV0,VOL0)
            call OUTTIME(MYLRANK(1),'MTZERO ......',TIME_I,ITER)
! =====================================================================
! ============================= ENERGY and FORCES =====================
! =====================================================================
          end if
        end do
        call OUTTIME(MYLRANK(1),'calculated pot ......',TIME_I,ITER)
! =====================================================================
! ======= I1 = 1,NAEZ ================================================
! =====================================================================


!****************************************************** MPI COLLECT DATA
        WORK1(1) = VAV0
        WORK1(2) = VOL0
        call MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM, &
        LCOMM(LMPIC),IERR)
        VAV0 = WORK2(1)
        VOL0 = WORK2(2)
!****************************************************** MPI COLLECT DATA


        VBC(1) = -VAV0/VOL0
        if (ISHIFT>0) VBC(1) = VBC(1) + E2SHIFT
        VBC(2) = VBC(1)

        if(MYRANK==0) then
          write (6,fmt=9103) VOL0,VAV0,VBC(1)
          write(6,'(79(1H=),/)')
        end if

9103    format ('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)
! =====================================================================

! ---------------------------------------------------------------------

! -->   shift potential to muffin tin zero and
!       convolute potential with shape function for next iteration

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
            ITER,RFPI,FPI,IPF, &
            MIX, &
            FCM,IRC,IRMIN,R,DRDI,VONS, &
            VISP,VINS, &
            naez, irmd, irnsd)

            I1BRYD=I1
          end if
        end do

! -->  potential mixing procedures: Broyden or Andersen updating schemes
        if (IMIX>=3) then
          call BRYDBM(VISP,VONS,VINS, &
          LMPOT,R,DRDI,MIX, &
          IRC,IRMIN,NSPIN,I1BRYD,NAEZ, &
          IMIX,IPF,ITER, &
          UI2,VI2,WIT,SM1S,FM1S, &
          LMPIC,MYLRANK, &
          LGROUP,LCOMM,LSIZE, &
          itdbryd, irmd, irnsd, nspind, &
          LMPID * SMPID * EMPID)
        endif

!----------------------------------------------------------------------
! -->    reset to start new iteration
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        do I1 = 1,NAEZ
          if(MYLRANK(LMPIC)== &
          MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then

            !initialise VINS
            do ISPIN = 1,2
              do LM = 1, LMPOTD
                do J = IRMIND, IRMD
                  VINS(J,LM,ISPIN) = 0.0D0
                enddo
              enddo
            enddo

            !initialise VISP
            do ISPIN = 1,2
              do J = 1, IRMD
                VISP(J,ISPIN) = 0.0D0
              enddo
            enddo

            do ISPIN = 1,NSPIN
              IPOT = (I1-1)*NSPIN + ISPIN

              IRC1 = IRC(I1)
              call DCOPY(IRC1,VONS(1,1,ISPIN),1,VISP(1,ISPIN),1)

              if (LPOT>0) then
                IRMIN1 = IRMIN(I1)
                do LM = 2,LMPOT
                  do J = IRMIN1,IRC1
                    VINS(J,LM,ISPIN) = VONS(J,LM,ISPIN)
                  end do
                end do
              end if
            enddo
! ----------------------------------------------------- output_potential
            write(66,rec=I1) VINS,VISP,ECORE
! ----------------------------------------------------- output_potential
          end if
        end do

        close(72)
! =====================================================================
! ============================= POTENTIAL MIXING OUTPUT ===============
! =====================================================================
! ====== write RMS convergency data - not parallelized, written by
! MYRANK=0   (RMSOUT) =================================================
        call RMSOUT(RMSAVQ,RMSAVM,ITER,E2,EFOLD, &
        SCFSTEPS,VBC,QBOUND,NSPIN,NAEZ, &
        KXC,LPOT,A,B,IRC, &
        VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT, &
        ECORE,LCORE,NCORE,ZAT,ITITLE, &
        LMPIC,MYLRANK, &
        LGROUP,LCOMM,LSIZE)

9020    format ('                old', &
        ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
9030    format ('                new', &
        ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)



! Wait here in order to guarantee regular and non-errorneous output
! in RESULTS

        call MPI_BARRIER(LCOMM(LMPIC),IERR)


! -----------------------------------------------------------------
! L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------
        if(MYLRANK(LMPIC)==0) then

          call RESULTS(LRECRES2,IELAST,ITER,LMAX,NAEZ,NPOL, &
          NSPIN,KPRE,KTE,LPOT,E1,E2,TK,EFERMI, &
          ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,LDAU)

! --> update energy contour

!          if(TEST('fix-EF  ')) rewind(41)
!          if(TEST('fix-EF  ')) read(41,*) E2

          call EMESHT(EZ,DEZ,IELAST,E1,E2,E2,TK, &
          NPOL,NPNT1,NPNT2,NPNT3,IEMXD)

          do IE = 1,IELAST
            WEZ(IE) = -2.D0/PI*DEZ(IE)
          end do

          write(6,'(79(1H=))')

! .. get info on MYACTVRANK of this processor: to be used in
!    subsequent reduce-commands
          MYBCRANK = MYACTVRANK
! ..
        endif
! -----------------------------------------------------------------
! L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
! -----------------------------------------------------------------

! ..
      endif
! -----------------------------------------------------------------
! L-MPI: only processes with LMPIC = 1 are working here
! -----------------------------------------------------------------
      close(66)  ! close 'vpotnew'
! -----------------------------------------------------------------

!      CALL MPI_BARRIER(ACTVCOMM,IERR)

      call MPI_ALLREDUCE(MYBCRANK,BCRANK,1,MPI_INTEGER,MPI_MAX, &
      ACTVCOMM,IERR)

      call MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION, &
      BCRANK,ACTVCOMM,IERR)

      call MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM, &
      ACTVCOMM,IERR)

!      CALL MPI_BARRIER(ACTVCOMM,IERR)

      if(MYLRANK(1)==0) then

        open (67,file='energy_mesh',form='unformatted')
        write (67) IELAST,EZ,WEZ,E1,E2
        write (67) NPOL,TK,NPNT1,NPNT2,NPNT3
        close (67)

        write(6,'(79(1H=))')
        write(6,'(19X,A,I3,A,I10)') '       ITERATION : ', &
        ITER,' SUM of QMR ',NOITER_ALL
        write(6,'(79(1H=),/)')
        call SYSTEM_CLOCK(SYSTEM_F,RATETIME,MAXTIME)
        WALLCLOCK_F=MPI_WTIME()
        write(6,'(79(1H=))')
        write(6,*) 'Wallclock, System and CPU-time compared:'
        write(6,*) 'MPI_WTIME=',(WALLCLOCK_F-WALLCLOCK_I),'sec'
        write(6,*) 'SYSTEM_CLOCK=',(SYSTEM_F-SYSTEM_I)/100,'sec'
        call OUTTIME(MYLRANK(1),'end .................', &
        TIME_I,ITER)
        write(6,'(79(1H=))')
        write(2,'(79(1H=))')
        call OUTTIME(MYLRANK(1),'finished in .........', &
        TIME_S,ITER)
        write(2,'(79(1H=))')
        write(6,'(79(1H=),/)')
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
    close(IPF)                       ! RMS
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

end program MAIN2
        
