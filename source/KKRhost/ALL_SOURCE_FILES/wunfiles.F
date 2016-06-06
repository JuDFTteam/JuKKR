      module mod_wunfiles
        
      implicit none
      
         
   !define type that replace wunfiles here, later define bcast routine
      type :: type_params
   
      integer :: Nscalars = 124

!     .. Scalars
      INTEGER IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,NEMBD1,       &
     &        LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,  &
     &        LM2D,LMAXD1,NRD,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,  &
     &        NSPIND,NSPINDD,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD,  &
     &        IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,LLY,          &
     &        NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,NOFGIJ,          &
     &        NCLS,IEND,ICC,IGF,NLBASIS,NRBASIS,NCPA,ITCPAMAX,ICST,      &
     &        KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC,               &
     &        INTERVX,INTERVY,INTERVZ,NMVECMAX,                          &
     &        LPOT,LMPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,                  &
     &        KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE,                        &
     &        IDOLDAU,ITRUNLDAU,NTLDAU,I1,ITAB,KPOIBZ,                   &
     &        NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,                  &
     &        ITMPDIR,ILTMP,NCHEB,NTOTD,NCHEBD
      DOUBLE PRECISION EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG,          &
     &                 E1,E2,TK,EFERMI,ALAT,CPATOL,LASTERR,              &
     &                 MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF,EFOLD,CHRGOLD &
     &                ,DENEF,CHRGSEMICORE
      DOUBLE COMPLEX DELTAE
      LOGICAL LINTERFACE,LRHOSYM
      CHARACTER*10 SOLVER
      CHARACTER*80 TMPDIR
      
!     .. Arrays
      DOUBLE COMPLEX, allocatable :: EZ(:),WEZ(:),DROTQ(:,:,:),          &
     &               DSYMLL(:,:,:),                                      &
     &               LEFTTINVLL(:,:,:,:,:),                              &
     &               RIGHTTINVLL(:,:,:,:,:),                             &
     &               CREL(:,:),RC(:,:),                                  &
     &               RREL(:,:),SRREL(:,:,:),                             &
     &               PHILDAU(:,:)
      DOUBLE PRECISION, allocatable :: VINS(:,:,:),VISP(:,:),            &
     &                 VBC(:),VOLBZ(:),                                  &
     &                 VTREL(:,:),BTREL(:,:),                            &
     &                 SOCSCALE(:),                                      &
     &                 DRDIREL(:,:),R2DRDIREL(:,:),                      &
     &                 RMREL(:,:),CMOMHOST(:,:),                         &
     &                 ECORE(:,:),QMTET(:),QMPHI(:),                     &
     &                 QMPHITAB(:,:),QMTETTAB(:,:),                      &
     &                 QMGAMTAB(:,:),ZAT(:),R(:,:),                      &
     &                 DRDI(:,:),RMTREF(:),                              &
     &                 VREF(:),CLEB(:,:),RCLS(:,:,:),                    &
     &                 SOCSCL(:,:),CSCL(:,:),                            &
     &                 RBASIS(:,:),RR(:,:),CONC(:),                      &
     &                 RROT(:,:,:),RATOM(:,:),                           &
     &                 A(:),B(:) ,THETAS(:,:,:),                         &
     &                 RMT(:),RMTNEW(:),RWS(:),GSH(:),                   &
     &                 EREFLDAU(:),UEFF(:),JEFF(:),                      &
     &                 ULDAU(:,:,:,:,:),                                 &
     &                 WLDAU(:,:,:,:),                                   &
     &                 RPAN_INTERVALL(:,:),                              &
     &                 RNEW(:,:),BZKP(:,:,:),VOLCUB(:,:),                &
     &                 THETASNEW(:,:,:),MVEVI(:,:,:),MVEVIEF(:,:,:),     &
     &                 RHO2NS(:,:,:,:),R2NEF(:,:,:,:),RHOC(:,:),
     &                 DENEFAT(:),ESPV(:,:),EDC(:),EU(:),RHOORB(:,:),
     &                 ECOREREL(:,:),theta(:),phi(:),RCLSIMP(:,:)
      INTEGER, allocatable :: LOPT(:),ITLDAU(:),                         &
     &        IRSHIFT(:),JWSREL(:),ZREL(:),                              &
     &        LCORE(:,:),NCORE(:),IPAN(:),                               &
     &        IRCUT(:,:),JEND(:,:,:),                                    &
     &        ICLEB(:,:),ATOM(:,:),CLS(:),                               &
     &        NACLS(:),NOFKS(:),                                         &
     &        LOFLM(:),EZOA(:,:),                                        &
     &        KAOEZ(:,:),IQAT(:),                                        &
     &        ICPA(:),NOQ(:),KMESH(:),                                   &
     &        NSHELL(:),NSH1(:),NSH2(:),                                 &
     &        IJTABCALC(:),IJTABCALC_I(:),IJTABSYM(:),IJTABSH(:),        &
     &        ISH(:,:),JSH(:,:),IQCALC(:),                               &
     &        ICHECK(:,:),ATOMIMP(:),REFPOT(:),                          &
     &        IRREL(:,:,:),NRREL(:,:),IFUNM1(:,:),                       &
     &        ITITLE(:,:),LMSP1(:,:),NTCELL(:),                          &
     &        IXIPOL(:),IRNS(:),IFUNM(:,:),                              &
     &        LLMSP(:,:),LMSP(:,:),IMT(:),                               &
     &        IRC(:),IRMIN(:),IRWS(:),NFU(:),                            &
     &        HOSTIMP(:),ILM(:,:),IMAXSH(:),                             &
     &        NPAN_LOG(:),NPAN_EQ(:),                                    &
     &        NPAN_TOT(:),                                               &
     &        IPAN_INTERVALL(:,:),NKCORE(:,:),KAPCORE(:,:)
      LOGICAL, allocatable :: SYMUNITARY(:),VACFLAG(:)
      CHARACTER*24, allocatable :: TXC(:)
      CHARACTER*8,  allocatable :: OPTC(:),TESTC(:)
      
      end type type_params


      type (type_params), save :: t_params

      
      
      contains


      SUBROUTINE WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ,
     &                    EFERMI,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,
     &                    IESEMICORE,TKSEMI,EBOTSEMI,EMUSEMI,
     &                    FSEMICORE,VINS,VISP,VBC,VTREL,BTREL,RMREL,
     &                    DRDIREL,R2DRDIREL,ZREL,JWSREL,IRSHIFT,
     &                    ITSCF,SCFSTEPS,CMOMHOST,ECORE,LCORE,NCORE,
     &                    QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,DROTQ,
     &                    NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,
     &                    NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI,
     &                    REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,
     &                    ATOM,CLS,RCLS,NACLS,LOFLM,SOLVER,SOCSCL,CSCL,
     &                    ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,
     &                    CPATOL,RBASIS,RR,EZOA,NSHELL,NSH1,NSH2,
     &                   IJTABCALC,IJTABCALC_I,ISH,JSH,IJTABSYM,IJTABSH,
     &                   NOFGIJ,NQCALC,IQCALC,KMROT,KAOEZ,IQAT,NOQ,CONC,
     &                    KMESH,MAXMESH,NSYMAT,SYMUNITARY,RROT,
     &                    DSYMLL,INVMOD,ICHECK,
     &                    NATOMIMP,RATOM,ATOMIMP,
     &                    RC,CREL,RREL,SRREL,NRREL,IRREL,
     &                    LEFTTINVLL,RIGHTTINVLL,VACFLAG,
     &                    A,B,IFUNM,IFUNM1,INTERVX,INTERVY,INTERVZ,
     &                    ITITLE,LMSP1,NTCELL,THETAS,
     &                    LPOT,LMPOT,NRIGHT,NLEFT,LINTERFACE,
     &                    IMIX,MIXING,QBOUND,FCM,ITDBRY,IRNS,KPRE,
     &                    KSHAPE,KTE,KVMAD,KXC,LAMBDA_XC,TXC,ISHIFT,
     &                    IXIPOL,LRHOSYM,KFORCE,LMSP,LLMSP,RMT,RMTNEW,
     &                    RWS,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,GSH,ILM,
     &                    IMAXSH,IDOLDAU,ITRUNLDAU,NTLDAU,LOPT,ITLDAU,
     &                    UEFF,JEFF,EREFLDAU,ULDAU,WLDAU,PHILDAU,
     &                    IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,
     &                    NEMBD1,LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,
     &                    NCLEB,NACLSD,NCLSD,LM2D,LMAXD1,MMAXD,NRD,
     &                    NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,
     &                    NSPIND,IRID,NFUND,NCELLD,LMXSPD,
     &                    NGSHD,KREL,NTOTD,NCHEBD,NPAN_LOG,NPAN_EQ,
     &                    NCHEB,R_LOG,NPAN_TOT,RNEW,RPAN_INTERVALL,
     &                    IPAN_INTERVALL,NSPINDD,THETASNEW,SOCSCALE,
     &                    TOLRDIF,LLY,DELTAE,RCLSIMP)
C **********************************************************************
C *                                                                    *
C *  This subroutine is part of the MAIN0 program in the tbkkr package *
C *  It writes out different unformatted files meant to provide the    *
C *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
C *  during an SCF cycle                                               *
C *  v.popescu, munich 2004                                            *
C *                                                                    *
C **********************************************************************

      use mod_types, only: t_inc, t_tgmat, t_lloyd, t_cpa

      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,NEMBD1,
     &        LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,
     &        LM2D,LMAXD1,NRD,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,
     &        NSPIND,NSPINDD,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD
C     .. nembd2 = naezd+nembd, lmaxd1=lmaxd+1, naezdpd=naezd/nprincd)
      INTEGER IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,LLY
      INTEGER NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,NOFGIJ,
     &        NCLS,ICST,IEND,ICC,IGF,NLBASIS,NRBASIS,NCPA,ITCPAMAX
      INTEGER KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC
      INTEGER INTERVX,INTERVY,INTERVZ
      INTEGER LPOT,LMPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,
     &        KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU
      INTEGER NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE
      DOUBLE PRECISION EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG
      DOUBLE PRECISION E1,E2,TK,EFERMI,ALAT,CPATOL
      DOUBLE PRECISION MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF
      LOGICAL LINTERFACE,LRHOSYM
      CHARACTER*10 SOLVER
C     ..
C     .. Array arguments
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD),DROTQ(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),
     &               LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &           RREL(LMMAXD,LMMAXD),SRREL(2,2,LMMAXD)
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
      DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD)
      DOUBLE PRECISION VBC(2)
      DOUBLE PRECISION VTREL(IRMD,NATYPD),BTREL(IRMD,NATYPD)
      DOUBLE PRECISION SOCSCALE(NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD,NATYPD),R2DRDIREL(IRMD,NATYPD),
     &                 RMREL(IRMD,NATYPD),CMOMHOST(LMPOTD,NEMBD1)
      DOUBLE PRECISION ECORE(20,NPOTD),QMTET(NAEZD),QMPHI(NAEZD)
      DOUBLE PRECISION QMPHITAB(NAEZD,3),QMTETTAB(NAEZD,3),
     &                 QMGAMTAB(NAEZD,3),ZAT(NATYPD),R(IRMD,NATYPD),
     &                 DRDI(IRMD,NATYPD),RMTREF(NREFD),
     &                 VREF(NREFD),CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD)
      DOUBLE PRECISION SOCSCL(LMAXD1,NATYPD),CSCL(LMAXD1,NATYPD)
      DOUBLE PRECISION RBASIS(3,NEMBD2),RR(3,0:NRD),CONC(NATYPD)
      DOUBLE PRECISION RROT(48,3,NSHELD),RATOM(3,NSHELD)
      DOUBLE PRECISION A(NATYPD),B(NATYPD) ,THETAS(IRID,NFUND,NCELLD)
      DOUBLE PRECISION RMT(NATYPD),RMTNEW(NATYPD),RWS(NATYPD),GSH(NGSHD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
C      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C     ..
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)
      INTEGER LCORE(20,NPOTD),NCORE(NPOTD),IPAN(NATYPD),
     &        IRCUT(0:IPAND,NATYPD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     &        ICLEB(NCLEB,4),ATOM(NACLSD,NEMBD2),CLS(NEMBD2),
     &        NACLS(NCLSD)
      INTEGER LOFLM(LM2D),EZOA(NACLSD,NEMBD2),
     &        KAOEZ(NATYPD,NEMBD2),IQAT(NATYPD),
     &        ICPA(NAEZD),NOQ(NAEZD),KMESH(IEMXD)
      INTEGER NSHELL(0:NSHELD),NSH1(NSHELD),NSH2(NSHELD),
     &        IJTABCALC(NOFGIJD),IJTABCALC_I(NOFGIJD),IJTABSYM(NOFGIJD),
     &        IJTABSH(NOFGIJD)
      INTEGER ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD),IQCALC(NAEZD)
      INTEGER ICHECK(NAEZDPD,NAEZDPD),ATOMIMP(NATOMIMPD),REFPOT(NEMBD2)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD),IFUNM1(LMXSPD,NATYPD),
     &        ITITLE(20,NPOTD),LMSP1(LMXSPD,NATYPD),NTCELL(NATYPD)
      INTEGER IXIPOL(NATYPD),IRNS(NATYPD),IFUNM(NATYPD,LMXSPD)
      INTEGER LLMSP(NATYPD,NFUND),LMSP(NATYPD,LMXSPD),IMT(NATYPD),
     &        IRC(NATYPD),IRMIN(NATYPD),IRWS(NATYPD),NFU(NATYPD),
     &        HOSTIMP(0:NATYPD),ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      INTEGER NTOTD,NCHEBD,NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),
     +        NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     &                 RNEW(NTOTD*(NCHEBD+1),NATYPD),
     &                 THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
      LOGICAL SYMUNITARY(NSYMAXD),VACFLAG(2)
      CHARACTER*24 TXC(5)
      CHARACTER*80 TMPDIR
      INTEGER ITMPDIR,ILTMP      
C     ..
C     .. Local scalars
      INTEGER I1,I2
      integer ic, naclsmax, nqdos  !variables for t_inc filling
      double precision phi(natyp),theta(natyp)
C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
C
      ITMPDIR=0
      ILTMP=0
      
      ! put information about scf steps into t_inc
      t_inc%i_iteration = ITSCF
      t_inc%N_iteration = SCFSTEPS

C ------------------------------------------------------------- itermdir
C
      IF (OPT('ITERMDIR')) THEN
         I1 = 0
         E1 = 0D0
         t_params%QMTET    = QMTET
         t_params%QMPHI    = QMPHI
         t_params%QMPHITAB = QMPHITAB
         t_params%QMTETTAB = QMTETTAB
         t_params%QMGAMTAB = QMGAMTAB
         t_params%I1       = I1
         t_params%E1       = E1
         t_params%DROTQ    = DROTQ
      END IF
C ---------------------------------------------------------------- lda+u
C                                               data in this file change
C
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
C
C ---------------------------------------------------------------- nonco_angle file
C
      call read_angles(t_params,NATYP,THETA,PHI)
C
C --------------------------------------------------------------- fill t_inc type 
C                                  meant for simpler passing of basic parameter to other routines

      ! find maximal cluster information (taken from main1b)
      NACLSMAX = 1
      DO IC = 1,NCLS
         IF (NACLS(IC).GT.NACLSMAX) NACLSMAX = NACLS(IC)
      ENDDO
      
      ! find NQDOS (taken from main1b)
      IF (OPT('qdos    ')) THEN
         OPEN(67,FILE='qvec.dat')
         READ(67,*) NQDOS
         CLOSE(67)
      ELSE
         NQDOS = 1
      END IF
   
      !fill t_inc
      t_inc%LMMAXD = LMMAXD
      t_inc%NSPIN  = NSPIN
      t_inc%IELAST = IELAST
      t_inc%NQDOS  = NQDOS
      t_inc%NATYP  = NATYP
      t_inc%LMGF0D = (LMAXD+1)**2  ! see main1b
      t_inc%NCLSD  = NCLS
      t_inc%NACLSD = NACLSMAX
      t_inc%NSHELL0 = NSHELL(0)
      IF(OPT('NEWSOSOL')) t_inc%NEWSOSOL = .true.
      IF(OPT('deci-out')) t_inc%deci_out = .true.
      
      !set logical switches in t_tgmat which control if tmat, gmat and gref are written to files or stored in memory
      if(TEST('tmatfile')) t_tgmat%tmat_to_file = .true.
      if(TEST('gmatfile')) t_tgmat%gmat_to_file = .true.
      if(TEST('greffile')) t_tgmat%gref_to_file = .true.
      if(TEST('projfile')) t_cpa%dmatproj_to_file = .true.

      !some special run options:
      if(OPT('KKRFLEX '))  t_tgmat%tmat_to_file = .true.  ! for KKRFLEX option tmat must be written to file 
      if(OPT('qdos    '))  t_tgmat%gmat_to_file = .true.  ! for qdos write gmat to file since it scales with NQDOS and can become huge
      
      !set logical switches in t_lloyd which control if files are written to files or stored in memory
      if(TEST('tmatfile').or.TEST('llyfiles')) 
     &                    t_lloyd%dtmat_to_file = .true.
      if(TEST('wrtdtmat').or.TEST('llyfiles'))
     &                    t_lloyd%tralpha_to_file = .true.
      if(TEST('wrtcdos ').or.TEST('llyfiles')) 
     &                    t_lloyd%cdos_diff_lly_to_file = .true.
      if(TEST('wrtdgref').or.TEST('llyfiles')) 
     &                    t_lloyd%dgref_to_file = .true.
      if(TEST('wrtgotr ').or.TEST('llyfiles')) 
     &                    t_lloyd%g0tr_to_file = .true.
     
      ! set verbosity level in t_inc%i_write = 0,1,2 for default, verbose1, verbose2
      t_inc%i_write = 0 !default: write only output.000.txt and reset file after each iteration
      if(TEST('verbose1')) t_inc%i_write = 1 !write on all processors but only the latest iteration
      if(TEST('verbose2')) t_inc%i_write = 2 !write everything
      ! and t_inc_i_time for timing writeout
      t_inc%i_time = 1  !default: only timings from master, all iterations
      if(TEST('timings0')) t_inc%i_time = 0  !only timings from master, only the last iteration
      if(TEST('timings2')) t_inc%i_time = 2  !all timing files, all iterations
      
      
      ! all parameters are stored in t_params fomr mod_wunfiles
      ! first fill scalar values
      call fill_t_params_scalars(
     &        IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,NEMBD1,
     &        LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,
     &        LM2D,LMAXD1,NRD,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,
     &        NSPIND,NSPINDD,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD,
     &        IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,LLY,
     &        NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,NOFGIJ,
     &        NCLS,ICST,IEND,ICC,IGF,NLBASIS,NRBASIS,NCPA,ITCPAMAX,
     &        KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC,INTERVX,
     &        INTERVY,INTERVZ,LPOT,LMPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,
     &        KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE,IDOLDAU,ITRUNLDAU,
     &        NTLDAU,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,
     &        EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG,E1,E2,TK,EFERMI,
     &        ALAT,CPATOL,MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF,
     &        LINTERFACE,LRHOSYM,SOLVER,TMPDIR,ITMPDIR,ILTMP,
     &        NTOTD,NCHEBD,NCHEB,DELTAE,0.d0,0.d0,t_params)
     
      ! initialize allocatable arrays
      call init_t_params(t_params)
           
      ! now fill arrays that have just been allocated
      call fill_t_params_arrays(t_params,IEMXD,LMMAXD,NAEZD,
     &       NSYMAXD,NEMBD1,NSPINDD,IRMIND,IRMD,
     &       LMPOTD,NSPOTD,NPOTD,NATYPD,NRD,NEMBD2,NREFD,NCLEB,NCLSD,
     &       NACLSD,NSHELD,NGSHD,NFUND,IRID,NCELLD,MMAXD,LM2D,LMXSPD,
     &       LMAXD1,NSPIND,NTOTD,NCHEBD,IPAND,LMAXD,NOFGIJD,NAEZDPD,
     &       NATOMIMPD,
     &       EZ,WEZ,DROTQ,DSYMLL,LEFTTINVLL,
     &       RIGHTTINVLL,CREL,RC,RREL,SRREL,PHILDAU,VINS,VISP,VBC,VTREL,
     &       BTREL,SOCSCALE,DRDIREL,R2DRDIREL,RMREL,CMOMHOST,ECORE,QMTET
     &      ,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ZAT,R,DRDI,RMTREF,VREF,
     &       CLEB,RCLS,SOCSCL,CSCL,RBASIS,RR,CONC,RROT,RATOM,A,B,THETAS,
     &       RMT,RMTNEW,RWS,GSH,EREFLDAU,UEFF,JEFF,ULDAU,WLDAU,
     &       RPAN_INTERVALL,RNEW,THETASNEW,LOPT,ITLDAU,IRSHIFT,JWSREL,
     &       ZREL,LCORE,NCORE,IPAN,IRCUT,JEND,ICLEB,ATOM,CLS,NACLS,LOFLM
     &      ,EZOA,KAOEZ,IQAT,ICPA,NOQ,KMESH,NSHELL,NSH1,NSH2,IJTABCALC,
     &       IJTABCALC_I,
     &       IJTABSYM,IJTABSH,ISH,JSH,IQCALC,ICHECK,ATOMIMP,REFPOT,IRREL
     &      ,NRREL,IFUNM1,ITITLE,LMSP1,NTCELL,IXIPOL,IRNS,IFUNM,LLMSP,
     &       LMSP,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,ILM,IMAXSH,NPAN_LOG,
     &       NPAN_EQ,NPAN_TOT,IPAN_INTERVALL,SYMUNITARY,VACFLAG,TXC,
     &       RCLSIMP)
     
      ! save information about the energy mesh
      call save_emesh(IELAST,EZ,WEZ,E1,E2,IESEMICORE,FSEMICORE,NPOL,TK,
     &                NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,
     &                NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD,t_params)

          
C ======================================================================
      END subroutine




      subroutine init_t_params(t_params)
      ! allocate arrays, has to be done after bcast t_params_scalars for myrank<>master
      ! otherwise are the parameters not set
      implicit none
      
      type(type_params), intent(inout) :: t_params
      
      ALLOCATE(t_params%EZ(t_params%IEMXD),   !DOUBLE COMPLEX
     &t_params%WEZ(t_params%IEMXD),
     &t_params%DROTQ(t_params%LMMAXD,t_params%LMMAXD,t_params%NAEZD),
     &t_params%DSYMLL(t_params%LMMAXD,t_params%LMMAXD,t_params%NSYMAXD),
     &t_params%LEFTTINVLL(t_params%LMMAXD,t_params%LMMAXD,
     &   t_params%NEMBD1,t_params%NSPINDD,t_params%IEMXD),
     &t_params%RIGHTTINVLL(t_params%LMMAXD,t_params%LMMAXD,
     &   t_params%NEMBD1,t_params%NSPINDD,t_params%IEMXD),
     &t_params%CREL(t_params%LMMAXD,t_params%LMMAXD),
     &t_params%RC(t_params%LMMAXD,t_params%LMMAXD),
     &t_params%RREL(t_params%LMMAXD,t_params%LMMAXD),
     &t_params%SRREL(2,2,t_params%LMMAXD),
     &t_params%PHILDAU(t_params%IRMD,t_params%NATYPD)) !DOUBLE COMPLEX
     
      ALLOCATE(t_params%VINS(t_params%IRMIND:t_params%IRMD,
     &             t_params%LMPOTD,t_params%NSPOTD), !DOUBLE PRECISION
     &t_params%VISP(t_params%IRMD,t_params%NPOTD), 
     &t_params%VBC(2),                                           
     &t_params%VTREL(t_params%IRMD,t_params%NATYPD),
     &t_params%BTREL(t_params%IRMD,t_params%NATYPD),            
     &t_params%SOCSCALE(t_params%NATYPD),                                 
     &t_params%DRDIREL(t_params%IRMD,t_params%NATYPD),
     &t_params%R2DRDIREL(t_params%IRMD,t_params%NATYPD),      
     &t_params%RMREL(t_params%IRMD,t_params%NATYPD),
     &t_params%CMOMHOST(t_params%LMPOTD,t_params%NEMBD1),       
     &t_params%ECORE(20,t_params%NPOTD),
     &t_params%QMTET(t_params%NAEZD),
     &t_params%QMPHI(t_params%NAEZD),        
     &t_params%QMPHITAB(t_params%NAEZD,3),
     &t_params%QMTETTAB(t_params%NAEZD,3),          
     &t_params%QMGAMTAB(t_params%NAEZD,3),
     &t_params%ZAT(t_params%NATYPD),
     &t_params%R(t_params%IRMD,t_params%NATYPD),     
     &t_params%DRDI(t_params%IRMD,t_params%NATYPD),
     &t_params%RMTREF(t_params%NREFD),                  
     &t_params%VREF(t_params%NREFD),
     &t_params%CLEB(t_params%NCLEB,2),
     &t_params%RCLS(3,t_params%NACLSD,t_params%NCLSD),   
     &t_params%SOCSCL(t_params%LMAXD1,t_params%NATYPD),
     &t_params%CSCL(t_params%LMAXD1,t_params%NATYPD),        
     &t_params%RBASIS(3,t_params%NEMBD2),
     &t_params%RR(3,0:t_params%NRD),
     &t_params%CONC(t_params%NATYPD),   
     &t_params%RROT(48,3,t_params%NSHELD),
     &t_params%RATOM(3,t_params%NSHELD),
     &t_params%A(t_params%NATYPD),
     &t_params%B(t_params%NATYPD),
     &t_params%THETAS(t_params%IRID,t_params%NFUND,
     &               t_params%NCELLD),
     &t_params%RMT(t_params%NATYPD),
     &t_params%RMTNEW(t_params%NATYPD),
     &t_params%RWS(t_params%NATYPD),
     &t_params%GSH(t_params%NGSHD),
     &t_params%EREFLDAU(t_params%NATYPD),
     &t_params%UEFF(t_params%NATYPD),
     &t_params%JEFF(t_params%NATYPD),
     &t_params%ULDAU(t_params%MMAXD,t_params%MMAXD,
     &         t_params%MMAXD,t_params%MMAXD,t_params%NATYPD),
     &t_params%WLDAU(t_params%MMAXD,t_params%MMAXD,
     &         t_params%NSPIND,t_params%NATYPD),              
     &t_params%RPAN_INTERVALL(0:t_params%NTOTD,t_params%NATYPD),
     &t_params%RNEW(t_params%NTOTD*(t_params%NCHEBD+1),
     &             t_params%NATYPD),             
     &t_params%MVEVI(t_params%NATYPD,3,t_params%NMVECMAX),
     &t_params%MVEVIEF(t_params%NATYPD,3,t_params%NMVECMAX),      
     &t_params%THETASNEW(t_params%NTOTD*(t_params%NCHEBD+1),
     &                  t_params%NFUND,t_params%NCELLD),
     &t_params%RHO2NS(t_params%IRMD,t_params%LMPOTD,t_params%NATYPD,2),
     &t_params%R2NEF(t_params%IRMD,t_params%LMPOTD,t_params%NATYPD,2),
     &t_params%RHOC(t_params%IRMD,t_params%NPOTD),
     &t_params%DENEFAT(t_params%NATYPD),
     &t_params%ESPV(0:t_params%LMAXD1,t_params%NPOTD),
     &t_params%EDC(t_params%NATYPD),
     &t_params%EU(t_params%NATYPD),
     &t_params%RHOORB(t_params%IRMD*t_params%KREL+(1-t_params%KREL),
     &                t_params%NATYPD),
     &t_params%ECOREREL(t_params%KREL*20+(1-t_params%KREL),
     &                  t_params%NPOTD),
     &t_params%RCLSIMP(3,t_params%NATOMIMPD))    !DOUBLE PRECISION
     
      ALLOCATE(t_params%LOPT(t_params%NATYPD),
     &t_params%ITLDAU(t_params%NATYPD),
     &t_params%IRSHIFT(t_params%NATYPD),
     &t_params%JWSREL(t_params%NATYPD),
     &t_params%ZREL(t_params%NATYPD),
     &t_params%LCORE(20,t_params%NPOTD),
     &t_params%NCORE(t_params%NPOTD),
     &t_params%IPAN(t_params%NATYPD),
     &t_params%IRCUT(0:t_params%IPAND,t_params%NATYPD),
     &t_params%JEND(t_params%LMPOTD,0:t_params%LMAXD,0:t_params%LMAXD),
     &t_params%ICLEB(t_params%NCLEB,4),
     &t_params%ATOM(t_params%NACLSD,t_params%NEMBD2),
     &t_params%CLS(t_params%NEMBD2),
     &t_params%NACLS(t_params%NCLSD),
     &t_params%LOFLM(t_params%LM2D),
     &t_params%EZOA(t_params%NACLSD,t_params%NEMBD2),
     &t_params%KAOEZ(t_params%NATYPD,t_params%NEMBD2),
     &t_params%IQAT(t_params%NATYPD),
     &t_params%ICPA(t_params%NAEZD),
     &t_params%NOQ(t_params%NAEZD),
     &t_params%KMESH(t_params%IEMXD),
     &t_params%NSHELL(0:t_params%NSHELD),
     &t_params%NSH1(t_params%NSHELD),
     &t_params%NSH2(t_params%NSHELD),
     &t_params%IJTABCALC(t_params%NOFGIJD),
     &t_params%IJTABCALC_I(t_params%NOFGIJD),
     &t_params%IJTABSYM(t_params%NOFGIJD),
     &t_params%IJTABSH(t_params%NOFGIJD),
     &t_params%ISH(t_params%NSHELD,t_params%NOFGIJD),
     &t_params%JSH(t_params%NSHELD,t_params%NOFGIJD),
     &t_params%IQCALC(t_params%NAEZD),
     &t_params%ICHECK(t_params%NAEZDPD,t_params%NAEZDPD),
     &t_params%ATOMIMP(t_params%NATOMIMPD),
     &t_params%REFPOT(t_params%NEMBD2),
     &t_params%IRREL(2,2,t_params%LMMAXD),
     &t_params%NRREL(2,t_params%LMMAXD),
     &t_params%IFUNM1(t_params%LMXSPD,t_params%NATYPD),
     &t_params%ITITLE(20,t_params%NPOTD),
     &t_params%LMSP1(t_params%LMXSPD,t_params%NATYPD),
     &t_params%NTCELL(t_params%NATYPD),
     &t_params%IXIPOL(t_params%NATYPD),
     &t_params%IRNS(t_params%NATYPD),
     &t_params%IFUNM(t_params%NATYPD,t_params%LMXSPD),
     &t_params%LLMSP(t_params%NATYPD,t_params%NFUND),
     &t_params%LMSP(t_params%NATYPD,t_params%LMXSPD),
     &t_params%IMT(t_params%NATYPD),
     &t_params%IRC(t_params%NATYPD),
     &t_params%IRMIN(t_params%NATYPD),
     &t_params%IRWS(t_params%NATYPD),
     &t_params%NFU(t_params%NATYPD),
     &t_params%HOSTIMP(0:t_params%NATYPD),
     &t_params%ILM(t_params%NGSHD,3),
     &t_params%IMAXSH(0:t_params%LMPOTD),
     &t_params%NPAN_LOG(t_params%NATYPD),
     &t_params%NPAN_EQ(t_params%NATYPD),
     &t_params%NPAN_TOT(t_params%NATYPD),
     &t_params%IPAN_INTERVALL(0:t_params%NTOTD,t_params%NATYPD),
     &t_params%NKCORE(20,t_params%NATYPD),
     &t_params%KAPCORE(20,t_params%NPOTD)) !INTEGER
     
      ALLOCATE(t_params%SYMUNITARY(t_params%NSYMAXD), !LOGICALS
     &         t_params%VACFLAG(2))
     
      ALLOCATE(t_params%TXC(5)) !CHARACTER*24
      
      if(.not.allocated(t_params%TESTC)) then
        ALLOCATE(t_params%TESTC(32),t_params%OPTC(32)) !CHARACTER*8
      end if
      
      if (.not.allocated(t_params%nofks)) then
        allocate(t_params%BZKP(3,t_params%KPOIBZ,t_params%MAXMESH), ! double precision
     &           t_params%VOLCUB(t_params%KPOIBZ,t_params%MAXMESH), ! double precision
     &           t_params%VOLBZ(t_params%MAXMESH), ! double precision
     &           t_params%NOFKS(t_params%MAXMESH)) ! integer
      end if
      
      if (.not.allocated(t_params%theta)) then
        allocate(t_params%THETA(t_params%natyp),    !double precision
     &           t_params%PHI(t_params%natyp))      !double precision
      end if
      
      end subroutine init_t_params


      
#ifdef CPP_MPI
      subroutine bcast_t_params_scalars(t_params)
      ! broadcast scalar parameters, deal with arrays later
      use mpi
      use mod_mympi, only: master
      implicit none

      type(type_params), intent(inout) :: t_params
      integer :: ierr
      integer :: blocklen1(t_params%Nscalars),etype1(t_params%Nscalars),
     &           myMPItype1
      integer(kind=MPI_ADDRESS_KIND) :: disp1(t_params%Nscalars), base
      
      !INTEGER 
      call MPI_Get_address(t_params%IEMXD,      disp1(1), ierr) 
      call MPI_Get_address(t_params%IRMIND,     disp1(2), ierr) 
      call MPI_Get_address(t_params%IRMD,       disp1(3), ierr) 
      call MPI_Get_address(t_params%LMPOTD,     disp1(4), ierr)  
      call MPI_Get_address(t_params%NSPOTD,     disp1(5), ierr)   
      call MPI_Get_address(t_params%NPOTD,      disp1(6), ierr)  
      call MPI_Get_address(t_params%NATYPD,     disp1(7), ierr)   
      call MPI_Get_address(t_params%NEMBD1 ,    disp1(8), ierr)          
      call MPI_Get_address(t_params%LMMAXD,     disp1(9), ierr)   
      call MPI_Get_address(t_params%NAEZD,      disp1(10), ierr)  
      call MPI_Get_address(t_params%IPAND,      disp1(11), ierr)  
      call MPI_Get_address(t_params%NEMBD2,     disp1(12), ierr)   
      call MPI_Get_address(t_params%NREFD,      disp1(13), ierr)  
      call MPI_Get_address(t_params%LMAXD,      disp1(14), ierr)  
      call MPI_Get_address(t_params%NCLEB,      disp1(15), ierr)  
      call MPI_Get_address(t_params%NACLSD,     disp1(16), ierr)   
      call MPI_Get_address(t_params%NCLSD  ,    disp1(17), ierr)    
      call MPI_Get_address(t_params%LM2D,       disp1(18), ierr) 
      call MPI_Get_address(t_params%LMAXD1,     disp1(19), ierr)   
      call MPI_Get_address(t_params%NRD,        disp1(20), ierr) 
      call MPI_Get_address(t_params%NSHELD,     disp1(21), ierr)   
      call MPI_Get_address(t_params%NSYMAXD,    disp1(22), ierr)    
      call MPI_Get_address(t_params%NAEZDPD,    disp1(23), ierr)    
      call MPI_Get_address(t_params%NATOMIMPD,  disp1(24), ierr)      
      call MPI_Get_address(t_params%NOFGIJD  ,  disp1(25), ierr)      
      call MPI_Get_address(t_params%NSPIND,     disp1(26), ierr)   
      call MPI_Get_address(t_params%NSPINDD,    disp1(27), ierr)    
      call MPI_Get_address(t_params%IRID,       disp1(28), ierr) 
      call MPI_Get_address(t_params%NFUND,      disp1(29), ierr)  
      call MPI_Get_address(t_params%NCELLD,     disp1(30), ierr)   
      call MPI_Get_address(t_params%LMXSPD,     disp1(31), ierr)   
      call MPI_Get_address(t_params%NGSHD,      disp1(32), ierr)  
      call MPI_Get_address(t_params%KREL,       disp1(33), ierr)  
      call MPI_Get_address(t_params%MMAXD,      disp1(34), ierr)      
      call MPI_Get_address(t_params%IELAST,     disp1(35), ierr)    
      call MPI_Get_address(t_params%NPOL,       disp1(36), ierr)  
      call MPI_Get_address(t_params%NPNT1,      disp1(37), ierr)  
      call MPI_Get_address(t_params%NPNT2,      disp1(38), ierr)  
      call MPI_Get_address(t_params%NPNT3,      disp1(39), ierr)  
      call MPI_Get_address(t_params%ITSCF,      disp1(40), ierr)  
      call MPI_Get_address(t_params%SCFSTEPS,   disp1(41), ierr)     
      call MPI_Get_address(t_params%LLY,        disp1(42), ierr)    
      call MPI_Get_address(t_params%NSRA,       disp1(43), ierr) 
      call MPI_Get_address(t_params%INS,        disp1(44), ierr) 
      call MPI_Get_address(t_params%NATYP,      disp1(45), ierr)  
      call MPI_Get_address(t_params%NAEZ,       disp1(46), ierr) 
      call MPI_Get_address(t_params%NINEQ,      disp1(47), ierr)  
      call MPI_Get_address(t_params%NREF,       disp1(48), ierr) 
      call MPI_Get_address(t_params%NSPIN,      disp1(49), ierr)  
      call MPI_Get_address(t_params%LMAX,       disp1(50), ierr) 
      call MPI_Get_address(t_params%NOFGIJ,     disp1(51), ierr)   
      call MPI_Get_address(t_params%NCLS,       disp1(52), ierr) 
      call MPI_Get_address(t_params%ICST,       disp1(53), ierr)  
      call MPI_Get_address(t_params%IEND,       disp1(54), ierr)  
      call MPI_Get_address(t_params%ICC,        disp1(55), ierr) 
      call MPI_Get_address(t_params%IGF,        disp1(56), ierr) 
      call MPI_Get_address(t_params%NLBASIS,    disp1(57), ierr)    
      call MPI_Get_address(t_params%NRBASIS,    disp1(58), ierr)    
      call MPI_Get_address(t_params%NCPA,       disp1(59), ierr) 
      call MPI_Get_address(t_params%ITCPAMAX,   disp1(60), ierr)     
      call MPI_Get_address(t_params%KMROT,      disp1(61), ierr)  
      call MPI_Get_address(t_params%MAXMESH,    disp1(62), ierr)    
      call MPI_Get_address(t_params%NSYMAT,     disp1(63), ierr)   
      call MPI_Get_address(t_params%NATOMIMP,   disp1(64), ierr)     
      call MPI_Get_address(t_params%INVMOD,     disp1(65), ierr)   
      call MPI_Get_address(t_params%NQCALC  ,   disp1(66), ierr)                  
      call MPI_Get_address(t_params%INTERVX,    disp1(67), ierr)    
      call MPI_Get_address(t_params%INTERVY,    disp1(68), ierr)    
      call MPI_Get_address(t_params%INTERVZ,    disp1(69), ierr)    
      call MPI_Get_address(t_params%LPOT,       disp1(70), ierr) 
      call MPI_Get_address(t_params%LMPOT,      disp1(71), ierr)  
      call MPI_Get_address(t_params%NRIGHT,     disp1(72), ierr)   
      call MPI_Get_address(t_params%NLEFT,      disp1(73), ierr)  
      call MPI_Get_address(t_params%IMIX,       disp1(74), ierr) 
      call MPI_Get_address(t_params%ITDBRY,     disp1(75), ierr)   
      call MPI_Get_address(t_params%KPRE   ,    disp1(76), ierr)                   
      call MPI_Get_address(t_params%KSHAPE,     disp1(77), ierr)   
      call MPI_Get_address(t_params%KTE,        disp1(78), ierr) 
      call MPI_Get_address(t_params%KVMAD,      disp1(79), ierr)  
      call MPI_Get_address(t_params%KXC,        disp1(80), ierr) 
      call MPI_Get_address(t_params%ISHIFT,     disp1(81), ierr)   
      call MPI_Get_address(t_params%KFORCE  ,   disp1(82), ierr)                           
      call MPI_Get_address(t_params%IDOLDAU,    disp1(83), ierr)    
      call MPI_Get_address(t_params%ITRUNLDAU,  disp1(84), ierr)      
      call MPI_Get_address(t_params%NTLDAU,     disp1(85), ierr)                                  
      call MPI_Get_address(t_params%NPOLSEMI,   disp1(86), ierr)     
      call MPI_Get_address(t_params%N1SEMI,     disp1(87), ierr)   
      call MPI_Get_address(t_params%N2SEMI,     disp1(88), ierr)   
      call MPI_Get_address(t_params%N3SEMI,     disp1(89), ierr)   
      call MPI_Get_address(t_params%IESEMICORE, disp1(90), ierr)       
      call MPI_Get_address(t_params%ITMPDIR,    disp1(91), ierr)    
      call MPI_Get_address(t_params%ILTMP,      disp1(92), ierr)  
      call MPI_Get_address(t_params%NCHEB,      disp1(93), ierr)  
      call MPI_Get_address(t_params%NTOTD,      disp1(94), ierr)  
      call MPI_Get_address(t_params%NCHEBD,     disp1(95), ierr)   
      !DOUBPLE PRECISION
      call MPI_Get_address(t_params%EBOTSEMI,   disp1(96), ierr)
      call MPI_Get_address(t_params%EMUSEMI,    disp1(97), ierr)
      call MPI_Get_address(t_params%TKSEMI,     disp1(98), ierr)
      call MPI_Get_address(t_params%FSEMICORE,  disp1(99), ierr)
      call MPI_Get_address(t_params%R_LOG,      disp1(100), ierr)
      call MPI_Get_address(t_params%E1,         disp1(101), ierr)
      call MPI_Get_address(t_params%E2,         disp1(102), ierr)
      call MPI_Get_address(t_params%TK,         disp1(103), ierr)
      call MPI_Get_address(t_params%EFERMI,     disp1(104), ierr)
      call MPI_Get_address(t_params%ALAT,       disp1(105), ierr)
      call MPI_Get_address(t_params%CPATOL,     disp1(106), ierr)               
      call MPI_Get_address(t_params%MIXING,     disp1(107), ierr)
      call MPI_Get_address(t_params%QBOUND,     disp1(108), ierr)
      call MPI_Get_address(t_params%FCM,        disp1(109), ierr)
      call MPI_Get_address(t_params%LAMBDA_XC,  disp1(110), ierr)
      call MPI_Get_address(t_params%TOLRDIF,    disp1(111), ierr)
      call MPI_Get_address(t_params%EFOLD,      disp1(112), ierr)
      call MPI_Get_address(t_params%CHRGOLD,    disp1(113), ierr)
      !DOUBLE COMPLEX
      call MPI_Get_address(t_params%DELTAE,     disp1(114), ierr)
      !LOGICAL
      call MPI_Get_address(t_params%LINTERFACE, disp1(115), ierr)
      call MPI_Get_address(t_params%LRHOSYM,    disp1(116), ierr)
      !CHARACTER*10
      call MPI_Get_address(t_params%SOLVER,     disp1(117), ierr)
      !CHARACTER*80
      call MPI_Get_address(t_params%TMPDIR,     disp1(118), ierr)
      !INTEGER
      call MPI_Get_address(t_params%I1,         disp1(119), ierr)
      call MPI_Get_address(t_params%NMVECMAX,   disp1(120), ierr)
      call MPI_Get_address(t_params%ITAB,       disp1(121), ierr)
      !DOUBLE PRECISION
      call MPI_Get_address(t_params%LASTERR,      disp1(122), ierr)
      call MPI_Get_address(t_params%DENEF,        disp1(123), ierr)
      call MPI_Get_address(t_params%CHRGSEMICORE, disp1(124), ierr)

      base  = disp1(1)
      disp1 = disp1 - base

      blocklen1(1:116)=1
      blocklen1(117)=10
      blocklen1(118)=80
      blocklen1(119:124)=1

      etype1(1:95) = MPI_INTEGER
      etype1(96:113) = MPI_DOUBLE_PRECISION
      etype1(114) = MPI_DOUBLE_COMPLEX
      etype1(115:116) = MPI_LOGICAL
      etype1(117:118) = MPI_CHARACTER
      etype1(119:121) = MPI_INTEGER
      etype1(122:124) = MPI_DOUBLE_PRECISION

      call MPI_Type_create_struct(t_params%Nscalars, blocklen1, disp1,
     &                            etype1, myMPItype1, ierr)
      if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_params'
  
      call MPI_Type_commit(myMPItype1, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error comiting create_mpimsk_t_params'
  
      call MPI_Bcast(t_params%Nscalars, 1, myMPItype1, master, 
     &                          MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_params'
  
      call MPI_Type_free(myMPItype1, ierr)
      
      ! somehow this parameter gets overlooked in the communication, possibly a but somewhere, but for now this workaround does the job
      call MPI_Bcast(t_params%NCHEBD, 1, MPI_INTEGER, master, 
     &               MPI_COMM_WORLD, ierr)
      
      end subroutine bcast_t_params_scalars
      
      
      
      subroutine bcast_t_params_arrays(t_params)
      ! broadcast arrays from t_params
      use mpi
      use mod_mympi, only: master, myrank
      implicit none
      
      type(type_params), intent(inout) :: t_params
      integer :: ierr
      
      !DOUBLE COMPLEX
      call MPI_Bcast(t_params%EZ,t_params%IEMXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%WEZ,t_params%IEMXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DROTQ,t_params%LMMAXD*
     &          t_params%LMMAXD*t_params%NAEZD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DSYMLL,t_params%LMMAXD*
     &          t_params%LMMAXD*t_params%NSYMAXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LEFTTINVLL,t_params%LMMAXD*t_params%LMMAXD
     &           *t_params%NEMBD1*t_params%NSPINDD*t_params%IEMXD,
     &            MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RIGHTTINVLL,t_params%LMMAXD*
     &          t_params%LMMAXD*
     &          t_params%NEMBD1*t_params%NSPINDD*t_params%IEMXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CREL,t_params%LMMAXD*t_params%LMMAXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RC,t_params%LMMAXD*t_params%LMMAXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RREL,t_params%LMMAXD*t_params%LMMAXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SRREL,2*2*t_params%LMMAXD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%PHILDAU,t_params%IRMD*t_params%NATYPD,
     &          MPI_DOUBLE_COMPLEX,master,MPI_COMM_WORLD,ierr)
      
      !DOUBLE PRECISION
      call MPI_Bcast(t_params%VINS,((t_params%IRMD-t_params%IRMIND+1)*
     &               t_params%LMPOTD*t_params%NSPOTD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VISP,(t_params%IRMD*t_params%NPOTD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VBC,2,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VTREL,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%BTREL,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SOCSCALE,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DRDIREL,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%R2DRDIREL,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMREL,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CMOMHOST,(t_params%LMPOTD*t_params%NEMBD1)
     &         ,MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ECORE,(20*t_params%NPOTD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMTET,(t_params%NAEZD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMPHI,(t_params%NAEZD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMPHITAB,(t_params%NAEZD*3),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMTETTAB,(t_params%NAEZD*3),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%QMGAMTAB,(t_params%NAEZD*3),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ZAT,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%R,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DRDI,(t_params%IRMD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMTREF,(t_params%NREFD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VREF,(t_params%NREFD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CLEB,(t_params%NCLEB*2),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RCLS,(3*t_params%NACLSD*t_params%NCLSD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%SOCSCL,(t_params%LMAXD1*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CSCL,(t_params%LMAXD1*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RBASIS,(3*t_params%NEMBD2),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RR,(3*(t_params%NRD+1)),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CONC,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RROT,(48*3*t_params%NSHELD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RATOM,(3*t_params%NSHELD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%A,t_params%NATYPD,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%B,t_params%NATYPD,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%THETAS,(t_params%IRID*t_params%NFUND*
     &          t_params%NCELLD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMT,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RMTNEW,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RWS,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%GSH,(t_params%NGSHD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EREFLDAU,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%UEFF,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JEFF,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ULDAU,(t_params%MMAXD*t_params%MMAXD*
     &          t_params%MMAXD*t_params%MMAXD*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%WLDAU,(t_params%MMAXD*t_params%MMAXD*
     &          t_params%NSPIND*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RPAN_INTERVALL,((t_params%NTOTD+1)*
     &          t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RNEW,(t_params%NTOTD*(t_params%NCHEBD+1)*
     &          t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%THETASNEW,(t_params%NTOTD*
     &         (t_params%NCHEBD+1)*
     &          t_params%NFUND*t_params%NCELLD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%MVEVI,t_params%NATYPD*3*t_params%NMVECMAX,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%MVEVIEF,
     &          t_params%NATYPD*3*t_params%NMVECMAX,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHO2NS,
     &         (t_params%IRMD*t_params%LMPOTD*t_params%NATYPD*2),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%R2NEF,
     &         (t_params%IRMD*t_params%LMPOTD*t_params%NATYPD*2),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHOC,(t_params%IRMD*t_params%NPOTD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%DENEFAT,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ESPV,(t_params%LMAXD1+1)*t_params%NPOTD,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EDC,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EU,(t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RHOORB,
     &  (t_params%IRMD*t_params%KREL+(1-t_params%KREL)*t_params%NATYPD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ECOREREL,
     &              (t_params%KREL*20+(1-t_params%KREL)*t_params%NPOTD),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%theta,t_params%natyp,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%phi,t_params%natyp,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%RCLSIMP,3*t_params%NATOMIMPD,
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr)
     
      !INTEGERS
      call MPI_Bcast(t_params%LOPT,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ITLDAU,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRSHIFT,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JWSREL,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ZREL,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LCORE,(20*t_params%NPOTD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NCORE,(t_params%NPOTD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IPAN,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRCUT,((t_params%IPAND+1)*
     &          t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JEND,(t_params%LMPOTD*(t_params%LMAXD+1)*
     &       (t_params%LMAXD+1)),MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICLEB,(t_params%NCLEB*4),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ATOM,(t_params%NACLSD*t_params%NEMBD2),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%CLS,(t_params%NEMBD2),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NACLS,(t_params%NCLSD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LOFLM,(t_params%LM2D),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%EZOA,(t_params%NACLSD*t_params%NEMBD2),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KAOEZ,(t_params%NATYPD*t_params%NEMBD2),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IQAT,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICPA,(t_params%NAEZD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NOQ,(t_params%NAEZD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KMESH,(t_params%IEMXD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSHELL,((t_params%NSHELD+1)),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSH1,(t_params%NSHELD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NSH2,(t_params%NSHELD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABCALC,(t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABCALC_I,(t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABSYM,(t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IJTABSH,(t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ISH,(t_params%NSHELD*t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%JSH,(t_params%NSHELD*t_params%NOFGIJD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IQCALC,(t_params%NAEZD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ICHECK,(t_params%NAEZDPD*
     &          t_params%NAEZDPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ATOMIMP,(t_params%NATOMIMPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%REFPOT,(t_params%NEMBD2),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRREL,(2*2*t_params%LMMAXD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NRREL,(2*t_params%LMMAXD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IFUNM1,(t_params%LMXSPD*t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ITITLE,(20*t_params%NPOTD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LMSP1,(t_params%LMXSPD*t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NTCELL,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IXIPOL,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRNS,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IFUNM,(t_params%NATYPD*t_params%LMXSPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LLMSP,(t_params%NATYPD*t_params%NFUND),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%LMSP,(t_params%NATYPD*t_params%LMXSPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IMT,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRC,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRMIN,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IRWS,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NFU,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%HOSTIMP,((t_params%NATYPD+1)),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%ILM,(t_params%NGSHD*3),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IMAXSH,((t_params%LMPOTD+1)),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_LOG,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_EQ,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NPAN_TOT,(t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%IPAN_INTERVALL,((t_params%NTOTD+1)*
     &          t_params%NATYPD)
     &         ,MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%NKCORE,(20*t_params%NATYPD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%KAPCORE,(20*t_params%NPOTD),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr)

      !LOGICAL
      call MPI_Bcast(t_params%SYMUNITARY,(t_params%NSYMAXD),
     &          MPI_LOGICAL,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t_params%VACFLAG,2,MPI_LOGICAL,master,
     &          MPI_COMM_WORLD,ierr)

      !CHARACTER
      call MPI_Bcast(t_params%TXC,5*24, !5 entries of length 24
     &          MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*24
      call MPI_Bcast(t_params%TESTC,32*8, !32 entries of length 8
     &          MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*8
      call MPI_Bcast(t_params%OPTC,32*8, !32 entries of length 8
     &          MPI_CHARACTER,master,MPI_COMM_WORLD,ierr) !CHARACTER*8
     
      ! kpoints
      call MPI_Bcast(t_params%BZKP,(3*t_params%KPOIBZ*t_params%MAXMESH),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! double precision
      call MPI_Bcast(t_params%VOLCUB,(t_params%KPOIBZ*t_params%MAXMESH),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! double precision
      call MPI_Bcast(t_params%VOLBZ,(t_params%MAXMESH),
     &          MPI_DOUBLE_PRECISION,master,MPI_COMM_WORLD,ierr) ! double precision
      call MPI_Bcast(t_params%NOFKS,(t_params%MAXMESH),
     &          MPI_INTEGER,master,MPI_COMM_WORLD,ierr) ! integer
      
      end subroutine bcast_t_params_arrays
#endif

      subroutine fill_t_params_scalars(
     &        IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,NEMBD1,
     &        LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,
     &        LM2D,LMAXD1,NRD,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,
     &        NSPIND,NSPINDD,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD,
     &        IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,LLY,
     &        NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,NOFGIJ,
     &        NCLS,ICST,IEND,ICC,IGF,NLBASIS,NRBASIS,NCPA,ITCPAMAX,
     &        KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC,INTERVX,
     &        INTERVY,INTERVZ,LPOT,LMPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,
     &        KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE,IDOLDAU,ITRUNLDAU,
     &        NTLDAU,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE,
     &        EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG,E1,E2,TK,EFERMI,
     &        ALAT,CPATOL,MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF,
     &        LINTERFACE,LRHOSYM,SOLVER,TMPDIR,ITMPDIR,ILTMP,
     &        NTOTD,NCHEBD,NCHEB,DELTAE,EFOLD,CHRGOLD,t_params)
      ! fill scalars into t_params
      IMPLICIT NONE
      type(type_params), intent(inout) :: t_params
C     ..
C     .. Scalar arguments
      INTEGER IEMXD,IRMIND,IRMD,LMPOTD,NSPOTD,NPOTD,NATYPD,NEMBD1,
     &        LMMAXD,NAEZD,IPAND,NEMBD2,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,
     &        LM2D,LMAXD1,NRD,NSHELD,NSYMAXD,NAEZDPD,NATOMIMPD,NOFGIJD,
     &        NSPIND,NSPINDD,IRID,NFUND,NCELLD,LMXSPD,NGSHD,KREL,MMAXD
      INTEGER IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,LLY
      INTEGER NSRA,INS,NATYP,NAEZ,NINEQ,NREF,NSPIN,LMAX,NOFGIJ,
     &        NCLS,ICST,IEND,ICC,IGF,NLBASIS,NRBASIS,NCPA,ITCPAMAX
      INTEGER KMROT,MAXMESH,NSYMAT,NATOMIMP,INVMOD,NQCALC
      INTEGER INTERVX,INTERVY,INTERVZ
      INTEGER LPOT,LMPOT,NRIGHT,NLEFT,IMIX,ITDBRY,KPRE,
     &        KSHAPE,KTE,KVMAD,KXC,ISHIFT,KFORCE
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU,NTOTD,NCHEBD,NCHEB
      INTEGER NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IESEMICORE
      DOUBLE PRECISION EBOTSEMI,EMUSEMI,TKSEMI,FSEMICORE,R_LOG
      DOUBLE PRECISION E1,E2,TK,EFERMI,ALAT,CPATOL,EFOLD,CHRGOLD
      DOUBLE PRECISION MIXING,QBOUND,FCM,LAMBDA_XC,TOLRDIF
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
      LOGICAL LINTERFACE,LRHOSYM
      CHARACTER*10 SOLVER
      CHARACTER*80 TMPDIR
      INTEGER ITMPDIR,ILTMP
C     ..
      ! fill scalars:
      t_params%IEMXD=      IEMXD       
      t_params%IRMIND=     IRMIND      
      t_params%IRMD=       IRMD        
      t_params%LMPOTD=     LMPOTD      
      t_params%NSPOTD=     NSPOTD      
      t_params%NPOTD=      NPOTD       
      t_params%NATYPD=     NATYPD      
      t_params%NEMBD1 =    NEMBD1      
      t_params%LMMAXD=     LMMAXD      
      t_params%NAEZD=      NAEZD       
      t_params%IPAND=      IPAND       
      t_params%NEMBD2=     NEMBD2      
      t_params%NREFD=      NREFD       
      t_params%LMAXD=      LMAXD       
      t_params%NCLEB=      NCLEB       
      t_params%NACLSD=     NACLSD      
      t_params%NCLSD  =    NCLSD       
      t_params%LM2D=       LM2D        
      t_params%LMAXD1=     LMAXD1      
      t_params%NRD=        NRD         
      t_params%NSHELD=     NSHELD      
      t_params%NSYMAXD=    NSYMAXD     
      t_params%NAEZDPD=    NAEZDPD     
      t_params%NATOMIMPD=  NATOMIMPD   
      t_params%NOFGIJD  =  NOFGIJD     
      t_params%NSPIND=     NSPIND      
      t_params%NSPINDD=    NSPINDD     
      t_params%IRID=       IRID        
      t_params%NFUND=      NFUND       
      t_params%NCELLD=     NCELLD      
      t_params%LMXSPD=     LMXSPD      
      t_params%NGSHD=      NGSHD       
      t_params%KREL=       KREL        
      t_params%MMAXD=      MMAXD       
      t_params%IELAST=     IELAST      
      t_params%NPOL=       NPOL        
      t_params%NPNT1=      NPNT1       
      t_params%NPNT2=      NPNT2       
      t_params%NPNT3=      NPNT3       
      t_params%ITSCF=      ITSCF       
      t_params%SCFSTEPS=   SCFSTEPS    
      t_params%LLY=        LLY         
      t_params%NSRA=       NSRA        
      t_params%INS=        INS         
      t_params%NATYP=      NATYP       
      t_params%NAEZ=       NAEZ        
      t_params%NINEQ=      NINEQ       
      t_params%NREF=       NREF        
      t_params%NSPIN=      NSPIN       
      t_params%LMAX=       LMAX        
      t_params%NOFGIJ=     NOFGIJ      
      t_params%NCLS=       NCLS        
      t_params%ICST=       ICST        
      t_params%IEND=       IEND        
      t_params%ICC=        ICC         
      t_params%IGF=        IGF         
      t_params%NLBASIS=    NLBASIS     
      t_params%NRBASIS=    NRBASIS     
      t_params%NCPA=       NCPA        
      t_params%ITCPAMAX=   ITCPAMAX    
      t_params%KMROT=      KMROT       
      t_params%MAXMESH=    MAXMESH     
      t_params%NSYMAT=     NSYMAT      
      t_params%NATOMIMP=   NATOMIMP    
      t_params%INVMOD=     INVMOD      
      t_params%NQCALC  =   NQCALC      
      t_params%INTERVX=    INTERVX     
      t_params%INTERVY=    INTERVY     
      t_params%INTERVZ=    INTERVZ     
      t_params%LPOT=       LPOT        
      t_params%LMPOT=      LMPOT       
      t_params%NRIGHT=     NRIGHT      
      t_params%NLEFT=      NLEFT       
      t_params%IMIX=       IMIX        
      t_params%ITDBRY=     ITDBRY      
      t_params%KPRE   =    KPRE        
      t_params%KSHAPE=     KSHAPE      
      t_params%KTE   =     KTE         
      t_params%KVMAD=      KVMAD       
      t_params%KXC=        KXC         
      t_params%ISHIFT=     ISHIFT      
      t_params%KFORCE  =   KFORCE      
      t_params%IDOLDAU=    IDOLDAU     
      t_params%ITRUNLDAU=  ITRUNLDAU   
      t_params%NTLDAU=     NTLDAU      
      t_params%NPOLSEMI=   NPOLSEMI    
      t_params%N1SEMI=     N1SEMI      
      t_params%N2SEMI=     N2SEMI      
      t_params%N3SEMI=     N3SEMI      
      t_params%IESEMICORE= IESEMICORE  
      t_params%ITMPDIR=    ITMPDIR     
      t_params%ILTMP=      ILTMP       
      t_params%NCHEB=      NCHEB       
      t_params%NTOTD=      NTOTD       
      t_params%NCHEBD=     NCHEBD      
      t_params%EBOTSEMI=   EBOTSEMI    
      t_params%EMUSEMI=    EMUSEMI     
      t_params%TKSEMI=     TKSEMI      
      t_params%FSEMICORE=  FSEMICORE   
      t_params%R_LOG=      R_LOG       
      t_params%E1=         E1          
      t_params%E2=         E2          
      t_params%TK=         TK          
      t_params%EFERMI=     EFERMI      
      t_params%ALAT=       ALAT        
      t_params%CPATOL=     CPATOL      
      t_params%MIXING=     MIXING      
      t_params%QBOUND=     QBOUND      
      t_params%FCM=        FCM         
      t_params%LAMBDA_XC=  LAMBDA_XC   
      t_params%TOLRDIF=    TOLRDIF     
      t_params%DELTAE=     DELTAE      
      t_params%LINTERFACE= LINTERFACE  
      t_params%LRHOSYM=    LRHOSYM     
      t_params%SOLVER=     SOLVER      
      t_params%TMPDIR=     TMPDIR      
      t_params%NMVECMAX=   4

      end subroutine fill_t_params_scalars
      
      
      
      subroutine fill_t_params_arrays(t_params,IEMXD,LMMAXD,NAEZD,
     &       NSYMAXD,NEMBD1,NSPINDD,IRMIND,IRMD,
     &       LMPOTD,NSPOTD,NPOTD,NATYPD,NRD,NEMBD2,NREFD,NCLEB,NCLSD,
     &       NACLSD,NSHELD,NGSHD,NFUND,IRID,NCELLD,MMAXD,LM2D,LMXSPD,
     &       LMAXD1,NSPIND,NTOTD,NCHEBD,IPAND,LMAXD,NOFGIJD,NAEZDPD,
     &       NATOMIMPD,
     &       EZ,WEZ,DROTQ,DSYMLL,LEFTTINVLL,
     &       RIGHTTINVLL,CREL,RC,RREL,SRREL,PHILDAU,VINS,VISP,VBC,VTREL,
     &       BTREL,SOCSCALE,DRDIREL,R2DRDIREL,RMREL,CMOMHOST,ECORE,QMTET
     &      ,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,ZAT,R,DRDI,RMTREF,VREF,
     &       CLEB,RCLS,SOCSCL,CSCL,RBASIS,RR,CONC,RROT,RATOM,A,B,THETAS,
     &       RMT,RMTNEW,RWS,GSH,EREFLDAU,UEFF,JEFF,ULDAU,WLDAU,
     &       RPAN_INTERVALL,RNEW,THETASNEW,LOPT,ITLDAU,IRSHIFT,JWSREL,
     &       ZREL,LCORE,NCORE,IPAN,IRCUT,JEND,ICLEB,ATOM,CLS,NACLS,LOFLM
     &      ,EZOA,KAOEZ,IQAT,ICPA,NOQ,KMESH,NSHELL,NSH1,NSH2,IJTABCALC,
     &       IJTABCALC_I,
     &       IJTABSYM,IJTABSH,ISH,JSH,IQCALC,ICHECK,ATOMIMP,REFPOT,IRREL
     &      ,NRREL,IFUNM1,ITITLE,LMSP1,NTCELL,IXIPOL,IRNS,IFUNM,LLMSP,
     &       LMSP,IMT,IRC,IRMIN,IRWS,NFU,HOSTIMP,ILM,IMAXSH,NPAN_LOG,
     &       NPAN_EQ,NPAN_TOT,IPAN_INTERVALL,SYMUNITARY,VACFLAG,TXC,
     &       RCLSIMP)   
      ! fill arrays after they have been allocated in init_t_params
C     ..
      IMPLICIT NONE
      type(type_params), intent(inout) :: t_params
C     ..
C     .. Scalars for array dimensions
      INTEGER IEMXD,LMMAXD,NAEZD,NSYMAXD,NEMBD1,NSPINDD,IRMIND,IRMD,
     &        LMPOTD,NSPOTD,NPOTD,NATYPD,NRD,NEMBD2,NREFD,NCLEB,NCLSD,
     &        NACLSD,NSHELD,NGSHD,NFUND,IRID,NCELLD,MMAXD,LM2D,LMXSPD,
     &        LMAXD1,NSPIND,NTOTD,NCHEBD,IPAND,LMAXD,NOFGIJD,NAEZDPD,
     &        NATOMIMPD
C     .. Array arguments
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD),DROTQ(LMMAXD,LMMAXD,NAEZD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),
     &               LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &               RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD)
      DOUBLE COMPLEX CREL(LMMAXD,LMMAXD),RC(LMMAXD,LMMAXD),
     &           RREL(LMMAXD,LMMAXD),SRREL(2,2,LMMAXD)
      DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD)
      DOUBLE PRECISION VBC(2)
      DOUBLE PRECISION VTREL(IRMD,NATYPD),BTREL(IRMD,NATYPD)
      DOUBLE PRECISION SOCSCALE(NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD,NATYPD),R2DRDIREL(IRMD,NATYPD),
     &                 RMREL(IRMD,NATYPD),CMOMHOST(LMPOTD,NEMBD1)
      DOUBLE PRECISION ECORE(20,NPOTD),QMTET(NAEZD),QMPHI(NAEZD)
      DOUBLE PRECISION QMPHITAB(NAEZD,3),QMTETTAB(NAEZD,3),
     &                 QMGAMTAB(NAEZD,3),ZAT(NATYPD),R(IRMD,NATYPD),
     &                 DRDI(IRMD,NATYPD),RMTREF(NREFD),
     &                 VREF(NREFD),CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD)
      DOUBLE PRECISION SOCSCL(LMAXD1,NATYPD),CSCL(LMAXD1,NATYPD)
      DOUBLE PRECISION RBASIS(3,NEMBD2),RR(3,0:NRD),CONC(NATYPD)
      DOUBLE PRECISION RROT(48,3,NSHELD),RATOM(3,NSHELD)
      DOUBLE PRECISION A(NATYPD),B(NATYPD) ,THETAS(IRID,NFUND,NCELLD)
      DOUBLE PRECISION RMT(NATYPD),RMTNEW(NATYPD),RWS(NATYPD),GSH(NGSHD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION RCLSIMP(3,NATOMIMPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C     ..
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)
      INTEGER LCORE(20,NPOTD),NCORE(NPOTD),IPAN(NATYPD),
     &        IRCUT(0:IPAND,NATYPD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     &        ICLEB(NCLEB,4),ATOM(NACLSD,NEMBD2),CLS(NEMBD2),
     &        NACLS(NCLSD)
      INTEGER LOFLM(LM2D),EZOA(NACLSD,NEMBD2),
     &        KAOEZ(NATYPD,NEMBD2),IQAT(NATYPD),
     &        ICPA(NAEZD),NOQ(NAEZD),KMESH(IEMXD)
      INTEGER NSHELL(0:NSHELD),NSH1(NSHELD),NSH2(NSHELD),
     &        IJTABCALC(NOFGIJD),IJTABCALC_I(NOFGIJD),IJTABSYM(NOFGIJD),
     &        IJTABSH(NOFGIJD)
      INTEGER ISH(NSHELD,NOFGIJD),JSH(NSHELD,NOFGIJD),IQCALC(NAEZD)
      INTEGER ICHECK(NAEZDPD,NAEZDPD),ATOMIMP(NATOMIMPD),REFPOT(NEMBD2)
      INTEGER IRREL(2,2,LMMAXD),NRREL(2,LMMAXD),IFUNM1(LMXSPD,NATYPD),
     &        ITITLE(20,NPOTD),LMSP1(LMXSPD,NATYPD),NTCELL(NATYPD)
      INTEGER IXIPOL(NATYPD),IRNS(NATYPD),IFUNM(NATYPD,LMXSPD)
      INTEGER LLMSP(NATYPD,NFUND),LMSP(NATYPD,LMXSPD),IMT(NATYPD),
     &        IRC(NATYPD),IRMIN(NATYPD),IRWS(NATYPD),NFU(NATYPD),
     &        HOSTIMP(0:NATYPD),ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),
     +        NPAN_TOT(NATYPD)
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     &                 RNEW(NTOTD*(NCHEBD+1),NATYPD),
     &                 THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
      LOGICAL SYMUNITARY(NSYMAXD),VACFLAG(2)
      CHARACTER*24 TXC(5)  
C     ..
      !fill arrays:
      t_params%EZ             = EZ             
      t_params%WEZ            = WEZ            
      t_params%DROTQ          = DROTQ          
      t_params%DSYMLL         = DSYMLL         
      t_params%LEFTTINVLL     = LEFTTINVLL     
      t_params%RIGHTTINVLL    = RIGHTTINVLL    
      t_params%CREL           = CREL           
      t_params%RC             = RC             
      t_params%RREL           = RREL           
      t_params%SRREL          = SRREL          
      t_params%PHILDAU        = PHILDAU        
      t_params%VINS           = VINS           
      t_params%VISP           = VISP           
      t_params%VBC            = VBC                                   
      if(t_params%KREL.gt.0) t_params%VTREL          = VTREL          
      if(t_params%KREL.gt.0) t_params%BTREL          = BTREL          
      t_params%SOCSCALE       = SOCSCALE       
      if(t_params%KREL.gt.0) t_params%DRDIREL        = DRDIREL        
      if(t_params%KREL.gt.0) t_params%R2DRDIREL      = R2DRDIREL      
      if(t_params%KREL.gt.0) t_params%RMREL          = RMREL          
      t_params%CMOMHOST       = CMOMHOST            
      t_params%ECORE          = ECORE          
      t_params%QMTET          = QMTET          
      t_params%QMPHI          = QMPHI          
      t_params%QMPHITAB       = QMPHITAB       
      t_params%QMTETTAB       = QMTETTAB         
      t_params%QMGAMTAB       = QMGAMTAB       
      t_params%ZAT            = ZAT            
      t_params%R              = R              
      t_params%DRDI           = DRDI           
      t_params%RMTREF         = RMTREF              
      t_params%VREF           = VREF           
      t_params%CLEB           = CLEB           
      t_params%RCLS           = RCLS           
      t_params%SOCSCL         = SOCSCL         
      t_params%CSCL           = CSCL                
      t_params%RBASIS         = RBASIS         
      t_params%RR             = RR             
      t_params%CONC           = CONC           
      t_params%RROT           = RROT           
      t_params%RATOM          = RATOM             
      t_params%A              = A              
      t_params%B              = B              
      t_params%THETAS         = THETAS         
      t_params%RMT            = RMT            
      t_params%RMTNEW         = RMTNEW         
      t_params%RWS            = RWS            
      t_params%GSH            = GSH            
      t_params%EREFLDAU       = EREFLDAU       
      t_params%UEFF           = UEFF           
      t_params%JEFF           = JEFF           
      if(t_params%IDOLDAU==1) t_params%ULDAU          = ULDAU             
      t_params%WLDAU          = WLDAU          
      t_params%RCLSIMP        = RCLSIMP
      t_params%RPAN_INTERVALL = RPAN_INTERVALL 
      t_params%RNEW           = RNEW                        
      t_params%THETASNEW      = THETASNEW      
      t_params%LOPT           = LOPT           
      t_params%ITLDAU         = ITLDAU         
      t_params%IRSHIFT        = IRSHIFT        
      t_params%JWSREL         = JWSREL         
      t_params%ZREL           = ZREL           
      t_params%LCORE          = LCORE          
      t_params%NCORE          = NCORE          
      t_params%IPAN           = IPAN                      
      t_params%IRCUT          = IRCUT          
      t_params%JEND           = JEND                       
      t_params%ICLEB          = ICLEB          
      t_params%ATOM           = ATOM           
      t_params%CLS            = CLS            
      t_params%NACLS          = NACLS                      
      t_params%LOFLM          = LOFLM          
      t_params%EZOA           = EZOA           
      t_params%KAOEZ          = KAOEZ          
      t_params%IQAT           = IQAT           
      t_params%ICPA           = ICPA           
      t_params%NOQ            = NOQ            
      t_params%KMESH          = KMESH          
      t_params%NSHELL         = NSHELL         
      t_params%NSH1           = NSH1           
      t_params%NSH2           = NSH2           
      t_params%IJTABCALC      = IJTABCALC      
      t_params%IJTABCALC_I    = IJTABCALC_I    
      t_params%IJTABSYM       = IJTABSYM       
      t_params%IJTABSH        = IJTABSH        
      t_params%ISH            = ISH            
      t_params%JSH            = JSH            
      t_params%IQCALC         = IQCALC             
      t_params%ICHECK         = ICHECK         
      t_params%ATOMIMP        = ATOMIMP        
      t_params%REFPOT         = REFPOT         
      t_params%IRREL          = IRREL          
      t_params%NRREL          = NRREL          
      t_params%IFUNM1         = IFUNM1         
      t_params%ITITLE         = ITITLE         
      t_params%LMSP1          = LMSP1          
      t_params%NTCELL         = NTCELL         
      t_params%IXIPOL         = IXIPOL         
      t_params%IRNS           = IRNS           
      t_params%IFUNM          = IFUNM          
      t_params%LLMSP          = LLMSP          
      t_params%LMSP           = LMSP           
      t_params%IMT            = IMT                       
      t_params%IRC            = IRC            
      t_params%IRMIN          = IRMIN          
      t_params%IRWS           = IRWS           
      t_params%NFU            = NFU            
      t_params%HOSTIMP        = HOSTIMP        
      t_params%ILM            = ILM            
      t_params%IMAXSH         = IMAXSH         
      t_params%NPAN_LOG       = NPAN_LOG       
      t_params%NPAN_EQ        = NPAN_EQ        
      t_params%NPAN_TOT       = NPAN_TOT       
      t_params%IPAN_INTERVALL = IPAN_INTERVALL 
      t_params%SYMUNITARY     = SYMUNITARY     
      t_params%VACFLAG        = VACFLAG        
      t_params%TXC            = TXC   
      t_params%NTOTD          = NTOTD          
      
      end subroutine fill_t_params_arrays



      subroutine get_params_1a(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,
     &       NCLSD,NREFD,NCLEB,NEMBD,NAEZD,LM2D,NSRA,INS,NAEZ,NATYP,
     &       NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY,
     &       KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,IEND,EZ,VINS,
     &       IRMIN,ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,VISP,
     &       RMTREF,VREF,CLEB,CSCL,SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,
     &       SOLVER,TMPDIR,DELTAE,tolrdif,NPAN_LOG,NPAN_EQ,
     &       NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,LMAXD,
     &       NTOTD,NRMAXD,R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,
     &       R2DRDIREL,RMREL,IRMIND,LMPOTD,NSPOTD,NPOTD,JWSREL,ZREL,
     &       ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
      ! get relevant parameters from t_params
!     ..
      use mod_types, only: t_inc
      IMPLICIT NONE
      type(type_params), intent(in) :: t_params 
      integer, intent(in)  :: IPAND,NATYPD,IRMD,NACLSD,NCLSD,NREFD,KREL,
     &                        NCLEB,NEMBD,NAEZD,LM2D,LMAXD,NTOTD,NRMAXD,
     &                        IEMXD,IRMIND,LMPOTD,NSPOTD,NPOTD,NATOMIMPD
      integer, intent(inout) :: NSRA,NAEZ,NATYP,NSPIN,ICST,
     &         IRCUT(0:IPAND,NATYPD),LMAX,NCLS,NINEQ,NREF,LLY,
     &         INS,ATOM(NACLSD,NAEZD+NEMBD),CLS(NAEZD+NEMBD),
     &         ICLEB(NCLEB,4),IPAN(NATYPD),IEND,LOPT(NATYPD),
     &         LOFLM(LM2D),NACLS(NCLSD),REFPOT(NAEZD+NEMBD),IELAST,
     &         IRWS(NATYPD),IRMIN(NATYPD),ITMPDIR,ILTMP,IDOLDAU,
     &         NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD),
     &         IPAN_INTERVALL(0:NTOTD,NATYPD),NTLDAU,ITLDAU(NATYPD),
     &         JWSREL(NATYPD),ZREL(NATYPD),ITSCF,NATOMIMP,
     &         ATOMIMP(NATOMIMPD),IQAT(NATYPD)
      double precision, intent(inout) :: ALAT,DRDI(IRMD,NATYPD),
     &         RMESH(IRMD,NATYPD),ZAT(NATYPD),RCLS(3,NACLSD,NCLSD),
     &         RMTREF(NREFD),VREF(NREFD),CLEB(NCLEB,2),
     &         CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL)),SOCSCALE(NATYPD),
     &         SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL)),
     &         EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD),tolrdif,
     &         RPAN_INTERVALL(0:NTOTD,NATYPD),RNEW(NRMAXD,NATYPD),R_LOG,
     &         VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD),
     &         VTREL(IRMD*KREL+(1-KREL),NATYPD),
     &         BTREL(IRMD*KREL+(1-KREL),NATYPD),
     &         DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &         R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &         RMREL(IRMD*KREL+(1-KREL),NATYPD)

      double complex, intent(inout) :: DELTAE,EZ(IEMXD)
      character*10, intent(inout) :: SOLVER
      character*80, intent(inout) :: TMPDIR
      
      NSRA    = t_params%NSRA   
      INS     = t_params%INS    
      NAEZ    = t_params%NAEZ   
      NATYP   = t_params%NATYP  
      NSPIN   = t_params%NSPIN  
      ICST    = t_params%ICST   
      IPAN    = t_params%IPAN   
      IRCUT   = t_params%IRCUT  
      LMAX    = t_params%LMAX   
      NCLS    = t_params%NCLS   
      NINEQ   = t_params%NINEQ  
      NREF    = t_params%NREF   
      IDOLDAU = t_params%IDOLDAU
      LLY     = t_params%LLY    
C ......................................................................
Consistency check 
C
      if( (KREL.EQ.1) .AND. (INS.NE.0) ) then
         write(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         stop ' set INS = 0 in the input'
      endif
C     
      if( NSRA.LE.2 ) then
         if( KREL.EQ.1 ) stop
     &        ' KVREL <= 1 in input, but relativistic program used'
      else
         if( KREL.EQ.0 ) stop
     &        ' KVREL > 1 in input, but non-relativistic program used'
      endif
C ......................................................................
      ALAT           = t_params%ALAT          
      ZAT            = t_params%ZAT           
      DRDI           = t_params%DRDI          
      RMESH          = t_params%R        
      RMTREF         = t_params%RMTREF        
      VREF           = t_params%VREF          
      IEND           = t_params%IEND          
      CLEB           = t_params%CLEB          
      RCLS           = t_params%RCLS          
      ATOM           = t_params%ATOM          
      CLS            = t_params%CLS           
      ICLEB          = t_params%ICLEB         
      LOFLM          = t_params%LOFLM         
      NACLS          = t_params%NACLS         
      REFPOT         = t_params%REFPOT        
      IRWS           = t_params%IRWS          
      IRMIN          = t_params%IRMIN         
      TOLRDIF        = t_params%TOLRDIF       
      DELTAE         = t_params%DELTAE          
      SOCSCALE       = t_params%SOCSCALE        
      TMPDIR         = t_params%TMPDIR        
      ITMPDIR        = t_params%ITMPDIR        
      ILTMP          = t_params%ILTMP         
      NPAN_LOG       = t_params%NPAN_LOG       
      NPAN_EQ        = t_params%NPAN_EQ       
      NCHEB          = t_params%NCHEB         
      R_LOG          = t_params%R_LOG         
      NPAN_TOT       = t_params%NPAN_TOT      
      RNEW           = t_params%RNEW          
      RPAN_INTERVALL = t_params%RPAN_INTERVALL 
      IPAN_INTERVALL = t_params%IPAN_INTERVALL 

      if(KREL.eq.1) then
         SOLVER = t_params%SOLVER
         SOCSCL = t_params%SOCSCL
         CSCL   = t_params%CSCL  
      endif
      
      if( IDOLDAU.eq.1 ) then
         NTLDAU   = t_params%NTLDAU  
         ITLDAU   = t_params%ITLDAU  
         LOPT     = t_params%LOPT    
         UEFF     = t_params%UEFF    
         JEFF     = t_params%JEFF    
         EREFLDAU = t_params%EREFLDAU
      endif
C ---------------------------------------------------------- energy_mesh
      IELAST = t_params%IELAST
      EZ     = t_params%EZ
C ------------------------------------------------------ input_potential
      VINS  = t_params%VINS 
      VISP  = t_params%VISP 
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL    
         DRDIREL   = t_params%DRDIREL  
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL     
         JWSREL    = t_params%JWSREL   
         VTREL     = t_params%VTREL    
         BTREL     = t_params%BTREL    
      END IF
      ITSCF    = t_params%ITSCF   
C ------------------------------------------------------ cluster atoms
      NATOMIMP   = t_params%NATOMIMP
      ATOMIMP    = t_params%ATOMIMP
      IQAT       = t_params%IQAT

      end subroutine get_params_1a


      subroutine get_params_1b(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,
     &       NCLSD,NREFD,NCLEB,NEMBD,NAEZD,NSRA,INS,NAEZ,NATYP,
     &       NSPIN,LMAX,NCLS,NREF,LLY,KREL,ATOM,CLS,NACLS,REFPOT,EZ,
     &       ITMPDIR,ILTMP,ALAT,RCLS,IEMXD,RMTREF,VREF,TMPDIR,
     &       LMAXD,NSHELD,NPRINCD,KPOIBZ,NTOTD,NSPOTD,ATOMIMP,NATOMIMPD,
     &       ICC,IGF,NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL,NRD,
     &       IDECI,RBASIS,RR,EZOA,NSHELL,KMROT,KAOEZ,ISH,JSH,NSH1,NSH2,
     &       NOQ,IQAT,NOFGIJD,NOFGIJ,NATOMIMP,CONC,KMESH,MAXMESH,NSYMAT,
     &       NQCALC,RATOM,RROT,DROTQ,IJTABCALC,IJTABCALC_I,IJTABSYM,
     &       IJTABSH,IQCALC,
     &       DSYMLL,INVMOD,ICHECK,SYMUNITARY,RC,CREL,RREL,SRREL,
     &       NRREL,IRREL,LEFTTINVLL,RIGHTTINVLL,VACFLAG,NOFKS,VOLBZ,
     &       BZKP,VOLCUB,WEZ,NEMBD1,LMMAXD,NSYMAXD,NSPINDD,MAXMSHD,
     &       RCLSIMP)
      ! get relevant parameters from t_params
!     ..
      use mod_types, only: t_inc
      IMPLICIT NONE
      type(type_params), intent(in) :: t_params 
      integer, intent(in)  :: IPAND,NATYPD,IRMD,NACLSD,NCLSD,NREFD,KREL,
     &         NCLEB,NEMBD,NAEZD,LMAXD,NTOTD,IEMXD,
     &         NSPOTD,NOFGIJD,NEMBD1,LMMAXD,NSYMAXD,NSPINDD
     &        ,MAXMSHD,NSHELD,NPRINCD,KPOIBZ,NATOMIMPD,NRD
      integer, intent(inout) :: NSRA,NAEZ,NATYP,NSPIN,MAXMESH,
     &         LMAX,NCLS,NREF,LLY,NSYMAT,
     &         INS,ATOM(NACLSD,NAEZD+NEMBD),CLS(NAEZD+NEMBD),NQCALC,
     &         NATOMIMP,ATOMIMP(NATOMIMPD),
     &         NACLS(NCLSD),REFPOT(NAEZD+NEMBD),IELAST,
     &         ITMPDIR,ILTMP,INVMOD,
     &         ICC,IGF,NLBASIS,NRBASIS
     &        ,NCPA,ICPA(NAEZD),ITCPAMAX,IDECI,EZOA(NACLSD,NAEZD+NEMBD),
     &         NSHELL(0:NSHELD),KMROT,KAOEZ(NATYPD,NAEZD+NEMBD),
     &         NOQ(NAEZD),IQAT(NATYPD),ISH(NSHELD,NOFGIJD),NOFGIJ,
     &         JSH(NSHELD,NOFGIJD),NSH1(NSHELD),NSH2(NSHELD),
     &         KMESH(IEMXD),IJTABCALC(NOFGIJD),IJTABCALC_I(NOFGIJD),
     &         IJTABSYM(NOFGIJD),
     &         IJTABSH(NOFGIJD),IQCALC(NAEZD),NOFKS(MAXMSHD),
     &         ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD),IRREL(2,2,LMMAXD),
     &         NRREL(2,LMMAXD)
      double precision, intent(inout) :: ALAT,
     &         RR(3,0:NRD),RCLS(3,NACLSD,NCLSD),
     &         RMTREF(NREFD),VREF(NREFD),
     &         CPATOL,VOLBZ(MAXMSHD),
     &         RATOM(3,NSHELD),
     &         RROT(48,3,NSHELD),
     &         CONC(NATYPD),
     &         RBASIS(3,NAEZD+NEMBD),
     &         BZKP(3,KPOIBZ,MAXMSHD),
     &         VOLCUB(KPOIBZ,MAXMSHD),
     &         RCLSIMP(3,NATOMIMPD)

      double complex, intent(inout) :: EZ(IEMXD),
     &         DROTQ(LMMAXD,LMMAXD,NAEZD),DSYMLL(LMMAXD,LMMAXD,NSYMAXD),
     &         RC(LMMAXD,LMMAXD),CREL(LMMAXD,LMMAXD),RREL(LMMAXD,LMMAXD)
     &        ,SRREL(2,2,LMMAXD),WEZ(IEMXD),
     &         LEFTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD),
     &         RIGHTTINVLL(LMMAXD,LMMAXD,NEMBD1,NSPINDD,IEMXD)
      character*80, intent(inout) :: TMPDIR
      logical, intent(inout) :: SYMUNITARY(NSYMAXD),VACFLAG(2)
      
      integer :: i,j,l,id
      
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST

      NSRA     = t_params%NSRA    
      INS      = t_params%INS   
      NATYP    = t_params%NATYP  
      NAEZ     = t_params%NAEZ  
      NSPIN    = t_params%NSPIN  
      LMAX     = t_params%LMAX   
      NREF     = t_params%NREF   
      ICC      = t_params%ICC    
      IGF      = t_params%IGF    
      NLBASIS  = t_params%NLBASIS 
      NRBASIS  = t_params%NRBASIS 
      NCPA     = t_params%NCPA    
      ICPA     = t_params%ICPA    
      ITCPAMAX = t_params%ITCPAMAX
      CPATOL   = t_params%CPATOL  
C ......................................................................
Consistency check 
C
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
C     
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
      IDECI = 0
C ......................................................................
      ALAT      = t_params%ALAT     
      RBASIS    = t_params%RBASIS   
      REFPOT    = t_params%REFPOT   
      RMTREF    = t_params%RMTREF   
      VREF      = t_params%VREF     
      RCLS      = t_params%RCLS     
      RR        = t_params%RR       
      ATOM      = t_params%ATOM     
      CLS       = t_params%CLS      
      NCLS      = t_params%NCLS     
      EZOA      = t_params%EZOA     
      NACLS     = t_params%NACLS    
      NSHELL    = t_params%NSHELL   
      KMROT     = t_params%KMROT    
      KAOEZ     = t_params%KAOEZ    
      IQAT      = t_params%IQAT      
      NOQ       = t_params%NOQ      
      CONC      = t_params%CONC     
      KMESH     = t_params%KMESH    
      MAXMESH   = t_params%MAXMESH  
      LLY       = t_params%LLY      
      NSYMAT    = t_params%NSYMAT   
      NATOMIMP  = t_params%NATOMIMP 
      NOFGIJ    = t_params%NOFGIJ   
      NQCALC    = t_params%NQCALC   
      RATOM     = t_params%RATOM    
      RROT      = t_params%RROT     
      NSH1      = t_params%NSH1     
      NSH2      = t_params%NSH2     
      DROTQ     = t_params%DROTQ    
      IJTABCALC = t_params%IJTABCALC
      IJTABCALC_I = t_params%IJTABCALC_I
      do I=1,NSHELL(0)
         do L=1,NOFGIJ
            ISH(I,L) = t_params%ISH(I,L)
            JSH(I,L) = t_params%JSH(I,L)
         end do
      end do
      IJTABSYM   = t_params%IJTABSYM  
      IJTABSH    = t_params%IJTABSH   
      IQCALC     = t_params%IQCALC    
      DSYMLL     = t_params%DSYMLL    
      INVMOD     = t_params%INVMOD    
      ICHECK     = t_params%ICHECK    
      ATOMIMP    = t_params%ATOMIMP   
      SYMUNITARY = t_params%SYMUNITARY
      TMPDIR     = t_params%TMPDIR    
      ITMPDIR    = t_params%ITMPDIR   
      ILTMP      = t_params%ILTMP     
      RCLSIMP    = t_params%RCLSIMP
      IF ( KREL.EQ.1 ) then
         RC    = t_params%RC   
         CREL  = t_params%CREL 
         RREL  = t_params%RREL 
         SRREL = t_params%SRREL
         NRREL = t_params%NRREL
         IRREL = t_params%IRREL
      endif
      IF ( OPT('DECIMATE') ) THEN
         LEFTTINVLL  = t_params%LEFTTINVLL 
         RIGHTTINVLL = t_params%RIGHTTINVLL
         VACFLAG     = t_params%VACFLAG    
         IDECI = 1
      END IF
C ------------------------------------------------------------- k-points
C
      IF(TEST('kptsfile')) THEN
      OPEN (52,FILE='kpoints',FORM='formatted')
      REWIND (52)
      write(1337,*) 'kpoints read from kpoints file due to test option 
     &              "kptsfile"'
      ENDIF
      DO L = 1,MAXMESH
         IF(TEST('kptsfile')) THEN 
            READ (52,FMT='(I8,f15.10)') NOFKS(L),VOLBZ(L)
            write(1337,*) 'kpts:',NOFKS(L),VOLBZ(L),t_params%NOFKS(L),
     &                 t_params%VOLBZ(L)
         ELSE
            NOFKS(L) = t_params%NOFKS(L)
            VOLBZ(L) = t_params%VOLBZ(L)
         ENDIF
         DO I=1,NOFKS(L)
            IF(TEST('kptsfile')) THEN 
               READ (52,FMT=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
            ELSE
               DO ID=1,3
!                   write(*,*)'bzkp', BZKP(ID,I,L), t_params%BZKP(ID,I,L)
                  BZKP(ID,I,L) = t_params%BZKP(ID,I,L)
               ENDDO
!                write(*,*)'volcub', VOLCUB(I,L), t_params%VOLCUB(I,L)
               VOLCUB(I,L) = t_params%VOLCUB(I,L)
            ENDIF
!             write(*,'(A,4F21.17,I5,F21.17)') 'bzkmesh input:',
!      &             (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L),NOFKS(L),VOLBZ(L)
         END DO
      END DO
      IF(TEST('kptsfile')) CLOSE (52)
C ---------------------------------------------------------- energy_mesh
C
      IELAST = t_params%IELAST
      EZ     = t_params%EZ
      WEZ    = t_params%WEZ
C -------------------------------------------------------------- itermdir
C
      IF (OPT('ITERMDIR')) THEN
         DROTQ = t_params%DROTQ
         IF ( KMROT.EQ.0 ) KMROT = 1
      END IF

      end subroutine get_params_1b


      subroutine get_params_1c(t_params,KREL,NAEZD,NATYPD,NCLEB,LM2D,
     &          NCHEBD,IPAND,LMPOTD,LMAXD,LMXSPD,NFUND,NPOTD,NCLSD,
     &          NTOTD,MMAXD,IEMXD,IRMD,NSRA,INS,NATYP,NAEZ,NSPIN,NACLS1,
     &          ICST,KMROT,IQAT,IDOLDAU,LMAX,IRWS,IPAN,IRCUT,IEND,ICLEB, 
     &          LOFLM,JEND,IFUNM1,LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,
     &          IRMIN,ITITLE,INTERVX,INTERVY,INTERVZ,LLY,ITMPDIR,
     &          ILTMP,NPAN_EQ,IPAN_INTERVALL,NPAN_LOG,NCHEB,NPAN_TOT,
     &          NTLDAU,LOPT,ITLDAU,IELAST,IESEMICORE,NPOL,IRSHIFT,JWSREL
     &         ,ZREL,ITRUNLDAU,QMTET,QMPHI,CONC,ALAT,ZAT,DRDI,RMESH,A,B,
     &          CLEB,THETAS,SOCSCALE,RPAN_INTERVALL,CSCL,RNEW,SOCSCL,
     &          THETASNEW,EFERMI,EREFLDAU,UEFF,JEFF,E1,E2,TK,VINS,VISP,
     &          ECORE,DRDIREL,R2DRDIREL,RMREL,VTREL,BTREL,WLDAU,ULDAU,EZ
     &         ,WEZ,PHILDAU,TMPDIR,SOLVER,NSPIND,NSPOTD,
     &          IRMIND,LMAXD1,NCELLD,IRID,R_LOG)
      ! get relevant parameters from t_params
!     ..
      use mod_types, only: t_inc
      IMPLICIT NONE
      type(type_params), intent(in) :: t_params 
      integer, intent(in)  :: KREL,NAEZD,NATYPD,NCLEB,LM2D,NTOTD,NCHEBD,
     &          IPAND,LMPOTD,LMAXD,LMXSPD,NFUND,NPOTD,NCLSD,MMAXD,IEMXD,
     &          IRMD,NSPIND,NSPOTD,IRMIND,LMAXD1,NCELLD,IRID
      integer, intent(inout) :: NSRA,INS,NATYP,NAEZ,NSPIN,ICST,KMROT,
     &          IQAT(NATYPD),IDOLDAU,LMAX,IRWS(NATYPD),IPAN(NATYPD),
     &          IRCUT(0:IPAND,NATYPD),IEND,ICLEB(NCLEB,4), LOFLM(LM2D),
     &          JEND(LMPOTD,0:LMAXD,0:LMAXD),IFUNM1(LMXSPD,NATYPD),
     &          LMSP1(LMXSPD,NATYPD),NFU(NATYPD),LLMSP(NATYPD,NFUND),
     &          LCORE(20,NPOTD),NCORE(NPOTD),NTCELL(NATYPD),
     &          IRMIN(NATYPD),ITITLE(20,NPOTD),INTERVX,INTERVY,INTERVZ,
     &          LLY,ITMPDIR,ILTMP,NPAN_EQ(NATYPD),
     &          IPAN_INTERVALL(0:NTOTD,NATYPD),NPAN_LOG(NATYPD),
     &          NCHEB,NPAN_TOT(NATYPD),NTLDAU,LOPT(NATYPD),
     &          ITLDAU(NATYPD),IELAST,IESEMICORE,NPOL,IRSHIFT(NATYPD),
     &          JWSREL(NATYPD),ZREL(NATYPD),ITRUNLDAU,NACLS1
      double precision, intent(inout) :: QMTET(NAEZD),QMPHI(NAEZD),
     &          CONC(NATYPD),ALAT,ZAT(NATYPD),DRDI(IRMD,NATYPD),
     &          RMESH(IRMD,NATYPD),A(NATYPD),B(NATYPD),CLEB(NCLEB,2),
     &          THETAS(IRID,NFUND,NCELLD),SOCSCALE(NATYPD),
     &          RPAN_INTERVALL(0:NTOTD,NATYPD),CSCL(LMAXD1,NATYPD),
     &          RNEW(NTOTD*(NCHEBD+1),NATYPD),SOCSCL(LMAXD1,NATYPD),
     &          THETASNEW(NTOTD*(NCHEBD+1),NFUND,NCELLD),EFERMI,
     &          EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD),E1,E2,TK,
     &          VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD),
     &          ECORE(20,NPOTD),DRDIREL(IRMD,NATYPD),
     &          R2DRDIREL(IRMD,NATYPD),RMREL(IRMD,NATYPD),
     &          VTREL(IRMD,NATYPD),BTREL(IRMD,NATYPD),
     &          WLDAU(MMAXD,MMAXD,NSPIND,NATYPD),
     &          ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD),R_LOG
      double complex, intent(inout) :: EZ(IEMXD),WEZ(IEMXD),
     &          PHILDAU(IRMD,NATYPD)
      character*80, intent(inout) :: TMPDIR
      character*10, intent(inout) :: SOLVER
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
          
C
      NSRA    = t_params%NSRA    
      INS     = t_params%INS    
      NATYP   = t_params%NATYP  
      NAEZ    = t_params%NAEZ   
      NSPIN   = t_params%NSPIN  
      ICST    = t_params%ICST   
      IPAN    = t_params%IPAN   
      IRCUT   = t_params%IRCUT  
      KMROT   = t_params%KMROT  
      IQAT    = t_params%IQAT   
      CONC    = t_params%CONC   
      QMTET   = t_params%QMTET  
      QMPHI   = t_params%QMPHI  
      IDOLDAU = t_params%IDOLDAU
      LMAX    = t_params%LMAX   
      IRWS    = t_params%IRWS   
C ......................................................................
Consistency check 
C
      IF ( ( KREL.EQ.1 ) .AND. ( INS.NE.0 ) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
C     
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
C ......................................................................
      ALAT           = t_params%ALAT          
      ZAT            = t_params%ZAT           
      DRDI           = t_params%DRDI          
      RMESH          = t_params%R         
      A              = t_params%A             
      B              = t_params%B             
      IEND           = t_params%IEND          
      CLEB           = t_params%CLEB          
      ICLEB          = t_params%ICLEB         
      LOFLM          = t_params%LOFLM         
      JEND           = t_params%JEND          
      THETAS         = t_params%THETAS         
      IFUNM1         = t_params%IFUNM1         
      LMSP1          = t_params%LMSP1          
      NFU            = t_params%NFU            
      LLMSP          = t_params%LLMSP          
      LCORE          = t_params%LCORE          
      NCORE          = t_params%NCORE          
      NTCELL         = t_params%NTCELL        
      IRMIN          = t_params%IRMIN         
      ITITLE         = t_params%ITITLE        
      INTERVX        = t_params%INTERVX       
      INTERVY        = t_params%INTERVY       
      INTERVZ        = t_params%INTERVZ       
      NACLS1         = t_params%NACLS(1)        
      LLY            = t_params%LLY           
      SOCSCALE       = t_params%SOCSCALE      
      TMPDIR         = t_params%TMPDIR        
      ITMPDIR        = t_params%ITMPDIR       
      ILTMP          = t_params%ILTMP         
      NPAN_LOG       = t_params%NPAN_LOG      
      NPAN_EQ        = t_params%NPAN_EQ       
      NCHEB          = t_params%NCHEB         
      R_LOG          = t_params%R_LOG         
      NPAN_TOT       = t_params%NPAN_TOT      
      RNEW           = t_params%RNEW          
      RPAN_INTERVALL = t_params%RPAN_INTERVALL
      IPAN_INTERVALL = t_params%IPAN_INTERVALL
      THETASNEW      = t_params%THETASNEW     
      IF (KREL.EQ.1) THEN
         SOLVER = t_params%SOLVER
         SOCSCL = t_params%SOCSCL 
         CSCL   = t_params%CSCL   
      ENDIF
      IF ( IDOLDAU.EQ.1 ) THEN
         NTLDAU   = t_params%NTLDAU  
         ITLDAU   = t_params%ITLDAU  
         LOPT     = t_params%LOPT    
         UEFF     = t_params%UEFF    
         JEFF     = t_params%JEFF    
         EREFLDAU = t_params%EREFLDAU
      ENDIF
C ---------------------------------------------------------- energy_mesh
C
      IELAST     = t_params%IELAST    
      EZ         = t_params%EZ        
      WEZ        = t_params%WEZ      
      E1         = t_params%E1        
      E2         = t_params%E2        
      IESEMICORE = t_params%IESEMICORE
      NPOL       = t_params%NPOL      
      TK         = t_params%TK        
      IF ( NPOL.EQ.0 ) EFERMI = t_params%EFERMI
C ------------------------------------------------------ input_potential
C
      VINS  = t_params%VINS 
      VISP  = t_params%VISP 
      ECORE = t_params%ECORE
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL    
         DRDIREL   = t_params%DRDIREL  
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL     
         JWSREL    = t_params%JWSREL   
         IRSHIFT   = t_params%IRSHIFT  
         VTREL     = t_params%VTREL    
         BTREL     = t_params%BTREL    
      END IF
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
C ------------------------------------------------------------- itermdir
C
      IF ( OPT('ITERMDIR') ) THEN
         QMTET = t_params%QMTET
         QMPHI = t_params%QMPHI
      END IF
C ---------------------------------------------------------------- lda+u   
C
      IF ( IDOLDAU.EQ.1 ) THEN
         ITRUNLDAU = t_params%ITRUNLDAU
         WLDAU     = t_params%WLDAU
         ULDAU     = t_params%ULDAU
         PHILDAU   = t_params%PHILDAU
      END IF

      end subroutine get_params_1c



      subroutine get_params_2(t_params,KREL,NATYPD,IPAND,NPOTD,NATOMIMPD
     &    ,LMXSPD,NFUND,LMPOTD,NCELLD,IRMD,NEMBD1,NEMBD,
     &     IRMIND,NSRA,INS,NATYP,NAEZ,NSPIN,IPAN,IRCUT,LCORE,NCORE,LMAX,
     &     NTCELL,LPOT,LMPOT,NLBASIS,NRBASIS,NRIGHT,NLEFT,NATOMIMP,
     &     ATOMIMP,IMIX,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,KVMAD,KXC
     &    ,ICC,ISHIFT,IXIPOL,KFORCE,IFUNM,LMSP,IMT,IRC,IRMIN,IRWS,LLMSP,
     &     ITITLE,NFU,HOSTIMP,ILM,IMAXSH,IELAST,NPOL,NPNT1,NPNT2,NPNT3,
     &     ITSCF,SCFSTEPS,IESEMICORE,KAOEZ,IQAT,NOQ,LLY,NPOLSEMI,N1SEMI,
     &     N2SEMI,N3SEMI,ZREL,JWSREL,IRSHIFT,MIXING,LAMBDA_XC,A,B,THETAS
     &    ,DRDI,R,ZAT,RMT,RMTNEW,RWS,E1,E2,TK,ALAT,EFOLD,CHRGOLD,
     &     CMOMHOST,CONC,GSH,EBOTSEMI,EMUSEMI,TKSEMI,VINS,VISP,RMREL,
     &     DRDIREL,VBC,FSOLD,R2DRDIREL,ECORE,EZ,WEZ,TXC,
     &     LINTERFACE,LRHOSYM,NGSHD,NAEZD,IRID,NSPOTD,IEMXD)
      ! get relevant parameters from t_params
!     ..
      use mod_types, only: t_inc
      IMPLICIT NONE
      type(type_params), intent(in) :: t_params 
      
      integer, intent(in)  :: KREL,NATYPD,IPAND,NPOTD,NATOMIMPD,LMXSPD,
     &     NFUND,LMPOTD,NCELLD,IRMD,NEMBD1,IRMIND,NGSHD
     &    ,NAEZD,IRID,NSPOTD,IEMXD,NEMBD
     
      integer, intent(inout) :: NSRA,INS,NATYP,NAEZ,NSPIN,IPAN(NATYPD),
     &     IRCUT(0:IPAND,NATYPD),LCORE(20,NPOTD),NCORE(NPOTD),LMAX,
     &     NTCELL(NATYPD),LPOT,LMPOT,NLBASIS,NRBASIS,NRIGHT,NLEFT,
     &     NATOMIMP,ATOMIMP(NATOMIMPD),IMIX,ITDBRY,
     &     IRNS(NATYPD),KPRE,KSHAPE,KTE,KVMAD,KXC,ICC,ISHIFT,
     &     IXIPOL(NATYPD),KFORCE,IFUNM(NATYPD,LMXSPD),
     &     LMSP(NATYPD,LMXSPD),IMT(NATYPD),IRC(NATYPD),IRMIN(NATYPD),
     &     IRWS(NATYPD),LLMSP(NATYPD,NFUND),ITITLE(20,NPOTD),
     &     NFU(NATYPD),HOSTIMP(0:NATYPD),ILM(NGSHD,3),IMAXSH(0:LMPOTD),
     &     IELAST,NPOL,NPNT1,NPNT2,NPNT3,ITSCF,SCFSTEPS,IESEMICORE,
     &     KAOEZ(NATYPD,NAEZD+NEMBD),IQAT(NATYPD),NOQ(NAEZD),LLY,
     &     N1SEMI,N2SEMI,N3SEMI,ZREL(NATYPD),JWSREL(NATYPD),NPOLSEMI,
     &     IRSHIFT(NATYPD)
     
      double precision, intent(inout) :: MIXING,LAMBDA_XC,A(NATYPD),
     &     B(NATYPD),THETAS(IRID,NFUND,NCELLD),DRDI(IRMD,NATYPD),
     &     R(IRMD,NATYPD),ZAT(NATYPD),RMT(NATYPD),RMTNEW(NATYPD),
     &     RWS(NATYPD),E1,E2,TK,ALAT,EFOLD,CHRGOLD,QBOUND,FCM,
     &     CMOMHOST(LMPOTD,NEMBD1),CONC(NATYPD),GSH(NGSHD),EBOTSEMI,
     &     EMUSEMI,TKSEMI,VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     &     VISP(IRMD,NPOTD),RMREL(IRMD*KREL+(1-KREL),NATYPD),
     &     DRDIREL(IRMD*KREL+(1-KREL),NATYPD),VBC(2),FSOLD,
     &     R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),ECORE(20,NPOTD)

      double complex, intent(inout) :: EZ(IEMXD),WEZ(IEMXD)
      character*24, intent(inout) :: TXC(5)
      logical, intent(inout) :: LINTERFACE,LRHOSYM
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
 
      NSRA       = t_params%NSRA      
      INS        = t_params%INS       
      NATYP      = t_params%NATYP     
      NAEZ       = t_params%NAEZ      
      NSPIN      = t_params%NSPIN     
      IPAN       = t_params%IPAN      
      IRCUT      = t_params%IRCUT     
      LCORE      = t_params%LCORE     
      NCORE      = t_params%NCORE     
      NTCELL     = t_params%NTCELL    
      LMAX       = t_params%LMAX          
      LPOT       = t_params%LPOT         
      LMPOT      = t_params%LMPOT     
      NLBASIS    = t_params%NLBASIS    
      NRBASIS    = t_params%NRBASIS    
      NRIGHT     = t_params%NRIGHT     
      NLEFT      = t_params%NLEFT      
      LINTERFACE = t_params%LINTERFACE
      ATOMIMP    = t_params%ATOMIMP   
      NATOMIMP   = t_params%NATOMIMP  

C ......................................................................
Consistency check 
C
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
C ......................................................................
      IMIX      = t_params%IMIX     
      MIXING    = t_params%MIXING   
      QBOUND    = t_params%QBOUND   
      FCM       = t_params%FCM      
      ITDBRY    = t_params%ITDBRY   
      IRNS      = t_params%IRNS     
      KPRE      = t_params%KPRE     
      KSHAPE    = t_params%KSHAPE   
      KTE       = t_params%KTE      
      KVMAD     = t_params%KVMAD    
      KXC       = t_params%KXC      
      LAMBDA_XC = t_params%LAMBDA_XC
      TXC       = t_params%TXC      
      ICC       = t_params%ICC      
      ISHIFT    = t_params%ISHIFT   
      IXIPOL    = t_params%IXIPOL   
      LRHOSYM   = t_params%LRHOSYM  
      KFORCE    = t_params%KFORCE   
      A         = t_params%A        
      B         = t_params%B        
      DRDI      = t_params%DRDI     
      R         = t_params%R        
      THETAS    = t_params%THETAS   
      ZAT       = t_params%ZAT       
      IFUNM     = t_params%IFUNM    
      LMSP      = t_params%LMSP     
      RMT       = t_params%RMT      
      RMTNEW    = t_params%RMTNEW   
      RWS       = t_params%RWS      
      IMT       = t_params%IMT      
      IRC       = t_params%IRC      
      IRMIN     = t_params%IRMIN    
      IRWS      = t_params%IRWS     
      ITITLE    = t_params%ITITLE   
      LLMSP     = t_params%LLMSP    
      NFU       = t_params%NFU      
      HOSTIMP   = t_params%HOSTIMP  
      ALAT      = t_params%ALAT     
      KAOEZ     = t_params%KAOEZ     
      IQAT      = t_params%IQAT     
      NOQ       = t_params%NOQ      
      CONC      = t_params%CONC     
      GSH       = t_params%GSH      
      ILM       = t_params%ILM      
      IMAXSH    = t_params%IMAXSH   
      LLY       = t_params%LLY      
C ---------------------------------------------------------- energy_mesh
C
      IELAST     = t_params%IELAST    
      EZ         = t_params%EZ        
      WEZ        = t_params%WEZ       
      E1         = t_params%E1        
      E2         = t_params%E2        
      IESEMICORE = t_params%IESEMICORE
      FSOLD      = t_params%FSEMICORE
      NPOL       = t_params%NPOL      
      TK         = t_params%TK        
      NPNT1      = t_params%NPNT1     
      NPNT2      = t_params%NPNT2     
      NPNT3      = t_params%NPNT3     
      EBOTSEMI   = t_params%EBOTSEMI  
      EMUSEMI    = t_params%EMUSEMI   
      TKSEMI     = t_params%TKSEMI    
      NPOLSEMI   = t_params%NPOLSEMI  
      N1SEMI     = t_params%N1SEMI    
      N2SEMI     = t_params%N2SEMI     
      N3SEMI     = t_params%N3SEMI    
C ------------------------------------------------------ input_potential
C
      VINS  = t_params%VINS 
      VISP  = t_params%VISP 
      ECORE = t_params%ECORE
      VBC   = t_params%VBC  
      IF (KREL.EQ.1) THEN
         RMREL     = t_params%RMREL    
         DRDIREL   = t_params%DRDIREL  
         R2DRDIREL = t_params%R2DRDIREL
         ZREL      = t_params%ZREL     
         JWSREL    = t_params%JWSREL   
         IRSHIFT   = t_params%IRSHIFT  
      END IF
      ITSCF    = t_params%ITSCF   
      SCFSTEPS = t_params%SCFSTEPS
      EFOLD    = t_params%EFOLD   
      CHRGOLD  = t_params%CHRGOLD 
      CMOMHOST = t_params%CMOMHOST
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0

      end subroutine get_params_2
      
      
      
      
      subroutine save_emesh(IELAST,EZ,WEZ,E1,E2,IESEMICORE,FSEMICORE,
     &                NPOL,TK,
     &                NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,
     &                NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD,t_params)
      ! save information of energy mesh in t_params
      IMPLICIT NONE
      type(type_params), intent(inout) :: t_params
      integer, intent(in) :: IEMXD
      integer, intent(in) :: IELAST,IESEMICORE,NPOL,NPNT1,NPNT2,NPNT3,
     &                       NPOLSEMI,N1SEMI,N2SEMI,N3SEMI
      double precision, intent(in) :: E1,E2,FSEMICORE,TK,EBOTSEMI,
     &                                EMUSEMI,TKSEMI
      double complex, intent(in) :: EZ(IEMXD),WEZ(IEMXD)
     
      t_params%IELAST     = IELAST    
      t_params%EZ         = EZ        
      t_params%WEZ        = WEZ
      t_params%E1         = E1       
      t_params%E2         = E2        
      t_params%IESEMICORE = IESEMICORE
      t_params%FSEMICORE  = FSEMICORE 
      t_params%NPOL       = NPOL      
      t_params%TK         = TK        
      t_params%NPNT1      = NPNT1     
      t_params%NPNT2      = NPNT2     
      t_params%NPNT3      = NPNT3     
      t_params%EBOTSEMI   = EBOTSEMI  
      t_params%EMUSEMI    = EMUSEMI   
      t_params%TKSEMI     = TKSEMI    
      t_params%NPOLSEMI   = NPOLSEMI  
      t_params%N1SEMI     = N1SEMI    
      t_params%N2SEMI     = N2SEMI    
      t_params%N3SEMI     = N3SEMI    

      end subroutine save_emesh
      
     
     
      subroutine save_scfinfo(t_params,VINS,VISP,ECORE,VBC,RMREL,DRDIREL
     &                       ,R2DRDIREL,ZREL,JWSREL,IRSHIFT,VTREL,BTREL,
     &                        ITSCF,SCFSTEPS,EFOLD,CHRGOLD,CMOMHOST,KREL
     &                       ,IRMIND,IRMD,LMPOTD,NSPOTD,NATYPD,NPOTD,
     &                        NEMBD1)
      ! save information that is needed in next iteration and that is changeing, i.e. potential etc.
      IMPLICIT NONE
      type(type_params) , intent(inout) :: t_params
      integer, intent(in) :: KREL,IRMIND,IRMD,LMPOTD,NSPOTD,NATYPD,NPOTD
     &                      ,NEMBD1
      integer, intent(in) :: ITSCF,ZREL(NATYPD),IRSHIFT(NATYPD),SCFSTEPS
     &                      ,JWSREL(NATYPD)
      double precision, intent(in) :: VINS(IRMIND:IRMD,LMPOTD,NSPOTD),
     &               VISP(IRMD,NPOTD),VTREL(IRMD*KREL+(1-KREL),NATYPD),
     &               BTREL(IRMD*KREL+(1-KREL),NATYPD),EFOLD,CHRGOLD,
     &               CMOMHOST(LMPOTD,NEMBD1),ECORE(20,NPOTD),VBC(2),
     &               RMREL(IRMD*KREL+(1-KREL),NATYPD),
     &               DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &               R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD)
      
      t_params%VINS  = VINS 
      t_params%VISP  = VISP 
      t_params%ECORE = ECORE
      t_params%VBC   = VBC  
      IF (KREL.EQ.1) THEN
         t_params%RMREL     = RMREL    
         t_params%DRDIREL   = DRDIREL  
         t_params%R2DRDIREL = R2DRDIREL
         t_params%ZREL      = ZREL     
         t_params%JWSREL    = JWSREL   
         t_params%IRSHIFT   = IRSHIFT  
         t_params%VTREL     = VTREL    
         t_params%BTREL     = BTREL    
      END IF
      t_params%ITSCF    = ITSCF   
      t_params%SCFSTEPS = SCFSTEPS
      t_params%EFOLD    = EFOLD   
      t_params%CHRGOLD  = CHRGOLD 
      t_params%CMOMHOST = CMOMHOST
      end subroutine save_scfinfo
      
      
      
      subroutine save_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,
     &                        ESPV,ECORE,IDOLDAU,LOPT,EU,EDC,
     &                        CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,
     &                        KAPCORE,KREL,NATYPD,NPOTD,IRMD,LMPOTD,
     &                        LMAXD1)
      ! save density after it has been calculated in main1c, is further processed in main2
      IMPLICIT NONE
      type(type_params) , intent(inout) :: t_params
      integer, intent(in) :: KREL,NATYPD,NPOTD,IRMD,LMPOTD,LMAXD1
      integer, intent(in) :: IDOLDAU,LOPT(NATYPD),NKCORE(20,NATYPD),
     &     KAPCORE(20,NPOTD)
      double precision, intent(in) :: RHO2NS(IRMD,LMPOTD,NATYPD,2),
     &     R2NEF(IRMD,LMPOTD,NATYPD,2),RHOC(IRMD,NPOTD),DENEF,
     &     DENEFAT(NATYPD),ESPV(0:LMAXD1,NPOTD),ECORE(20,NPOTD),
     &     EDC(NATYPD),EU(NATYPD),CHRGSEMICORE,
     &     RHOORB(IRMD*KREL+(1-KREL),NATYPD),
     &     ECOREREL(KREL*20+(1-KREL),NPOTD)
      
      t_params%RHO2NS       = RHO2NS      
      t_params%R2NEF        = R2NEF       
      t_params%RHOC         = RHOC        
      t_params%DENEF        = DENEF       
      t_params%DENEFAT      = DENEFAT     
      t_params%ESPV         = ESPV        
      t_params%ECORE        = ECORE       
      t_params%IDOLDAU      = IDOLDAU     
      t_params%LOPT         = LOPT        
      t_params%EU           = EU          
      t_params%EDC          = EDC         
      t_params%CHRGSEMICORE = CHRGSEMICORE
      IF (KREL.EQ.1) THEN
         t_params%RHOORB   = RHOORB  
         t_params%ECOREREL = ECOREREL
         t_params%NKCORE   = NKCORE  
         t_params%KAPCORE  = KAPCORE 
      ENDIF
      
      end subroutine save_density
      
      
      
      subroutine read_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,
     &                        ESPV,ECORE,IDOLDAU,LOPT,EU,EDC,
     &                        CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,
     &                        KAPCORE,KREL,NATYPD,NPOTD,IRMD,LMPOTD,
     &                        LMAXD1)
      ! read density in main2 
      IMPLICIT NONE
      type(type_params) , intent(inout) :: t_params
      integer, intent(in) :: KREL,NATYPD,NPOTD,IRMD,LMPOTD,LMAXD1
      integer, intent(out) :: IDOLDAU,LOPT(NATYPD),NKCORE(20,NATYPD),
     &     KAPCORE(20,NPOTD)
      double precision, intent(out) :: RHO2NS(IRMD,LMPOTD,NATYPD,2),
     &     R2NEF(IRMD,LMPOTD,NATYPD,2),RHOC(IRMD,NPOTD),DENEF,
     &     DENEFAT(NATYPD),ESPV(0:LMAXD1,NPOTD),ECORE(20,NPOTD),
     &     EDC(NATYPD),EU(NATYPD),CHRGSEMICORE,
     &     RHOORB(IRMD*KREL+(1-KREL),NATYPD),
     &     ECOREREL(KREL*20+(1-KREL),NPOTD)
      
      RHO2NS       = t_params%RHO2NS      
      R2NEF        = t_params%R2NEF       
      RHOC         = t_params%RHOC        
      DENEF        = t_params%DENEF       
      DENEFAT      = t_params%DENEFAT     
      ESPV         = t_params%ESPV        
      ECORE        = t_params%ECORE       
      IDOLDAU      = t_params%IDOLDAU     
      LOPT         = t_params%LOPT        
      EU           = t_params%EU          
      EDC          = t_params%EDC         
      CHRGSEMICORE = t_params%CHRGSEMICORE
      
      IF (KREL.EQ.1) THEN
         RHOORB   = t_params%RHOORB  
         ECOREREL = t_params%ECOREREL
         NKCORE   = t_params%NKCORE  
         KAPCORE  = t_params%KAPCORE 
      ENDIF
      
      end subroutine read_density
      
      
      
      subroutine read_angles(t_params,NATYP,THETA,PHI)
      ! read nonco_angles
      use mod_types, only: t_inc
      use mod_mympi, only: myrank, master
      IMPLICIT NONE
      type(type_params), intent(inout) :: t_params
      integer, intent(in) :: NATYP
      double precision, intent(out) :: theta(natyp),phi(natyp)
      
      logical :: lread, lcheckangles
      integer :: i1
      double precision :: th1, ph1
      double precision, parameter :: PI=4.d0*datan(1.d0), eps=1d-5
      
      ! if executed first in wunfiles theta is not allocated, thus read angles from file
      if(.not.allocated(t_params%theta)) then
            
       THETA(:) = 0.D0
       PHI(:) = 0.D0
       LREAD = .FALSE.
       LCHECKANGLES = .false.
       INQUIRE(file='nonco_angle.dat',EXIST=LREAD)
       IF (LREAD) THEN
         OPEN(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
         DO I1 = 1,NATYP
           READ(10,*) th1, ph1
           IF((ABS(th1).LT.(PI+eps)   .AND. ABS(th1).GT.eps).OR.
     +        (ABS(ph1).LT.(2*PI+eps) .AND. ABS(ph1).GT.eps))
     &          LCHECKANGLES = .true.
           THETA(I1)=th1*(PI/180.0D0)
           PHI(I1)  =ph1*(PI/180.0D0)
         ENDDO
         CLOSE(10)
         IF ( LCHECKANGLES.and.((t_inc%i_write>0).or.(myrank==master)) )
     & write(1337,*) 'WARNING: Check if your
     & nonco_angels file is correct! Found only values that are
     & smaller than pi for theta and 2pi for phi, respectively. But
     & angles are given in degree (0...360)'
       write(1337,'(A)') '      I1  THETA[deg]  PHI[deg]'
       do i1=1,NATYP
         write(1337,'(I8,2F12.6)') i1,THETA(I1)*180D0/PI,
     +                                 PHI(I1)*180D0/PI
       end do!i1
       ENDIF ! LREAD
       
       !now save this also to t_params
       allocate(t_params%theta(natyp),t_params%phi(natyp))
       t_params%THETA = THETA
       t_params%PHI   = PHI
      
      else ! not first run: information saved in t_params
      
       THETA = t_params%THETA
       PHI   = t_params%PHI
       
      endif
      
      end subroutine read_angles
  

      end module
