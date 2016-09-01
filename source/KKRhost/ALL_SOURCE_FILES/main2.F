      MODULE MOD_MAIN2
      
            
      !use modulename, only rountine1
      
      implicit none
      
      contains
      
      subroutine main2()
      
      use mod_types, only: t_inc
      use mod_wunfiles
#ifdef CPP_TIMING
      use mod_timing
#endif
      
      implicit none
      
C     .. Parameters ..
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (KREL+KORBIT+1)*(LMAXD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD=(2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C parameter nembd1 avoids zero sized arrays.(2.1.01 R.Zeller)
      INTEGER NEMBD1
      PARAMETER (NEMBD1=NEMBD+1)
C     ..
      INTEGER IOBROY
      PARAMETER ( IOBROY = 20 )
      INTEGER NMVECMAX
      PARAMETER (NMVECMAX = 4)
C     ..
C     .. Local Scalars ..
      INTEGER IELAST,ICC,INS,LMAX,NATYP,
     &        NPNT1,NPNT2,NPNT3,NPOL,NSPIN,NSRA,IE,ICONT
      INTEGER ITSCF,SCFSTEPS
      INTEGER IPF,IMIX,ISHIFT,ITDBRY,KSHAPE  
      INTEGER KPRE,KTE,KVMAD,KXC,KFORCE,LSMEAR
      INTEGER LPOT,LMPOT,NAEZ
      INTEGER I,J,IPOT,ISPIN,I1,IH,IRC1,IRMIN1,IT,IO,LM,IR
      INTEGER NLBASIS,NRBASIS,NLEFT,NRIGHT
      INTEGER ATOMIMP(NATOMIMPD),NATOMIMP

C     ..
C     .. old/new FERMI energy
      DOUBLE PRECISION EFOLD,EFNEW
      DOUBLE PRECISION DENEF,E1,E2,TK
      DOUBLE PRECISION CHRGOLD,CHRGNT,E2SHIFT
      DOUBLE PRECISION DF
      DOUBLE PRECISION ALAT,FCM,MIX,MIXING,QBOUND,
     &                 FPI,PI,RFPI,RMSAVM,RMSAVQ,RMSAV0
      DOUBLE PRECISION SUM,RV
C     .. ITERMDIR variables
      INTEGER NK,NMVEC
      DOUBLE PRECISION ERRAVANG
C     .. semicore ..
      INTEGER IDOSEMICORE,IESEMICORE,NPOLSEMI,N1SEMI,N2SEMI,N3SEMI
      DOUBLE PRECISION FSEMICORE,EBOTSEMI,EMUSEMI,TKSEMI
      DOUBLE PRECISION FSOLD,CHRGSEMICORE
C     ..
      LOGICAL LRHOSYM,OPT,TEST,LINTERFACE
      CHARACTER*24 TXC(5)
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX EZ(IEMXD),DEZ(IEMXD),WEZ(IEMXD)
C ----------------------------------------------------------------------
C     ITERMDIR variables
      DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX)
      DOUBLE COMPLEX MVEVIEF(NATYPD,3,NMVECMAX)
C     ITERMDIR variables
C ----------------------------------------------------------------------
      DOUBLE PRECISION VBC(2)
C     .. input potential (spherical VISP, nonspherical VINS)
!       DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,NSPOTD),       
      DOUBLE PRECISION VISP(IRMD,NPOTD)
C     .. output potential (nonspherical VONS)
!       DOUBLE PRECISION VONS(IRMD,LMPOTD,NPOTD)
      DOUBLE PRECISION, allocatable :: VINS(:,:,:)
! C     .. output potential (nonspherical VONS)
      DOUBLE PRECISION, allocatable :: VONS(:,:,:)
C     ..
C ----------------------------------------------------------------------
C   ECOU(0:LPOTD,NATYPD)    ! Coulomb energy                    
C   EPOTIN(NATYPD),         ! energy of input potential (EPOTINB
C   ESPC(0:3,NPOTD),        ! energy single particle core       
C   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence    
C                           ! both changed for the relativistic 
C                           ! case
C   EXC(0:LPOTD,NATYPD),    ! E_xc
C ----------------------------------------------------------------------
      DOUBLE PRECISION ECOU(0:LPOTD,NATYPD),EPOTIN(NATYPD),      
     +                 ESPC(0:3,NPOTD),ESPV(0:LMAXD1,NPOTD),
     +                 EXC(0:LPOTD,NATYPD)
      DOUBLE PRECISION A(NATYPD),B(NATYPD),DRDI(IRMD,NATYPD)
      DOUBLE PRECISION RMT(NATYPD),RMTNEW(NATYPD),RWS(NATYPD)
      DOUBLE PRECISION R(IRMD,NATYPD),ECORE(20,NPOTD)
!       DOUBLE PRECISION THETAS(IRID,NFUND,NCELLD),ZAT(NATYPD)
      DOUBLE PRECISION ZAT(NATYPD)
      !allocatables:
      double precision, allocatable :: THETAS(:,:,:), thesme(:,:,:)
C     .. Dummy variables needed only in IMPURITY program 
      DOUBLE PRECISION VSPSMDUM(IRMD,NPOTD)
C ----------------------------------------------------------------------
C  R2NEF (IRMD,LMPOTD,NATYPD,2)  ! rho at FERMI energy
C  RHO2NS(IRMD,LMPOTD,NATYPD,2)  ! radial density
C   nspin=1            : (*,*,*,1) radial charge density
C   nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
C                               (*,*,*,2) rho(2) - rho(1) -> mag. moment
C  RHOC(IRMD,NPOTD)              ! core charge density
C ----------------------------------------------------------------------
      DOUBLE PRECISION RHOC(IRMD,NPOTD)
      double precision, allocatable :: R2NEF(:,:,:,:),
     +                                 RHO2NS(:,:,:,:)
C     ..
      DOUBLE PRECISION GSH(NGSHD)
      DOUBLE PRECISION RHOORB(IRMD*KREL + (1-KREL),NATYPD)
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 RMREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION QMTET(NAEZD),QMPHI(NAEZD)
      DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),NPOTD) 
C ----------------------------------------------------------------------
C  CMINST(LMPOTD,NATYPD)            ! charge moment of interstitial
C  CMOM(LMPOTD,NATYPD)              ! LM moment of total charge
C  CHRGATOM(NATYPD,
C           2*KREL+(1-KREL)*NSPIND) ! total charge per atom
C ----------------------------------------------------------------------
      DOUBLE PRECISION CMINST(LMPOTD,NATYPD),CMOM(LMPOTD,NATYPD),
     +                 CHRGATOM(NATYPD,2*KREL+(1-KREL)*NSPIND)
C     ,,
      DOUBLE PRECISION C00(LMPOTD),CMOMHOST(LMPOTD,NEMBD1)
      DOUBLE PRECISION CONC(NATYPD),DENEFAT(NATYPD)
C     .. FORCES
      DOUBLE PRECISION FLM(-1:1,NATYPD),FLMC(-1:1,NATYPD)
C ----------------------------------------------------------------------
C     ITERMDIR variables
      DOUBLE PRECISION QMGAM(NAEZD)
      DOUBLE PRECISION FACT(0:100)
      DOUBLE PRECISION MVGAM(NATYPD,NMVECMAX),MVPHI(NATYPD,NMVECMAX),
     &                 MVTET(NATYPD,NMVECMAX)
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U
      INTEGER IDOLDAU,LOPT(NATYPD)
      DOUBLE PRECISION EDC(NATYPD),EU(NATYPD)
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

!--------------
! Scale magn. part of xc-potential:
      DOUBLE PRECISION LAMBDA_XC, EXCDIFF   
      DOUBLE PRECISION EXCNM(0:LPOTD,NATYPD)   
      DOUBLE PRECISION, allocatable :: VXCM(:,:,:),VXCNM(:,:,:),
     +                                 RHO2NSNM(:,:,:,:)   
!------------ 


      INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)
      INTEGER NKCORE(20,NATYPD),KAPCORE(20,NPOTD)
      INTEGER IFUNM(NATYPD,LMXSPD),
     +        KAOEZ(NATYPD,NAEZD+NEMBD)
      INTEGER IXIPOL(NATYPD), ITOQ(NATYPD,NAEZD)
      INTEGER IMT(NATYPD),IPAN(NATYPD),IRC(NATYPD),
     +        IRCUT(0:IPAND,NATYPD),IRMIN(NATYPD),IRNS(NATYPD),
     +        IRWS(NATYPD),ITITLE(20,NPOTD)
      INTEGER LCORE(20,NPOTD),LCOREMAX(NATYPD),LLMSP(NATYPD,NFUND)
      INTEGER LMSP(NATYPD,LMXSPD),NCORE(NPOTD),NFU(NATYPD),
     +        NSHELL(0:NSHELD),NTCELL(NATYPD)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      INTEGER NOQ(NAEZD),IQAT(NATYPD)
C     ..
C     changes for impurity 20/02/2004 -- v.popescu according to 
C     ..                                 n.papanikolaou 
C     ..
      INTEGER HOSTIMP(0:NATYPD)
      DOUBLE PRECISION VINTERS(LMPOTD,NAEZD)
C     ..
      INTEGER LLY ! LLY <> 0 : apply Lloyd's formula
C     ..
C     .. External Subroutines ..
      EXTERNAL BRYDBM,CONVOL,DAXPY,DCOPY,ECOUB,EPATHTB,EPOTINB,ESPCB,
     &         ETOTB1,FORCE,FORCEH,FORCXC,MDIRNEWANG,MIXSTR,MTZERO,OPT,
     &         RELPOTCVT,RHOSYMM,RHOTOTB,RINIT,RITES,SCFITERANG,
     &         TEST,VINTERFACE,VINTRAS,VMADELBLK,VXCDRV
c
C     .. Intrinsic Functions ..
      INTRINSIC DABS,ATAN,DMIN1,DSIGN,SQRT,MAX,DBLE
C     ..
C declarations needed:
       DOUBLE PRECISION  AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
       INTEGER LRECABMAD,I2,IREC
       LOGICAL LPOTSYMM(NATYPD,LMPOTD)
c New variables for fixed Fermi energy calculations; noted as fxf
       DOUBLE PRECISION VMT_INIT(2),VSHIFT  ! fxf
C     ..
! for SIMULASA:
      INTEGER*4 IPOS,ILMP,IAS
! allocatable arrays
      integer :: ierr
C     ..

! allocations:
      allocate(VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VONS(IRMD,LMPOTD,NPOTD),
     &         VXCM(IRMD,LMPOTD,NPOTD),VXCNM(IRMD,LMPOTD,NPOTD),
     &         R2NEF(IRMD,LMPOTD,NATYPD,2),
     &         RHO2NS(IRMD,LMPOTD,NATYPD,2),
     &         RHO2NSNM(IRMD,LMPOTD,NATYPD,2),
     &         THETAS(IRID,NFUND,NCELLD), 
     &         THESME(IRID,NFUND,NCELLD), stat=ierr)
      if(ierr/=0) stop '[main2] Error allocating arrays'

C     ..................................................................
Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &     STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C --------------------------------------------------------------- input2
      call get_params_2(t_params,KREL,NATYPD,IPAND,NPOTD,NATOMIMPD
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
C -------------------------------------------------------------- density
C
      call read_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,
     &                  ESPV,ECORE,IDOLDAU,LOPT,EU,EDC,
     &                  CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,KAPCORE,
     &                  KREL,NATYPD,NPOTD,IRMD,LMPOTD,LMAXD1)
C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
C============================================================= CONSTANTS
      PI = 4.0D0*ATAN(1.0D0)
      FPI = 4.0D0*PI
      RFPI = SQRT(FPI)
      RMSAV0 = 1.0D10
C
C  Setting dummy argument LSMEAR to allow compatibility with IMPURITY
C
      LSMEAR=0
C
      ICONT = 1
      IPF = 1337
      NSPIN = 2*KREL + (1-KREL)*NSPIN
      IDOSEMICORE = 0
      IF ( OPT('SEMICORE') ) IDOSEMICORE = 1
C=======================================================================
C **********************************************************************
C **************************  ITERATION BEGIN  *************************
C **********************************************************************
C
      ITSCF = ITSCF + 1         ! initialised to 0 in main0
      t_inc%i_iteration = ITSCF
C          
      WRITE(1337,'(/,79(1H*))')
      WRITE(1337,'(19X,A,I3,A,I3,A)') '****** ITERATION : ',
     &                             ITSCF,' OUT OF ',SCFSTEPS,' ******'
      WRITE(1337,'(79(1H*),/)')
C
      IF (IMIX.GE.3) THEN
         OPEN (IOBROY,FILE='broy_io.unformatted',FORM='unformatted',
     &         STATUS='unknown')
         OPEN (IOBROY+2,FILE='broy_io2.unformatted',FORM='unformatted',
     &         STATUS='unknown')
      ENDIF  
C ----------------------------------------------------------------------
C
C the next four lines may not always work
C
      NSHELL(0) = NATYP
      DO I1 = 1,NATYP
         NSHELL(I1) = 1
      END DO
C ----------------------------------------------------------------------
C
C -->   determine total charge density expanded in spherical harmonics
C
      IF(TEST('flow    ')) write(1337,*) '>>> RHOTOTB'
      CALL RHOTOTB(IPF,NATYP,NAEZ,NSPIN,RHO2NS,RHOC,RHOORB,
     +             ZAT,DRDI,IRWS,IRCUT,
     +             LPOT,NFU,LLMSP,THETAS,NTCELL,KSHAPE,IPAN,CHRGNT,
     +             ITSCF,NSHELL,NOQ,CONC,KAOEZ,CHRGATOM)

      IF(TEST('flow    ')) write(1337,*) '<<< RHOTOTB'

      IF ( TEST('RHOVALTW') ) THEN !Bauer
        DO I1 = 1,NATYP
          open(unit=324234,file='out_rhotot')
          WRITE(324234,*) '#IATOM',I1
          write(324234,'(50000F)') RHO2NS(:,:,I1,1)
          IF (NSPIN==2) write(324234,'(50000F)') RHO2NS(:,:,I1,2)
        END DO
      END IF


C ----------------------------------------------------------------------
C
C --> determine new Fermi level due to valence charge up to 
C     old Fermi level E2 and density of states DENEF
C
      IF ( ITSCF.GT.1 .AND. CHRGNT*CHRGOLD .LT. 0.D0 .AND.
     +     ABS(CHRGNT) .GT. 5.D-2) THEN
         E2SHIFT = CHRGNT/(CHRGNT-CHRGOLD)*(E2-EFOLD)
      ELSE
         E2SHIFT = CHRGNT/DENEF
      END IF
C
      E2SHIFT = DMIN1(DABS(E2SHIFT),0.05D0)*DSIGN(1.0D0,E2SHIFT)
      EFOLD = E2
      CHRGOLD = CHRGNT
      IF (TEST('no-neutr').OR.OPT('no-neutr')) THEN
         WRITE(1337,*) 
     &        'test-opt no-neutr: Setting FERMI level shift to zero'
         E2SHIFT = 0.d0
      ENDIF
      IF (TEST('slow-neu').OR.OPT('slow-neu')) THEN
         WRITE(1337,*) 
     &        'test-opt slow-neu: FERMI level shift * STRMIX'
         E2SHIFT = E2SHIFT * MIXING
      ENDIF
      IF (ISHIFT.EQ.0) E2 = E2 - E2SHIFT
C ----------------------------------------------------------------------
      FSEMICORE=0d0
      IF ( IDOSEMICORE.EQ.1 ) THEN
C
C --> semicore treatment, recalculate the normalisation factor
C
         IF ( CHRGSEMICORE.LT.1D-10 ) CHRGSEMICORE = 1D-10
C
C    -> number of semicore bands
         I1 = NINT(CHRGSEMICORE)
         FSEMICORE = DBLE(I1)/CHRGSEMICORE * FSOLD
         WRITE(1337,'(6X,"< SEMICORE > : ",/,
     &        21X,"charge found in semicore :",F10.6,/,
     &        21X,"new normalisation factor :",F20.16,/)')
     &        CHRGSEMICORE,FSEMICORE
      END IF
C ----------------------------------------------------------------------
C
      !WRITE (6,FMT=9020) EFOLD,E2SHIFT
      WRITE (1337,FMT=9020) EFOLD,E2SHIFT
C
C --> divided by NAEZ because the weight of each atom has been already
C     taken into account in 1c
C
      WRITE (1337,FMT=9030) E2,DENEF/DBLE(NAEZ) 
      WRITE (6,FMT=9030) E2,DENEF/DBLE(NAEZ) 
      WRITE(1337,'(79(1H+),/)')
C ----------------------------------------------------------------------
      DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIN)
C ----------------------------------------------------------------------
C
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS ISPIN LOOP
      DO ISPIN = 1,NSPIN
C ----------------------------------------------------------------------
         IF (KTE.EQ.1) THEN
            DO I1 = 1,NATYP
               IPOT = (I1-1)*NSPIN + ISPIN
               ESPV(0,IPOT) = ESPV(0,IPOT) -
     +                        EFOLD*CHRGNT/DBLE(NSPIN*NAEZ)
            END DO
         END IF                 ! (KTE.EQ.1)
C ----------------------------------------------------------------------
C
C -->     get correct density
C
         IF (.NOT.(OPT('DECIMATE'))) THEN 
            DO I1 = 1,NATYP
               DO LM = 1,LMPOT
                  CALL DAXPY(IRC(I1),DF,R2NEF(1,LM,I1,ISPIN),1,
     &                       RHO2NS(1,LM,I1,ISPIN),1)
               END DO
            END DO
         END IF
C ----------------------------------------------------------------------
      END DO
C SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM ITERMDIR
C
      IF ((KREL.EQ.1).AND.(OPT('ITERMDIR'))) THEN
         MVEVI   = t_params%MVEVI
         MVEVIEF = t_params%MVEVIEF
C     
         CALL RINIT(NAEZD,QMGAM)
         DO I1 = 1,NAEZ
            ITOQ(1,I1) = KAOEZ(1,I1)
         END DO
         NK = 2 * LMAXD + 1
         NMVEC = 2
C
         FACT(0) = 1.0D0
         DO I = 1,100
            FACT(I) = FACT(I-1)*DBLE(I)
         END DO
C ----------------------------------------------------------------------
         IF (.NOT.(OPT('DECIMATE'))) THEN
            DO I1 = 1,NATYP
               DO LM = 1, NMVEC
                  DO IT = 1,3
                     MVEVI(I1,IT,LM) = MVEVI(I1,IT,LM)
     &                                 + E2SHIFT*MVEVIEF(I1,IT,LM)
                  END DO
               END DO
            END DO
         END IF
C ----------------------------------------------------------------------
         DO I1 = 1,NATYP
            CALL MDIRNEWANG(I1,NMVEC,MVEVI,MVPHI,MVTET,MVGAM,
     &                      NATYPD,LMAXD,NMVECMAX)
         END DO
C
         OPEN(67,FILE='itermdir.unformatted',FORM='unformatted')
         CALL SCFITERANG(ITSCF,ITOQ,FACT,
     &                   MVPHI,MVTET,MVGAM,QMPHI,QMTET,QMGAM,
     &                   NAEZ,NK,ERRAVANG,NAEZD,NATYPD,NMVECMAX,LMMAXD)
         t_params%MVEVI   = MVEVI
         t_params%MVEVIEF = MVEVIEF
      END IF
C MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
C
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C                           POTENTIAL PART
C
      IF (LRHOSYM) CALL RHOSYMM(LMPOT,NSPIN,1,NATYP,RHO2NS,IXIPOL,
     +                          IRWS,IRCUT,IPAN,KSHAPE)


C
      CMINST(:,:) = 0.d0
      CMOM(:,:) =0.d0
      VONS(:,:,:) = 0.d0
      CALL VINTRAS(CMOM,CMINST,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,
     +     R,DRDI,IRWS,IRCUT,IPAN,KSHAPE,NTCELL,ILM,IFUNM,IMAXSH,GSH,
     +     THETAS,LMSP)

      if ( TEST('vintrasp') ) then !Bauer
        open(unit=786785,file='test_vintraspot')
        do i1=1,nspin*natyp
          write(786785,*) '# atom/spin index ',i1
          write(786785,'(50000E)') vons(:,:,i1)
        end do !iatom
        close(786785)
      end if !



C
C =====================================================================
cfivos     IF ( .NOT.TEST('NoMadel ').AND. ( SCFSTEPS.GT.1 ) 
cfivos     &     .OR. (ICC .GT. 0 ) )THEN
C ---------------------------------------------------------------------
         IF ( LINTERFACE ) THEN  
            CALL VINTERFACE(CMOM,CMINST,LPOT,NSPIN,NAEZ,
     &                      NATYP,VONS,ZAT,R,IRWS,IRCUT,IPAN,KSHAPE,
     &                      NOQ,KAOEZ,IQAT,CONC,CHRGATOM(1,1),
     &                      ICC,HOSTIMP,NLBASIS,NLEFT,NRBASIS,NRIGHT,
     &                      CMOMHOST,CHRGNT,VINTERS)
C ---------------------------------------------------------------------
         ELSE
C ---------------------------------------------------------------------
            CALL VMADELBLK(CMOM,CMINST,LPOT,NSPIN,NAEZ,
     &                     VONS,ZAT,R,IRWS,IRCUT,IPAN,KSHAPE,
     &                     NOQ,KAOEZ,CONC,CHRGATOM(1,1),
     &                     ICC,HOSTIMP,VINTERS)
        END IF

      IF (OPT('KKRFLEX ')) THEN
        CALL WRITEKKRFLEX(NATOMIMP,NSPIN,IELAST,(LPOT+1)**2,LMMAXD, 
     +                  ALAT,NATYP,KSHAPE,VBC,ATOMIMP,HOSTIMP,NOQ,ZAT,
     +                  KAOEZ,CONC,CMOM,CMINST,VINTERS)
      END IF


C ---------------------------------------------------------------------
cfivos      END IF
C =====================================================================
        IF ( TEST('Vspher  ') ) VONS(1:IRMD,2:LMPOTD,1:NPOTD) = 0.D0
        if ( TEST('vpotout ') ) then !bauer
          open(unit=54633163,file='test_vpotout_inter')
          do i1=1,natyp*nspin
              write(54633163,*) '# atom ',i1
              write(54633163,'(50000E)') vons(:,:,i1)
          end do !iatom
          close(54633163)
        end if ! config_testflag('write_gmatonsite')



C
C =====================================================================
C
C --> write the CMOMS to a file
C
C In case of DECIMATION output, we store ONLY information connected
C with the effective CPA-medium (for the case of NO-CPA NAEZ=NATYP)
C hence the CMOMS are calculated site-dependent. In the same format
C are read in by <MAIN0> -- < CMOMSREAD >     v.popescu 01/02/2002
C 
      IF (OPT('deci-out').AND.(ITSCF.EQ.1)) THEN
         open(37, file='decifile', form='formatted', position='append')
         WRITE(37,1080) NAEZ,LMPOT
         DO IH=1,NAEZ
            WRITE(37,*) IH
            DO LM=1,LMPOT
               C00(LM) = 0.0D0
C ------------------------------------------ store the charge on SITE IH
               DO IO=1,NOQ(IH)
                  IT = KAOEZ(IO,IH)
                  C00(LM) = C00(LM) + CMOM(LM,IT) * CONC(IT)
                  IF (INS.NE.0) 
     &                 C00(LM) = C00(LM)+CMINST(LM,IT) * CONC(IT)
                  IF (LM.EQ.1) 
     &                 C00(1) = C00(1) - ZAT(IT)/RFPI*CONC(IT)
               END DO
C ----------------------------------------------------------------------
            END DO
            WRITE(37,1090) (C00(LM),LM=1,LMPOT)
         END DO  
      CLOSE(37)
      END IF
C =====================================================================
C
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES
C
      IF ( (KFORCE.EQ.1).AND.(KREL.NE.1) ) THEN
C ---------------------------------------------------------------------
         IF (INS.EQ.0) THEN
            CALL FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,
     &                  R,DRDI,IRWS,ZAT)
            CALL FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     &                 DRDI,IRWS)
C ---------------------------------------------------------------------
         ELSE
C ---------------------------------------------------------------------
            CALL FORCEH(CMOM,FLM,LPOT,NSPIN,1,NATYP,RHO2NS,VONS,
     +                  R,DRDI,IMT,ZAT)
            CALL FORCE(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                 DRDI,IMT)
         END IF
C ---------------------------------------------------------------------
      END IF
C
C Force Calculation stops here look after VXCDRV
C
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C       
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
C
      IF (KTE.EQ.1) THEN
         CALL ESPCB(ESPC,NSPIN,NATYP,ECORE,LCORE,LCOREMAX,NCORE)    ! single-particle core energy
C
         CALL EPOTINB(EPOTIN,NSPIN,NATYP,RHO2NS,VISP,R,DRDI,        ! "energy of the input potential"
     &                INS,IRMIN,IRWS,LPOT,VINS,IRCUT,IPAN,ZAT)      ! Int V(r) rho(r) d^3r
C
         CALL ECOUB(CMOM,ECOU,LPOT,NSPIN,NATYP,RHO2NS,VONS,ZAT,R,   ! Coulomb hartree energy
     &              DRDI,IRWS,KVMAD,KSHAPE,IRCUT,IPAN,IMAXSH,IFUNM,
     &              ILM,NTCELL,GSH,THETAS,LMSP)

      END IF
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C =====================================================================
      VXCM(:,:,:) = 0.D0 
      CALL VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NS,VXCM,R,DRDI,
     +            A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,ILM,IMAXSH,IFUNM,
     +            THETAS,LMSP)

        IF ( TEST('Vspher  ') ) VONS(1:IRMD,2:LMPOTD,1:NPOTD) = 0.D0

! Recalculate XC-potential with zero spin density for magn. moment scaling 
      VXCNM(:,:,:) = 0.D0                 ! Initialize
      EXCNM(:,:) = 0.D0
      IF (LAMBDA_XC.NE.1.D0.AND.NSPIN.EQ.2) THEN 
         RHO2NSNM(:,:,:,1) = RHO2NS(:,:,:,1) ! Copy charge density
         RHO2NSNM(:,:,:,2) = 0.D0            ! Set spin density to zero
         CALL VXCDRV(EXCNM,KTE,KXC,LPOT,NSPIN,1,NATYP,RHO2NSNM,VXCNM,
     +               R,DRDI,A,IRWS,IRCUT,IPAN,NTCELL,KSHAPE,GSH,
     +               ILM,IMAXSH,IFUNM,THETAS,LMSP)
!        compute the EXC-difference 
         EXCDIFF = 0.D0
         DO I1 = 1,NATYP
            DO LM = 0,LPOT
               EXCDIFF = EXCDIFF + EXC(LM,I1) - EXCNM(LM,I1)
            END DO
         END DO
         WRITE(1337,*) 'LAMBDA_XC=',LAMBDA_XC,'EXCDIF=',EXCDIFF
      ENDIF

      VONS(:,:,:) = VONS(:,:,:) + ! Add xc-potential with magn. part weighted by lambda_xc
     &           LAMBDA_XC*VXCM(:,:,:)  + (1.D0-LAMBDA_XC)*VXCNM(:,:,:)
      EXC(:,:) = LAMBDA_XC*EXC(:,:)     + (1.D0-LAMBDA_XC)*EXCNM(:,:)


        if ( TEST('vpotout ') ) then !bauer
          open(unit=57633263,file='test_vpotout_xc')
          do i1=1,natyp*nspin
              write(57633263,*) '# atom ',i1
              write(57633263,'(50000E)') vons(:,:,i1)
          end do !iatom
          close(57633263)
        end if ! config_testflag('write_gmatonsite')



C =====================================================================
C
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES
C
C Force calculation continues here 
C
      IF ( (KFORCE.EQ.1).AND.(KREL.NE.1) ) THEN
C ---------------------------------------------------------------------
         IF (KSHAPE.EQ.0) THEN
            CALL FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                  ALAT,DRDI,IRWS,0)
C ---------------------------------------------------------------------
         ELSE
C ---------------------------------------------------------------------
            CALL FORCXC(FLM,FLMC,LPOT,NSPIN,1,NATYP,RHOC,VONS,R,
     +                  ALAT,DRDI,IMT,0)
         END IF
C ---------------------------------------------------------------------
      END IF
C
C Force calculation ends
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C 
      IF (ISHIFT.EQ.2) THEN                                           ! fxf
c        OPEN (67,FILE='vmtzero_init',FORM='formatted')               ! fxf
c        READ (67,*) VMT_INIT(1)                                      ! fxf
c        CLOSE(67)                                                    ! fxf
         VMT_INIT(1) = 0.d0
         VMT_INIT(2) = VMT_INIT(1) ! read initial muffin-tin zero     ! fxf
c vmt_init is needed as a common reference for potential mixing if more 
c iterations are to be remembered, e.g., in Anderson mixing.
C =====================================================================
         ! Shift new potential to initial muffin-tin zero             ! fxf
         CALL MTZERO(LMPOT,NATYP,CONC,NSPIN,VONS,VMT_INIT,ZAT,R,DRDI, ! fxf
     &        IMT,IRCUT,IPAN,NTCELL,LMSP,IFUNM,THETAS,IRWS,E2SHIFT,   ! fxf
     &        ISHIFT,NSHELL,LINTERFACE)                               ! fxf
         OPEN (67,FILE='vmtzero',FORM='formatted')                    ! fxf
         WRITE (67,*) VMT_INIT(1)                                     ! fxf
         CLOSE(67)                                                    ! fxf
         ! Shift old potential to initial muffin-tin zero for correct mixing !fxf
         VSHIFT = - VBC(1)                                            ! fxf
         CALL POTENSHIFT(                                             ! fxf
     +     VISP,VINS,NATYP,NSPIN,IRCUT,IRC,IRMIN,NTCELL,              ! fxf
     +     IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,THETAS,THESME,             ! fxf
     +     RFPI,R,KSHAPE,VSHIFT)                                      ! fxf
      ELSE IF (ISHIFT.EQ.1) THEN
         ! Shift new potential to old MT-zero for correct mixing 
         ! (convolution with shapes is done later)
         DO ISPIN = 1,NSPIN
            DO IH = 1,NATYP
               IPOT = NSPIN * (IH-1) + ISPIN
               IRC1 = IRC(IH)
               VSHIFT = RFPI * VBC(ISPIN)
               VONS(1:IRC1,1,IPOT) = VONS(1:IRC1,1,IPOT) + VSHIFT
            ENDDO
         ENDDO

      ELSE  ! fxf
            ! Before fxf, only the following call was present.
         CALL MTZERO(LMPOT,NATYP,CONC,NSPIN,VONS,VBC,ZAT,R,DRDI,
     &        IMT,IRCUT,IPAN,NTCELL,LMSP,IFUNM,THETAS,IRWS,E2SHIFT
     &        ,ISHIFT,NSHELL,LINTERFACE)
      END IF                    ! fxf
C =====================================================================
C
      WRITE(1337,'(79(1H=),/)')
C ---------------------------------------------------------------------
C
C -->   convolute potential with shape function for next iteration
C


        if ( TEST('vpotout ') ) then !bauer
          open(unit=12633269,file='test_vpotout_shift')
          do i1=1,natyp*nspin
              write(12633269,*) '# atom ',i1
              write(12633269,'(50000E)') vons(:,:,i1)
          end do !iatom
          close(12633269)
        end if ! config_testflag('write_gmatonsite')


      IF (KSHAPE.NE.0) THEN
         DO ISPIN = 1,NSPIN
            DO I1 = 1,NATYP
               IPOT = NSPIN* (I1-1) + ISPIN

               if ( TEST('vpotout ') ) then !bauer
                  open(unit=12642269,file='test_convol')
                  write(12642269,*) '# atom ',i1

                  WRITE(12642269,*) IRCUT(1,I1),IRC(I1),
     &                 IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     &                 THETAS,ZAT(I1),RFPI,
     &                 R(:,I1),VONS(:,:,IPOT),LMSP
                  close(12642269)
               end if           ! config_testflag('write_gmatonsite')



               CALL CONVOL(IRCUT(1,I1),IRC(I1),NTCELL(I1),
     &                     IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     &                     THETAS,THESME,ZAT(I1),RFPI,
     &                     R(:,I1),VONS(:,:,IPOT),VSPSMDUM(1,1),LMSP)
            END DO
         END DO
      END IF


        if ( TEST('vpotout ') ) then !bauer
          open(unit=57633269,file='test_vpotout_conv')
          do i1=1,natyp*nspin
              write(57633269,*) '# atom ',i1
              write(57633269,'(50000E)') vons(:,:,i1)
          end do !iatom
          close(57633269)
        end if ! config_testflag('write_gmatonsite')


     
C ---------------------------------------------------------------------
C
C -->   symmetrisation of the potentials
C
C Keep only symmetric part of the potential
      if (test('potcubic')) then
      write(1337,*) 'Keeping only symmetric part of potential:'
      write(1337,*) 'Components L = 1, 11, 21, 25, 43, 47.'
      do ipot = 1,npotd
         do lm = 1, lmpotd
            if (lm.ne.1.and.lm.ne.11.and.lm.ne.21
     &                 .and.lm.ne.25.and.lm.ne.43.and.lm.ne.47) then
               do i = 1,irmd
                  VONS(i,lm,ipot) = 0.d0
               enddo
            endif
         enddo
      enddo
      endif

       IF(TEST('potsymm ')) THEN
C declarations needed:
C     DOUBLE PRECISION AVMAD(LMPOTD,LMPOTD),BVMAD(LMPOTD)
C     INTEGER LRECABMAD,I2,IREC
C     LOGICAL LPOTSYMM(NATYPD,LMPOTD)
C ..................
      LRECABMAD = WLENGTH*2*LMPOTD*LMPOTD + WLENGTH*2*LMPOTD
       OPEN (69,ACCESS='direct',RECL=LRECABMAD,FILE='abvmad.unformatted'
     &     ,FORM='unformatted')
       DO I1 = 1,NATYP
         DO LM = 1,LMPOTD
            LPOTSYMM(I1,LM)=.FALSE.
         END DO
         DO I2 = 1,NATYP
            IREC = I2 + NAEZ*(I1-1)
            READ (69,REC=IREC) AVMAD,BVMAD
            DO LM = 1,LMPOTD
               IF(ABS(BVMAD(LM)).GT.1D-10) LPOTSYMM(I1,LM)=.TRUE.
            END DO
         END DO
         DO LM = 1,LMPOTD
           IF(LPOTSYMM(I1,LM)) THEN
              WRITE(1337,*) 'atom ',I1,'lm = ',LM,' contribution used'
              ELSE
              DO ISPIN=1,NSPIN
                 IPOT = NSPIN* (I1-1) + ISPIN
                 DO IR = 1,IRMD
                    VONS(IR,LM,IPOT) = 0.0D0
                 END DO
              END DO
           ENDIF
         END DO
        END DO
         CLOSE (69)
          END IF


        if ( TEST('vpotout ') ) then !bauer
          open(unit=54633563,file='test_vpotout')
          do i1=1,natyp*nspin
              write(54633563,*) '# atom ',i1
              write(54633563,'(50000E)') vons(:,:,i1)
          end do !iatom
          close(54633563)
        end if ! config_testflag('write_gmatonsite')

! for simulasa:
      if(opt('simulasa')) then
         DO IAS = 1, NPOTD
          DO ILMP = 1, LMPOTD
            DO IPOS = 1, IRMD
              IF (ILMP .NE. 1) THEN
                VONS(IPOS, ILMP, IAS) = 0.D0
              ENDIF
            ENDDO
          ENDDO
         ENDDO
      endif

C ---------------------------------------------------------------------
C
C -->   final construction of the potentials (straight mixing)
C
      MIX = MIXING
      IF (TEST('alt mix ')) MIX = MIXING/DBLE(1+MOD(ITSCF,2))
      IF (TEST('spec mix')) 
     +       MIX = MIXING/
     +       (1.0D0 + 1.0D+3 * ABS(CHRGNT)/DBLE(NAEZ*NSPIN))
      write(1337,*) 'MIXSTR',MIX
      CALL MIXSTR(RMSAVQ,RMSAVM,INS,LPOT,LMPOT,0,NSHELL,
     +            1,NATYP,CONC,NSPIN,
     +            ITSCF,RFPI,FPI,IPF,
     +            MIX,
     +            FCM,IRC,IRMIN,R,DRDI,VONS,
     +            VISP,VINS,
     +            VSPSMDUM,VSPSMDUM,LSMEAR)
C
C                           POTENTIAL PART
C PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
C      
C ======================================================================
      IF (ITSCF.NE.1) RMSAV0 = 1.0d2*MAX(RMSAVQ,RMSAVM)
C
      WRITE(1337,FMT=9160) MIX
      WRITE(1337,'(79(1H=),/)')
C
C ======================================================================
      IF (MAX(RMSAVQ,RMSAVM).LT.QBOUND) THEN
         t_inc%i_iteration = t_inc%N_iteration
      ELSE
C ----------------------------------------------------------------------
C
C -->  potential mixing procedures: Broyden or Andersen updating schemes
C
         IF (IMIX.GE.3) 
     +       CALL BRYDBM(VISP,VONS,VINS,VSPSMDUM,VSPSMDUM,INS,
     +                   LMPOT,R,DRDI,MIX,CONC,
     +                   IRC,IRMIN,NSPIN,1,NATYP,ITDBRY,
     +                   IMIX,IOBROY,IPF,LSMEAR)
C----------------------------------------------------------------------
C
C -->    reset to start new iteration
C
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         DO I = 1,NSPIN*NATYP
            IT = I
            IF (NSPIN.EQ.2) IT = (I+1)/2
C
            IRC1 = IRC(IT)
            CALL DCOPY(IRC1,VONS(1,1,I),1,VISP(1,I),1)
C
            IF ( ( INS.NE.0 ).AND.( LPOT.GT.0 ) ) THEN
               IRMIN1 = IRMIN(IT)
               DO LM = 2,LMPOT
                  DO J = IRMIN1,IRC1
                     VINS(J,LM,I) = VONS(J,LM,I)
                  END DO
                  SUM = 0.0D0
                  DO IR = IRMIN1,IRC1
                     RV = VINS(IR,LM,I)*R(IR,IT)
                     SUM = SUM + RV*RV*DRDI(IR,IT)
                  END DO
                  IF ( SQRT(SUM).LT.QBOUND ) THEN
                     DO J = IRMIN1,IRC1
                        VINS(J,LM,I) = 0.0D0
                     END DO
                  END IF
               END DO
            END IF
         END DO
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END IF
C ======================================================================
C
C=======================================================================
C
      REWIND 11
C     
      EFNEW = E2
      IF (OPT('rigid-ef').OR.OPT('DECIMATE')) EFNEW = EFOLD
C
      IF (ISHIFT.EQ.2) THEN ! Shift mixed potential to new muffin-tin zero ! fxf
         VBC(1) = VBC(1) + E2SHIFT                          ! fxf
         VBC(2) = VBC(1)                                    ! fxf
         VSHIFT = VBC(1)                                    ! fxf
         CALL POTENSHIFT(                                   ! fxf
     +        VISP,VINS,NATYP,NSPIN,IRCUT,IRC,IRMIN,NTCELL, ! fxf
     +        IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,THETAS,THESME,! fxf
     +        RFPI,R,KSHAPE,VSHIFT)                         ! fxf
         WRITE(1337,*) 'New VMT ZERO:',VBC(1)                  ! fxf
      END IF                                                ! fxf
C
      CALL RITES(11,1,NATYP,NSPIN,ZAT,ALAT,RMT,RMTNEW,RWS,ITITLE,
     +           R,DRDI,VISP,IRWS,A,B,TXC,KXC,INS,IRNS,LPOT,VINS,
     +           QBOUND,IRC,KSHAPE,EFNEW,VBC,ECORE,LCORE,NCORE,
     +           ECOREREL,NKCORE,KAPCORE)
      CLOSE (11)
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE ENERGIES
C
      IF ((KTE.EQ.1 .AND. ICC.EQ.0) .or. OPT('KKRFLEX ')) 
     &     CALL ETOTB1(ECOU,EPOTIN,ESPC,ESPV,EXC,KPRE,LMAX,LPOT,
     &                 LCOREMAX,NSPIN,NATYP,NSHELL(1),CONC,IDOLDAU,
     &                 LOPT,EU,EDC)
C
C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC CONVERGENCY TESTS
C
      IF ( OPT('SEARCHEF') .AND. (DABS(E2SHIFT).LT.1D-8)) THEN
         t_inc%i_iteration = t_inc%N_iteration
         ICONT = 0
         GOTO 260
      END IF
C ----------------------------------------------------------------------
      IF (MAX(RMSAVQ,RMSAVM).LT.QBOUND) THEN
         WRITE(6,'(17X,A)') '++++++ SCF ITERATION CONVERGED ++++++'
         WRITE(6,'(79(1H*))')
         ICONT = 0
         GO TO 260
C ----------------------------------------------------------------------
      ELSE
C ----------------------------------------------------------------------
         IF (MAX(RMSAVQ,RMSAVM).GT.RMSAV0) THEN
            WRITE(6,*) 'ITERATION DIVERGED ---'
            ICONT = 0
            GO TO 260
         END IF
C ----------------------------------------------------------------------
         IF (ITSCF.GE.SCFSTEPS) THEN
         t_inc%i_iteration = t_inc%N_iteration
            WRITE(6,'(12X,A)')
     &           '++++++ NUMBER OF SCF STEPS EXHAUSTED ++++++'
            WRITE(6,'(79(1H*))')
            ICONT = 0
            GOTO 260
         END IF
      END IF
C ----------------------------------------------------------------------
C
 260  CONTINUE                  ! jump mark
C     
C **********************************************************************
C ************************    ITERATION END    *************************
C **********************************************************************
C
C ======================================================================
C
C --> update energy contour
C
      IF ( ICONT.EQ.1 ) THEN
         CALL EPATHTB(EZ,DEZ,E2,IELAST,IESEMICORE,IDOSEMICORE,E1,E2,TK,
     &                NPOL,NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,
     &                NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD)
         DO IE = 1,IELAST
            WEZ(IE) = -2.D0/PI*DEZ(IE)
            IF ( IE.LE.IESEMICORE ) WEZ(IE) = WEZ(IE)*FSEMICORE
         END DO
         WRITE(1337,'(79(1H=))')
      END IF
C ======================================================================
C
C --> convert VISP potential to the relativistic form VTREL,BTREL. 
C
      IF ( KREL.EQ.1 ) CALL RELPOTCVT(2,VISP,ZAT,R,DRDI,IRCUT,
     &                 VTREL,BTREL,ZREL,RMREL,JWSREL,DRDIREL,R2DRDIREL,
     &                 IRSHIFT,IPAND,IRMD,NPOTD,NATYPD)
C
C ======================================================================
C =             write out information for the next iteration           =
C ======================================================================
C
C ------------------------------------------------------ new_energy_mesh
C
      call save_emesh(IELAST,EZ,WEZ,E1,E2,IESEMICORE,FSEMICORE,NPOL,TK,
     &                NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,
     &                NPOLSEMI,N1SEMI,N2SEMI,N3SEMI,IEMXD,t_params) 
C ----------------------------------------------------- output_potential
C
      call save_scfinfo(t_params,VINS,VISP,ECORE,VBC,RMREL,DRDIREL
     &                 ,R2DRDIREL,ZREL,JWSREL,IRSHIFT,VTREL,BTREL,
     &                  ITSCF,SCFSTEPS,EFOLD,CHRGOLD,CMOMHOST,KREL
     &                 ,IRMIND,IRMD,LMPOTD,NSPOTD,NATYPD,NPOTD,
     &                  NEMBD1)
C ======================================================================
 9020 FORMAT ('                old',
     &     ' E Fermi ',F14.10,' Delta E_F = ',E16.8)
 9030 FORMAT ('                new',
     &     ' E FERMI ',F14.10,'  DOS(E_F) = ',F12.6)
 9100 FORMAT (19X,' TIME IN ITERATION : ',f9.2,/,79('*'))
 9160 FORMAT(20X,'mixing factor used : ',1P,D12.2)
 1080       FORMAT('CMOMC',2I6)
 1090       FORMAT(4D22.14)
 
      deallocate(thetas, thesme, vins, vons, vxcm, vxcnm, stat=ierr)
      if(ierr/=0) stop '[main2] Error deallocating thetas etc.'
      deallocate(R2NEF, RHO2NS, RHO2NSNM, stat=ierr)
      if(ierr/=0) stop '[main2] Error deallocating r2nef etc.'
 
 
      END SUBROUTINE !main2




c********************************************************************
      SUBROUTINE POTENSHIFT(VISP,VINS,NATYP,NSPIN,
     +     IRCUT,IRC,IRMIN,NTCELL,IMAXSH,ILM,IFUNM,LMSP,LMPOT,GSH,
     +     THETAS,THESME,RFPI,RMESH,KSHAPE,VSHIFT)
      implicit none
c Adds a constant (=VSHIFT) to the potentials of atoms
c
c Parameters:
      include 'inc.p'
      INTEGER NPOTD,LMPOTD,LMXSPD,IRMIND
      PARAMETER (NPOTD=NSPIND*NATYPD,LMPOTD= (LPOTD+1)**2
     &     ,LMXSPD= (2*LPOTD+1)**2,IRMIND=IRMD-IRNSD)
c Input
      INTEGER KSHAPE,LMPOT,NATYP,NSPIN
      INTEGER IRCUT(0:IPAND,NATYPD),IRC(NATYPD),NTCELL(NATYPD)
     &     ,IMAXSH(0:LMPOTD),ILM(NGSHD,3),IFUNM(NATYPD,LMXSPD)
     &     ,LMSP(NATYPD,LMXSPD),IRMIN(NATYPD)
      REAL*8 GSH(NGSHD),THETAS(IRID,NFUND,NCELLD),RFPI
     &     ,RMESH(IRMD,NATYPD)
      REAL*8 THESME(IRID,NFUND,NCELLD)
      REAL*8 VSHIFT


c Input/Output:
      REAL*8 VISP(IRMD,NPOTD),VINS(IRMIND:IRMD,LMPOTD,NSPOTD)

c Inside
      REAL*8 PSHIFTLMR(IRMD,LMPOTD),PSHIFTR(IRMD)
      INTEGER ISPIN,IH,IPOT,IR,LM,IMT1,IRC1,IRMIN1


      DO IH = 1,NATYP

         IMT1 = IRCUT(1,IH)
         IRC1 = IRC(IH)
         IRMIN1 = IRMIN(IH)

         DO ISPIN = 1,NSPIN

            WRITE (1337,*) 'SHIFTING OF THE POTENTIALS OF ATOM',IH,
     &           ' BY', VSHIFT, 'RY.'
            IPOT = NSPIN * (IH-1) + ISPIN

            CALL RINIT(IRMD*LMPOTD,PSHIFTLMR)
            CALL RINIT(IRMD,PSHIFTR)
            DO IR = 1,IRC1
               PSHIFTLMR(IR,1) = VSHIFT
            ENDDO

            IF (KSHAPE.EQ.0) THEN ! ASA

               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               END DO

            ELSE                ! Full-potential


               CALL CONVOL(IMT1,IRC1,NTCELL(IH),
     &              IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
     &              THETAS,THESME,0.d0,RFPI,
     &              RMESH(1,IH),PSHIFTLMR,PSHIFTR,LMSP)


               DO IR = 1,IRC1
                  VISP(IR,IPOT) = VISP(IR,IPOT) + PSHIFTLMR(IR,1)
               ENDDO


               DO LM = 2,LMPOT
                  DO IR = IRMIN1,IRC1
                VINS(IR,LM,IPOT)=VINS(IR,LM,IPOT)+PSHIFTLMR(IR,LM)*RFPI
                  ENDDO
               ENDDO

            END IF              ! (KSHAPE.EQ.0)


         END DO

      END DO



      END SUBROUTINE !potenshift
      
      END MODULE

