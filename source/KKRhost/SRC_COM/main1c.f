      MODULE MOD_MAIN1C

      
      !use modulename, only rountine1
      
      implicit none
      
      contains
      
      subroutine main1c()


CMPI  include 'mpif.h'
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
C     .. Parameters ..
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER MMAXD
      PARAMETER (MMAXD=2*LMAXD+1)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NRMAXD
      PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
C     .. itermdir parameter
      INTEGER NMVECMAX
      PARAMETER (NMVECMAX = 4)
C     .. 
      INTEGER LRECTMT,LRECTMT2
      PARAMETER (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
!     ..
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.D0,0.D0))
!     ..
      INTEGER LLY ! LLY <> 0 : apply Lloyd's formula

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,DENEF,E1,E2,PI,TK,EFERMI,CHRGSEMICORE
      INTEGER I1,ICELL,ICST,IE,IELAST,IEND,INS,INTERVX,INTERVY,INTERVZ,
     +        IPOT,IPOT1,IR,IS,ISPIN,L,LM,NACLS1,
     +        NATYP,NAEZ,NPOL,NSPIN,NSRA,LMAX,LMAXP1
      INTEGER KMROT,IQ,NSPINPOT,IHOST,ITMPDIR,ILTMP,IESEMICORE
      CHARACTER*5 TEXTNS
      CHARACTER*10 SOLVER
      CHARACTER*80 TMPDIR
      LOGICAL ITERMVDIR,LDORHOEF,LMOMVEC,LREAD
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD),
     &                 EDC(NATYPD),EU(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION WLDAUOLD(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
      DOUBLE COMPLEX DENMATC(MMAXD,MMAXD,NPOTD)
C for new spin-orbit solver
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION R_LOG,THETA(NATYPD),PHI(NATYPD)
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     +                 RNEW(NRMAXD,NATYPD),
     +                 VINSNEW(NRMAXD,LMPOTD,NSPOTD),
     +                 THETASNEW(NRMAXD,NFUND,NCELLD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C     ..
C     .. Local Arrays ..
CMPI  INTEGER MYRANK,NROFNODES,IERR,IDIM
CMPI  COMMON /MPI/MYRANK,NROFNODES
CMPI  DOUBLE COMPLEX WORK(IRMD,LMPOTD,NATYPD,2)
CMPI  EXTERNAL MPI_ALLREDUCE,MPI_COMM_RANK,MPI_COMM_SIZE,
CMPI +         MPI_FINALIZE,MPI_INIT,ZCOPY,DCOPY
CT3E  INTEGER ICLKTCK,ICLOCKS,MEND,MSTART
CT3E  DOUBLE PRECISION TIME1,TIME2,TIME3,TIME4
      DOUBLE COMPLEX DEN(0:LMAXD1,IEMXD,NPOTD),EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX DEN1(0:LMAXD1,IEMXD,2)
      DOUBLE COMPLEX CDOSAT0(IEMXD),CDOSAT1(IEMXD),CSUM                       ! LLY Lloyd
      DOUBLE COMPLEX CDOS0(IEMXD),CDOS_LLY(IEMXD,NSPIND)                      ! LLY Lloyd
      REAL*8 EREAD,CHARGE_LLY(NSPIND)                                         ! LLY Lloyd
      REAL*8 RENORM_AT(2,NATYPD)  ! 1: charge renormalization per atom        ! LLY Lloyd
                                  ! 2: spin moment renormalization per atom   ! LLY Lloyd
      DOUBLE PRECISION A(NATYPD),B(NATYPD),CHARGE(0:LMAXD1,NATYPD,2),
     +               CLEB(NCLEB,2),DRDI(IRMD,NATYPD),RMESH(IRMD,NATYPD),
     +                 DOSTOT(0:LMAXD1,2),ECORE(20,NPOTD),
     +                 ESPV(0:LMAXD1,NPOTD),DENEFAT(NATYPD)
      DOUBLE PRECISION ESPV1(0:LMAXD1,2)
C ----------------------------------------------------------------------
C  R2NEF (IRMD,LMPOTD,NATYPD,2)  ! rho at FERMI energy
C  RHO2NS(IRMD,LMPOTD,NATYPD,2)  ! radial density
C   nspin=1            : (*,*,*,1) radial charge density
C   nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
C                               (*,*,*,2) rho(2) - rho(1) -> mag. moment
C  RHOC(IRMD,NPOTD)              ! core charge density
C ----------------------------------------------------------------------
      DOUBLE PRECISION R2NEF(IRMD,LMPOTD,NATYPD,2),
     +                 RHO2NS(IRMD,LMPOTD,NATYPD,2),RHOC(IRMD,NPOTD),
     +                 RHO2N1(IRMD,LMPOTD,NPOTD),
     +                 RHO2N2(IRMD,LMPOTD,NPOTD),
     +                 RHO2M1(IRMD,LMPOTD,4),
     +                 RHO2M2(IRMD,LMPOTD,4)
      DOUBLE PRECISION VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD),
     +                 THETAS(IRID,NFUND,NCELLD),
     +                 ZAT(NATYPD),CONC(NATYPD),
     &                 SOCSCALE(NATYPD)
C----------------------------------- orbital magnetic moment
C     attention: muorb second index means both spins and total
C----------------------------------- orbital density
      DOUBLE PRECISION  MUORB(0:LMAXD1+1,3,NATYPD)
      DOUBLE PRECISION  RHOORB(IRMD*KREL + (1-KREL),NATYPD)
C---------------------------------------------------------------      
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 RMREL(IRMD*KREL+(1-KREL),NATYPD)
      INTEGER JWSREL(NATYPD),ZREL(NATYPD),IRSHIFT(NATYPD)
C ===================================================================
C   RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
C   SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
C
      DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),NPOTD) 
      INTEGER NKCORE(20,NATYPD),KAPCORE(20,NPOTD)
C ===================================================================
      INTEGER ICLEB(NCLEB,4),IFUNM1(LMXSPD,NATYPD),IPAN(NATYPD),
     +        IRCUT(0:IPAND,NATYPD),ITITLE(20,NPOTD),IRMIN(NATYPD),
     +        LMSP1(LMXSPD,NATYPD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LCORE(20,NPOTD),NCORE(NPOTD),LOFLM(LM2D),NTCELL(NATYPD),
     +        IRWS(NATYPD),NFU(NATYPD),LLMSP(NATYPD,NFUND)
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*7 TEXTS(3)
C     .. itermdir variables 
      DOUBLE PRECISION QMTET(NAEZD),QMPHI(NAEZD)
      INTEGER IQAT(NATYPD),NQDOS,LMMAXSO
      PARAMETER (LMMAXSO = 2*LMMAXD)
      DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX),         ! OUTPUT
     &               MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX) 
      DOUBLE COMPLEX MVEVIEF(NATYPD,3,NMVECMAX)        ! OUTPUT
      DOUBLE COMPLEX MVEVIL1(0:LMAXD,3,NMVECMAX)
      DOUBLE COMPLEX MVEVIL2(0:LMAXD,3,NMVECMAX) ! WORK ARRAYS
C     .. qdos and lmlm-dos ..
      CHARACTER*8 QDOSOPT
      DOUBLE COMPLEX DF(IEMXD),GFLLE(LMMAXSO,LMMAXSO,IEMXD,100)
      INTEGER IREC,LM1,LM2
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(32),TESTC(32)
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC
      COMMON /TESTC/TESTC
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,RHOCORE,RHOVAL,WMATLDAU,WRLDAUPOT,WRLDOS,WRMOMS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,DREAL
C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
C     ..
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/'spin dn','spin up','       '/
      DATA TEXTNS/' ns ='/
      DATA LDORHOEF/.TRUE./
      DATA IHOST / 1 /          ! this is the host program
C     ..
C
Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C
CT3E  CALL SYSTEM_CLOCK(MSTART)
CMPI  CALL MPI_INIT(IERR)
CMPI  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
CMPI  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)
CMPI  IF(MYRANK.NE.0) OPEN (6,STATUS='scratch',FORM='formatted')
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C -------------------------------------------------------------- input1c
C
      OPEN (67,FILE='input1c.unformatted',FORM='unformatted')
      READ (67) NSRA,INS,NATYP,NAEZ,NSPIN,ICST,IPAN,IRCUT,KMROT,IQAT,
     &          CONC,QMTET,QMPHI,IDOLDAU,LMAX,IRWS
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
      READ (67) ALAT,ZAT,DRDI,RMESH,A,B,IEND,CLEB,ICLEB,LOFLM,JEND,
     &    THETAS,IFUNM1,LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,IRMIN,ITITLE,
     &          INTERVX,INTERVY,INTERVZ,NACLS1,TESTC,OPTC,LLY,SOCSCALE
      READ (67) TMPDIR,ITMPDIR,ILTMP
      READ (67) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      READ (67) RNEW,RPAN_INTERVALL,IPAN_INTERVALL,THETASNEW
      IF (KREL.EQ.1) READ (67) SOLVER,SOCSCL,CSCL
      IF ( IDOLDAU.EQ.1 ) READ(67) NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU
      CLOSE (67)
C ---------------------------------------------------------- energy_mesh
C
      OPEN (67,FILE='energy_mesh',FORM='unformatted')
      READ (67) IELAST,EZ,WEZ,E1,E2,IESEMICORE
      READ (67) NPOL,TK
      IF ( NPOL.EQ.0 ) READ(67) EFERMI
      CLOSE (67)
C ------------------------------------------------------ input_potential
C
      OPEN (67,FILE='input_potential',FORM='unformatted')
      READ (67) VINS,VISP,ECORE
      IF (KREL.EQ.1) THEN
         READ (67) RMREL,DRDIREL,R2DRDIREL
         READ (67) ZREL,JWSREL,IRSHIFT
         READ (67) VTREL,BTREL
      END IF
      CLOSE (67)
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
C ------------------------------------------------------------- itermdir
C
      IF ( OPT('ITERMDIR') ) THEN
         OPEN (67,FILE='itermdir.unformatted',FORM='unformatted')
         READ (67) QMTET,QMPHI
         CLOSE (67)
      END IF
C ---------------------------------------------------------------- lda+u   
C
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         READ (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C     
CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME1 = (MEND-MSTART)
CT3E  TIME2 = (MEND-MSTART)
CT3E  CALL SYSTEM_CLOCK(MSTART)
CMPI  CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME3 = (MEND-MSTART)
CT3E  CALL SYSTEM_CLOCK(MSTART)
C
C ======================================================================
C


      CALL CINIT(IEMXD*(LMAXD+2)*NPOTD,DEN)
      DENEF = 0.0D0
      CALL RINIT(NATYPD,DENEFAT)
C
      ITERMVDIR = OPT('ITERMDIR') 
      LMOMVEC = ( ITERMVDIR .OR. ( KMROT.NE.0 ) )
      PI = 4.0D0*ATAN(1.0D0)
      NSPINPOT = KREL*2 + (1-KREL)*NSPIN
C
C --> no need to calculate charge correction if no host program, if 
C     decimation or if no energy contour
C
      LDORHOEF = (IHOST.EQ.1) .AND. (.NOT.OPT('DECIMATE')) 
     &           .AND. (NPOL.NE.0)
C
C LDA+U 
      IF ( IDOLDAU.EQ.1 ) CALL CINIT(MMAXD*MMAXD*NPOTD,DENMATC(1,1,1))
C LDA+U 
C
      CALL OPENDAFILE(69,'gmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)

!    write parameters file that contains passed parameters for further treatment of gflle
      IF (OPT('lmlm-dos')) THEN                                          ! lmlm-dos
        QDOSOPT = 'n'                                                    ! lmlm-dos
        IF (OPT('qdos    ')) THEN                                        ! lmlm-dos qdos
           QDOSOPT = 'y'                                                 ! lmlm-dos qdos
        ENDIF                                                            ! lmlm-dos qdos
        OPEN(67,FORM='formatted',FILE='parameters.gflle')                ! lmlm-dos
        DF(:)=WEZ(:)/DBLE(NSPIN)                                         ! lmlm-dos
        WRITE(67,*) IELAST,IEMXD,NATYP,NSPIN,LMAX,QDOSOPT,DF(1:IELAST),  ! lmlm-dos
     &              EZ(1:IELAST),KORBIT                                  ! lmlm-dos
        CLOSE(67)                                                        ! lmlm-dos
      ENDIF  ! OPT('lmlm-dos')                                           ! lmlm-dos

! -------------------------------------------------------------------------! LLY Lloyd
      IF (LLY.NE.0) THEN                                                   ! LLY Lloyd
                                                                           ! LLY Lloyd
! Calculate free-space contribution to dos                                 ! LLY Lloyd
         CDOS0(1:IEMXD) = CZERO                                            ! LLY Lloyd
         DO I1 = 1,NAEZ                                                    ! LLY Lloyd
            CDOSAT0(1:IEMXD) = CZERO                                       ! LLY Lloyd
            ICELL = NTCELL(I1)                                             ! LLY Lloyd
            DO IE = 1,IELAST                                               ! LLY Lloyd
               CALL RHOVAL0(                                               ! LLY Lloyd
     &           EZ(IE),WEZ(IE),                                           ! LLY Lloyd
     &           DRDI(1,I1),RMESH(1,I1),                                   ! LLY Lloyd
     &           IPAN(I1),IRCUT(0,I1),                                     ! LLY Lloyd
     &           THETAS(1,1,ICELL),LMAX,                                   ! LLY Lloyd
     &           CDOSAT0(IE),CDOSAT1(IE))                                  ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
            CDOS0(1:IEMXD) = CDOS0(1:IEMXD) + CDOSAT0(1:IEMXD)             ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CDOS0(:) = -CDOS0(:) / PI                                         ! LLY Lloyd
                                                                           ! LLY Lloyd
         OPEN(701,FILE='freedos.dat',FORM='FORMATTED')                     ! LLY Lloyd
         DO IE = 1,IELAST                                                  ! LLY Lloyd
            WRITE(701,FMT='(10E16.8)') EZ(IE),CDOS0(IE)                    ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
                                                                           ! LLY Lloyd
         CDOS_LLY(1:IEMXD,1:NSPIND) = CZERO                                ! LLY Lloyd
         OPEN (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')               ! LLY Lloyd
         DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
            DO IE = 1,IELAST                                               ! LLY Lloyd
               READ(701,FMT='(10E25.16)') EREAD,CDOS_LLY(IE,ISPIN)         ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
                                                                           ! LLY Lloyd
         ! Add free-space contribution cdos0                               ! LLY Lloyd
         DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
            CDOS_LLY(1:IEMXD,ISPIN) = CDOS_LLY(1:IEMXD,ISPIN) +            ! LLY Lloyd
     &                                CDOS0(1:IEMXD)                       ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
                                                                           ! LLY Lloyd
         CHARGE_LLY(1:NSPIND) = 0.D0                                       ! LLY Lloyd
         DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
            CSUM = CZERO                                                   ! LLY Lloyd
            DO IE = 1,IELAST                                               ! LLY Lloyd
               CSUM = CSUM + CDOS_LLY(IE,ISPIN) * WEZ(IE)                  ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
            CHARGE_LLY(ISPIN) = -DIMAG(CSUM) * PI / NSPINPOT               ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
                                                                           ! LLY Lloyd
                                                                           ! LLY Lloyd
         OPEN (701,FILE='cdos_lloyd.dat',FORM='FORMATTED')                 ! LLY Lloyd
         DO ISPIN=1,NSPIN                                                  ! LLY Lloyd
            DO IE=1,IELAST                                                 ! LLY Lloyd
               WRITE(701,FMT='(10E16.8)') EZ(IE),CDOS_LLY(IE,ISPIN)        ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
         WRITE(*,*) 'Valence charge from Lloyds formula:',                 ! LLY Lloyd
     &        (CHARGE_LLY(ISPIN),ISPIN=1,NSPIN)                            ! LLY Lloyd
         ENDIF                                                             ! LLY Lloyd
                                                                           ! LLY Lloyd
! -------------------------------------------------------------------------! LLY Lloyd


      IF (.NOT.OPT('NEWSOSOL')) THEN
C
C ================================================================ NATYP
      DO I1 = 1,NATYP
C ----------------------------------------------------------------- SPIN
         IQ = IQAT(I1)
         DO ISPIN = 1,NSPIN
            ICELL = NTCELL(I1)
            IPOT = (I1-1) * NSPINPOT + ISPIN
            IPOT1 = (I1-1) * NSPINPOT + 1

            CALL RHOVAL(IHOST,LDORHOEF,ICST,INS,IELAST,
     &           NSRA,ISPIN,NSPIN,NSPINPOT,I1,EZ,WEZ,
     &           DRDI(1,I1),RMESH(1,I1),
     &           VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),VISP(1,IPOT),
     &           ZAT(I1),IPAN(I1),IRCUT(0,I1),IRMIN(I1),
     &           THETAS(1,1,ICELL),IFUNM1(1,ICELL),LMSP1(1,ICELL),
     &           RHO2N1(1,1,IPOT1),RHO2N2(1,1,IPOT1),
     &           RHOORB(1,I1),DEN(0,1,IPOT),
     &           MUORB(0,1,I1),ESPV(0,IPOT1),
     &           CLEB,LOFLM,ICLEB,IEND,JEND,SOLVER,
     &           SOCSCL(1,KREL*I1+(1-KREL)),
     &           CSCL(1,KREL*I1+(1-KREL)),
     &           VTREL(1,I1),BTREL(1,I1),
     &           RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
     &           ZREL(I1),JWSREL(I1),IRSHIFT(I1),
     &           LMOMVEC,QMTET(IQ),QMPHI(IQ),MVEVIL1,MVEVIL2,NMVECMAX,
     &           IDOLDAU,LOPT(I1),PHILDAU(1,I1),WLDAU(1,1,1,I1),
     &           DENMATC(1,1,IPOT),
     &           LLY,NATYP)              ! LLY Lloyd
         END DO
C ----------------------------------------------------------------- SPIN
C
        IPOT1 = (I1-1)*NSPINPOT + 1
C
        DO LM = 1,LMPOTD
           DO IR = 1,IRMD
              RHO2NS(IR,LM,I1,1) = RHO2N1(IR,LM,IPOT1)
              R2NEF(IR,LM,I1,1)  = RHO2N2(IR,LM,IPOT1)
           END DO
        END DO

C     
        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +           DIMAG(DEN(L,IELAST,IPOT1))/PI/DBLE(NSPINPOT)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &           DIMAG(DEN(L,IELAST,IPOT1))/PI/DBLE(NSPINPOT)
        END DO
C     
        IF (NSPINPOT.EQ.2) THEN
           DO LM = 1,LMPOTD
              DO IR = 1,IRMD
                 RHO2NS(IR,LM,I1,2) = RHO2N1(IR,LM,IPOT1+1)
                 R2NEF(IR,LM,I1,2)  = RHO2N2(IR,LM,IPOT1+1)
              END DO
           END DO
C     
           DO L = 0,LMAXD1
              DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +              DIMAG(DEN(L,IELAST,IPOT1+1))/PI/DBLE(NSPINPOT)
              DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &              DIMAG(DEN(L,IELAST,IPOT1+1))/PI/DBLE(NSPINPOT) 
           END DO
        END IF

        IF ( TEST('RHOVALW ') ) THEN !Bauer
        open(unit=324234,file='out_rhoval')
        WRITE(324234,*) '#IATOM',I1
        write(324234,'(50000F)') RHO2NS(:,:,I1,1)
        IF (NSPIN==2) write(324234,'(50000F)') RHO2NS(:,:,I1,2)
        END IF


C ----------------------------------------------------------------------
C        
C --> itermdir/kmrot <> 0 
C
        IF (LMOMVEC) THEN
           DO IS = 1, NMVECMAX
              DO LM=1, 3
                 MVEVI(I1,LM,IS) = (0.0D0,0.0D0)
                 MVEVIEF(I1,LM,IS) = (0.0D0,0.0D0) 
                 DO L = 0, LMAXD
                    MVEVIL(L,I1,LM,IS) = MVEVIL1(L,LM,IS)
                    MVEVI(I1,LM,IS) = MVEVI(I1,LM,IS) +
     &                   MVEVIL1(L,LM,IS)
C     
                    MVEVIEF(I1,LM,IS) = MVEVIEF(I1,LM,IS) + 
     &                   MVEVIL2(L,LM,IS)
                 END DO
              END DO
           END DO
        END IF
C ----------------------------------------------------------------------
      END DO
C ================================================================ NATYP
      CLOSE (69)

      ELSE ! new spin-orbit solver

      THETA(:) = 0.D0
      PHI(:) = 0.D0
      LREAD = .FALSE.
      INQUIRE(file='nonco_angle.dat',EXIST=LREAD)
      IF (LREAD) THEN 
         OPEN(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
         DO I1 = 1,NATYP
            READ(10,*) THETA(I1),PHI(I1)
         ENDDO
      ENDIF

      IF (.NOT.TEST('FIXMOM  ')) THEN
       OPEN(UNIT=11,file='nonco_angle_new.dat',form='formatted')
       OPEN(UNIT=12,file='nonco_angle_out.dat',form='formatted')
      ENDIF

c interpolate potential

       CALL INTERPOLATE_POTEN(NSPIN,RMESH,IRMIN,IRWS,IPAN,IRCUT,VINS,
     +                        VISP,NPAN_LOG,NPAN_EQ,NCHEB,NPAN_TOT,
     +                        RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                        VINSNEW)

      DO I1 = 1,NATYP
       ICELL = NTCELL(I1)
       IPOT = (I1-1) * NSPIN + 1
       
        CALL RHOVALNEW(LDORHOEF,IELAST,NSRA,NSPIN,LMAX,EZ,WEZ,
     &           ZAT(I1),SOCSCALE(I1),CLEB(1,1),ICLEB,IEND,
     &           IFUNM1(1,ICELL),LMSP1(1,ICELL),
     &           NCHEB,NPAN_TOT(I1),NPAN_LOG(I1),
     &           NPAN_EQ(I1),RMESH(1,I1),IRWS(I1),RPAN_INTERVALL(0,I1),
     &           IPAN_INTERVALL(0,I1),RNEW(1,I1),VINSNEW,
     &           THETASNEW(1,1,ICELL),THETA(I1),PHI(I1),I1,IPOT,
     &           DEN1(0,1,1),ESPV1(0,1),RHO2M1,RHO2M2,MUORB(0,1,I1))

        DO L = 0,LMAXD1
         ESPV(L,IPOT)=ESPV1(L,1)
         ESPV(L,IPOT+NSPIN-1)=ESPV1(L,NSPIN)
         DO IE=1,IELAST
          DEN(L,IE,IPOT)=DEN1(L,IE,1)
          DEN(L,IE,IPOT+NSPIN-1)=DEN1(L,IE,NSPIN)
         ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
         DO LM = 1,LMPOTD
          DO IR = 1,IRMD
           RHO2NS(IR,LM,I1,ISPIN) = RHO2M1(IR,LM,ISPIN)
           R2NEF(IR,LM,I1,ISPIN)  = RHO2M2(IR,LM,ISPIN)
          END DO
         END DO
        END DO

        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +           DIMAG(DEN(L,IELAST,IPOT))/PI/DBLE(NSPIN)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &           DIMAG(DEN(L,IELAST,IPOT))/PI/DBLE(NSPIN)
        END DO
        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +            DIMAG(DEN(L,IELAST,IPOT+1))/PI/DBLE(NSPIN)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &            DIMAG(DEN(L,IELAST,IPOT+1))/PI/DBLE(NSPIN) 
        END DO
       DO ISPIN=1,NSPIN
        DO L=0,LMAXD1
         MUORB(L,3,I1)=MUORB(L,3,I1)+MUORB(L,ISPIN,I1)
        ENDDO
       ENDDO
       DO ISPIN=1,3
        DO L=0,LMAXD1
         MUORB(LMAXD1+1,ISPIN,I1)=MUORB(LMAXD1+1,ISPIN,I1)+
     &                            MUORB(L,ISPIN,I1)
        ENDDO
       ENDDO
      END DO !I1
      CLOSE(10)
      CLOSE(11)

      CLOSE (69)

c rewrite new theta and phi to nonco_angle.dat
      IF (.NOT.TEST('FIXMOM  ')) THEN
       OPEN(UNIT=10,file='nonco_angle.dat',form='formatted')
       CLOSE(10,status='delete')
       OPEN(UNIT=11,file='nonco_angle_new.dat',form='formatted')
       OPEN(UNIT=13,file='nonco_angle.dat',form='formatted')
       DO I1=1,NATYP
        READ(11,*) THETA(I1),PHI(I1)
        WRITE(13,*) THETA(I1),PHI(I1)
       ENDDO
       CLOSE(11,status='delete')
       CLOSE(13)
      ENDIF
      
      ENDIF


! In case of Lloyds formula renormalize valence charge                   ! LLY Lloyd
      IF (LLY.GT.0) THEN                                                 ! LLY Lloyd
         LMAXP1 = LMAX                                                   ! LLY Lloyd
         IF (INS.NE.0) LMAXP1 = LMAX + 1                                 ! LLY Lloyd
         CALL RENORM_LLY(
     >     CDOS_LLY,IELAST,NSPIN,NATYP,DEN,LMAXP1,CONC,
     >     1,IELAST,WEZ,IRCUT,IPAN,
     X     RHO2NS)


      ENDIF                                                              ! LLY Lloyd


C
C****************************************************** MPI COLLECT DATA
C
CMPI    IDIM = IRMD*LMPOTD*NATYPD*2
CMPI    CALL MPI_ALLREDUCE(RHO2NS,WORK,IDIM,
CMPI +                     MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                     MPI_COMM_WORLD,IERR)
CMPI    CALL DCOPY(IDIM,WORK,1,RHO2NS,1)
C
CMPI    IDIM = IRMD*LMPOTD*NATYPD*2
CMPI    CALL MPI_ALLREDUCE(R2NEF,WORK,IDIM,
CMPI +                     MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                     MPI_COMM_WORLD,IERR)
CMPI    CALL DCOPY(IDIM,WORK,1,R2NEF,1)
C
CMPI    IDIM = (LMAXD+2)*NPOTD
CMPI    CALL MPI_ALLREDUCE(ESPV,WORK,IDIM,
CMPI +                     MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                     MPI_COMM_WORLD,IERR)
CMPI    CALL DCOPY(IDIM,WORK,1,ESPV,1)
C  
CMPI    IDIM = IEMXD*(LMAXD+2)*NPOTD
CMPI    CALL MPI_ALLREDUCE(DEN,WORK,IDIM,
CMPI +                     MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,
CMPI +                     IERR)
CMPI    CALL ZCOPY(IDIM,WORK,1,DEN,1)
C  
CMPI    IF (IDOLDAU.EQ.1) THEN 
CMPI       IDIM = MMAXD*MMAXD*NPOTD
CMPI       CALL MPI_ALLREDUCE(DENMATC,WORK,IDIM,
CMPI +                        MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,
CMPI +                        IERR)
CMPI       CALL ZCOPY(IDIM,WORK,1,DENMATC,1)
CMPI    END IF
C
CMPI    IDIM = 1
CMPI    CALL MPI_ALLREDUCE(DENEF,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                     MPI_COMM_WORLD,IERR)
CMPI    CALL DCOPY(IDIM,WORK,1,DENEF,1)
C
CMPI    IF (KREL.EQ.1) THEN 
CMPI      IDIM = IRMD*NATYPD
CMPI      CALL MPI_ALLREDUCE(RHOORB,WORK,IDIM,
CMPI +                       MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                       MPI_COMM_WORLD,IERR)
CMPI      CALL DCOPY(IDIM,WORK,1,RHOORB,1)
C
CMPI      IDIM = (LMAXD+3)*NATYPD*3
CMPI      CALL MPI_ALLREDUCE(MUORB,WORK,IDIM,
CMPI +                       MPI_DOUBLE_PRECISION,MPI_SUM,
CMPI +                       MPI_COMM_WORLD,IERR)
CMPI      CALL DCOPY(IDIM,WORK,1,MUORB,1)
C
CMPI      IF (LMOMVEC) THEN
CMPI         IDIM = NATYPD*3*NMVECMAX
CMPI         CALL MPI_ALLREDUCE(MVEVI,WORK,IDIM,
CMPI +                          MPI_DOUBLE_COMPLEX,MPI_SUM,
CMPI +                          MPI_COMM_WORLD,IERR)
CMPI         CALL ZCOPY(IDIM,WORK,1,MVEVI,1)
C
CMPI         IDIM = (LMAXD+1)*NATYPD*3*NMVECMAX
CMPI         CALL MPI_ALLREDUCE(MVEVIL,WORK,IDIM,
CMPI +                          MPI_DOUBLE_COMPLEX,MPI_SUM,
CMPI +                          MPI_COMM_WORLD,IERR)
CMPI         CALL ZCOPY(IDIM,WORK,1,MVEVIL,1)
C
CMPI         IDIM = NATYPD*3*NMVECMAX
CMPI         CALL MPI_ALLREDUCE(MVEVIEF,WORK,IDIM,
CMPI +                          MPI_DOUBLE_COMPLEX,MPI_SUM,
CMPI +                          MPI_COMM_WORLD,IERR)
CMPI         CALL ZCOPY(IDIM,WORK,1,MVEVIEF,1)
C
CMPI      END IF          ! KREL.EQ.1
CMPI    END IF
C
C****************************************************** MPI COLLECT DATA
C
CMPI  IF(MYRANK.EQ.0) THEN
C
C================================================================ NATYP
      CHRGSEMICORE = 0D0
      DO I1 = 1,NATYP
C
C ---> l/m_s/atom-resolved charges 
C
          DO ISPIN = 1,NSPINPOT
             IPOT = (I1-1)*NSPINPOT + ISPIN
             DO L = 0,LMAXD1
                CHARGE(L,I1,ISPIN) = 0.0D0
C
                DO IE = 1,IELAST
                   CHARGE(L,I1,ISPIN) = CHARGE(L,I1,ISPIN) +
     +                  DIMAG(WEZ(IE)*DEN(L,IE,IPOT))/
     +                  DBLE(NSPINPOT)
                   IF ( IE.EQ.IESEMICORE ) CHRGSEMICORE = 
     +                  CHRGSEMICORE + CONC(I1)*CHARGE(L,I1,ISPIN)
                END DO
C
             END DO
          END DO
          EU(I1) = 0D0
          EDC(I1) = 0D0
C
C ---> orbital magnetic moments (array initialised to 0.0D0 in rhoval)
C
          IF (KREL.EQ.1) THEN
             DO ISPIN = 1,3
                DO L = 0,LMAXD + 1
                   MUORB(LMAXD1+1,ISPIN,I1) = MUORB(LMAXD1+1,ISPIN,I1)
     &                                       + MUORB(L,ISPIN,I1)
                END DO
             END DO
          END IF
      END DO
C================================================================ NATYP
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LDA+U
      IF ( IDOLDAU.EQ.1 ) THEN
C -> Save old LDA+U interaction matrix for mixing
         CALL DCOPY(MMAXD*MMAXD*NSPIND*NATYPD,WLDAU,1,WLDAUOLD,1)
C
C -> Construct LDA+U interaction matrix for next iteration
C
         CALL WMATLDAU(NTLDAU,ITLDAU,NSPINPOT,DENMATC,LOPT,
     &                 UEFF,JEFF,ULDAU,WLDAU,EU,EDC,MMAXD,NPOTD)
C -> Mix old and new LDA+U interaction matrices
         CALL MIXLDAU(
     >        MMAXD,NSPIND,NATYPD,NATYP,NSPIN,LOPT,WLDAUOLD,
     X        WLDAU)
C
C -> update variables-file
C
         ITRUNLDAU = ITRUNLDAU + 1
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
C
C -> write full lda+u information in ascii file ldaupot_new
C
         CALL WRLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,EREFLDAU,NATYP,WLDAU,
     &                  ULDAU,PHILDAU,IRMD,NATYPD,NSPIND,MMAXD)
      END IF
C
C LDA+U
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       CALL WRMOMS(KREL+KORBIT,NATYP,NSPINPOT,TEXTS,TEXTL,TEXTNS,CHARGE,
     &             MUORB,LMAXD,LMAXD1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
      IF ( ( KREL.EQ.1 ) .AND. LMOMVEC ) THEN
         DO I1=1,NATYP
            IQ = IQAT(I1)
            CALL MVECGLOBAL(I1,IQ,NATYP,QMPHI(IQ),QMTET(IQ),
     &                      MVEVI,MVEVIL,MVEVIEF,
     &                      NATYPD,LMAXD,NMVECMAX)
         END DO
      END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c    TEST BRAHIM
      IF (NPOL.EQ.0 .OR.TEST('DOS     ')) CALL WRLDOS(DEN,EZ,WEZ,
     +    LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,E1,E2,ALAT,TK,NACLS1,
     +    NSPINPOT,NATYP,CONC,IELAST,INTERVX,INTERVY,INTERVZ,DOSTOT)


C
C ---------------------------------------------------------- CORE STATES
C   RHO_core is calculated only if also RHO_valence was
C
      IF (NPOL.NE.0) THEN 
         WRITE (6,*)
         WRITE (6,'(78(1H#))')
         WRITE (6,'(33X,A)') 'CORE  STATES'
         WRITE (6,'(78(1H#))')
         DO I1 = 1,NATYP
            DO ISPIN = 1,NSPIN
               IPOT = (I1-1) * NSPINPOT + ISPIN
               IPOT1 = (I1-1) * NSPINPOT + 1
C
               CALL RHOCORE(NSRA,ISPIN,NSPIN,I1,
     &              DRDI(1,I1),RMESH(1,I1),VISP(1,IPOT),
     &              A(I1),B(I1),ZAT(I1),
     &              IRCUT(0,I1),RHOC(1,IPOT1),
     &              ECORE(1,IPOT),NCORE(IPOT),LCORE(1,IPOT),
     &              CSCL(1,KREL*I1+(1-KREL)),
     &              VTREL(1,I1),BTREL(1,I1),
     &              RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
     &              ZREL(I1),JWSREL(I1),IRSHIFT(I1),
     &              ECOREREL(1,IPOT1),NKCORE(1,I1),
     &              KAPCORE(1,IPOT1))
C
            END DO
         END DO
         WRITE (6,*)
         WRITE (6,'(78(1H#))')
         WRITE (6,*)
      END IF
C
C ---------------------------------------------------------- CORE STATES
C
      OPEN (67,FILE='density',FORM='unformatted')
      WRITE (67) RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,ESPV,ECORE,
     &           IDOLDAU,LOPT,EU,EDC,CHRGSEMICORE
      IF (KREL.EQ.1) WRITE(67) RHOORB,ECOREREL,NKCORE,KAPCORE
      CLOSE (67)
      IF (TEST('den-asci')) THEN
         OPEN (67,FILE='densitydn.ascii',FORM='formatted')
         DO I1 = 1,NATYP
            DO LM = 1,LMPOTD
               DO IR = 1,IRMD
                  WRITE(67,FMT='(I6,2I5,2E25.16)') 
     &               I1,LM,IR,RHO2NS(IR,LM,I1,1),RHO2NS(IR,LM,I1,2)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(67)
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C      ITERMDIR
C
      IF (ITERMVDIR) THEN
         OPEN (67,FILE='itermdir.unformatted',FORM='unformatted')
         READ (67) 
         READ (67) 
         WRITE (67) MVEVI,MVEVIEF
         CLOSE (67) 
      END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 


C
CMPI  ENDIF
C
CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME4 = (MEND-MSTART)
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1,' seconds',
CT3E +  TIME2,' seconds',
CT3E +  TIME3,' seconds',
CT3E +  TIME4,' seconds'
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1+TIME2+TIME3+TIME4,' seconds'
CT3E  TIME1=TIME1*NROFNODES
CT3E  TIME2=TIME2*NROFNODES
CT3E  TIME3=TIME3*NROFNODES
CT3E  TIME4=TIME4*NROFNODES
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1,' seconds',
CT3E +  TIME2,' seconds',
CT3E +  TIME3,' seconds',
CT3E +  TIME4,' seconds'
CT3E  IF (MYRANK.EQ.0) WRITE (6,FMT=*) 'elapsed time ',
CT3E +  TIME1+TIME2+TIME3+TIME4,' seconds'
CMPI  CALL MPI_FINALIZE(IERR)
      !STOP

      WRITE (6,'(79(1H=),/,30X,"< KKR1c finished >",/,79(1H=),/)')
      
      END SUBROUTINE !MAIN1c
      
      END MODULE