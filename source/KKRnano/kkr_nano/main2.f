C
C KKRnano
C massive parallel KKR for nanoscaled systems
C
C
      PROGRAM MAIN2
C
      IMPLICIT NONE
C
      INCLUDE 'mpif.h'
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
C     .. Parameters ..
      INTEGER   LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
      INTEGER   NPOTD
      PARAMETER (NPOTD=NSPIND*NAEZD)
      INTEGER   LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER    MMAXD
      PARAMETER (MMAXD  = 2*LMAXD + 1)
      INTEGER   LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER   LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER   LASSLD
      PARAMETER (LASSLD=4*LMAXD)
      INTEGER   LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER   IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER   LRECPOT
      PARAMETER (LRECPOT=8*(LMPOTD*(IRNSD+1)+IRMD+20))
      INTEGER   LRECRES1,LRECRES2
      PARAMETER (LRECRES2=4+8*(NSPIND*(LMAXD+7)+2*LPOTD+4+2))
      INTEGER   NTIRD
      PARAMETER (NTIRD=(IRMD+(IRNSD+1)*(LMPOTD-1))*NSPIND)
      INTEGER   MAXMSHD
      PARAMETER (MAXMSHD=8)
      INTEGER   LLYALM 
      PARAMETER (LLYALM=LLY*(NAEZD*LMMAXD-1)+1)
      INTEGER   LLYNGD 
      PARAMETER (LLYNGD=LLY*(NACLSD*LMMAXD-1)+1)
      INTEGER   NSYMAXD
      PARAMETER (NSYMAXD=48)
      DOUBLE COMPLEX CZERO
      PARAMETER      (CZERO=(0.0D0,0.0D0))
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX   JSCAL   ! scaling factor for Jij calculation
      DOUBLE PRECISION DENEF,E1,E2,TK,EFERMI,EBOT,EFOLD,
     +                 CHRGNT,E2SHIFT,RCUTJIJ,DF,
     +                 ALAT,FCM,MIX,MIXING,QBOUND,VOLUME0,
     +                 RMAX,GMAX,FPI,PI,RFPI,RMSAVM,RMSAVQ,
     +                 EREFLDAU,                                 !LDA+U
     +                 WALLCLOCK_I,WALLCLOCK_F   ! time management
      REAL             TIME_I,TIME_S,            ! time management
     +                 TIME_E,TIME_EX
      INTEGER          ICELL,IR,LMAX,
     +                 NPNT1,NPNT2,NPNT3,NPOL,
     +                 NSPIN,
     +                 ITER,SCFSTEPS,IMIX,NOITER,NOITER_ALL,
     +                 IPF,ISHIFT,
     +                 KPRE,KTE,KVMAD,KXC,KFORCE,
     +                 L,LPOT,LMPOT,NAEZ,IEND,NCLEBD,LM1,LM2,
     +                 IEND1,
     +                 I,J,IPOT,ISPIN,I1,I1BRYD,
     +                 IH,IRC1,IRMIN1,LM,NR,EKM,
     +                 SYSTEM_I,SYSTEM_F,
     +                 RATETIME,MAXTIME 
      LOGICAL          TEST,LCORDENS,XCCPL,JIJ,STOPIT,LDORHOEF,LDAU,
     +                 ERESJIJ
C     ..
C
C     .. Local Arrays ..
      DOUBLE COMPLEX   EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX   DEZ(IEMXD),
     +                 WEZRN(IEMXD,2),
     +                 DEN(0:LMAXD1,IEMXD,NSPIND),
     +                 DSYMLL(LMMAXD,LMMAXD,48),
     +                 PHILDAU(IRMD,LMAXD1),                     !LDA+U
     +                 DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
      DOUBLE PRECISION WG(LASSLD),YRG(LASSLD,0:LASSLD,0:LASSLD),
     +                 BRAVAIS(3,3),RBASIS(3,NAEZD),RECBV(3,3),
     +                 VBC(2),VAV0,VOL0,
     +                 SMAT(LMXSPD,NAEZD),CLEB(LMXSPD*LMPOTD),
     +                 CLEB1C(NCLEB,2),RNORM(IEMXD,2),
     +                 DFAC(0:LPOTD,0:LPOTD),RWS(NAEZD),RMT(NAEZD),
     +                 GN(3,NMAXD),RM(3,NMAXD),
     +                 BZKP(3,KPOIBZ,MAXMSHD),
     +                 VOLCUB(KPOIBZ,MAXMSHD),VOLBZ(MAXMSHD),
     +                 RR(3,0:NRD),
C     .. input potential (spherical VISP, nonspherical VINS)
     +                 VINS(IRMIND:IRMD,LMPOTD,2),
     +                 VISP(IRMD,2),
C     .. output potential (nonspherical VONS)
     +                 VONS(IRMD,LMPOTD,2),
     +                 ULDAU(LMAXD1),JLDAU(LMAXD1),              !LDA+U
     +                 UMLDAU(MMAXD,MMAXD,MMAXD,MMAXD,LMAXD1),   !LDA+U
     +                 WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)         !LDA+U
      INTEGER          NMESH,KMESH(IEMXD),NOFKS(MAXMSHD),
     +                 EZOA(NACLSD,NAEZD),NSYMAT,
     +                 NUMN0(NAEZD),INDN0(NAEZD,NACLSD),
     +                 ID,MAXMESH,
     +                 NSG(ISHLD),NSR(ISHLD),
     +                 LLDAU(LMAXD1),                            !LDA+U
     +                 NGMAX,NRMAX,NSHLG,NSHLR
C     ..
C ----------------------------------------------------------------------
C   ECOU(0:LPOTD,NAEZD)    ! Coulomb energy                    
C   EPOTIN(NAEZD),         ! energy of input potential (EPOTINB
C   ESPC(0:3,NPOTD),        ! energy single particle core       
C   ESPV(0:LMAXD1,NPOTD)    ! energy single particle valence    
C   EXC(0:LPOTD,NAEZD),    ! E_xc
C ----------------------------------------------------------------------
      DOUBLE COMPLEX
     +        TMATN(LMMAXD,LMMAXD,NSPIND),
     +        DTDE(LMMAXD,LMMAXD,NSPIND),
     +        TREFLL(LMMAXD,LMMAXD,NREFD),
     +        DTREFLL(LMMAXD,LMMAXD,NREFD),
     +        DGREFN(LMMAXD,LMMAXD,NACLSD,NCLSD),
     +        GREFN(LMMAXD,LMMAXD,NACLSD,NCLSD),
     +        GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
     +        GMATN_ALL(LMMAXD,LMMAXD,IEMXD,NSPIND),
     +        DTIXIJ(LMMAXD,LMMAXD),
     +        GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND),
     +        GXIJ_ALL(LMMAXD,LMMAXD,NXIJD,NSPIND)
C----- preconditioning ---------------------------------------------------
      COMPLEX          PRSC(NGUESSD*LMMAXD,EKMD,NSPIND-SMPID+1)
      DOUBLE PRECISION QMRBOUND
      DOUBLE PRECISION CNVFAC(EKMD,NSPIND-SMPID+1)
      INTEGER          IGUESS,BCP,PRSPIN
      INTEGER          SPRS(NGUESSD*LMMAXD+1,EKMD+1,NSPIND-SMPID+1)
C----- Lloyd ----------------------------------------------------------- 
      DOUBLE COMPLEX LLY_G0TR(IEMXD,NCLSD),
     +               LLY_GRDT(IEMXD,NSPIND),
     +               LLY_GRDT_ALL(IEMXD,NSPIND),
     +               TR_ALPH(NSPIND),
     +               JXCIJINT(NXIJD)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      DOUBLE PRECISION ECOU(0:LPOTD),EPOTIN,
     +                 ESPC(0:3,NSPIND),ESPV(0:LMAXD1,NSPIND),
     +                 EXC(0:LPOTD),EULDAU,EDCLDAU
      DOUBLE PRECISION A(NAEZD),B(NAEZD),DRDI(IRMD,NAEZD)
      DOUBLE PRECISION R(IRMD,NAEZD),ECORE(20,2)
      DOUBLE PRECISION THETAS(IRID,NFUND,NCELLD),ZAT(NAEZD)
C ----------------------------------------------------------------------
      DOUBLE PRECISION RHOCAT(IRMD,2),
     +                 R2NEF(IRMD,LMPOTD,2)
      DOUBLE PRECISION RHO2NS(IRMD,LMPOTD,2)
      DOUBLE PRECISION CHARGE(0:LMAXD1,2)
C     ..
      DOUBLE PRECISION GSH(NGSHD)
C ----------------------------------------------------------------------
C  CMINST(LMPOTD,NAEZD)            ! charge moment of interstitial
C  CMOM(LMPOTD,NAEZD)              ! LM moment of total charge
C  CATOM                            ! total charge per atom
C  QC                               ! core charge 
C ----------------------------------------------------------------------
      DOUBLE PRECISION CMINST(LMPOTD),CMOM(LMPOTD),CATOM(NSPIND),QC
      DOUBLE PRECISION VMAD
C     ,,
C     .. FORCES
      DOUBLE PRECISION FLM(-1:1,NAEZD),FLMC(-1:1,NAEZD)
C     .. MIXING
      DOUBLE PRECISION SM1S(NTIRD),FM1S(NTIRD)
      DOUBLE PRECISION UI2(NTIRD,2:ITDBRYD),VI2(NTIRD,2:ITDBRYD),
     +                 WIT(2:ITDBRYD) 
C ----------------------------------------------------------------------
      INTEGER IMT(NAEZD),IPAN(NAEZD),IRC(NAEZD),
     +        IRNS(NAEZD),
     +        IRCUT(0:IPAND,NAEZD),IRMIN(NAEZD),
     +        IRWS(NAEZD),ITITLE(20,NPOTD)
      INTEGER LCORE(20,NPOTD),LCOREMAX,LLMSP(NFUND,NAEZD)
      INTEGER NCORE(NPOTD),NFU(NAEZD),
     +        NTCELL(NAEZD),ISYMINDEX(48)
      INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      INTEGER ICLEB(LMXSPD*LMPOTD,3),LOFLM(LMXSPD)
      INTEGER ICLEB1C(NCLEB,3),LOFLM1C(LM2D)
      INTEGER IFUNM(LMXSPD,NAEZD),ICST,NSRA,
     +        LMSP(LMXSPD,NAEZD),JEND(LMPOTD,0:LMAXD,0:LMAXD)
      INTEGER NCLS,NREF,RF,NLDAU
      DOUBLE PRECISION RCLS(3,NACLSD,NCLSD),RMTREF(NREFD),VREF(NAEZD)
      DOUBLE PRECISION RXIJ(NXIJD),       ! interatomic distance Ri-Rj
     +                 RXCCLS(3,NXIJD),   ! position relative of j rel. to i (sorted) 
     +                 ZKRXIJ(48,3,NXIJD) ! set up in clsjij, used in kkrmat01
      INTEGER          IXCP(NXIJD),       ! index to atom in elem/cell at site in cluster
     +                 NXCP(NXIJD),       ! index to bravais lattice at site in cluster
     +                 XIJ,NXIJ
      INTEGER          ATOM(NACLSD,NAEZD),
     +                 CLS(NAEZD),
     +                 NACLS(NCLSD),
     +                 REFPOT(NAEZD)
      INTEGER          NUTRC,              ! number of inequivalent atoms in the cluster
     +                 INTRC(NATRCD),      ! pointer to atoms in the unit cell
     +                 NATRC,              ! number of atoms in cluster
     +                 ATTRC(NATRCD),      ! index to atom in elem/cell at site in cluster
     +                 EZTRC(NATRCD)       ! index to bravais lattice  at site in cluster
C-----------------------------------------------------------------------
C     ..
C-----------------------------------------------------------------------
C     .. MPI ..
C     .. N-MPI
      INTEGER      MYRANK,NROFNODES,IERR
      INTEGER      MAPBLOCK
      EXTERNAL     MAPBLOCK
      COMMON       /MPI/MYRANK,NROFNODES
C     .. L-MPI
      INTEGER      MYLRANK(LMPID*SMPID*EMPID),
     +             LCOMM(LMPID*SMPID*EMPID),
     +             LGROUP(LMPID*SMPID*EMPID),
     +             LSIZE(LMPID*SMPID*EMPID),
     +             LMPIC
C     .. LS-MPI
      INTEGER      LSRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSMYRANK(LMPID,NAEZD*SMPID*EMPID),
     +             LSMPIC,LSMPIB
C     .. S-MPI
      INTEGER      SRANK(SMPID,NAEZD*LMPID*EMPID),
     +             SMYRANK(SMPID,NAEZD*LMPID*EMPID),
     +             SMPIC,SMPIB,MAPSPIN
C     .. E-MPI
      REAL         ETIME(IEMXD)
      INTEGER      EMYRANK(EMPID,NAEZD*LMPID*SMPID),
     +             ERANK(EMPID,NAEZD*LMPID*SMPID),
     +             EMPIC,EMPIB,
     +             IE,IELAST,
     +             EPROC(IEMXD),EPROCO(IEMXD)
C     .. ACTV-MPI
      INTEGER      MYACTVRANK,
     +             ACTVCOMM,ACTVGROUP,ACTVSIZE,
     +             MYBCRANK,BCRANK
C     .. MPI_OP_CREATE
!       INTEGER      CLS_OP
!       EXTERNAL     CLS_SCAT

      DOUBLE PRECISION WORK1(2),WORK2(2)

      EXTERNAL MPI_ALLREDUCE,
     +         MPI_FINALIZE,MPI_SEND
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(8),TESTC(16)
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC
      COMMON /TESTC/TESTC
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,CONVOL,ECOUB,EPOTINB,ESPCB,
     &         ETOTB1,FORCE,FORCEH,FORCXC,MTZERO,KLOOPZ1,
     &         TEST,VMADELBLK,VXCDRV,RMSOUT,
     &         OUTTIME,TIME
c
C     .. Intrinsic Functions ..
      INTRINSIC DABS,ATAN,DMIN1,DSIGN,SQRT,MAX,DBLE
C     ..
C============================================================= CONSTANTS
      LCORDENS=.TRUE.
      PI = 4.0D0*ATAN(1.0D0)
      FPI = 4.0D0*PI
      RFPI = SQRT(FPI)
      IPF = 74


C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
      OPEN (67,FILE='inp.unf',FORM='unformatted')
      READ (67) KMESH,MAXMESH,RR,EZOA,NUMN0,INDN0,NSYMAT,DSYMLL
      READ (67) NAEZ,NSPIN,IPAN,IRNS,IRCUT,LCORE,NCORE,NTCELL,
     &          LMAX,LPOT,LMPOT
      READ (67) IMIX,MIXING,QBOUND,FCM,KPRE,KTE,
     &          KVMAD,KXC,ISHIFT,KFORCE,IGUESS,BCP,QMRBOUND
      READ (67) A,B,DRDI,R,THETAS,ZAT,IMT,IRC,
     &          IRMIN,IRWS,RWS,RMT,ITITLE,LLMSP,NFU
      READ (67) ALAT,GSH,ILM,IMAXSH,TESTC,OPTC
      READ (67) BRAVAIS,RBASIS,RECBV,VOLUME0,RMAX,GMAX
      READ (67) IEND1,CLEB1C,ICLEB1C,LOFLM1C,JEND,IFUNM,LMSP,NSRA,ICST
      READ (67) NCLS,NREF,RMTREF,VREF,RCLS,
     &          ATOM,CLS,NACLS,REFPOT
      READ (67) NR,RCUTJIJ,JIJ,LDAU
      READ (67) ISYMINDEX,SCFSTEPS
      CLOSE (67)

C ---------------------------------------------------------- k_mesh
C
      OPEN (52,FILE='kpoints',FORM='formatted')
      REWIND (52)
      DO L = 1,MAXMESH
         READ (52,FMT='(I8,f15.10)') NOFKS(L),VOLBZ(L)
            READ (52,FMT=*) (BZKP(ID,1,L),ID=1,3),VOLCUB(1,L)
         DO I=2,NOFKS(L)
            READ (52,FMT=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
         END DO
      END DO
      CLOSE (52)

      OPEN (67,FILE='energy_mesh',FORM='unformatted')
      READ (67) IELAST,EZ,WEZ,E1,E2
      READ (67) NPOL,TK,NPNT1,NPNT2,NPNT3
      IF ( NPOL.EQ.0 ) READ(67) EFERMI
      CLOSE (67)

      IF (KFORCE.EQ.1) OPEN (54,FILE='force',FORM='formatted')

C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
      CALL IMPI(
     >          NAEZ,
     <          MYRANK,NROFNODES,
     <          LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     <          LSMPIB,LSMPIC,LSRANK,LSMYRANK,
     <          SMPIB,SMPIC,SRANK,SMYRANK,
     <          EMPIB,EMPIC,ERANK,EMYRANK,
     <          MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE)
C
C
C
C====================================================================
C
C=====================================================================
C     processors not fitting in NAEZ*LMPID*SMPID*EMPID do nothing ...
C ... and wait after SC-ITER loop
C=====================================================================
C
      IF (LMPIC.NE.0.OR.LSMPIC.NE.0) THEN ! ACTVGROUP
C
      MYBCRANK = 0
C
C
      OPEN (IPF,ACCESS='direct',RECL=NSPIN*8,FILE='RMS.unf',
     +     FORM='unformatted')
C
C
C ======================================================================
C========= TIMING ======================================================
      IF (MYLRANK(1).EQ.0) THEN
        RATETIME = 100
        MAXTIME  = 100000
        CALL SYSTEM_CLOCK(SYSTEM_I,RATETIME,MAXTIME)
        WALLCLOCK_I = MPI_WTIME()
        CALL CPU_TIME(TIME_I)
        OPEN (2,FILE='time-info',FORM='formatted')
      ENDIF
C========= TIMING ======================================================
C -------------------------------------------------------- not.converged 
C
      OPEN (28,FILE='not.converged',FORM='formatted')
      READ (28,'(1P,4D17.10)') EFOLD,VBC
      CLOSE (28)
C
      DO IH=1,EKMD
      DO ISPIN=1,NSPIND-SMPID+1
        CNVFAC(IH,ISPIN) = 1000.0D0
      ENDDO
      ENDDO
C
      DO IH=1,NTIRD
      DO LM1=2,ITDBRYD
           UI2(IH,LM1)=0.00
           VI2(IH,LM1)=0.00
      ENDDO
           SM1S(IH)=0.00
           FM1S(IH)=0.00
      ENDDO
C
C
C ######################################################################
C ######################################################################
      DO ITER = 1, SCFSTEPS 
C ######################################################################
C ######################################################################
C
C
        CALL CPU_TIME(TIME_S)
C
        EKM    = 0
        NOITER = 0
C
C
        IF (MYLRANK(1).EQ.0) THEN
          WRITE(2,'(79(1H=))')
      CALL OUTTIME(MYLRANK(1),'started at ..........',TIME_I,ITER)
          WRITE(2,'(79(1H=))')
        ENDIF
C
C
        LDORHOEF = NPOL.NE.0
C
        CALL GAUNT2(WG,YRG)
C
        CALL MADELUNG3D(LPOT,YRG,WG,ALAT,
     &                  RMAX,GMAX,BRAVAIS,RECBV,
     &                  LMXSPD,LASSLD,LPOTD,LMPOTD,
     &                  NMAXD,ISHLD,
     &                  LMPOT,CLEB,ICLEB,IEND,
     &                  NCLEBD,LOFLM,DFAC,
     &                  NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM)
C
        DO LM = 1,LMPOTD
          CMOM(LM) = 0.0D0
          CMINST(LM) = 0.0D0
        END DO
C
        CHRGNT = 0.0D0
C
        LRECRES1 = 8*43 + 16*(LMAXD+2)
        IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN
          LRECRES1 = LRECRES1 + 32*(LMAXD+2)*IEMXD
        END IF
C
      OPEN (71,ACCESS='direct',RECL=LRECRES1,FILE='results1',
     +     FORM='unformatted')
      OPEN (66,ACCESS='direct',RECL=LRECPOT*2,FILE='vpotnew',
     +     FORM='unformatted')
      IF (TRC.EQ.1) OPEN (37,ACCESS='direct',RECL=LRECTRC,
     +     FILE='trnc.unf',FORM='unformatted')
C
C
CN ====================================================================
C     BEGIN do loop over atoms (NMPID-parallel)
CN ====================================================================
C
      DO I1 = 1,NAEZ
      IF(MYLRANK(LMPIC).EQ.MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN
c
c=======================================================================
cxccpl
c
      XCCPL = .FALSE.
      IF ((ITER.EQ.SCFSTEPS).AND.JIJ) XCCPL = .TRUE.
c
      IF (XCCPL) THEN
c
        INQUIRE(FILE='ERESJIJ',EXIST=ERESJIJ)
c
        CALL CLSJIJ(
     >              I1,NAEZ,RR,NR,RBASIS,RCUTJIJ,LMPIC,NSYMAT,ISYMINDEX,
     <              IXCP,NXCP,NXIJ,RXIJ,RXCCLS,ZKRXIJ)
c
         DO ISPIN = 1, NSPIN
           DO XIJ = 1, NXIJ
             DO LM1 = 1,LMMAXD
               DO LM2 = 1,LMMAXD
                 GMATXIJ(LM1,LM2,XIJ,ISPIN) = CZERO
               ENDDO
             ENDDO
           ENDDO
         ENDDO
c
      ENDIF
c
cxccpl
c=======================================================================
c
C
        READ(66,REC=I1) VINS,VISP,ECORE
        IF (TRC.EQ.1) READ(37,REC=I1) NATRC,ATTRC,EZTRC,NUTRC,INTRC
C
C
C LDA+U
        IF (LDAU) THEN
C
          EREFLDAU = EFERMI
          EREFLDAU = 0.48
C
!          IF (ITER.EQ.1) THEN
            CALL LDAUINIT(
     >                    I1,ITER,NSRA,NLDAU,LLDAU,ULDAU,JLDAU,EREFLDAU,
     >                    VISP,NSPIN,R(1,I1),DRDI(1,I1),
     >                    ZAT(I1),IPAN(I1),IRCUT(0,I1),
     <                    PHILDAU,UMLDAU,WMLDAU)
C
!          ENDIF
C
        ENDIF
C LDA+U
C
C
C
CTIME
      CALL OUTTIME(MYLRANK(1),'initialized .........',TIME_I,ITER)
CTIME
C
C
C
C
CIE ====================================================================
C     BEGIN do loop over energies (EMPID-parallel)
CIE ====================================================================
C
        DO IE = 1,IELAST

C
          CALL CPU_TIME(TIME_E)
C
          ETIME(IE) = 0.0D0
C
          IF (ITER.EQ.1.AND.IE.EQ.1) THEN
            CALL EBALANCE(
     >                    'I',ITER,SCFSTEPS,
     >                    IELAST,NPNT1,
     >                    MYACTVRANK,ACTVCOMM,
     >                    ETIME,EPROC,EPROCO)
          ENDIF
C
CIE ====================================================================
        IF (EMPIB.EQ.EPROC(IE)) THEN
CIE ====================================================================
C
 
C

          DO RF = 1,NREF
C
            CALL TREF(EZ(IE),VREF(RF),LMAX,RMTREF(RF),
     +                TREFLL(1,1,RF),DTREFLL(1,1,RF))
C
          END DO
C
C
          CALL GREF(EZ(IE),ALAT,IEND1,NCLS,NAEZ,
     +              CLEB1C,RCLS,ATOM,CLS,ICLEB1C,LOFLM1C,NACLS,
     +              REFPOT,
     <              TREFLL(1,1,1),DTREFLL(1,1,1),GREFN,DGREFN,
     +              IE,
     +              LLY_G0TR,I1,
     +              LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE)
C
        DO 370 ISPIN = 1,NSPIN
 
           CALL CALCTMAT(LDAU,NLDAU,ICST,
     +                   NSRA,EZ(IE),
     +                   DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN),
     +                   VISP(1,ISPIN),ZAT(I1),IPAN(I1),
     +                   IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1,
     <                   TMATN(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN,
     >                   LLDAU,WMLDAU)
         IF(LLY.EQ.1) THEN
           CALL CALCDTMAT(LDAU,NLDAU,ICST,
     +                   NSRA,EZ(IE),IE,NPNT1,NPNT2,NPNT3,PI,TK,
     +                   DRDI(1,I1),R(1,I1),VINS(IRMIND,1,ISPIN),
     +                   VISP(1,ISPIN),ZAT(I1),IPAN(I1),
     +                   IRCUT(0,I1),CLEB1C,LOFLM1C,ICLEB1C,IEND1,
     <                   DTDE(1,1,ISPIN),TR_ALPH(ISPIN),LMAX,ISPIN,
     >                   LLDAU,WMLDAU)
         END IF
C
C
         IF (XCCPL) THEN
           IF (ISPIN.EQ.1) THEN
             DO LM1 = 1,LMMAXD
               DO LM2 = 1,LMMAXD
                 DTIXIJ(LM1,LM2) = TMATN(LM1,LM2,ISPIN)
               ENDDO
             ENDDO
           ELSE
             DO LM1 = 1,LMMAXD
               DO LM2 = 1,LMMAXD
                 DTIXIJ(LM1,LM2) = DTIXIJ(LM1,LM2)-TMATN(LM1,LM2,ISPIN)
               ENDDO
             ENDDO
           ENDIF
         ENDIF
C
C
          RF = REFPOT(I1)
c         DO ISPIN = 1,NSPIN
            DO LM1 = 1,LMMAXD
               TMATN(LM1,LM1,ISPIN) =  TMATN(LM1,LM1,ISPIN)
     +                      - TREFLL(LM1,LM1,RF)
               DTDE(LM1,LM1,ISPIN) =  DTDE(LM1,LM1,ISPIN)
     +                      - DTREFLL(LM1,LM1,RF)
            END DO
c         END DO
C
C
C
CSPIN ==================================================================
C     BEGIN do loop over spins (SMPID-parallel)
CSPIN===================================================================
C
C         initialize LLY_GRDT
C
          TR_ALPH(ISPIN) = TR_ALPH(ISPIN) - LLY_G0TR(IE,CLS(I1))
          LLY_GRDT(IE,ISPIN) = CZERO
C
          IF (SMPID.EQ.1) THEN
            MAPSPIN = 0
            PRSPIN   = ISPIN
          ELSE
            MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
            PRSPIN   = 1
          ENDIF
C
C         true beginning of SMPID-parallel section
C
          IF(SRANK(SMPIB,SMPIC).EQ.MAPSPIN) THEN
C
            NMESH = KMESH(IE)
C
            IF( MYLRANK(LMPIC).EQ.0 ) THEN
            WRITE (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)') 
     &             ' ** IE = ',IE,' ENERGY =',EZ(IE),
     &                ' KMESH = ', NMESH,' ISPIN = ',ISPIN
            END IF
C
C
C <<>>
C
C
            CALL KLOOPZ1(
     <                   GMATN(1,1,1,ISPIN),
     >                   ALAT,IE,IELAST,ITER,NAEZ,
     >                   NOFKS(NMESH),VOLBZ(NMESH),
     >                   BZKP(1,1,NMESH),VOLCUB(1,NMESH),
     >                   CLS,NACLS,RR,
     >                   EZOA,ATOM,GREFN,DGREFN,
     >                   NSYMAT,DSYMLL,
     >                   TMATN(:,:,ISPIN),DTDE(:,:,ISPIN),
     >                   NUMN0,INDN0,I1,
     >                   NATRC,ATTRC,EZTRC,NUTRC,INTRC,
     <                   SPRS(1,1,PRSPIN),PRSC(1,1,PRSPIN),
     >                   EKM,NOITER,
     >                   EZ,QMRBOUND,IGUESS,BCP,CNVFAC(1,PRSPIN),
     >                   NXIJ,XCCPL,IXCP,ZKRXIJ,
     >                   LLY_GRDT(IE,ISPIN),TR_ALPH(ISPIN),
     <                   GMATXIJ(1,1,1,ISPIN),
     >                   LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE,
     >                   LSMYRANK,LSRANK,LSMPIB,LSMPIC)
C
C
C

        ENDIF

 370    CONTINUE                    !  ISPIN = 1,NSPIN
C
C
C
CSPIN ==================================================================
C     END do loop over spins (SMPID-parallel)
CSPIN===================================================================
C
C
C
c
c
c
c =====================================================================
c Calculate Jij for the in CLSJIJ predefined atom pairs i,j
cxccpl
c
        IF (XCCPL) THEN
c
            CALL SREDGX(
     >                  ISPIN,NSPIN,
     >                  MYRANK,NROFNODES,
     >                  SMPIB,SMPIC,SMYRANK,SRANK,
     >                  GMATXIJ,
     <                  GXIJ_ALL)
c
            JSCAL = WEZ(IE)/DBLE(NSPIN)
c
            CALL XCCPLJIJ(
     >                    'R',I1,IE,JSCAL,
     >                    RXIJ,NXIJ,IXCP,RXCCLS,
     >                    GXIJ_ALL,DTIXIJ,
     >                    LMPIC,LCOMM,
     >                    MYRANK,EMPIC,EMYRANK,
     <                    JXCIJINT,ERESJIJ)
c
        END IF
c
cxccpl
c End of Jij calculation
c =====================================================================
C
          CALL CPU_TIME(TIME_EX)
C
          ETIME(IE) = TIME_EX-TIME_E
C
C
CIE ====================================================================
       ENDIF
CIE ====================================================================
C
C for preconditioning purposes calculate sparse indices combining IE.KPT

       EKM = EKM + NOFKS(KMESH(IE))
C
C
       END DO                       ! IE = 1,IELAST
C
       IF (ERESJIJ) CLOSE(75)
C
CIE ====================================================================
C     END do loop over energies (EMPID-parallel) to be implemented
CIE ====================================================================
C
C
C
C
C=======================================================================
C     "allreduce" information of 1 .. EMPID and 1 .. SMPID processors
C=======================================================================
      CALL SREDGM(
     >            NSPIN,IELAST,
     >            MYRANK,
     >            SMPIC,SMYRANK,SRANK,
     >            EMPIC,EMYRANK,ERANK,EPROC,
     >            GMATN,LLY_GRDT,
     <            GMATN_ALL,LLY_GRDT_ALL)
C=======================================================================
C=======================================================================
C
C
C
C
CTIME
      CALL OUTTIME(MYLRANK(1),'G obtained ..........',TIME_I,ITER)
CTIME
C
C
C
C
C=======================================================================
C     output of Jij's .. calling xccpljij with flag 'F'
C=======================================================================
      IF (XCCPL) THEN
      CALL XCCPLJIJ(
     >              'F',I1,IE,JSCAL,
     >              RXIJ,NXIJ,IXCP,RXCCLS,
     >              GXIJ_ALL,DTIXIJ,
     >              LMPIC,LCOMM,
     >              MYRANK,EMPIC,EMYRANK,
     <              JXCIJINT,ERESJIJ)
      ENDIF
C=======================================================================
C=======================================================================
C
C
C
C
C=======================================================================
C     on the basis of new timings determine now new distribution of
C     work to 1 .. EMPID processors
C=======================================================================
      CALL EBALANCE(
     >              'R',ITER,SCFSTEPS,
     >              IELAST,NPNT1,
     >              MYACTVRANK,ACTVCOMM,
     >              ETIME,EPROC,EPROCO)

C=======================================================================
C=======================================================================
C
C
C=======================================================================
C     in case of IGUESS and EMPID > 1 preconditioning arrays might
C     have to be adjusted to new distributions
C=======================================================================
      IF ((IGUESS.EQ.1).AND.(EMPID.GT.1)) THEN
C
        DO ISPIN = 1,NSPIN
C
        IF (SMPID.EQ.1) THEN
          MAPSPIN = 0
          PRSPIN   = ISPIN
        ELSE
          MAPSPIN = MAPBLOCK(ISPIN,1,SMPID,1,0,SMPID-1)
          PRSPIN   = 1
        ENDIF
C
C       true beginning of SMPID-parallel section
C
        IF(SRANK(SMPIB,SMPIC).EQ.MAPSPIN) THEN
C
          CALL EPRDIST(
     >                 IELAST,KMESH,NOFKS,
     >                 PRSC(1,1,PRSPIN),
     >                 SPRS(1,1,PRSPIN),
     >                 CNVFAC(1,PRSPIN),
     >                 MYRANK,EMPIC,EMYRANK,
     >                 EPROC,EPROCO)
C
        ENDIF
        ENDDO
C
      ENDIF
C=======================================================================
C=======================================================================
C
C
C
C
C
C
C----------------------------------------------------------------------
C BEGIN only processes with LMPIC = 1 are working
C----------------------------------------------------------------------
      IF (LMPIC.EQ.1) THEN
C
        IF (LLY.EQ.1) THEN

        CALL LLOYD0(EZ,WEZ,CLEB1C,DRDI,R,IRMIN,VINS,VISP,
     +              THETAS,ZAT,ICLEB1C,
     +              IFUNM,IPAN,IRCUT,LMSP,JEND,LOFLM1C,
     +              NTCELL,ICST,
     +              IELAST,IEND1,NAEZ,NSPIN,NSRA,
     <              WEZRN,RNORM,
     >              GMATN_ALL,
     >              LLY_GRDT_ALL,
     >              LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <              DMATLDAU,
     >              LMPIC,MYLRANK,
     >              LGROUP,LCOMM,LSIZE)

C           IF (MYLRANK(LMPIC).EQ.0) THEN
C           DO IE=1,IELAST
C        WRITE(6,*) 'WEZRN1',IE,DREAL(WEZRN(IE,1)),DIMAG(WEZRN(IE,1))
C        WRITE(6,*) 'WEZRN2',IE,DREAL(WEZRN(IE,2)),DIMAG(WEZRN(IE,2))
C           ENDDO
C           ENDIF
C
CTIME
      CALL OUTTIME(MYLRANK(1),'Lloyd processed......',TIME_I,ITER)
CTIME
C
        ELSE
          DO IE=1,IELAST
            WEZRN(IE,1) = WEZ(IE)
            WEZRN(IE,2) = WEZ(IE)
          ENDDO
        ENDIF
C
        CALL CINIT(IEMXD*(LMAXD+2)*NSPIND,DEN)
        DENEF = 0.0D0
C
C
        IF (LDAU) THEN
          CALL CINIT(MMAXD*MMAXD*NSPIND*LMAXD1,DMATLDAU(1,1,1,1))
        ENDIF
C
CSPIN ==================================================================
C     BEGIN do loop over spins
CSPIN===================================================================
C
        DO ISPIN = 1,NSPIN
          ICELL = NTCELL(I1)
          IPOT = (I1-1) * NSPIN + ISPIN
C
          CALL RHOVAL(LDORHOEF,ICST,IELAST,
     &                NSRA,ISPIN,NSPIN,EZ,WEZRN(1,ISPIN),
     &                DRDI(1,I1),R(1,I1),IRMIN(I1),
     +                VINS(IRMIND,1,ISPIN),VISP(1,ISPIN),
     &                ZAT(I1),IPAN(I1),IRCUT(0,I1),
     &                THETAS(1,1,ICELL),IFUNM(1,ICELL),LMSP(1,ICELL),
     &                RHO2NS,R2NEF,
     &                DEN(0,1,ISPIN),ESPV(0,ISPIN),
     &                CLEB1C,LOFLM1C,ICLEB1C,IEND1,JEND,
     >                GMATN_ALL,
     >                LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU,
     <                DMATLDAU)



          IF (LLY.EQ.1) THEN
            DO IE=1,IELAST
              DO L=0,LMAXD1
                DEN(L,IE,ISPIN)=DEN(L,IE,ISPIN)*RNORM(IE,ISPIN)
              END DO
            END DO
          END IF
C
          EBOT = E1

          CALL RHOCORE(EBOT,NSRA,ISPIN,NSPIN,I1,
     &                 DRDI(1,I1),R(1,I1),VISP(1,ISPIN),
     &                 A(I1),B(I1),ZAT(I1),
     &                 IRCUT(0,I1),RHOCAT,QC,
     &                 ECORE(1,ISPIN),NCORE(IPOT),LCORE(1,IPOT))

        END DO

C
CSPIN ==================================================================
C     END do loop over spins
CSPIN===================================================================
C
         DO ISPIN = 1,NSPIN
            DO L = 0,LMAXD1
               DENEF = DENEF - 2.0D0 *
     +                 DIMAG(DEN(L,IELAST,ISPIN))/PI/DBLE(NSPIN)
            END DO
         END DO

C
C ---> l/m_s/atom-resolved charges 
C
        DO ISPIN = 1,NSPIN
          DO L = 0,LMAXD1
             CHARGE(L,ISPIN) = 0.0D0
C
             DO IE = 1,IELAST
                 CHARGE(L,ISPIN) = CHARGE(L,ISPIN) +
     +                             DIMAG(WEZ(IE)*DEN(L,IE,ISPIN))/
     +                             DBLE(NSPIN)

             END DO
C
          END DO
        END DO
C
C
CLDAU
C
        EULDAU = 0.0D0
        EDCLDAU = 0.0D0
C
        IF (LDAU.AND.NLDAU.GE.1) THEN
C
          CALL LDAUWMAT(I1,NSPIN,ITER,MIXING,DMATLDAU,NLDAU,LLDAU,
     &                  ULDAU,JLDAU,UMLDAU,WMLDAU,EULDAU,EDCLDAU)
C
        ENDIF
C
CLDAU
C
C

C ----------------------------------------------------------------------
C
C -->   determine total charge density expanded in spherical harmonics
C
C -------------------------------------------------------------- density
C
        CALL RHOTOTB(NAEZ,I1,NSPIN,RHO2NS,RHOCAT,
     +               ZAT,DRDI,IRCUT,
     +               LPOT,NFU,LLMSP(1,ICELL),THETAS,ICELL,IPAN,
     +               CATOM)
C
        CHRGNT = CHRGNT + CATOM(1) - ZAT(I1)

        IF (NPOL.EQ.0 .OR. TEST('DOS     ')) THEN 
          WRITE(71,REC=I1) QC,CATOM,CHARGE,ECORE,DEN
        ELSE
          WRITE(71,REC=I1) QC,CATOM,CHARGE,ECORE
        END IF
C
      ENDIF
C----------------------------------------------------------------------
C END L-MPI: only processes with LMPIC = 1 are working
C----------------------------------------------------------------------
C
C
C
      END IF
      END DO
C
CN ====================================================================
C     END do loop over atoms (NMPID-parallel)
CN ====================================================================
C
      IF (TRC.EQ.1) CLOSE(37)
      CLOSE(66)
      CLOSE(71)
C
C
      CALL OUTTIME(MYLRANK(1),'density calculated ..',TIME_I,ITER)
C
C
C----------------------------------------------------------------------
C BEGIN L-MPI: only processes with LMPIC = 1 are working
C----------------------------------------------------------------------
      IF (LMPIC.EQ.1) THEN
C
C****************************************************** MPI COLLECT DATA
C
      CALL MPI_ALLREDUCE(CHRGNT,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                   LCOMM(LMPIC),IERR)
      CALL DCOPY(1,WORK1,1,CHRGNT,1)
C
      CALL MPI_ALLREDUCE(DENEF,WORK1,1,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                   LCOMM(LMPIC),IERR)
      CALL DCOPY(1,WORK1,1,DENEF,1)

C****************************************************** MPI COLLECT DATA
C
C **************************  ITERATION NUMBER  ************************
      OPEN (28,FILE='not.converged',FORM='formatted')
      READ (28,'(1P,4D17.10)') EFOLD,VBC
      CLOSE (28)
C **********************************************************************

C ----------------------------------------------------------------------
C
C --> determine new Fermi level due to valence charge up to 
C     old Fermi level E2 and density of states DENEF
C
      E2SHIFT = CHRGNT/DENEF
C
      E2SHIFT = DMIN1(DABS(E2SHIFT),0.03D0)*DSIGN(1.0D0,E2SHIFT)

      EFOLD = E2
      IF (ISHIFT.LT.2) E2 = E2 - E2SHIFT
C
      IF( MYLRANK(LMPIC).EQ.0 ) THEN
      WRITE (6,FMT=9020) EFOLD,E2SHIFT
C
C --> divided by NAEZ because the weight of each atom has been already
C     taken into account in 1c
C
      WRITE (6,FMT=9030) E2,DENEF/DBLE(NAEZ) 
      WRITE(6,'(79(1H+),/)')
      END IF
C ----------------------------------------------------------------------
      DF = 2.0D0/PI*E2SHIFT/DBLE(NSPIN)
C ----------------------------------------------------------------------

      OPEN (66,ACCESS='direct',RECL=LRECPOT*2,FILE='vpotnew',
     +     FORM='unformatted')
      OPEN (72,ACCESS='direct',RECL=LRECRES2,FILE='results2',
     +     FORM='unformatted')

C =====================================================================
C ======= I1 = 1,NAEZ ================================================
C =====================================================================

      DO I1 = 1,NAEZ
      IF(MYLRANK(LMPIC).EQ.
     +   MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN

      ICELL = NTCELL(I1)

      DO ISPIN = 1,NSPIN
C
C -->     get correct density and valence band energies
C
         ESPV(0,ISPIN) = ESPV(0,ISPIN) -
     +                        EFOLD*CHRGNT/DBLE(NSPIN*NAEZ)
         IF (LCORDENS) THEN 
               DO LM = 1,LMPOT
                  CALL DAXPY(IRC(I1),DF,R2NEF(1,LM,ISPIN),1,
     &                       RHO2NS(1,LM,ISPIN),1)
               END DO
         END IF

C ----------------------------------------------------------------------
      END DO

      CALL RHOMOM(CMOM,CMINST,LPOT,I1,RHO2NS,
     +     R,DRDI,IRWS,IRCUT,IPAN,ICELL,ILM,IFUNM(1,ICELL),IMAXSH,GSH,
     +     THETAS,LMSP(1,ICELL))
      CALL OUTTIME(MYLRANK(1),'RHOMOM ......',TIME_I,ITER)

C =====================================================================
C ============================= ENERGY and FORCES =====================
C =====================================================================

      CALL VINTRAS(LPOT,NSPIN,I1,RHO2NS,VONS,
     +     R,DRDI,IRWS,IRCUT,IPAN,ICELL,ILM,IFUNM(1,ICELL),IMAXSH,GSH,
     +     THETAS,LMSP(1,ICELL))
      CALL OUTTIME(MYLRANK(1),'VINTRAS ......',TIME_I,ITER)

      CALL STRMAT(ALAT,LPOT,NAEZ,NGMAX,NRMAX,NSG,NSR,NSHLG,NSHLR,GN,RM,
     &            RBASIS,SMAT,VOLUME0,LASSLD,LMXSPD,NAEZD,I1)
      CALL OUTTIME(MYLRANK(1),'STRMAT ......',TIME_I,ITER)

      CALL VMADELBLK(CMOM,CMINST,LPOT,NSPIN,
     &               NAEZ,VONS,ZAT,R,IRCUT,IPAN,
     &               VMAD,
     &               LMPOT,SMAT,CLEB,ICLEB,IEND,
     &               LMXSPD,NCLEBD,LOFLM,DFAC,I1,
     >               LMPIC,MYLRANK,
     >               LGROUP,LCOMM,LSIZE)
      CALL OUTTIME(MYLRANK(1),'VMADELBLK ......',TIME_I,ITER)

C =====================================================================
C
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES
C
      IF (KFORCE.EQ.1.AND.ITER.EQ.SCFSTEPS) THEN
C ---------------------------------------------------------------------
         CALL FORCEH(CMOM,FLM,LPOT,NSPIN,I1,RHO2NS,VONS,
     +               R,DRDI,IMT,ZAT)
         CALL FORCE(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R,
     +              DRDI,IMT)
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
         CALL ESPCB(ESPC,NSPIN,I1,ECORE,LCORE,LCOREMAX,NCORE)
C
         CALL EPOTINB(EPOTIN,NSPIN,I1,RHO2NS,VISP,R,DRDI,
     &                IRMIN,IRWS,LPOT,VINS,IRCUT,IPAN,ZAT)
C

         CALL ECOUB(CMOM,ECOU,LPOT,NSPIN,I1,RHO2NS,
     &              VONS,ZAT,R,
     &              DRDI,IRWS,KVMAD,IRCUT,IPAN,IMAXSH,IFUNM(1,ICELL),
     &              ILM,ICELL,GSH,THETAS,LMSP(1,ICELL))
      END IF
      CALL OUTTIME(MYLRANK(1),'KTE ......',TIME_I,ITER)

C EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
C
C =====================================================================
      CALL VXCDRV(EXC,KTE,KXC,LPOT,NSPIN,I1,RHO2NS,
     +            VONS,R,DRDI,A,
     +            IRWS,IRCUT,IPAN,ICELL,GSH,ILM,IMAXSH,IFUNM(1,ICELL),
     +            THETAS,LMSP(1,ICELL))
      CALL OUTTIME(MYLRANK(1),'VXCDRV ......',TIME_I,ITER)

C =====================================================================
C
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES
C
C Force calculation continues here 
C
      IF (KFORCE.EQ.1.AND.ITER.EQ.SCFSTEPS) THEN
C ---------------------------------------------------------------------
         CALL FORCXC(FLM,FLMC,LPOT,NSPIN,I1,RHOCAT,VONS,R,
     >               ALAT,DRDI,IMT,ZAT,
     >               LMPIC,MYLRANK,
     >               LGROUP,LCOMM,LSIZE)
C ---------------------------------------------------------------------
      END IF

      WRITE(72,REC=I1) CATOM,VMAD,ECOU,EPOTIN,ESPC,ESPV,EXC,LCOREMAX,
     +                 EULDAU,EDCLDAU
C
C Force calculation ends
C FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
C =====================================================================

      VAV0 = 0.0D0
      VOL0 = 0.0D0

      CALL MTZERO(LMPOT,I1,NSPIN,VONS,ZAT,R,DRDI,IMT,IRCUT,
     &            IPAN,ICELL,LMSP(1,ICELL),IFUNM(1,ICELL),
     &            THETAS,IRWS,VAV0,VOL0)
      CALL OUTTIME(MYLRANK(1),'MTZERO ......',TIME_I,ITER)


C =====================================================================
C ============================= ENERGY and FORCES =====================
C =====================================================================

      END IF
      END DO

      CALL OUTTIME(MYLRANK(1),'calculated pot ......',TIME_I,ITER)

C =====================================================================
C ======= I1 = 1,NAEZ ================================================
C =====================================================================


C****************************************************** MPI COLLECT DATA
      WORK1(1) = VAV0 
      WORK1(2) = VOL0 
      CALL MPI_ALLREDUCE(WORK1,WORK2,2,MPI_DOUBLE_PRECISION,MPI_SUM,
     +                   LCOMM(LMPIC),IERR)
      VAV0 = WORK2(1)
      VOL0 = WORK2(2)
C****************************************************** MPI COLLECT DATA

      VBC(1) = -VAV0/VOL0
      IF (ISHIFT.GT.0) VBC(1) = VBC(1) + E2SHIFT
      VBC(2) = VBC(1) 
      IF(MYRANK.EQ.0) THEN
      WRITE (6,FMT=9103) VOL0,VAV0,VBC(1)
      WRITE(6,'(79(1H=),/)')
      END IF
 9103 FORMAT ('  VOL INT.',F16.9,'  VAV INT.',F16.9,'  VMT ZERO',F16.9)
C =====================================================================
C
C ---------------------------------------------------------------------
C
c -->   shift potential to muffin tin zero and
C       convolute potential with shape function for next iteration
C


       DO I1 = 1,NAEZ
       IF(MYLRANK(LMPIC).EQ.
     +    MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN

       ICELL = NTCELL(I1)

C =====================================================================
C ============================= POTENTIAL MIXING OUTPUT ===============
C =====================================================================

          DO ISPIN = 1,NSPIN
             IPOT = NSPIN* (I1-1) + ISPIN
             DO IR = 1,IRCUT(IPAN(I1),I1)
                VONS(IR,1,ISPIN) = VONS(IR,1,ISPIN) + RFPI*VBC(ISPIN)
             END DO
             CALL CONVOL(IRCUT(1,I1),IRC(I1),ICELL,
     &                   IMAXSH(LMPOT),ILM,IFUNM(1,ICELL),LMPOT,GSH,
     &                   THETAS,ZAT(I1),RFPI,
     &                   R(1,I1),VONS(1,1,ISPIN),LMSP(1,ICELL))
          END DO

C -->   final construction of the potentials (straight mixing)


      MIX = MIXING
      RMSAVQ = 0.0D0
      RMSAVM = 0.0D0

      CALL MIXSTR(RMSAVQ,RMSAVM,LPOT,LMPOT,
     +            I1,NSPIN,
     +            ITER,RFPI,FPI,IPF,
     +            MIX,
     +            FCM,IRC,IRMIN,R,DRDI,VONS,
     +            VISP,VINS)

      I1BRYD=I1
      END IF
      END DO

C -->  potential mixing procedures: Broyden or Andersen updating schemes

        IF (IMIX.GE.3) THEN

           CALL BRYDBM(VISP,VONS,VINS,
     +                 LMPOT,R,DRDI,MIX,
     +                 IRC,IRMIN,NSPIN,I1BRYD,NAEZ,
     +                 IMIX,IPF,ITER,
     +                 UI2,VI2,WIT,SM1S,FM1S,
     >                 LMPIC,MYLRANK,
     >                 LGROUP,LCOMM,LSIZE)

        ENDIF

C----------------------------------------------------------------------
C -->    reset to start new iteration
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

       DO I1 = 1,NAEZ
       IF(MYLRANK(LMPIC).EQ.
     +    MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) THEN

         DO ISPIN = 1,NSPIN
         IPOT = (I1-1)*NSPIN + ISPIN
C
            IRC1 = IRC(I1)
            CALL DCOPY(IRC1,VONS(1,1,ISPIN),1,VISP(1,ISPIN),1)
C
            IF (LPOT.GT.0) THEN
               IRMIN1 = IRMIN(I1)
               DO LM = 2,LMPOT
                  DO J = IRMIN1,IRC1
                     VINS(J,LM,ISPIN) = VONS(J,LM,ISPIN)
                  END DO
               END DO
            END IF
         ENDDO

C ----------------------------------------------------- output_potential
      WRITE(66,REC=I1) VINS,VISP,ECORE
C ----------------------------------------------------- output_potential

      END IF
      END DO

!       CLOSE(66)
      CLOSE(72)

C =====================================================================
C ============================= POTENTIAL MIXING OUTPUT ===============
C =====================================================================

C ====== write RMS convergency data - not parallelized, written by 
C MYRANK=0   (RMSOUT) ================================================= 

      CALL RMSOUT(RMSAVQ,RMSAVM,ITER,E2,EFOLD,
     +            SCFSTEPS,VBC,QBOUND,NSPIN,NAEZ,
     >            KXC,LPOT,A,B,IRC,
     >            VINS,VISP,DRDI,IRNS,R,RWS,RMT,ALAT,
     >            ECORE,LCORE,NCORE,ZAT,ITITLE,
     >            LMPIC,MYLRANK,
     >            LGROUP,LCOMM,LSIZE)


 9020 FORMAT ('                old',
     &     ' E FERMI ',F12.6,' Delta E_F = ',f12.6)
 9030 FORMAT ('                new',
     &     ' E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)
C
C
C
C Wait here in order to guarantee regular and non-errorneous output
C in RESULTS
C
        CALL MPI_BARRIER(LCOMM(LMPIC),IERR)
C
C
C ----------------------------------------------------------------- 
C L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
C ----------------------------------------------------------------- 
        IF(MYLRANK(LMPIC).EQ.0) THEN
C
          CALL RESULTS(LRECRES2,IELAST,ITER,LMAX,NAEZ,NPOL,
     +                 NSPIN,KPRE,KTE,LPOT,E1,E2,TK,EFERMI,
     +                 ALAT,ITITLE,CHRGNT,ZAT,EZ,WEZ,LDAU)
C
C --> update energy contour
C
C
          IF(TEST('fix-EF  ')) REWIND(41)
          IF(TEST('fix-EF  ')) READ(41,*) E2
C
          CALL EMESHT(EZ,DEZ,IELAST,E1,E2,E2,TK,
     &                NPOL,NPNT1,NPNT2,NPNT3,IEMXD)

C
          DO IE = 1,IELAST
            WEZ(IE) = -2.D0/PI*DEZ(IE)
          END DO
C
          WRITE(6,'(79(1H=))')
C
C .. get info on MYACTVRANK of this processor: to be used in
C    subsequent reduce-commands
          MYBCRANK = MYACTVRANK
C ..
        ENDIF
C ----------------------------------------------------------------- 
C L-MPI: only process with MYLRANK(LMPIC = 1) = 0 is working here
C ----------------------------------------------------------------- 
C
C
C ..
      ENDIF
C -----------------------------------------------------------------
C L-MPI: only processes with LMPIC = 1 are working here
C ----------------------------------------------------------------- 

      CLOSE(66)
!      CLOSE(72)
C -----------------------------------------------------------------
C
C      CALL MPI_BARRIER(ACTVCOMM,IERR)
C
      CALL MPI_ALLREDUCE(MYBCRANK,BCRANK,1,MPI_INTEGER,MPI_MAX,
     +                   ACTVCOMM,IERR)
C
      CALL MPI_BCAST(EZ,IEMXD,MPI_DOUBLE_COMPLEX,
     +               BCRANK,ACTVCOMM,IERR)
C
      CALL MPI_BCAST(WEZ,IEMXD,MPI_DOUBLE_COMPLEX,
     +               BCRANK,ACTVCOMM,IERR)
C
      CALL MPI_BCAST(E1,1,MPI_DOUBLE_PRECISION,
     +               BCRANK,ACTVCOMM,IERR)
C
      CALL MPI_BCAST(E2,1,MPI_DOUBLE_PRECISION,
     +               BCRANK,ACTVCOMM,IERR)
C
      CALL MPI_ALLREDUCE(NOITER,NOITER_ALL,1,MPI_INTEGER,MPI_SUM,
     +                   ACTVCOMM,IERR)
C
C      CALL MPI_BARRIER(ACTVCOMM,IERR)
C
      IF(MYLRANK(1).EQ.0) THEN
C
        OPEN (67,FILE='energy_mesh',FORM='unformatted')
        WRITE (67) IELAST,EZ,WEZ,E1,E2
        WRITE (67) NPOL,TK,NPNT1,NPNT2,NPNT3
        CLOSE (67)
C
        WRITE(6,'(79(1H=))')
        WRITE(6,'(19X,A,I3,A,I10)') '       ITERATION : ',
     &                             ITER,' SUM of QMR ',NOITER_ALL
        WRITE(6,'(79(1H=),/)')

        CALL SYSTEM_CLOCK(SYSTEM_F,RATETIME,MAXTIME)
        WALLCLOCK_F=MPI_WTIME()
        WRITE(6,'(79(1H=))')
        WRITE(6,*) 'Wallclock, System and CPU-time compared:'
        WRITE(6,*) 'MPI_WTIME=',(WALLCLOCK_F-WALLCLOCK_I),'sec'
        WRITE(6,*) 'SYSTEM_CLOCK=',(SYSTEM_F-SYSTEM_I)/100,'sec'
        CALL OUTTIME(MYLRANK(1),'end .................',
     &               TIME_I,ITER)
        WRITE(6,'(79(1H=))')
        WRITE(2,'(79(1H=))')
        CALL OUTTIME(MYLRANK(1),'finished in .........',
     &               TIME_S,ITER)
        WRITE(2,'(79(1H=))')
        WRITE(6,'(79(1H=),/)')
      ENDIF
C
C manual exit possible by creation of file 'STOP' in home directory
C
      INQUIRE(FILE='STOP',EXIST=STOPIT)
      IF (STOPIT) GOTO 200
C
C
C ######################################################################
C ######################################################################
      ENDDO          ! SC ITERATION LOOP ITER=1, SCFSTEPS
C ######################################################################
C ######################################################################
C
C
  200 CONTINUE
C
C
      IF (MYLRANK(1).EQ.0) CLOSE(2)            ! TIME
      CLOSE(IPF)                                ! RMS
      IF (KFORCE.EQ.1) CLOSE(54)

C ======================================================================
C ======================================================================
C
C Free communicators and groups ..
C ..
      IF (MYLRANK(LMPIC).GE.0) THEN
      CALL MPI_COMM_FREE(LCOMM(LMPIC),IERR)
      ENDIF
      CALL MPI_GROUP_FREE(LGROUP(LMPIC),IERR)
C
      CALL MPI_COMM_FREE(ACTVCOMM,IERR)
      CALL MPI_GROUP_FREE(ACTVGROUP,IERR)
C .. .
C
C
C
      ENDIF ! ACTVGROUP
C
C=====================================================================
C     processors not fitting in NAEZ*LMPID do nothing ...
C ... and wait here
C=====================================================================
C
C
C      WRITE(6,*) 'BARRIER i:',MYRANK

      CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

C      WRITE(6,*) 'BARRIER f:',MYRANK

      CALL MPI_FINALIZE(IERR)

      STOP

      END

