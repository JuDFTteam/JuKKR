      module mod_wunfiles
        
      implicit none
      
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
     &                    IJTABCALC,ISH,JSH,IJTABSYM,IJTABSH,NOFGIJ,
     &                    NQCALC,IQCALC,KMROT,KAOEZ,IQAT,NOQ,CONC,
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
     &                    TOLRDIF,LLY,DELTAE)
C **********************************************************************
C *                                                                    *
C *  This subroutine is part of the MAIN0 program in the tbkkr package *
C *  It writes out different unformatted files meant to provide the    *
C *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
C *  during an SCF cycle                                               *
C *  v.popescu, munich 2004                                            *
C *                                                                    *
C **********************************************************************

      use mod_types

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
      DOUBLE PRECISION A(NATYPD),B(NATYPD),THETAS(IRID,NFUND,NCELLD)
      DOUBLE PRECISION RMT(NATYPD),RMTNEW(NATYPD),RWS(NATYPD),GSH(NGSHD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
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
     &        IJTABCALC(NOFGIJD),IJTABSYM(NOFGIJD),IJTABSH(NOFGIJD)
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
      CHARACTER*24 TXC(4)
      CHARACTER*80 TMPDIR
      INTEGER ITMPDIR,ILTMP
C     ..
C     .. Arrays in common
      CHARACTER*8 TESTC(32),OPTC(32)
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
C     ..
C     .. Local scalars
      INTEGER I1,I2
C     ..
C     .. External Functions ..
      LOGICAL OPT
      EXTERNAL OPT
C -------------------------------------------------------------- SCRATCH
C Looking for SCRATCH system variable to store files: gmat, tmat 
C and gref on the local file system (this is necessary if you want
C MPI to run properly when paralellising over energies)
C
c call to SCRATCHDIR commented out 15.09.2006 by fivos in iff820c.. cluster with gfortran
c     CALL SCRATCHDIR(TMPDIR,ITMPDIR,ILTMP)
C
C ---------------------------------------------------------- energy_mesh
C                                    some data in this file might change
C
      ITMPDIR=0
      ILTMP=0

      OPEN (67,FILE='energy_mesh',FORM='unformatted')
      WRITE (67) IELAST,EZ,WEZ,E1,E2,IESEMICORE,FSEMICORE
      WRITE (67) NPOL,TK,NPNT1,NPNT2,NPNT3,EBOTSEMI,EMUSEMI,TKSEMI,
     &           NPOLSEMI,N1SEMI,N2SEMI,N3SEMI
      IF ( NPOL.EQ.0 ) WRITE(67) EFERMI
      CLOSE (67)
C ------------------------------------------------------ input_potential
C                                          some data in this file change
C
      E1 = 0D0
!       OPEN (67,FILE='input_potential',FORM='unformatted')
      OPEN (67,FILE='input_scf.unformatted',FORM='unformatted')
      WRITE (67) VINS,VISP,ECORE,VBC
      IF (KREL.EQ.1) THEN
         WRITE (67) RMREL,DRDIREL,R2DRDIREL
         WRITE (67) ZREL,JWSREL,IRSHIFT
         WRITE (67) VTREL,BTREL
      END IF
      WRITE (67) ITSCF,SCFSTEPS,E1,E1,CMOMHOST
      CLOSE (67)
      
      type0%i_iteration = ITSCF
      type0%N_iteration = SCFSTEPS

C ------------------------------------------------------------- itermdir
C                                               data in this file change
C
      IF (OPT('ITERMDIR')) THEN
         OPEN (67,FILE='itermdir.unformatted',FORM='unformatted')
         I1 = 0
         E1 = 0D0
         WRITE (67) QMTET,QMPHI,QMPHITAB,QMTETTAB,QMGAMTAB,I1,E1
         WRITE (67) DROTQ
         CLOSE(67)
      END IF
C ---------------------------------------------------------------- lda+u
C                                               data in this file change
C
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
C -------------------------------------------------------------- input1a
C                                 meant for MAIN1a, data does not change
C
      OPEN (67,FILE='input1a.unformatted',FORM='unformatted')
      WRITE(67) NSRA,INS,NAEZ,NATYP,NSPIN,ICST,IPAN,IRCUT,   ! LLY added NAEZ
     &          LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY
      WRITE(67) ALAT,ZAT,DRDI,R,RMTREF,VREF,IEND,CLEB,RCLS,
     &          ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,TESTC,OPTC,
     &          IRWS,IRMIN,TOLRDIF,DELTAE,SOCSCALE
      WRITE(67) TMPDIR,ITMPDIR,ILTMP
      WRITE(67) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      WRITE(67) RNEW,RPAN_INTERVALL,IPAN_INTERVALL
      IF ( KREL.EQ.1 ) WRITE(67) SOLVER,SOCSCL,CSCL
      IF ( IDOLDAU.EQ.1 ) WRITE(67) NTLDAU,ITLDAU,LOPT,UEFF,JEFF,
     &                              EREFLDAU
      CLOSE(67)
      ! test fivos begins
      OPEN (67,FILE='input1a.formatted')
      WRITE(67,*) NSRA,INS,NAEZ,NATYP,NSPIN,ICST,IPAN,IRCUT,   ! LLY added NAEZ
     &          LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY
      WRITE(67,*) ALAT,ZAT,DRDI,R,RMTREF,VREF,IEND,CLEB,RCLS,
     &          ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,TESTC,OPTC,
     &          IRWS,IRMIN,TOLRDIF,DELTAE,SOCSCALE
      WRITE(67,*) TMPDIR,ITMPDIR,ILTMP
      WRITE(67,*) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      WRITE(67,*) RNEW,RPAN_INTERVALL,IPAN_INTERVALL
      CLOSE(67)
      ! test fivos ends

C -------------------------------------------------------------- input1b
C                                 meant for MAIN1b, data does not change
C
      OPEN (67,FILE='input1b.unformatted',FORM='unformatted')
      WRITE(67) NSRA,INS,NATYP,NAEZ,NSPIN,LMAX,NREF,ICC,IGF,
     &          NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL
      WRITE(67) ALAT,RBASIS,REFPOT,RMTREF,VREF,RCLS,RR,
     &          ATOM,CLS,NCLS,EZOA,NACLS,NSHELL,KMROT,KAOEZ,IQAT,NOQ,
     &          CONC,KMESH,MAXMESH,TESTC,OPTC,LLY
      WRITE(67) NSYMAT,NATOMIMP,NOFGIJ,NQCALC,RATOM,RROT,NSH1,NSH2,DROTQ
      WRITE(67) IJTABCALC
      WRITE(67) ((ISH(I1,I2),I2=1,NOFGIJ),I1=1,NSHELL(0))
      WRITE(67) ((JSH(I1,I2),I2=1,NOFGIJ),I1=1,NSHELL(0))
      WRITE(67) IJTABSYM,IJTABSH,IQCALC
      WRITE(67) DSYMLL,INVMOD,ICHECK,ATOMIMP,SYMUNITARY
      WRITE(67) TMPDIR,ITMPDIR,ILTMP
      IF ( KREL.EQ.1 ) WRITE(67) RC,CREL,RREL,SRREL,NRREL,IRREL
      IF ( OPT('DECIMATE') ) WRITE(67) LEFTTINVLL,RIGHTTINVLL,VACFLAG
      CLOSE(67)

      ! test fivos begins
      OPEN (67,FILE='input1b.formatted')
      WRITE(67,*) NSRA,INS,NATYP,NAEZ,NSPIN,LMAX,NREF,ICC,IGF,
     &          NLBASIS,NRBASIS,NCPA,ICPA,ITCPAMAX,CPATOL
      WRITE(67,*) ALAT,RBASIS,REFPOT,RMTREF,VREF,RCLS,RR,
     &          ATOM,CLS,NCLS,EZOA,NACLS,NSHELL,KMROT,KAOEZ,IQAT,NOQ,
     &          CONC,KMESH,MAXMESH,TESTC,OPTC,LLY
      WRITE(67,*) NSYMAT,NATOMIMP,NOFGIJ,NQCALC,RATOM,RROT,NSH1,NSH2,
     &                                                          DROTQ
      WRITE(67,*) IJTABCALC
      WRITE(67,*) ((ISH(I1,I2),I2=1,NOFGIJ),I1=1,NSHELL(0))
      WRITE(67,*) ((JSH(I1,I2),I2=1,NOFGIJ),I1=1,NSHELL(0))
      WRITE(67,*) IJTABSYM,IJTABSH,IQCALC
      WRITE(67,*) DSYMLL,INVMOD,ICHECK,ATOMIMP,SYMUNITARY
      WRITE(67,*) TMPDIR,ITMPDIR,ILTMP
      CLOSE(67)
      ! test fivos ends
C -------------------------------------------------------------- input1c
C                                 meant for MAIN1c, data does not change
C
      OPEN (67,FILE='input1c.unformatted',FORM='unformatted')
      WRITE(67) NSRA,INS,NATYP,NAEZ,NSPIN,ICST,IPAN,IRCUT,KMROT,IQAT,
     &          CONC,QMTET,QMPHI,IDOLDAU,LMAX,IRWS
      WRITE(67) ALAT,ZAT,DRDI,R,A,B,IEND,CLEB,ICLEB,LOFLM,JEND,THETAS,
     &   IFUNM1,LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,IRMIN,ITITLE,INTERVX,
     &          INTERVY,INTERVZ,NACLS(1),TESTC,OPTC,LLY,SOCSCALE
      WRITE(67) TMPDIR,ITMPDIR,ILTMP
      WRITE(67) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      WRITE(67) RNEW,RPAN_INTERVALL,IPAN_INTERVALL,THETASNEW
      IF (KREL.EQ.1) WRITE(67) SOLVER,SOCSCL,CSCL
      IF ( IDOLDAU.EQ.1 ) WRITE(67) NTLDAU,ITLDAU,LOPT,UEFF,JEFF,
     &                              EREFLDAU
      CLOSE (67)


      ! test fivos begins
      OPEN (67,FILE='input1c.formatted')
      WRITE(67,*) NSRA,INS,NATYP,NAEZ,NSPIN,ICST,IPAN,IRCUT,KMROT,IQAT,
     &          CONC,QMTET,QMPHI,IDOLDAU,LMAX,IRWS
      WRITE(67,*) ALAT,ZAT,DRDI,R,A,B,IEND,CLEB,ICLEB,LOFLM,JEND,THETAS,
     &   IFUNM1,LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,IRMIN,ITITLE,INTERVX,
     &          INTERVY,INTERVZ,NACLS(1),TESTC,OPTC,LLY,SOCSCALE
      WRITE(67,*) TMPDIR,ITMPDIR,ILTMP
      WRITE(67,*) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      WRITE(67,*) RNEW,RPAN_INTERVALL,IPAN_INTERVALL,THETASNEW
      CLOSE (67)
      ! test fivos ends

C --------------------------------------------------------------- input2
C                                  meant for MAIN2, data does not change
C
      OPEN (67,FILE='input2.unformatted',FORM='unformatted')
      WRITE(67) NSRA,INS,NATYP,NAEZ,NSPIN,IPAN,IRCUT,LCORE,NCORE,NTCELL,
     &          LMAX,LPOT,LMPOT,NLBASIS,NRBASIS,NRIGHT,NLEFT,LINTERFACE
      WRITE(67) ATOMIMP,NATOMIMP
      WRITE(67) IMIX,MIXING,QBOUND,FCM,ITDBRY,IRNS,KPRE,KSHAPE,KTE,
     &          KVMAD,KXC,LAMBDA_XC,TXC,ICC,ISHIFT,IXIPOL,LRHOSYM,KFORCE
      WRITE(67) A,B,DRDI,R,THETAS,ZAT,IFUNM,LMSP,RMT,RMTNEW,RWS,IMT,IRC,
     &          IRMIN,IRWS,ITITLE,LLMSP,NFU,HOSTIMP
      WRITE(67) ALAT,KAOEZ,IQAT,NOQ,CONC,GSH,ILM,IMAXSH,TESTC,OPTC,LLY
      CLOSE(67)
C ======================================================================
      END subroutine
      
      end module
