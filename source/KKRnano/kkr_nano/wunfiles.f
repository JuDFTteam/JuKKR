      SUBROUTINE WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ,
     &                    BRAVAIS,RECBV,VOLUME0,RMAX,GMAX,
     &                    EFERMI,VBC,
     &                    SCFSTEPS,LCORE,NCORE,
     &                    NSRA,NAEZ,NREF,NSPIN,LMAX,
     &                    NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI,
     &                    REFPOT,RMTREF,VREF,IEND,JEND,CLEB,ICLEB,
     &                    ATOM,CLS,RCLS,NACLS,LOFLM,
     &                    RBASIS,RR,EZOA,
     &                    KMESH,MAXMESH,NSYMAT,
     &                    DSYMLL,
     &                    A,B,IFUNM1,
     &                    ITITLE,LMSP1,NTCELL,THETAS,
     &                    LPOT,LMPOT,
     &                    IMIX,MIXING,QBOUND,FCM,KPRE,
     &                    KTE,KVMAD,KXC,ISHIFT,
     &                    KFORCE,LLMSP1,IMT,
     &                    IRC,IRMIN,IRNS,IRWS,RWS,RMT,NFU,GSH,ILM,
     &                    IMAXSH,IEMXD,IRMD,LMPOTD,NPOTD,NAEZD,
     &                    LMMAXD,IPAND,NREFD,LMAXD,
     &                    NCLEB,NACLSD,NCLSD,LM2D,NRD,
     &                    NSYMAXD,IRID,NFUND,
     &                    NCELLD,LMXSPD,NGSHD,NUMN0,INDN0,
     &                    IGUESS,BCP,QMRBOUND,
     &                    NR,RCUTJIJ,JIJ,LDAU,ISYMINDEX)
C **********************************************************************
C *                                                                    *
C *  This subroutine is part of the MAIN0 program in the tbkkr package *
C *  It writes out different unformatted files meant to provide the    *
C *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
C *  during an SCF cycle                                               *
C *  v.popescu, munich 2004                                            *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar arguments
      INTEGER IEMXD,IRMD,LMPOTD,NPOTD,NAEZD,
     &        LMMAXD,IPAND,NREFD,LMAXD,NCLEB,NACLSD,NCLSD,
     &        LM2D,NRD,NSYMAXD,
     &        IRID,NFUND,NCELLD,LMXSPD,NGSHD,IGUESS,BCP,NR
      INTEGER IELAST,NPOL,NPNT1,NPNT2,NPNT3,SCFSTEPS
      INTEGER NSRA,NAEZ,NREF,NSPIN,LMAX,
     &        NCLS,ICST,IEND
      INTEGER MAXMESH,NSYMAT
      INTEGER LPOT,LMPOT,IMIX,KPRE,
     &        KTE,KVMAD,KXC,ISHIFT,KFORCE
      DOUBLE PRECISION E1,E2,TK,EFERMI,ALAT,VOLUME0,RMAX,GMAX
      DOUBLE PRECISION MIXING,QBOUND,FCM,QMRBOUND,RCUTJIJ
      LOGICAL JIJ,LDAU
C     ..
C     .. Array arguments
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
      DOUBLE PRECISION VBC(2)
      DOUBLE PRECISION BRAVAIS(3,3),RECBV(3,3)
      DOUBLE PRECISION ZAT(NAEZD),R(IRMD,NAEZD),
     &                 DRDI(IRMD,NAEZD),RMTREF(NREFD),
     &                 VREF(NAEZD),CLEB(NCLEB,2),RCLS(3,NACLSD,NCLSD)
      DOUBLE PRECISION RBASIS(3,NAEZD),RR(3,0:NRD)
      DOUBLE PRECISION A(NAEZD),B(NAEZD),THETAS(IRID,NFUND,NCELLD)
      DOUBLE PRECISION GSH(NGSHD),RWS(NAEZD),RMT(NAEZD)
C     ..
      INTEGER LCORE(20,NPOTD),NCORE(NPOTD),IPAN(NAEZD),
     &        IRCUT(0:IPAND,NAEZD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     &        ICLEB(NCLEB,3),ATOM(NACLSD,NAEZD),CLS(NAEZD),
     &        NACLS(NCLSD),ISYMINDEX(48),
     &        IRNS(NAEZD)
      INTEGER LOFLM(LM2D),EZOA(NACLSD,NAEZD),
     &        KMESH(IEMXD)
      INTEGER REFPOT(NAEZD)
      INTEGER IFUNM1(LMXSPD,NAEZD),
     &        ITITLE(20,NPOTD),LMSP1(LMXSPD,NAEZD),NTCELL(NAEZD)
      INTEGER LLMSP1(NFUND,NAEZD),IMT(NAEZD),
     &        IRC(NAEZD),IRMIN(NAEZD),IRWS(NAEZD),NFU(NAEZD),
     &        ILM(NGSHD,3),IMAXSH(0:LMPOTD)
      INTEGER NUMN0(NAEZD),INDN0(NAEZD,NACLSD)
C     ..
C     .. Arrays in common
      CHARACTER*8 TESTC(16),OPTC(8)
      COMMON /TESTC/TESTC
      COMMON /OPTC/OPTC
C     ..
C ---------------------------------------------------------- energy_mesh
C                                    some data in this file might change
C
      OPEN (67,FILE='energy_mesh',FORM='unformatted')
      WRITE (67) IELAST,EZ,WEZ,E1,E2
      WRITE (67) NPOL,TK,NPNT1,NPNT2,NPNT3
      IF ( NPOL.EQ.0 ) WRITE(67) EFERMI
      CLOSE (67)
C ------------------------------------------------------ input_potential
C                                          some data in this file change
C
C
C
      E1 = 0D0
C
      OPEN (28,FILE='not.converged',FORM='formatted')
      WRITE (28,'(1P,3D17.10)') E1,VBC
      CLOSE (28)
C -------------------------------------------------------------- input1b
C                                 meant for MAIN1b, data does not change
C
      OPEN (67,FILE='inp.unf',FORM='unformatted')
      WRITE(67) KMESH,MAXMESH,RR,EZOA,NUMN0,INDN0,NSYMAT,DSYMLL
      WRITE(67) NAEZ,NSPIN,IPAN,IRNS,IRCUT,LCORE,NCORE,NTCELL,
     &          LMAX,LPOT,LMPOT
      WRITE(67) IMIX,MIXING,QBOUND,FCM,KPRE,KTE,
     &          KVMAD,KXC,ISHIFT,KFORCE,IGUESS,BCP,QMRBOUND
      WRITE(67) A,B,DRDI,R,THETAS,ZAT,IMT,IRC,
     &          IRMIN,IRWS,RWS,RMT,ITITLE,LLMSP1,NFU
      WRITE(67) ALAT,GSH,ILM,IMAXSH,TESTC,OPTC
      WRITE(67) BRAVAIS,RBASIS,RECBV,VOLUME0,RMAX,GMAX
      WRITE(67) IEND,CLEB,ICLEB,LOFLM,JEND,IFUNM1,LMSP1,NSRA,ICST
      WRITE(67) NCLS,NREF,RMTREF,VREF,RCLS,
     &          ATOM,CLS,NACLS,REFPOT
      WRITE(67) NR,RCUTJIJ,JIJ,LDAU
      WRITE(67) ISYMINDEX,SCFSTEPS
      CLOSE(67)

C ======================================================================
      END
