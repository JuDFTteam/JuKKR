      SUBROUTINE WUNFILES(NPOL,NPNT1,NPNT2,NPNT3,IELAST,TK,E1,E2,EZ,WEZ,
     &                    BRAVAIS,RMAX,GMAX,
     &                    EFERMI,
     &                    SCFSTEPS,LCORE,NCORE,
     &                    NSRA,NREF,
     &                    NCLS,ICST,IPAN,IRCUT,ALAT,ZAT,R,DRDI,
     &                    REFPOT,RMTREF,VREF,
     &                    ATOM,CLS,RCLS,NACLS,
     &                    RBASIS,RR,EZOA,
     &                    KMESH,MAXMESH,NSYMAT,
     &                    DSYMLL,
     &                    A,B,IFUNM1,
     &                    ITITLE,LMSP1,NTCELL,THETAS,
     &                    IMIX,MIXING,FCM,KPRE,
     &                    KTE,KVMAD,KXC,ISHIFT,
     &                    KFORCE,LLMSP1,IMT,
     &                    IRC,IRMIN,IRNS,IRWS,RWS,RMT,NFU,
     &                    IEMXD,IRMD,LMPOTD,NPOTD,NAEZD,
     &                    LMMAXD,IPAND,NREFD,LMAXD,
     &                    NACLSD,NCLSD,NRD,
     &                    NSYMAXD,IRID,NFUND,
     &                    NCELLD,LMXSPD,NUMN0,INDN0,
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
     &        LMMAXD,IPAND,NREFD,LMAXD,NACLSD,NCLSD,
     &        NRD,NSYMAXD,
     &        IRID,NFUND,NCELLD,LMXSPD,IGUESS,BCP,NR
      INTEGER IELAST,NPOL,NPNT1,NPNT2,NPNT3,SCFSTEPS
      INTEGER NSRA,NREF,NSPIN,LMAX,
     &        NCLS,ICST
      INTEGER MAXMESH,NSYMAT
      INTEGER LPOT,LMPOT,IMIX,KPRE,
     &        KTE,KVMAD,KXC,ISHIFT,KFORCE
      DOUBLE PRECISION E1,E2,TK,EFERMI,ALAT,RMAX,GMAX
      DOUBLE PRECISION MIXING,QBOUND,FCM,QMRBOUND,RCUTJIJ
      LOGICAL JIJ,LDAU
C     ..
C     .. Array arguments
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
      DOUBLE PRECISION BRAVAIS(3,3)
      DOUBLE PRECISION ZAT(NAEZD),R(IRMD,NAEZD),
     &                 DRDI(IRMD,NAEZD),RMTREF(NREFD),
     &                 VREF(NAEZD),RCLS(3,NACLSD,NCLSD)
      DOUBLE PRECISION RBASIS(3,NAEZD),RR(3,0:NRD)
      DOUBLE PRECISION A(NAEZD),B(NAEZD),THETAS(IRID,NFUND,NCELLD)
      DOUBLE PRECISION RWS(NAEZD),RMT(NAEZD)
C     ..
      INTEGER LCORE(20,NPOTD),NCORE(NPOTD),IPAN(NAEZD),
     &        IRCUT(0:IPAND,NAEZD),
     &        ATOM(NACLSD,NAEZD),CLS(NAEZD),
     &        NACLS(NCLSD),ISYMINDEX(48),
     &        IRNS(NAEZD)
      INTEGER EZOA(NACLSD,NAEZD),
     &        KMESH(IEMXD)
      INTEGER REFPOT(NAEZD)
      INTEGER IFUNM1(LMXSPD,NAEZD),
     &        ITITLE(20,NPOTD),LMSP1(LMXSPD,NAEZD),NTCELL(NAEZD)
      INTEGER LLMSP1(NFUND,NAEZD),IMT(NAEZD),
     &        IRC(NAEZD),IRMIN(NAEZD),IRWS(NAEZD),NFU(NAEZD)
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
C     IF ( NPOL.EQ.0 ) WRITE(67) EFERMI
      WRITE(67) EFERMI
      CLOSE (67)
C ------------------------------------------------------ input_potential
C                                          some data in this file change
C
C
C
C     E1 = 0D0
C
C     OPEN (28,FILE='not.converged',FORM='formatted')
C     WRITE (28,'(1P,3D17.10)') E1,VBC
C     CLOSE (28)
C -------------------------------------------------------------- input1b
C                                 meant for MAIN1b, data does not change
C
      OPEN (67,FILE='inp.unf',FORM='unformatted')
      write (67) KMESH
      write (67) MAXMESH
      write (67) RR
      write (67) EZOA
      write (67) NUMN0
      write (67) INDN0
      write (67) NSYMAT
      write (67) DSYMLL
      write (67) IPAN
      write (67) IRNS
      write (67) IRCUT
      write (67) LCORE
      write (67) NCORE
      write (67) NTCELL
      write (67) IMIX
      write (67) MIXING
      write (67) FCM
      write (67) KPRE
      write (67) KTE
      write (67) KVMAD
      write (67) KXC
      write (67)  ISHIFT
      write (67) KFORCE
      write (67) IGUESS
      write (67) BCP
      write (67) QMRBOUND
      write (67) A
      write (67) B
      write (67) DRDI
      write (67) R
      write (67) THETAS
      write (67) ZAT
      write (67) IMT
      write (67) IRC
      write (67) IRMIN
      write (67) IRWS
      write (67) RWS
      write (67)  RMT
      write (67) ITITLE
      write (67) LLMSP1
      write (67) NFU
      write (67)  ALAT
      write (67) TESTC
      write (67) OPTC
      write (67) BRAVAIS
      write (67) RBASIS
      write (67) RMAX
      write (67) GMAX
      write (67) IFUNM1
      write (67) LMSP1
      write (67) NSRA
      write (67) ICST
      write (67) NCLS
      write (67) NREF
      write (67) RMTREF
      write (67) VREF
      write (67) RCLS
      write (67) ATOM
      write (67) CLS
      write (67) NACLS
      write (67) REFPOT
      write (67) NR
      write (67) RCUTJIJ
      write (67) JIJ
      write (67) LDAU
      write (67) ISYMINDEX
      write (67) SCFSTEPS
      CLOSE(67)

C ======================================================================
      END
