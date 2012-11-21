SUBROUTINE WUNFILES_NEW(NPOL,NPNT1,NPNT2, &
NPNT3,IELAST,TK,E1,E2,EZ,WEZ, &
BRAVAIS,RMAX,GMAX, &
EFERMI, &
SCFSTEPS,LCORE,NCORE, &
NSRA,NREF, &
NCLS,ICST,ALAT,ZAT, &
REFPOT,RMTREF,VREF, &
ATOM,CLS,RCLS,NACLS, &
RBASIS,RR,EZOA, &
KMESH,MAXMESH,NSYMAT, &
DSYMLL, &
ITITLE, &
IMIX,MIXING,FCM,KPRE, &
KTE,KVMAD,KXC,ISHIFT, &
KFORCE, &
IEMXD,NPOTD,NAEZD, &
LMMAXD,NREFD, &
NACLSD,NCLSD,NRD, &
NSYMAXD, &
NUMN0,INDN0, &
IGUESS,BCP,QMRBOUND, &
NR,RCUTJIJ,JIJ,LDAU,ISYMINDEX)
  ! **********************************************************************
  ! *                                                                    *
  ! *  This subroutine is part of the MAIN0 program in the tbkkr package *
  ! *  It writes out different unformatted files meant to provide the    *
  ! *  communication between the other parts (MAIN1a, 1b, 1c and 2)      *
  ! *  during an SCF cycle                                               *
  ! *  v.popescu, munich 2004                                            *
  ! *                                                                    *
  ! **********************************************************************
  IMPLICIT NONE
  !     ..
  !     .. Scalar arguments
  INTEGER IEMXD,NPOTD,NAEZD, &
  LMMAXD,NREFD,NACLSD,NCLSD, &
  NRD,NSYMAXD, &
  IGUESS,BCP,NR
  INTEGER IELAST,NPOL,NPNT1,NPNT2,NPNT3,SCFSTEPS
  INTEGER NSRA,NREF, &
  NCLS,ICST
  INTEGER MAXMESH,NSYMAT
  INTEGER IMIX,KPRE, &
  KTE,KVMAD,KXC,ISHIFT,KFORCE
  DOUBLE PRECISION E1,E2,TK,EFERMI,ALAT,RMAX,GMAX
  DOUBLE PRECISION MIXING,FCM,QMRBOUND,RCUTJIJ
  LOGICAL JIJ,LDAU
  !     ..
  !     .. Array arguments
  DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
  DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD)
  DOUBLE PRECISION BRAVAIS(3,3)
  DOUBLE PRECISION ZAT(NAEZD), &
  RMTREF(NREFD), &
  VREF(NAEZD),RCLS(3,NACLSD,NCLSD)
  DOUBLE PRECISION RBASIS(3,NAEZD),RR(3,0:NRD)
  !     ..
  INTEGER LCORE(20,NPOTD),NCORE(NPOTD), &
  ATOM(NACLSD,NAEZD),CLS(NAEZD), &
  NACLS(NCLSD),ISYMINDEX(48)

  INTEGER EZOA(NACLSD,NAEZD), &
  KMESH(IEMXD)
  INTEGER REFPOT(NAEZD)
  INTEGER ITITLE(20,NPOTD)

  INTEGER NUMN0(NAEZD),INDN0(NAEZD,NACLSD)
  !     ..
  !     .. Arrays in common
  CHARACTER*8 TESTC(16),OPTC(8)
  COMMON /TESTC/TESTC
  COMMON /OPTC/OPTC
  !     ..
  ! ---------------------------------------------------------- energy_mesh
  !                                    some data in this file might change

  OPEN (67,FILE='energy_mesh',FORM='unformatted')
  WRITE (67) IELAST,EZ,WEZ,E1,E2
  WRITE (67) NPOL,TK,NPNT1,NPNT2,NPNT3
  !     IF ( NPOL.EQ.0 ) WRITE(67) EFERMI
  WRITE(67) EFERMI
  CLOSE (67)


  OPEN (67,FILE='inpn.unf',FORM='unformatted')
  write (67) KMESH
  write (67) MAXMESH
  write (67) RR
  write (67) EZOA
  write (67) NUMN0
  write (67) INDN0
  write (67) NSYMAT
  write (67) DSYMLL
  write (67) LCORE
  write (67) NCORE
  write (67) IMIX
  write (67) MIXING
  write (67) FCM
  write (67) KPRE
  write (67) KTE
  write (67) KVMAD
  write (67) KXC
  write (67) ISHIFT
  write (67) KFORCE
  write (67) IGUESS
  write (67) BCP
  write (67) QMRBOUND
  write (67) ZAT
  write (67) ITITLE
  write (67) ALAT
  write (67) TESTC
  write (67) OPTC
  write (67) BRAVAIS
  write (67) RBASIS
  write (67) RMAX
  write (67) GMAX
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

! ======================================================================
END
