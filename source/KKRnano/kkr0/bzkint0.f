      SUBROUTINE BZKINT0(NAEZ,
     +                   RBASIS,BRAVAIS,RECBV,
     +                   NSYMAT,ISYMINDEX,
     +                   DSYMLL,
     +                   INTERVX,INTERVY,INTERVZ,
     +                   IELAST,EZ,KMESH,MAXMESH,MAXMSHD)
C
      IMPLICIT NONE
C     .. Parameters ..
      INCLUDE 'inc.p'
      INCLUDE 'inc.cls'
      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMMAXD
      PARAMETER (LMMAXD= (LMAXD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER NAEZ,NSYMAT
      INTEGER INTERVX,INTERVY,INTERVZ,MAXMESH,MAXMSHD,IELAST
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),EZ(IEMXD)
      DOUBLE PRECISION BRAVAIS(3,3),
     +                 RBASIS(3,NAEZD),
     +                 RECBV(3,3),
     +                 RSYMAT(64,3,3)
      INTEGER ISYMINDEX(NSYMAXD)
      INTEGER KMESH(IEMXD)
C
C     changes for impurity 20/02/2004 -- v.popescu according to 
C                                        n.papanikolaou 
C
C     ..
C     .. Local Scalars ..
      INTEGER IPRINT
      LOGICAL LIRR
C     ..
C     .. Local Arrays ..
      CHARACTER*10 ROTNAME(64)
C     ..
C     .. External Functions ..
      LOGICAL TEST
      EXTERNAL TEST
C     ..
C     .. External Subroutines ..
      EXTERNAL BZKMESH,FINDGROUP,POINTGRP,SYMTAUMAT
C     ..
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE(6,'(79(1H=),/,15X,A)') 
     &     'BZKINT0: finding symmetry, setting BZ integration'
      WRITE (6,'(79(1H=),/)')
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      CALL POINTGRP(RSYMAT,ROTNAME)
      CALL FINDGROUP(BRAVAIS,RECBV,RBASIS,NAEZ,RSYMAT,ROTNAME,ISYMINDEX,
     +               NSYMAT)

      LIRR = .TRUE.
      IPRINT = 0    
      IF (TEST('TAUSTRUC')) IPRINT = 2
C
C --> test: full BZ integration
C
      IF ( TEST('fullBZ  ') ) THEN
         NSYMAT = 1
         LIRR = .FALSE.
         WRITE(6,'(8X,2A,/)') 
     &        'Test option < fullBZ > : overriding NSYMAT,',
     &        ' generate full BZ k-mesh'
      END IF
C
C --> generate BZ k-mesh
C
      CALL BZKMESH(INTERVX,INTERVY,INTERVZ,MAXMESH,LIRR,BRAVAIS,RECBV,
     &             NSYMAT,RSYMAT,ISYMINDEX,
     &             IELAST,EZ,KMESH,IPRINT,MAXMSHD)
C
      CALL SYMTAUMAT(ROTNAME,RSYMAT,DSYMLL,NSYMAT,ISYMINDEX,
     &               NAEZD,LMMAXD,NAEZ,LMAX+1,KREL,
     &               IPRINT,NSYMAXD)
C
C Now DSYMLL hold NSYMAT symmetrization matrices
C

C ----------------------------------------------------------------------
      END
