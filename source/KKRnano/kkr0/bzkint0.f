      SUBROUTINE BZKINT0(NAEZ,
     +                   RBASIS,BRAVAIS,RECBV,
     +                   NSYMAT,ISYMINDEX,
     +                   DSYMLL,
     +                   INTERVX,INTERVY,INTERVZ,
     +                   IELAST,EZ,KMESH,MAXMESH,MAXMSHD,
C                        new after inc.p replacement
     &                   LMAX, IEMXD, KREL, KPOIBZ,
C                        output: (required EKMD value)
     &                   EKMD)
C
      IMPLICIT NONE

C     new after inc.p replace
      INTEGER LMAX
      INTEGER IEMXD
      INTEGER KREL
      INTEGER KPOIBZ
      INTEGER EKMD

      INTEGER NSYMAXD
      PARAMETER (NSYMAXD=48)
C     ..
C     .. Scalar Arguments ..
      INTEGER NAEZ,NSYMAT
      INTEGER INTERVX,INTERVY,INTERVZ,MAXMESH,MAXMSHD,IELAST
C     ..
C     .. Array Arguments ..
C     DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),EZ(IEMXD)
      DOUBLE COMPLEX DSYMLL((LMAX+1)**2,(LMAX+1)**2,NSYMAXD),EZ(IEMXD)
      DOUBLE PRECISION BRAVAIS(3,3),
     +                 RBASIS(3,NAEZ),
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

      INTEGER LMMAXD
      INTEGER NAEZD

      LMMAXD= (LMAX+1)**2
      NAEZD = NAEZ
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
     &             IELAST,EZ,KMESH,IPRINT,MAXMSHD,
     &             IEMXD, KPOIBZ, EKMD)
C
      CALL SYMTAUMAT(ROTNAME,RSYMAT,DSYMLL,NSYMAT,ISYMINDEX,
     &               NAEZD,LMMAXD,NAEZ,LMAX+1,KREL,
     &               IPRINT,NSYMAXD)
C
C Now DSYMLL hold NSYMAT symmetrization matrices
C

C ----------------------------------------------------------------------
      END
