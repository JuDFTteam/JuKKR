      SUBROUTINE BZKINT0(NSHELL,NAEZ,NATYP,NOQ,
     +                   RBASIS,KAOEZ,ICC,BRAVAIS,RECBV,ATOMIMP,
     +                   RSYMAT,ISYMINDEX,NSYMAT,IFILIMP,NATOMIMP,
     +                   NSH1,NSH2,RCLSIMP,RATOM,
     +                   IJTABSYM,IJTABSH,IJTABCALC,
     +                   IOFGIJ,JOFGIJ,NOFGIJ,ISH,JSH,
     +                   RROT,DSYMLL,PARA,QMTET,QMPHI,SYMUNITARY,
     +                   HOSTIMP,INTERVX,INTERVY,INTERVZ,
     +                   IELAST,EZ,KMESH,MAXMESH,MAXMSHD,
     +                   NSYMAXD,KREL,LMAXD,LMMAXD,KPOIBZ,NAEZD,NATYPD,
     +                   NATOMIMPD,NSHELD,NEMBD)
      IMPLICIT NONE
C     .. Parameters ..
      INTEGER NSYMAXD,KREL,LMAXD,LMMAXD
      INTEGER KPOIBZ,NAEZD,NATYPD,NATOMIMPD,NSHELD,NEMBD
C     ..
C     .. Scalar Arguments ..
      INTEGER ICC,NAEZ,NATOMIMP,NATYP,NSYMAT,NOFGIJ
      INTEGER INTERVX,INTERVY,INTERVZ,MAXMESH,MAXMSHD,IELAST
      CHARACTER*40 IFILIMP
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),EZ(*)
      DOUBLE PRECISION BRAVAIS(3,3),RATOM(3,NSHELD),
     +                 RBASIS(3,NAEZD+NEMBD),
     +                 RCLSIMP(3,NATOMIMPD),RECBV(3,3),
     +                 RROT(48,3,NSHELD),RSYMAT(64,3,3)
      INTEGER ATOMIMP(NATOMIMPD),ISYMINDEX(NSYMAXD),
     +        KAOEZ(NATYPD,NAEZD+NEMBD),NOQ(NAEZD),
     +        KMESH(*),NSH1(*),NSH2(*),NSHELL(0:NSHELD),
     +        IJTABSYM(*),IJTABSH(*),IJTABCALC(*),IOFGIJ(*),JOFGIJ(*),
     +        ISH(NSHELD,*),JSH(NSHELD,*)
C
C     changes for impurity 20/02/2004 -- v.popescu according to 
C                                        n.papanikolaou 
C
      INTEGER HOSTIMP(0:NATYPD)
C     ..
C     .. Local Scalars ..
      INTEGER I,ISHELL,IU,IPRINT
      LOGICAL LIRR
C     ..
C     .. Local Arrays ..
      CHARACTER*10 ROTNAME(64)
C     .. magnetisation angles ..
      REAL*8 QMTET(NAEZD),QMPHI(NAEZD)
C     .. unitary/antiunitary symmetry flag
      LOGICAL SYMUNITARY(NSYMAXD),PARA
C     ..
C     .. External Functions ..
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT
C     ..
C     .. External Subroutines ..
      EXTERNAL BZKMESH,CRTSTAR,FINDGROUP,GFSHELLS,POINTGRP,SYMTAUMAT
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
     +               NSYMAT,PARA,QMTET,QMPHI,SYMUNITARY,
     +               KREL,NAEZD,NEMBD,NSYMAXD)

      LIRR = .TRUE.
      IPRINT = 0    
      IF (TEST('TAUSTRUC')) IPRINT = 2
C
C --> test: full BZ integration
C
      IF ( TEST('fullBZ  ').OR.OPT('NEWSOSOL') ) THEN
         NSYMAT = 1
         LIRR = .FALSE.
         WRITE(6,'(8X,2A,/)') 
     &        'Test option < fullBZ > or Run option < NEWSOSOL >: ',
     &        ' overriding NSYMAT, generate full BZ k-mesh'
      END IF
C
C --> generate BZ k-mesh
C
      CALL BZKMESH(INTERVX,INTERVY,INTERVZ,MAXMESH,LIRR,BRAVAIS,RECBV,
     &             NSYMAT,RSYMAT,ISYMINDEX,SYMUNITARY,
     &             IELAST,EZ,KMESH,IPRINT,KREL,KPOIBZ,MAXMSHD)
C
      CALL SYMTAUMAT(ROTNAME,RSYMAT,DSYMLL,NSYMAT,ISYMINDEX,
     &               SYMUNITARY,NAEZD,LMMAXD,NAEZ,LMAXD+1,KREL,
     &               IPRINT,NSYMAXD)
C
C Now DSYMLL hold NSYMAT symmetrization matrices
C
      CALL GFSHELLS(ICC,NATOMIMP,NSH1,NSH2,
     +              IJTABSYM,IJTABSH,IJTABCALC,
     +              IOFGIJ,JOFGIJ,NOFGIJ,ISH,JSH,
     +              NSHELL,NAEZ,NATYP,NOQ,RBASIS,BRAVAIS,
     +              IFILIMP,RATOM,RCLSIMP,
     +              NSYMAT,ISYMINDEX,RSYMAT,KAOEZ,ATOMIMP,
     +              ROTNAME,IPRINT,HOSTIMP,                 ! 20.02.2004
     +              LMAXD,LMMAXD,NAEZD,NATYPD,NATOMIMPD,NEMBD,NSHELD)
C
C -->  creates difference vectors RROT for BZ integration in KKRMAT01
C
      CALL CRTSTAR(RATOM,NSHELL(0),RSYMAT,NSYMAT,ISYMINDEX,RROT)
C ----------------------------------------------------------------------
      IF ( IPRINT.GT.2 ) THEN
         DO ISHELL = 1,NSHELL(0)
            WRITE (6,FMT='(I4)') ISHELL
            WRITE (6,FMT='((I4,3F10.1))') 
     &           (IU, (RROT(IU,I,ISHELL),I=1,3),IU=1,NSYMAT)
         END DO
      END IF
C ----------------------------------------------------------------------
      END
