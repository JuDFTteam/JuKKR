SUBROUTINE bzkint0(nshell,naez,natyp,noq,  &
        rbasis,kaoez,icc,bravais,recbv,atomimp,  &
        rsymat,isymindex,nsymat,ifilimp,natomimp,  &
        nsh1,nsh2,rclsimp,ratom,  &
        ijtabsym,ijtabsh,ijtabcalc,  &
        iofgij,jofgij,nofgij,ish,jsh,  &
        rrot,dsymll,para,qmtet,qmphi,symunitary,  &
        hostimp,intervx,intervy,intervz,  &
        ielast,ez,kmesh,maxmesh,maxmshd,  &
        nsymaxd,krel,lmaxd,lmmaxd,kpoibz,naezd,natypd,  &
        natomimpd,nsheld,nembd)
      IMPLICIT NONE
!.. Parameters ..
      INTEGER NSYMAXD,KREL,LMAXD,LMMAXD
      INTEGER KPOIBZ,NAEZD,NATYPD,NATOMIMPD,NSHELD,NEMBD
!..
!.. Scalar Arguments ..
      INTEGER ICC,NAEZ,NATOMIMP,NATYP,NSYMAT,NOFGIJ
      INTEGER INTERVX,INTERVY,INTERVZ,MAXMESH,MAXMSHD,IELAST
      CHARACTER*40 IFILIMP
!..
!.. Array Arguments ..
DOUBLE COMPLEX DSYMLL(LMMAXD,LMMAXD,NSYMAXD),EZ(*)
DOUBLE PRECISION BRAVAIS(3,3),RATOM(3,NSHELD), &
                 RBASIS(3,NAEZD+NEMBD), &
                 RCLSIMP(3,NATOMIMPD),RECBV(3,3), &
                 RROT(48,3,NSHELD),RSYMAT(64,3,3)
INTEGER ATOMIMP(NATOMIMPD),ISYMINDEX(NSYMAXD), &
        KAOEZ(NATYPD,NAEZD+NEMBD),NOQ(NAEZD), &
        KMESH(*),NSH1(*),NSH2(*),NSHELL(0:NSHELD), &
        IJTABSYM(*),IJTABSH(*),IJTABCALC(*),IOFGIJ(*),JOFGIJ(*), &
        ISH(NSHELD,*),JSH(NSHELD,*)

!..anges for impurity 20/02/2004 -- v.popescu according to 
!..                                 n.papanikolaou 

      INTEGER HOSTIMP(0:NATYPD)
!..
!.. Local Scalars ..
      INTEGER I,ISHELL,IU,IPRINT
      LOGICAL LIRR
!..
!.. Local Arrays ..
      CHARACTER*10 ROTNAME(64)
!.. magnetisation angles ..
      REAL*8 QMTET(NAEZD),QMPHI(NAEZD)
!.. unitary/antiunitary symmetry flag
      LOGICAL SYMUNITARY(NSYMAXD),PARA
!..
!.. External Functions ..
      LOGICAL TEST,OPT
      EXTERNAL TEST,OPT
!..
!.. External Subroutines ..
      EXTERNAL BZKMESH,CRTSTAR,FINDGROUP,GFSHELLS,POINTGRP,SYMTAUMAT
!..

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE(1337,'(79(1H=),/,15X,A)')  &
    'BZKINT0: finding symmetry, setting BZ integration'
WRITE (1337,'(79(1H=),/)')
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

CALL pointgrp(rsymat,rotname)
CALL findgroup(bravais,recbv,rbasis,naez,rsymat,rotname,isymindex,  &
    nsymat,para,qmtet,qmphi,symunitary, krel,naezd,nembd,nsymaxd)

lirr = .true.
iprint = 0
IF (test('TAUSTRUC')) iprint = 2

! --> test: full BZ integration

IF ( test('fullBZ  ').OR.opt('NEWSOSOL') ) THEN
  nsymat = 1
  lirr = .false.
  WRITE(1337,'(8X,2A,/)')  &
      'Test option < fullBZ > or Run option < NEWSOSOL >: ',  &
      ' overriding NSYMAT, generate full BZ k-mesh'
END IF

! --> generate BZ k-mesh

CALL bzkmesh(intervx,intervy,intervz,maxmesh,lirr,bravais,recbv,  &
    nsymat,rsymat,isymindex,symunitary,  &
    ielast,ez,kmesh,iprint,krel,kpoibz,maxmshd)

CALL symtaumat(rotname,rsymat,dsymll,nsymat,isymindex,  &
    symunitary,naezd,lmmaxd,naez,lmaxd+1,krel, iprint,nsymaxd)

! Now DSYMLL hold NSYMAT symmetrization matrices

CALL gfshells(icc,natomimp,nsh1,nsh2, ijtabsym,ijtabsh,ijtabcalc,  &
    iofgij,jofgij,nofgij,ish,jsh, nshell,naez,natyp,noq,rbasis,bravais,  &
    ifilimp,ratom,rclsimp, nsymat,isymindex,rsymat,kaoez,atomimp,  &
    rotname,hostimp, &               ! 20.02.2004
    lmaxd,lmmaxd,naezd,natypd,natomimpd,nembd,nsheld)

! -->  creates difference vectors RROT for BZ integration in KKRMAT01

CALL crtstar(ratom,nshell(0),rsymat,nsymat,isymindex,rrot)
! ----------------------------------------------------------------------
IF ( iprint > 2 ) THEN
  DO ishell = 1,nshell(0)
    WRITE (1337,FMT='(I4)') ishell
    WRITE (1337,FMT='((I4,3F10.1))')  &
        (iu, (rrot(iu,i,ishell),i=1,3),iu=1,nsymat)
  END DO
END IF
! ----------------------------------------------------------------------
END SUBROUTINE bzkint0
