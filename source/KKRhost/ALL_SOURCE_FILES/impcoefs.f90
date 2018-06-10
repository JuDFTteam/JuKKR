SUBROUTINE impcoefs(natomimp,naez,atomimp,rclsimp,nshell,nsh1,  &
    nsh2,ratom,nsymat,isymindex,rotname,hostimp, natypd,lmaxd,nsheld,nsize)
! **********************************************************************
! *                                                                    *
! * Writes out the auxiliary file impurity.coefs which is needed for   *
! * impurity calculations                                              *
! * Sets up the array HOSTIMP -- also needed for impurity case         *
! *                                adopted from N. Papanikolaou        *
! **********************************************************************

      IMPLICIT NONE
!.. 
!.. Scalar arguments
      INTEGER LMAXD,NAEZ,NATOMIMP,NATYPD,NSHELD,NSIZE,NSYMAT
!..
!.. Array arguments
      INTEGER ATOMIMP(*),HOSTIMP(0:NATYPD),ISYMINDEX(*),NSH1(*),NSH2(*), &
     &        NSHELL(0:NSHELD)
      DOUBLE PRECISION RATOM(3,NSHELD),RCLSIMP(3,*)
      CHARACTER*10 ROTNAME(*)
!..
!.. Local scalars
      INTEGER AI,I,II,J,NB,NDIM,NHOST,NREP,NS
      DOUBLE PRECISION R1
!..
!.. Local arrays
      INTEGER IMPHOST(NAEZ),NSHOUT(NATOMIMP)
      LOGICAL EXIST(NAEZ)
!     ..

! -->  shells around atom icc are prepared for storing the
!      cluster-gf in subroutine kkrmat in GMATLL(LMMAXD,LMMAXD,*)

DO i = 1,naez
  EXIST(i) = .false.
END DO
DO i = 1,natomimp
  EXIST(atomimp(i)) = .true.
END DO

nhost = 0
DO i = 1,naez
  imphost(i) = 0
  IF ( EXIST(i) ) THEN
    nhost = nhost + 1
    imphost(i) = nhost
    hostimp(nhost) = i
  END IF
END DO
hostimp(0) = nhost
IF ( nhost /= naez ) WRITE (6,99001)

nrep = 1
ndim = 1

DO i = 1,natomimp
  nshout(i) = 1
END DO

OPEN (58,FILE='impurity.coefs',FORM='FORMATTED')
WRITE (58,99002) nrep,natomimp,lmaxd,natomimp, (nshout(i),i=1,natomimp)
WRITE (58,99003)
!-----------------------------------------------------------------------
DO i = 1,natomimp
  
  r1 = SQRT(rclsimp(1,i)**2+rclsimp(2,i)**2+rclsimp(3,i)**2)
  
  IF ( naez == nhost ) THEN
    ai = atomimp(i)
  ELSE
    ai = imphost(atomimp(i))
  END IF
  
  WRITE (58,99004) (rclsimp(j,i),j=1,3),ai,i,i,r1,atomimp(i)
END DO
!-----------------------------------------------------------------------

nb = 0
DO ns = 1,nshell(0)
  nb = nb + nshell(ns)
END DO

WRITE (58,99011) nsize,nb
WRITE (58,99002) ndim
WRITE (58,99005) nhost
WRITE (58,99006) (hostimp(i),i=1,nhost)
WRITE (58,99007) nsymat
WRITE (58,99008) (rotname(isymindex(i)),i=1,nsymat)
WRITE (58,99009)
WRITE (58,99002) nshell(0)
WRITE (58,99010) (ns,nsh1(ns),nsh2(ns),(ratom(ii,ns),ii=1,3), nshell(ns),  &
    SQRT(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns) **2),ns=1,nshell(0))
CLOSE (58)
! ======================================================================
99001 FORMAT (8X,'WARNING: Some host atoms are missing in the ',  &
    'impurity cluster',/,8X, '         Indexing will be changed. Check ',  &
    'impurity.coefs file?',/)
99002 FORMAT (11I5)
99003 FORMAT ('     Position of Impurity            Host Imp Shell',  &
    '   Dist     Host id in Bulk')
99004 FORMAT (3F12.8,i4,i4,i5,f10.6,i5)
99005 FORMAT ('Host order, no of host sites: ',i5)
99006 FORMAT (12I4)
99007 FORMAT (i5,'    Symmetries for the Bulk')
99008 FORMAT (5A10)
99009 FORMAT ('Shells in the reduced format')
99010 FORMAT (3I5,3F12.8,i8,f10.5)
99011 FORMAT (11I20)
END SUBROUTINE impcoefs
