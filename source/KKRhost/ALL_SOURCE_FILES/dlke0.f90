! 04.10.95 *************************************************************
SUBROUTINE dlke0(gllke,alat,naez,cls,nacls,naclsmax,rr,ezoa,atom,  &
    bzkp,rcls,ginp)
! **********************************************************************
IMPLICIT NONE
!     .. Parameters ..
INCLUDE 'inc.p'
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
! *          function, set up in the spin-independent non-relativstic *
! *          (l,m_l)-representation                                   *
! *                                                                   *
! *********************************************************************
!     ..
      INTEGER LMAX
      PARAMETER (LMAX=LMAXD)
      INTEGER LMGF0D
      PARAMETER (LMGF0D= (LMAX+1)**2)
      INTEGER ALMGF0
      PARAMETER (ALMGF0=LMGF0D*NAEZD)
!..
!.. Scalar Arguments ..
      DOUBLE PRECISION ALAT
      INTEGER NAEZ,NACLSMAX
!..
!.. Array Arguments ..


      DOUBLE COMPLEX GINP(LMGF0D*NACLSMAX,LMGF0D,*),GLLKE(ALMGF0,*)
      DOUBLE PRECISION BZKP(*),RCLS(3,NACLSD,*),RR(3,0:NRD)
      INTEGER ATOM(NACLSD,*),CLS(*),EZOA(NACLSD,*),NACLS(*)
!..
!.. Local Scalars ..
      INTEGER I,IC,IM,J,JN,M,N
!..
!.. Local Arrays ..
      DOUBLE COMPLEX GLLKE1(ALMGF0,LMGF0D)
      DOUBLE PRECISION KP(6)
!..
!.. External Subroutines ..
      EXTERNAL CINIT,DLKE1
!..
!.. Save statement ..
      SAVE
!      write(6,*) '>>> DLKE0 : Fourier-transforms the ',
!     +           'GF of reference system'
! ----------------------------------------------------------------------

!     .. External Functions ..
LOGICAL :: opt
EXTERNAL opt
!     ..
CALL cinit(almgf0*almgf0,gllke(1,1))

DO  i = 1,naez
  
  
  kp(1) = bzkp(1)
  kp(2) = bzkp(2)
  kp(3) = bzkp(3)
  IF (opt('COMPLEX ')) THEN
    kp(4) = bzkp(4)
    kp(5) = bzkp(5)
    kp(6) = bzkp(6)
  END IF
  
  ic = cls(i)
  CALL dlke1(gllke1,alat,nacls,naclsmax,rr,ezoa(1,i),atom(1,i),kp,  &
      ic,ginp(1,1,ic),rcls(1,1,ic))
  
  DO  m = 1,lmgf0d
    im = (i-1)*lmgf0d + m
    DO  jn = 1,lmgf0d*naez
      gllke(jn,im) = gllke(jn,im) + gllke1(jn,m)
    END DO
  END DO
  
  
  
END DO

! ----------------------------------------------------------------------
IF (opt('symG(k) ')) THEN
  
! -->   symmetrization
  
  DO  i = 1,naez
    
    kp(1) = -bzkp(1)
    kp(2) = -bzkp(2)
    kp(3) = -bzkp(3)
    IF (opt('COMPLEX ')) THEN
      kp(4) = -bzkp(4)
      kp(5) = -bzkp(5)
      kp(6) = -bzkp(6)
    END IF
    
    ic = cls(i)
    CALL dlke1(gllke1,alat,nacls,naclsmax,rr,ezoa(1,i),atom(1,i),  &
        kp,ic,ginp(1,1,ic),rcls(1,1,ic))
    
    DO  j = 1,naez
      DO  m = 1,lmgf0d
        im = (i-1)*lmgf0d + m
        DO  n = 1,lmgf0d
          jn = (j-1)*lmgf0d + n
          gllke(im,jn) = (gllke(im,jn)+gllke1(jn,m))/2.0D0
        END DO
      END DO
    END DO
    
  END DO
  
END IF
! ----------------------------------------------------------------------

RETURN

END SUBROUTINE dlke0
