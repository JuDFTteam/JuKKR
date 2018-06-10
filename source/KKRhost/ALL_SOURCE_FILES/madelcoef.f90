SUBROUTINE madelcoef(linterface,lpot,a,b,smat,cleb,icleb,iend,  &
    lpotd,lmpotd,lmxspd,nclebd)

      IMPLICIT NONE
!..
!.. Scalar arguments
      INTEGER LPOT,IEND,LPOTD,LMPOTD,LMXSPD,NCLEBD
      LOGICAL LINTERFACE
!..
!.. Array arguments
      DOUBLE PRECISION A(LMPOTD,LMPOTD),B(LMPOTD)
      DOUBLE PRECISION SMAT(LMXSPD),CLEB(NCLEBD)
      INTEGER ICLEB(NCLEBD,3)
!..
!.. Local scalars
      DOUBLE PRECISION, parameter :: PI=4.0D0*ATAN(1.0D0)
      DOUBLE PRECISION, parameter :: FPI=16.0D0*ATAN(1.0D0)
      INTEGER I,L,L1,L2,LM1,LM2,LM3,LMPOT,LOFLM(LMXSPD),M
!      INTEGER ICALL_madelcoef
!..
!.. Local arrays
      DOUBLE PRECISION DFAC(0:LPOTD,0:LPOTD)
!..
!.. Data statements
!      DATA ICALL_madelcoef /0/
!       integer, save :: icall_madelcoef=0
!..
!.. Intrinsic functions
      INTRINSIC ABS,DBLE
!     ..................................................................

lmpot = (lpot+1)**2

i = 1

! --> determine the l-value for given lm

DO l = 0,2*lpot
  DO m = -l,l
    loflm(i) = l
    i = i + 1
  END DO
END DO

! --> calculate:                             (2*(l+l')-1)!!
!                 dfac(l,l') = 4pi**2 *  ----------------------
!                                        (2*l+1)!! * (2*l'+1)!!

dfac(0,0) = fpi*fpi
DO l1 = 1,lpot
  dfac(l1,0) = dfac(l1-1,0)*DBLE(2*l1-1)/DBLE(2*l1+1)
  dfac(0,l1) = dfac(l1,0)
  DO l2 = 1,l1
    dfac(l1,l2) = dfac(l1,l2-1)*DBLE(2*(l1+l2)-1)/DBLE(2*l2+1)
    dfac(l2,l1) = dfac(l1,l2)
  END DO
END DO

! --> initialize

DO lm1 = 1,lmpot
  DO lm2 = 1,lmpot
!            write(*,*) 'test',LM1,LM2,LMPOT,LMPOTD
    a(lm1,lm2) = 0.0D0
  END DO
END DO

! --> calculate a(lm1,lm2)

DO i = 1,iend
  lm1 = icleb(i,1)
  lm2 = icleb(i,2)
  lm3 = icleb(i,3)
  l1 = loflm(lm1)
  l2 = loflm(lm2)
  
! --> this loop has to be calculated only for l1+l2=l3
  
  a(lm1,lm2) = a(lm1,lm2) + 2.0D0*dfac(l1,l2)*smat(lm3)*cleb(i)
END DO

IF ( linterface ) RETURN

! --> initialize

DO lm1 = 1,lmpot
  b(lm1) = 0.0D0
END DO

! --> calculate b(lm1)

DO lm1 = 1,lmpot
  l1 = loflm(lm1)
  b(lm1) = b(lm1) - 2.0D0*fpi/DBLE(2*l1+1)*smat(lm1)
END DO
END SUBROUTINE madelcoef
