subroutine madelcoef(linterface, lpot, a, b, smat, cleb, icleb, iend, lpotd, &
  lmpotd, lmxspd, nclebd)

  implicit none
!..
!.. Scalar arguments
  integer :: lpot, iend, lpotd, lmpotd, lmxspd, nclebd
  logical :: linterface
!..
!.. Array arguments
  double precision :: a(lmpotd, lmpotd), b(lmpotd)
  double precision :: smat(lmxspd), cleb(nclebd)
  integer :: icleb(nclebd, 3)
!..
!.. Local scalars
  double precision, parameter :: pi = 4.0d0*atan(1.0d0)
  double precision, parameter :: fpi = 16.0d0*atan(1.0d0)
  integer :: i, l, l1, l2, lm1, lm2, lm3, lmpot, loflm(lmxspd), m
!      INTEGER ICALL_madelcoef
!..
!.. Local arrays
  double precision :: dfac(0:lpotd, 0:lpotd)
!..
!.. Data statements
!      DATA ICALL_madelcoef /0/
!       integer, save :: icall_madelcoef=0
!..
!.. Intrinsic functions
  intrinsic :: abs, dble
!     ..................................................................

  lmpot = (lpot+1)**2

  i = 1

! --> determine the l-value for given lm

  do l = 0, 2*lpot
    do m = -l, l
      loflm(i) = l
      i = i + 1
    end do
  end do

! --> calculate:                             (2*(l+l')-1)!!
!                 dfac(l,l') = 4pi**2 *  ----------------------
!                                        (2*l+1)!! * (2*l'+1)!!

  dfac(0, 0) = fpi*fpi
  do l1 = 1, lpot
    dfac(l1, 0) = dfac(l1-1, 0)*dble(2*l1-1)/dble(2*l1+1)
    dfac(0, l1) = dfac(l1, 0)
    do l2 = 1, l1
      dfac(l1, l2) = dfac(l1, l2-1)*dble(2*(l1+l2)-1)/dble(2*l2+1)
      dfac(l2, l1) = dfac(l1, l2)
    end do
  end do

! --> initialize

  do lm1 = 1, lmpot
    do lm2 = 1, lmpot
!            write(*,*) 'test',LM1,LM2,LMPOT,LMPOTD
      a(lm1, lm2) = 0.0d0
    end do
  end do

! --> calculate a(lm1,lm2)

  do i = 1, iend
    lm1 = icleb(i, 1)
    lm2 = icleb(i, 2)
    lm3 = icleb(i, 3)
    l1 = loflm(lm1)
    l2 = loflm(lm2)

! --> this loop has to be calculated only for l1+l2=l3

    a(lm1, lm2) = a(lm1, lm2) + 2.0d0*dfac(l1, l2)*smat(lm3)*cleb(i)
  end do

  if (linterface) return

! --> initialize

  do lm1 = 1, lmpot
    b(lm1) = 0.0d0
  end do

! --> calculate b(lm1)

  do lm1 = 1, lmpot
    l1 = loflm(lm1)
    b(lm1) = b(lm1) - 2.0d0*fpi/dble(2*l1+1)*smat(lm1)
  end do
end subroutine
