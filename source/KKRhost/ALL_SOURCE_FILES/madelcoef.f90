    Subroutine madelcoef(linterface, lpot, a, b, smat, cleb, icleb, iend, &
      lpotd, lmpotd, lmxspd, nclebd)
      Use mod_datatypes, Only: dp

      Implicit None
!..
!.. Scalar arguments
      Integer :: lpot, iend, lpotd, lmpotd, lmxspd, nclebd
      Logical :: linterface
!..
!.. Array arguments
      Real (Kind=dp) :: a(lmpotd, lmpotd), b(lmpotd)
      Real (Kind=dp) :: smat(lmxspd), cleb(nclebd)
      Integer :: icleb(nclebd, 3)
!..
!.. Local scalars
      Real (Kind=dp), Parameter :: pi = 4.0E0_dp*atan(1.0E0_dp)
      Real (Kind=dp), Parameter :: fpi = 16.0E0_dp*atan(1.0E0_dp)
      Integer :: i, l, l1, l2, lm1, lm2, lm3, lmpot, loflm(lmxspd), m
!      INTEGER ICALL_madelcoef
!..
!.. Local arrays
      Real (Kind=dp) :: dfac(0:lpotd, 0:lpotd)
!..
!.. Data statements
!      DATA ICALL_madelcoef /0/
!       integer, save :: icall_madelcoef=0
!..
!.. Intrinsic functions
      Intrinsic :: abs, real
!     ..................................................................

      lmpot = (lpot+1)**2

      i = 1

! --> determine the l-value for given lm

      Do l = 0, 2*lpot
        Do m = -l, l
          loflm(i) = l
          i = i + 1
        End Do
      End Do

! --> calculate:                             (2*(l+l')-1)!!
!                 dfac(l,l') = 4pi**2 *  ----------------------
!                                        (2*l+1)!! * (2*l'+1)!!

      dfac(0, 0) = fpi*fpi
      Do l1 = 1, lpot
        dfac(l1, 0) = dfac(l1-1, 0)*real(2*l1-1, kind=dp)/ &
          real(2*l1+1, kind=dp)
        dfac(0, l1) = dfac(l1, 0)
        Do l2 = 1, l1
          dfac(l1, l2) = dfac(l1, l2-1)*real(2*(l1+l2)-1, kind=dp)/ &
            real(2*l2+1, kind=dp)
          dfac(l2, l1) = dfac(l1, l2)
        End Do
      End Do

! --> initialize

      Do lm1 = 1, lmpot
        Do lm2 = 1, lmpot
!            write(*,*) 'test',LM1,LM2,LMPOT,LMPOTD
          a(lm1, lm2) = 0.0E0_dp
        End Do
      End Do

! --> calculate a(lm1,lm2)

      Do i = 1, iend
        lm1 = icleb(i, 1)
        lm2 = icleb(i, 2)
        lm3 = icleb(i, 3)
        l1 = loflm(lm1)
        l2 = loflm(lm2)

! --> this loop has to be calculated only for l1+l2=l3

        a(lm1, lm2) = a(lm1, lm2) + 2.0E0_dp*dfac(l1, l2)*smat(lm3)*cleb(i)
      End Do

      If (linterface) Return

! --> initialize

      Do lm1 = 1, lmpot
        b(lm1) = 0.0E0_dp
      End Do

! --> calculate b(lm1)

      Do lm1 = 1, lmpot
        l1 = loflm(lm1)
        b(lm1) = b(lm1) - 2.0E0_dp*fpi/real(2*l1+1, kind=dp)*smat(lm1)
      End Do
    End Subroutine
