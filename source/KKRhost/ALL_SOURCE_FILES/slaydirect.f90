!-------------------------------------------------------------------------------
! SUBROUTINE: SLAYDIRECT
!> @brief This subroutine returns
!> \f$ SUM_L^{i_1,i_2}=\sum_n \frac{Y_L\left(R^{i_1}-R_{n}^{i_2}\right)}{\left|R^{i_1}-R^{i_2}_n \right|^{l+1}}\f$
!
!> @details  \f$i_1\f$, \f$i_2\f$ are plane indeces and the sum is done in the
!> \f$i_2\f$ plane
!>
!> In case \f$i_1=i_2\f$ the \f$r=0\f$ term is excluded!
!>
!> Error check is performed!
!> @note Jonathan Chico Apr. 2018: Removed inc.p dependencies, rewrote to Fortran90
!> and renamed the SUM array to SUML to avoid confusion with intrinsic SUM function
!-------------------------------------------------------------------------------
subroutine slaydirect(lpot, vec1, vec2, alat, br, suml)

  use :: constants
  use :: global_variables

  implicit none

! .. Parameters
  integer :: lmpotd
  integer :: l2potd
  integer :: lm2potd
  parameter (l2potd=2*lpotd)
  parameter (lmpotd=(lpotd+1)**2)
  parameter (lm2potd=(2*lpotd+1)**2)
! .. Input variables
  double precision, intent (in) :: alat !< Lattice constant in a.u.
  double precision, dimension (3), intent (in) :: vec1
  double precision, dimension (3), intent (in) :: vec2
  double precision, dimension (3, 3), intent (in) :: br
! .. Output variables
  double precision, dimension (lm2potd), intent (out) :: suml
! .. Local variables
  integer :: i1, i2
  integer :: l, lm, m, k, lpot, n
  integer :: i, j, lm2pot, n1, maxn, ia, n2
  integer :: ipar
  double precision :: zzz, r, cc, zz, cr, x0, y0, r2
  double precision :: zoffx, zoffy, rmax
  double precision :: r0, a1, a2, b1, b2, rtest, ccc
  logical :: ltest, test
  double precision, dimension (lm2potd) :: ylm
  double precision, dimension (lm2potd) :: sumt
!----------------------------------------------------------------------------
  do lm = 1, lm2potd
    suml(lm) = 0.d0
    sumt(lm) = 0.d0
  end do
! Run this sub only in the test option "electro"
  if (.not. test('electro ')) return

  zoffx = (vec2(1)-vec1(1))*alat
  zoffy = (vec2(2)-vec1(2))*alat
  zz = (vec2(3)-vec1(3))*alat
  zzz = zz*zz
  maxn = 130
  ia = 0
  ipar = 0
  ltest = .true.
  a1 = br(1, 1)
  a2 = br(2, 1)
  b1 = br(1, 2)
  b2 = br(2, 2)
  rmax = maxn*min(sqrt(a1*a1+a2*a2), sqrt(b1*b1+b2*b2))
  rmax = rmax*1.00001d0
  rtest = (maxn-20)*min(sqrt(a1*a1+a2*a2), sqrt(b1*b1+b2*b2))
  rtest = rtest*1.00001d0

  do n1 = -maxn, maxn
    n2 = -maxn
    r2 = 1e9
    r0 = 0.d0
    do while (n2<=maxn .and. (r2<=1.d-6) .and. (sqrt(r0)>rmax))
      n2 = n2 + 1
!do 10 n2=-MAXN,MAXN
      x0 = zoffx - (a1*n1+b1*n2)
      y0 = zoffy - (a2*n1+b2*n2)
      r2 = x0*x0 + y0*y0 + zzz
      r0 = (a1*n1+b1*n2)**2 + (a2*n1+b2*n2)**2
!IF (r2.LE.1.D-6) GO TO 10
!if (sqrt(r0).gt.rmax) goto 10
      ipar = ipar + n1 + n2
      ia = ia + 1
      call ymy(x0, y0, zz, r, ylm, l2potd)
      cr = 1.d0/r
      do l = 1, 2*lpot
        cc = cr**(l+1)
        do m = -l, l
          lm = l*(l+1) + m + 1
          ccc = cc*ylm(lm)
          suml(lm) = suml(lm) + ccc
        end do
      end do
!     test the convergence
      if (sqrt(r0)>rtest) then
        do l = 1, 2*lpot
          cc = cr**(l+1)
          do m = -l, l
            lm = l*(l+1) + m + 1
            ccc = cc*ylm(lm)
            sumt(lm) = sumt(lm) + ccc
          end do
        end do
      end if
    end do ! Do while
  end do ! all atoms in plane

  if (ipar/=0) write (1337, *) 'SUMLAYER  Asymetric sum '
  lm2pot = (2*lpot+1)**2
! Now test the sum
  if (test('electro ')) then
    do lm = 2, lm2pot ! test only l=3 terms!
      if (abs(suml(lm))>1.d-8) then ! test
        write (1337, 100) vec1(1), vec1(2), vec2(1), vec2(2), lm, suml(lm), &
          sumt(lm)
      end if
    end do
  end if
100 format ('SUMLAYER ', 4f8.4, ' LM=', i3, ' SUM=', 2d12.5)
110 format (i5, 5f12.6)
!         do lm=2,lm2pot
!            if (dabs(sum(lm)).gt.1.d-8) then
!               write(6,4005) i1,i2,lm,sum(lm)
!            end if
!         end do
!     -------------------------------------------------------------
120 format (1x, ' I1, I2 LM, DIRECT SUM=', 3i4, (d18.9,d18.9))

end subroutine
