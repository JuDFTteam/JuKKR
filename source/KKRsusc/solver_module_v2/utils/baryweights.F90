  subroutine baryweights(n,d,x,w)
! Weights for barycentric interpolation
! Rescaled by the factorial of the interpolation order

  implicit none

! Number of points, order of interpolation
  integer(kind=i4b), intent(in)  :: n, d
! Abcissas
  complex(kind=c8b), intent(in)  :: x(n)
! Weights
  complex(kind=c8b), intent(out) :: w(n)
! -----------------------------------------
  integer(kind=i4b) :: i, j, k, imin, imax, iset(d+1)
  real(kind=r8b)    :: minus, fact, hmin, h
  complex(kind=c8b) :: prod

  if (d > n) stop 'baryweights: d > n'
! factorial of the interpolation order
  fact = 1.d0
  do i=2,d
    fact = fact*i
  end do
  hmin = 1.d99
! smallest step size
  do k=1,n
    do i=k+1,n
      h = abs(x(k) - x(i))
      if (h < hmin) hmin = h
    end do
  end do
  hmin = hmin**d
! loop over interpolation points
  do k=1,n
!   form weights from neighbours
    imin = max(k-d,1)
    imax = min(k,n-d)
!   sum over products
    w(k) = 0.d0
    write(iodb,'("k imin imax:",3i6)') k, imin, imax
    do i=imin,imax
      minus = (-1.d0)**(i-1)
!     product
      prod = 1.d0
      do j=i,i+d
        if (j == k) cycle  ! skip 1/0
        write(iodb,'("k i j:",3i6)') k, i, j
        prod = prod/(x(k) - x(j))
      end do
      w(k) = w(k) + minus*prod
    end do
    w(k) = w(k)*fact*hmin
!    w(k) = (-1.d0)**(k-1)
  end do
!  w(1) = 0.5d0
!  w(n) = 0.5d0*(-1.d0)**(n-1)
! All done!
  end subroutine baryweights

