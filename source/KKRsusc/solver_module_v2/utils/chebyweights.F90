  subroutine chebyweights(a,b,n,np,z,w,w0,w1,w00,w10,w01,w11)
! Weights for integration on line segment [a,b]
! z, w: Chebyshev points and weights
! w0: weights to extrapolate f(b) from calculated f(x)
! w1: weights to extrapolate df/dx(b) from calculated f(x)
! w00: weights to integrate product f(x) g(x)
! w10: weights to integrate product df/dx g(x)
! w01: weights to integrate product f(x) dg/dx
! w11: weights to integrate product df/dx dg/dx
  use global

  implicit none

! Mesh start and end
  complex(kind=c8b), intent(in)  :: a, b
! Number of points, size of arrays
  integer(kind=i4b), intent(in)  :: n, np
! Points at which f(x) and g(x) have to be calculated
  complex(kind=c8b), intent(out) :: z(np)
! Weights for direct integration
  complex(kind=c8b), intent(out) :: w(np)
! Weights for extrapolation
  complex(kind=c8b), intent(out) :: w0(np), w1(np)
! Weights for integration of products
  complex(kind=c8b), intent(out) :: w00(np,np), w10(np,np), w01(np,np), w11(np,np)
! ----------------------------------------------------------------------
  real(kind=r8b),    parameter :: pi = 4.d0*atan(1.d0)
  integer(kind=i4b), parameter :: m = 2
  integer(kind=i4b) :: i, j, k, l
  real(kind=r8b)    :: tk, tl, fac, xc(m*n), wc(m*n)
  complex(kind=c8b) :: plus, minus


! Basic Chebyshev points and weights
  plus  = 0.5d0*(b+a)
  minus = 0.5d0*(b-a)
  do k=1,n
    tk = pi*(n-k+0.5d0)/n
    z(k) = plus + minus*cos(tk)
    w(k) = 2.d0
    do i=2,n-1,2
      w(k) = w(k) - 4.d0*cos(i*tk)/((i+1)*(i-1))
    end do
    w(k) = minus*w(k)/n
  end do
!  write(*,'(2es16.8)') sum(w)
! ----------------------------------------------------------------------
! Weights to extrapolate to f(b) and to df/dx(b)
  do k=1,n
    tk = pi*(n-k+0.5d0)/n
    w0(k) = 0.d0
    w1(k) = 0.d0
    do i=1,n-1
      w0(k) = w0(k) + cos(i*tk)
      w1(k) = w1(k) + cos(i*tk)*i*i
    end do
    w0(k) = (2.d0*w0(k) + 1.d0)/n
    w1(k) = 2.d0*w1(k)/(minus*n)
  end do
! ----------------------------------------------------------------------
! Weights to integrate f(z) g(z)
  do l=1,n
    tl = pi*(n-l+0.5d0)/n
    do k=1,n
      tk = pi*(n-k+0.5d0)/n
      w00(k,l) = 1.d0
      do i=2,n-1,2
        w00(k,l) = w00(k,l) - 2.d0*(cos(i*tk)+cos(i*tl))/((i+1)*(i-1))
      end do
      do j=1,n-1
        do i=1,n-1
          if (mod(i+j,2) == 1) cycle
          w00(k,l) = w00(k,l) - 2.d0*cos(i*tk)*cos(j*tl)*(1.d0/((i+j)**2-1) + 1.d0/((i-j)**2-1))
        end do
      end do
      w00(k,l) = minus*2.d0*w00(k,l)/(n*n)
    end do
  end do
! ----------------------------------------------------------------------
! Weights to integrate f'(z)g(z) and f(z)g'(z)
  do l=1,n
    tl = pi*(n-l+0.5d0)/n
    do k=1,n
      tk = pi*(n-k+0.5d0)/n
      w10(k,l) = 0.d0
      w01(k,l) = 0.d0
      do j=1,n-1
        do i=1,n-1
          if (mod(i+j,2) == 0) cycle
          w10(k,l) = w10(k,l) + 2.d0*cos(i*tk)*cos(j*tl)*i*(1.d0/(i+j) + 1.d0/(i-j))
          w01(k,l) = w01(k,l) + 2.d0*cos(i*tk)*cos(j*tl)*j*(1.d0/(i+j) + 1.d0/(j-i))
        end do
      end do
      do i=1,n-1,2
        w10(k,l) = w10(k,l) + 2.d0*cos(i*tk)
        w01(k,l) = w01(k,l) + 2.d0*cos(i*tl)
      end do
      w10(k,l) = 2.d0*w10(k,l)/(n*n)
      w01(k,l) = 2.d0*w01(k,l)/(n*n)
    end do
  end do
! ----------------------------------------------------------------------
! Weights to integrate f'(z) g'(z)
! mesh for numerical integration
  do k=1,m*n
    tk = pi*(m*n-k+0.5d0)/(m*n)
    xc(k) = 0.5d0*pi*(1.d0 + cos(tk))
    wc(k) = 2.d0
    do i=2,m*n-1,2
      wc(k) = wc(k) - 4.d0*cos(i*tk)/((i+1)*(i-1))
    end do
    wc(k) = 0.5d0*pi*wc(k)/(m*n)
  end do
! fill in the weights
  do l=1,n
    tl = pi*(n-l+0.5d0)/n
    do k=1,n
      tk = pi*(n-k+0.5d0)/n
      w11(k,l) = 0.d0
      do j=1,n-1
        do i=1,n-1
          fac = sum(wc*sin(i*xc)*sin(j*xc)/sin(xc))
          w11(k,l) = w11(k,l) + cos(i*tk)*cos(j*tl)*i*j*fac
        end do
      end do
      w11(k,l) = 4.d0*w11(k,l)/(minus*n*n)
    end do
  end do
! All done!
  end subroutine chebyweights
