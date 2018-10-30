  subroutine zratint(n,xa,ya,x,y,dy)
! Rational function interpolation from NR
! P_n(x)/Q_n(x), degree (n/2,n/2+remainder)
  use global, only: i4b, r8b, c8b

  implicit none

  integer(kind=i4b), intent(in)  :: n
  complex(kind=c8b), intent(in)  :: xa(n), ya(n), x
  complex(kind=c8b), intent(out) :: y, dy
! -----------------------------------------------------------------
  real(kind=r8b), parameter :: small = 1.d-16, tol = 1.d-8
  integer(kind=i4b) :: i, m, ns
  complex(kind=c8b) :: dd, t, w, c(n), d(n), hp
  real(kind=r8b)    :: h, hh

! Find point closest to desired one
  ns = 1
  hh = abs(x - xa(1))
  do i=1,n
    h = abs(x - xa(i))
! Return value at mesh point
    if (h < tol) then
      y = ya(i)
      dy = 0.d0
      return
    else if (h < hh) then
      ns = i
      hh = h
    end if
  end do
! Initialize
  c = ya; d = ya + small
! Zeroth approximation
  y = ya(ns); dy = 0.d0
!  write(*,'(i4,2es16.8)') ns, y
!  if (abs(y) < small) return
! Remaining corrections
  ns = ns - 1
  do m=1,n-1
    do i=1,n-m
      w  = c(i+1) - d(i)
      hp = xa(i+m) - x
      t  = (xa(i) - x)*d(i)/hp
      dd = w/(t - c(i+1))
      d(i) = c(i+1)*dd
      c(i) = t*dd
!      write(*,'(2i4,6es16.8)') n, m, c(i+1), d(i), (xa(i+m) - x)/(xa(i) - x)
    end do
    if (2*ns < n-m) then
      dy = c(ns+1)
    else
      dy = d(ns)
      ns = ns - 1
    end if
    y = y + dy
  end do
! All done!
  end subroutine zratint
