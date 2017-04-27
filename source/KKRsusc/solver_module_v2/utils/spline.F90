  subroutine spline(x,y,n,yp1,ypn,y2)
! copied from 'Numerical recipes in Fortran 77', page 109
! subroutine calculates second derivative of input data
! input: data points x(n),y(n)
!       first derivatives at border yp1,ypn
! output: second derivative y2
  use global, only: i4b,r8b,c8b
  
  implicit none
  
  integer(kind=i4b), intent(in)       :: n
  complex(kind=c8b), intent(in)       :: yp1,ypn,y(n)
  real(kind=r8b),    intent(in)       :: x(n)
  complex(kind=c8b), intent(out)      :: y2(n)
! #############################################################
  integer(kind=i4b), parameter        :: nmax=500
  integer(kind=i4b)                   :: i,k
  complex(kind=c8b)                   :: p,qn,un,u(nmax)
  real(kind=r8b)                      :: sig
  
! if derivative at the left border is too large take second derivative as zero
  if(abs(yp1) .gt. 0.99e30) then
    y2(1)=0.d0
    u(1)=0.d0
  else
    y2(1)=-0.5d0
    u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  end if
  
  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  
  if (abs(ypn) .gt. 0.99e30) then
    qn=0.d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  end if
  
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
  
  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do
  
  end subroutine spline
