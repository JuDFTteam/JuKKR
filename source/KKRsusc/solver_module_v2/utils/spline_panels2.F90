  module mod_spline_panels2

  interface spline_panels2
    module procedure spline_panels2_real, spline_panels2_complex
  end interface spline_panels2

  contains

  subroutine spline_panels2_complex(x,y,n0,n1,ypn0,ypn1,y2)
! copied from 'Numerical recipes in Fortran 77', page 109
! subroutine calculates second derivative of input data
! input: data points x(n),y(n)
!       first derivatives at border ypn0,ypn1
! output: second derivative y2
  use global, only: i4b,r8b,c8b
  
  implicit none
  
  integer(kind=i4b), intent(in)       :: n0,n1
  complex(kind=c8b), intent(in)       :: ypn0,ypn1,y(n0:n1)
  real(kind=r8b),    intent(in)       :: x(n0:n1)
  complex(kind=c8b), intent(out)      :: y2(n0:n1)
! #############################################################
  integer(kind=i4b), parameter        :: nmax=1000
  integer(kind=i4b)                   :: i,k
  complex(kind=c8b)                   :: p,qn,un,u(nmax)
  real(kind=r8b)                      :: sig
  
! if derivative at the left border is too large take second derivative as zero
  if(abs(ypn0) .gt. 0.99e30) then
    y2(n0)=0.d0
    u(n0)=0.d0
  else
    y2(n0)=-0.5d0
    u(n0)=(3.d0/(x(n0+1)-x(n0)))*((y(n0+1)-y(n0))/(x(n0+1)-x(n0))-ypn0)
  end if
  
  do i=n0+1,n1-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  
  if (abs(ypn1) .gt. 0.99e30) then
    qn=0.d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.d0/(x(n1)-x(n1-1)))*(ypn1-(y(n1)-y(n1-1))/(x(n1)-x(n1-1)))
  end if
  
  y2(n1)=(un-qn*u(n1-1))/(qn*y2(n1-1)+1.d0)
  
  do k=n1-1,n0,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do
  
  end subroutine spline_panels2_complex


  subroutine spline_panels2_real(x,y,n0,n1,ypn0,ypn1,y2)
! copied from 'Numerical recipes in Fortran 77', page 109
! subroutine calculates second derivative of input data
! input: data points x(n),y(n)
!       first derivatives at border ypn0,ypn1
! output: second derivative y2
  use global, only: i4b,r8b,c8b
  
  implicit none
  
  integer(kind=i4b), intent(in)       :: n0,n1
  real(kind=r8b),    intent(in)       :: ypn0,ypn1,y(n0:n1)
  real(kind=r8b),    intent(in)       :: x(n0:n1)
  real(kind=r8b),    intent(out)      :: y2(n0:n1)
! #############################################################
  integer(kind=i4b), parameter        :: nmax=1000
  integer(kind=i4b)                   :: i,k
  complex(kind=c8b)                   :: p,qn,un,u(nmax)
  real(kind=r8b)                      :: sig
  
! if derivative at the left border is too large take second derivative as zero
  if(abs(ypn0) .gt. 0.99e30) then
    y2(n0)=0.d0
    u(n0)=0.d0
  else
    y2(n0)=-0.5d0
    u(n0)=(3.d0/(x(n0+1)-x(n0)))*((y(n0+1)-y(n0))/(x(n0+1)-x(n0))-ypn0)
  end if
  
  do i=n0+1,n1-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  
  if (abs(ypn1) .gt. 0.99e30) then
    qn=0.d0
    un=0.d0
  else
    qn=0.5d0
    un=(3.d0/(x(n1)-x(n1-1)))*(ypn1-(y(n1)-y(n1-1))/(x(n1)-x(n1-1)))
  end if
  
  y2(n1)=(un-qn*u(n1-1))/(qn*y2(n1-1)+1.d0)
  
  do k=n1-1,n0,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do
  
  end subroutine spline_panels2_complex
