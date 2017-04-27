  subroutine splint(xa,ya,y2a,n,x,y)
! calculates cubic spline interpolation
! input: input data xa(n),ya(n) and second derivative y2a(n) from spline routine
!       point at which function should be interpolated
! output: interpolated y
  use global, only: i4b,r8b,c8b
  
  implicit none
  
  integer(kind=i4b), intent(in)    :: n
  real(kind=r8b),    intent(in)    :: x,xa(n)
  complex(kind=c8b), intent(in)    :: y2a(n),ya(n)
  complex(kind=c8b), intent(out)   :: y
! ###############################################################
  integer(kind=i4b)   :: k,khi,klo
  real(kind=r8b)      :: a,b,h
  
  klo=1
  khi=n
  
  1  if (khi-klo .gt. 1) then
       k=(khi+klo)/2
       if (xa(k) .gt. x) then
         khi=k
       else
         klo=k
       end if
     goto 1
     end if
  
  h=xa(khi)-xa(klo)
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
  end subroutine splint
