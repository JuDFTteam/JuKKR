  module bessel_functions
! Quick and dirty implementation of spherical bessel functions
  use global, only: r8b, c8b

  implicit none

!  private
!  public :: bessj, dbessj, bessn, dbessn, bessh1, dbessh1, bessh2, dbessh2

  complex(kind=c8b), parameter :: iu = (0.d0,1.d0)


! Bessel functions of the 1st kind, and derivatives
! integer index, real(kind=r8b) or complex(kind=c8b) x
  interface bessj
    module procedure rbess, zbess
  end interface

  interface dbessj
    module procedure drbess, dzbess
  end interface

! Bessel functions of the 2nd kind, and derivatives
! integer index, real(kind=r8b) or complex(kind=c8b) x
  interface bessn
    module procedure rneum, zneum
  end interface

  interface dbessn
    module procedure drneum, dzneum
  end interface

! Bessel functions of the 3rd kind, and derivatives
! integer index, real(kind=r8b) or complex(kind=c8b) x
  interface bessh1
    module procedure rhank1, zhank1
  end interface

  interface dbessh1
    module procedure drhank1, dzhank1
  end interface

  interface bessh2
    module procedure rhank2, zhank2
  end interface

  interface dbessh2
    module procedure drhank2, dzhank2
  end interface


  contains


! ----------------------------------------------------------------------
! bessel function of the 1st kind, real and complex arguments
! includes the analytical derivatives
! ----------------------------------------------------------------------

  real(kind=r8b) function rbess(n,x)
! Bessel function
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x

  rbess = real(rhank1(n,x))
! All done
  end function rbess


  real(kind=r8b) function drbess(n,x)
! derivative of Bessel function
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x

  drbess = (n*rbess(n-1,x)-(n+1)*rbess(n+1,x))/(2*n+1)
! All done!
  end function drbess


  complex(kind=c8b) function zbess(n,x)
! Bessel function
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x
  zbess = (zhank1(n,x)+zhank2(n,x))/2
! All done!
  end function zbess

  
  complex(kind=c8b) function dzbess(n,x)
! derivative of Bessel function
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  dzbess = (n*zbess(n-1,x)-(n+1)*zbess(n+1,x))/(2*n+1)
! All done!
  end function dzbess


! ----------------------------------------------------------------------
! bessel function of the 2nd kind, real and complex arguments
! includes the analytical derivatives
! ----------------------------------------------------------------------

  real(kind=r8b) function rneum(n,x)
! Neumann function
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x

  rneum = aimag(rhank1(n,x))
! All done!
  end function rneum


  real(kind=r8b) function drneum(n,x)
! derivative of Neumann function
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x
  drneum = (n*rneum(n-1,x)-(n+1)*rneum(n+1,x))/(2*n+1)
! All done!
  end function drneum


  complex(kind=c8b) function zneum(n,x)
! Neumann function
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  zneum = (zhank1(n,x)-zhank2(n,x))/(2*iu)
! All done!
  end function zneum


  complex(kind=c8b) function dzneum(n,x)
! derivative of Neumann function
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  dzneum = (n*zneum(n-1,x)-(n+1)*zneum(n+1,x))/(2*n+1)
! All done!
  end function dzneum


! ----------------------------------------------------------------------
! bessel functions of the 3rd kind, real and complex arguments
! includes the analytical derivatives
! ----------------------------------------------------------------------

  complex(kind=c8b) function rhank1(n,x)
! Hankel function of the first kind
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x
  complex(kind=c8b) :: z

  z = cmplx(x,kind=c8b)
  rhank1 = zhank1(n,z)
! All done!
  end function rhank1


  complex(kind=c8b) function rhank2(n,x)
! Hankel function of the 2nd kind
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x
  complex(kind=c8b) :: z

  z = cmplx(x,kind=c8b)
  rhank2 = zhank2(n,z)
! All done!
  end function rhank2


  complex(kind=c8b) function drhank1(n,x)
! derivative of Hankel function of the 1st kind
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x

  drhank1 = (n*rhank1(n-1,x)-(n+1)*rhank1(n+1,x))/(2*n+1)
! All done!
  end function drhank1


  complex(kind=c8b) function drhank2(n,x)
! derivative of Hankel function of the 2nd kind
! real argument

  implicit none

  integer,        intent(in) :: n
  real(kind=r8b), intent(in) :: x

  drhank2 = (n*rhank2(n-1,x)-(n+1)*rhank2(n+1,x))/(2*n+1)
! All done!
  end function drhank2


  recursive function zhank1(n,x) result(hank)
! Hankel function of the 1st kind
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x
  complex(kind=c8b) ::  hank

  if (n < 0 .or. abs(x) < 1.e-12) then
    hank = 0.0_r8b
    write(*,'("hank n=",i12)') n
  else if (abs(x) < 0.01*(2*n+1)) then
    hank = (1-x**2/(2*(2*n+3))) * x**n / dfac(2*n+1)              &
 &          - iu*(1-x**2/(2*(1-2*n)))*dfac(2*n-1)/x**(n+1)
    write(*,'("hank n=",i12)') n
  else if (n == 0) then
    hank = exp(x*iu)/(x*iu)
    write(*,'("hank n=",i12)') n
  else if (n == 1) then
    hank = -(1/x+iu/x**2)*exp(x*iu)
    write(*,'("hank n=",i12)') n
  else
    hank = (2*n-1)*zhank1(n-1,x)/x - zhank1(n-2,x)
    write(*,'("hank n=",i12)') n
  end if
! All done!
  end function zhank1


  recursive function zhank2(n,x) result(hank)
! Hankel function of the 2nd kind
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x
  complex(kind=c8b) ::  hank

  if (n < 0 .or. abs(x) < 1.e-12) then
    hank = 0.0_r8b
  else if (abs(x) < 0.01*(2*n+1)) then
    hank = (1-x**2/(2*(2*n+3))) * x**n / dfac(2*n+1)              &
 &          + iu*(1-x**2/(2*(1-2*n)))*dfac(2*n-1)/x**(n+1)
  else if (n == 0) then
    hank = exp(-x*iu)/(-x*iu)
  else if (n == 1) then
    hank = (-1/x+iu/x**2)*exp(-x*iu)
  else
    hank = (2*n-1)*zhank2(n-1,x)/x - zhank2(n-2,x)
  end if
! All done!
  end function zhank2


  complex(kind=c8b) function dzhank1(n,x)
! Derivative of Hankel function of the 1st kind
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  dzhank1 = (n*zhank1(n-1,x)-(n+1)*zhank1(n+1,x))/(2*n+1)
! All done!
  end function dzhank1


  complex(kind=c8b) function dzhank2(n,x)
! Derivative of Hankel function of the 2nd kind
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  dzhank2 = (n*zhank2(n-1,x)-(n+1)*zhank2(n+1,x))/(2*n+1)
! All done!
  end function dzhank2


  real(kind=r8b) function dfac(i)
! Double factorial

  integer, intent(in) :: i
! ----------------------------------------------------------------------
  integer :: k

  dfac = 1.d0
  k = i
  do
    dfac = dfac*k
    k = k - 2
    if (k == 1 .or. k == 0) exit
  end do
! All done!
  end function dfac


! All done!
  end module bessel_functions
