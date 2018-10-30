  module bessel_new
! Quick and dirty implementation of spherical bessel functions
  use global, only: r8b, c8b

  implicit none

  private
  public :: bessj, dbessj, bessn, dbessn, bessh1, dbessh1, bessh2, dbessh2

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


  complex(kind=c8b) function zhank1(n,x)
! Hankel function of the 1st kind
! check branch cut for complex x
! H_n^{(1)}(x) = (-i)^{n+1} \frac{e^{i x}}{x} P_n(x)
! P_n(x) = \sum_{m=0}^n \frac{i^m}{m!(2x)^m} \frac{(n+m)!}{(n-m)!}
! This polynomial in 1/x that can be written in nested form:
! P_2(x) = ( 1 + \frac{3i}{x} ( 1 + \frac{i}{x} ) )
! and the general form is implemented here

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x
! ----------------------------------------------------------------------
  complex(kind=c8b) :: coeff, ratio, phase
  integer           :: m

  ratio = 1.d0; phase = -iu
  do m=n,1,-1
    phase = -phase*iu
    coeff = iu*(n-m+1.d0)*((n+m)/(2.d0*m))
    ratio = 1.d0 + coeff*ratio/x
  end do
  zhank1 = phase*ratio*(exp(iu*x)/x)
! All done!
  end function zhank1


  complex(c8b) function zhank2(n,x)
! Hankel function of the 2nd kind
! check branch cut for complex x

  implicit none

  integer,           intent(in) :: n
  complex(kind=c8b), intent(in) :: x

  zhank2 = conjg(zhank1(n,conjg(x)))
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


! All done!
  end module bessel_new
