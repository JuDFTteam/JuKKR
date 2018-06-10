subroutine csimpk(cf, cfint, ipan, ircut, drdi)
!-----------------------------------------------------------------------
!     this subroutine does an integration up to rcut of an
!     complex function cf with an extended 3-point-simpson :

!                             rcut
!                      cfint = { cf(r') dr'
!                              0

!     modified for functions with kinks - at each kink the
!     integration is restarted .

!     attention : input cf is destroyed !

!-----------------------------------------------------------------------
!.. Scalar Arguments ..
  double complex :: cfint
  integer :: ipan
!..
!.. Array Arguments ..
  double complex :: cf(*)
  double precision :: drdi(*)
  integer :: ircut(0:ipan)
!..
!.. Local Scalars ..
  double precision :: a1, a2
  integer :: i, ien, ip, ist, n
!..
!.. External Functions ..
  double complex :: csum
  external :: csum
!..
!.. Intrinsic Functions ..
  intrinsic :: mod
!     ..
  a1 = 4.0d0/3.0d0
  a2 = 2.0d0/3.0d0
  cfint = 0.0d0

  do ip = 1, ipan

!---> loop over kinks

    ist = ircut(ip-1) + 1
    ien = ircut(ip)

    do i = ist, ien
      cf(i) = cf(i)*drdi(i)
    end do

    if (mod(ien-ist,2)==0) then
      cfint = cfint + (cf(ist)-cf(ien))/3.0d0
      ist = ist + 1
      n = (ien-ist+1)/2

    else
!---> four point lagrange integration for the first step
      cfint = cfint + (9.0d0*cf(ist)+19.0d0*cf(ist+1)-5.0d0*cf(ist+2)+cf(ist+3 &
        ))/24.0d0 + (cf(ist+1)-cf(ien))/3.0d0
      ist = ist + 2
      n = (ien-ist+1)/2
    end if

!---> calculate with an extended 3-point-simpson

    cfint = cfint + a1*csum(n, cf(ist), 2) + a2*csum(n, cf(ist+1), 2)
  end do

end subroutine
