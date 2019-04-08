! 14.10.95 ***************************************************************
      subroutine sinwk(f, fint, ipan, ircut)
! ************************************************************************
!    this subroutine does an inwards integration of a function
!    with kinks
!
!
!                             rc
!                   fint(r) = s f(r') dr'
!                             r
!
!    at each kink the integration is restarted
!    the starting value for this integration is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)
!
!    the weights drdi have to be multiplied before calling this
!    subroutine .
!
!                                     b. drittler oct. 1989
!-----------------------------------------------------------------------
      implicit none
!     ..
!     .. scalar arguments ..
      integer, intent(in) :: ipan
!     ..
!     .. array arguments ..
      double precision, intent(in) :: f(*)
      double precision, intent(out) :: fint(*)
      integer, intent(in) :: ircut(0:ipan)
!     ..
!     .. local scalars ..
      double precision, parameter :: a1=1.d0/3.d0, a2=4*a1
      integer :: i, ien, ip, ist
!
!---> loop over kinks
!
      do ip = ipan, 1, -1
        ist = ircut(ip)
        ien = ircut(ip-1) + 1

        if (ip == ipan) then
          fint(ist) = 0.d0
!---> integrate fint(ist-1) with a 4 point lagrangian
          fint(ist-1) =               (f(ist-3) - 5.d0*f(ist-2) + 19.d0*f(ist-1) + 9.d0*f(ist))/24.d0
        else
          fint(ist) = fint(ist+1)
!---> integrate fint(ist-1) with a 4 point lagrangian
          fint(ist-1) = fint(ist+1) + (f(ist-3) - 5.d0*f(ist-2) + 19.d0*f(ist-1) + 9.d0*f(ist))/24.d0
        endif
!
!---> calculate fint with an extended 3-point-simpson
!
        do i = ist - 2, ien, -1
          fint(i) = ((fint(i+2) + a1*f(i+2)) + a2*f(i+1)) + a1*f(i)
        enddo ! i
      enddo ! ip

      endsubroutine ! sinwk
