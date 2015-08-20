      subroutine csinwd(f, fint, lmmsqd, irmind, irmd, ipan, ircut)
!-----------------------------------------------------------------------
!     this subroutine does an inwards integration of llmax
!     functions f with an extended 3-point-simpson :
!
!
!                               irmax
!                   fint(:,i) = { f(:,i') di'
!                                ir
!
!  the starting value for this integration at is - 1 is determined by
!    a 4 point lagrangian integration  , coefficients given by
!    m. abramowitz and i.a. stegun, handbook of mathematical functions,
!    nbs applied mathematics series 55 (1968)
!
!  attention in case of radial integration :
!       the weights drdi have to be multiplied before calling this
!       subroutine .
!
!                                     b. drittler mar. 1989
!
!    modified for functions with kinks - at each kink the integration
!      is restarted
!
!    attention : it is supposed that irmin + 3 is less than imt !
!
!
!                                     b. drittler july 1989
!    modified by m. ogura, june 2015
!-----------------------------------------------------------------------
      implicit none
      integer, intent(in) :: ipan, irmd, irmind, lmmsqd
      double complex, intent(in) :: f(lmmsqd,irmind:irmd)
      double complex, intent(out) :: fint(lmmsqd,irmind:irmd)
      integer, intent(in) :: ircut(0:ipan)
      
      double precision, parameter :: a1=5.d0/12.d0, a2=8.d0/12.d0, a3=-1.d0/12.d0
      integer :: ir, ie, ip, is
!
!---> loop over kinks
!
      do ip = ipan, 1, -1
        is = ircut(ip)
        ie = ircut(ip-1) + 1
        if (ip == 1) ie = irmind ! first panel
!
        if (ip == ipan) then
          fint(:,is) = 0.0d0 ! last panel
        else
          fint(:,is) = fint(:,is+1)
        endif ! last panel
!
!---> calculate fint with an extended 3-point-simpson
!
        do ir = is, ie+2, -2
          fint(:,ir-1) = fint(:,ir-0) + f(:,ir)*a1 + f(:,ir-1)*a2 + f(:,ir-2)*a3
          fint(:,ir-2) = fint(:,ir-1) + f(:,ir)*a3 + f(:,ir-1)*a2 + f(:,ir-2)*a1
        enddo ! ir
        if (mod(is-ie,2) == 1) then
          fint(:,ie  ) = fint(:,ie+1) + f(:,ie)*a1 + f(:,ie+1)*a2 + f(:,ie+2)*a3
        endif
      enddo ! ip
     
      endsubroutine ! csinwd
