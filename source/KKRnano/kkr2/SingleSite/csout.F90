      subroutine csout(f,fint,lmmsqd,irmind,irmd,ipan,ircut)
!-----------------------------------------------------------------------
!     this subroutine does an outwards integration of llmax
!     functions f with an extended 3-point-simpson :
!
!
!                                ir
!                   fint(ll,i) = { f(ll,i') di'
!                               irmin
!
!  the starting value for this integration at irmin+1 is determined by
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
      integer, intent(in) :: ipan, irmd, irmind, lmmsqd
      double complex, intent(in) ::  f(lmmsqd,irmind:irmd)
      double complex, intent(inout) :: fint(lmmsqd,irmind:irmd)
      integer, intent(in) :: ircut(0:ipan)

      double precision, parameter :: a1=5.d0/12.d0, a2=8.d0/12.d0, a3=-1.d0/12.d0
      integer :: i, ien, ip, ist
!
!---> loop over kinks
!
      do ip = 1, ipan
        ien = ircut(ip)
        ist = ircut(ip-1) + 1

        if (ip == 1) then
          ist = irmind
          fint(:,ist) = 0.0d0
        else
          fint(:,ist) = fint(:,ist-1)
        endif
!
!---> calculate fint with an extended 3-point-simpson
!
        do i = ist, ien-2, 2
          fint(:,i+1) = fint(:,i+0) + f(:,i)*a1 + f(:,i+1)*a2 + f(:,i+2)*a3
          fint(:,i+2) = fint(:,i+1) + f(:,i)*a3 + f(:,i+1)*a2 + f(:,i+2)*a1
        enddo ! i
          
        if (mod(ien-ist, 2) == 1) then
          fint(:,ien) = fint(:,ien-1) + f(:,ien)*a1 + f(:,ien-1)*a2 + f(:,ien-2)*a3
        endif
      enddo ! ip

      endsubroutine ! csout
