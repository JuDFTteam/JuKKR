module mod_soutk

contains

  !-------------------------------------------------------------------------------
  !> Summary: Outward integration with kinks
  !> Author: B. Drittler
  !> Date: Oct 1989
  !> Category: KKRimp, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> this subroutine does an outwards integration of a function
  !> with kinks
  !>
  !>
  !>                          r
  !>                fint(r) = s f(r') dr'
  !>                          0
  !>
  !> at each kink the integration is restarted
  !> the starting value for this integration is determined by
  !> a 4 point lagrangian integration  , coefficients given by
  !> m. abramowitz and i.a. stegun, handbook of mathematical functions,
  !> nbs applied mathematics series 55 (1968)
  !>
  !> @warning the weights drdi have to be multiplied before calling this
  !> subroutine. @endwarning
  !-------------------------------------------------------------------------------
  subroutine soutk(f,fint,ipan,ircut)
    implicit none
    !     .. scalar arguments ..
    integer ipan
    !     .. array arguments ..
    real*8 f(*),fint(*)
    integer ircut(0:ipan)
    !     .. local scalars ..
    real*8 a1,a2
    integer i,ien,ip,ist


    a1 = 1.0d0/3.0d0
    a2 = 4.0d0/3.0d0

    !---> loop over kinks
    do ip = 1,ipan
      ien = ircut(ip)
      ist = ircut(ip-1) + 1

      if (ip.eq.1) then
      fint(ist) = 0.0d0
      !---> integrate fint(ist+1) with a 4 point lagrangian
      fint(ist+1) = (f(ist+3)-5.0d0*f(ist+2)+19.0d0*f(ist+1)+ 9.0d0*f(ist))/24.0d0

      else
        fint(ist) = fint(ist-1)
        !---> integrate fint(ist+1) with a 4 point lagrangian
        fint(ist+1) = fint(ist-1) + (f(ist+3)-5.0d0*f(ist+2)+ 19.0d0*f(ist+1)+9.0d0*f(ist))/24.0d0
      end if

      !---> calculate fint with an extended 3-point-simpson
      do i = ist + 2,ien
        fint(i) = ((fint(i-2)+f(i-2)*a1)+f(i-1)*a2) + f(i)*a1
      end do
    end do

  end subroutine

end module mod_soutk
