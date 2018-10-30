!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_hankel
  
  private
  public hankel

contains

  !-------------------------------------------------------------------------------
  !> Summary: Analysical expression of Hankel functions
  !> Author: 
  !> Category: KKRhost, special-functions, core-electrons
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine uses the explicit formulas for the Hankel
  !> functions. For higher l-values these formulas may lead to
  !> loss of significant figures. 
  !> @warning 
  !> This subroutine should be used
  !> only for core states.
  !> @endwarning
  !-------------------------------------------------------------------------------
  subroutine hankel(h, l, arg)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: ci
    implicit none
    ! .. Scalar Arguments ..
    complex (kind=dp) :: arg
    integer :: l
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: h(*)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: a1, a2, a3, a4
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: exp
    

    h(1) = -exp(arg*ci)/arg
    if (l/=1) then
      a1 = (1.e0_dp, 0.e0_dp) - arg*ci
      h(2) = h(1)*a1/arg
      if (l/=2) then
        a1 = 3.e0_dp*a1
        a2 = arg*arg
        h(3) = h(1)*(a1-a2)/a2
        if (l/=3) then
          a1 = 5.e0_dp*a1
          a3 = a2*arg*ci
          a4 = a2*arg
          a2 = 6.e0_dp*a2
          h(4) = h(1)*(a1-a2+a3)/a4
          if (l/=4) then
            a1 = 7.e0_dp*a1
            a2 = 7.5e0_dp*a2
            a3 = 10.e0_dp*a3
            a4 = a4*arg
            h(5) = h(1)*(a1-a2+a3+a4)/a4
            if (l/=5) then
              h(6) = (9.0e0_dp, 0.0e0_dp)*h(5)/arg - h(4)
              if (l/=6) then
                write (6, fmt=100) l
                stop 'HANKEL'

              end if

            end if

          end if

        end if

      end if

    end if

    return


100 format (2x, ' hankel :  l=', i2, ' is too large')
  end subroutine hankel

end module mod_hankel
