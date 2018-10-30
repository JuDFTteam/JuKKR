!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_gamfc
  
  private
  public :: gamfc

contains

  !-------------------------------------------------------------------------------
  !> Summary: Convergence function of Ewald sum
  !> Author: 
  !> Category: KKRhost, geometry, special-functions
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculation of convergence function
  !>
  !> glh = i(alpha,l)/r**(l+1)*sqrt(pi)
  !>
  !> with
  !> alpha = r times the splitting paramter lamda
  !> and
  !> i(x,l) = erfc(x) + exp(-x*x)/sqrt(pi) *
  !>
  !> sum ( 2**i * x**(2i-1) / (2i-1)!! )
  !> 1..i..l
  !-------------------------------------------------------------------------------
  subroutine gamfc(alpha, glh, lmax, r)

    use :: mod_datatypes, only: dp
    use :: mod_erfcex, only: erfcex
    implicit none

    real (kind=dp) :: alpha, r
    integer :: lmax
    real (kind=dp) :: glh(0:lmax)
    real (kind=dp) :: arg, facl, fex
    integer :: l

    arg = alpha*alpha
    glh(0) = erfcex(alpha)
    facl = 2.0e0_dp*alpha

    ! ---> recursion

    do l = 1, lmax
      glh(l) = glh(l-1) + facl
      facl = facl*arg/(real(l,kind=dp)+0.5e0_dp)
    end do

    ! changed 21/10/99
    ! if arg is to big then cannot calculate 1/exp(arg) !!

    fex = exp(-arg)

    do l = 0, lmax
      fex = fex/r
      glh(l) = glh(l)*fex
    end do

  end subroutine gamfc

end module mod_gamfc
