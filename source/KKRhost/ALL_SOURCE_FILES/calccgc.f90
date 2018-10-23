module mod_calccgc

  private
  public :: calccgc

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculate Clebsh-Gordon coefficients in kappa-mue representation
  !> Author: 
  !> Category: KKRhost, special-functions, dirac
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)  in kappa-mue representation
  !>                                              
  !> IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE
  !> IKM  = L*2*(J+1/2) + J + MUE + 1             
  !> IS= 1/2  SPIN DOWN/UP
  !-------------------------------------------------------------------------------
  subroutine calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)
    use :: mod_datatypes, only: dp

    implicit none

    ! Dummy arguments
    integer :: nkmax, nkmpmax, nmuemax
    real (kind=dp) :: cgc(nkmpmax, 2)
    integer :: kaptab(nmuemax), ltab(nmuemax), nmuetab(nmuemax)

    ! Local variables
    integer :: ikm, k, kappa, m
    real (kind=dp) :: j, l, mue, twolp1

    ikm = 0
    do k = 1, (nkmax+1)
      l = ltab(k)
      kappa = kaptab(k)
      j = abs(kappa) - 0.5e0_dp
      mue = -j - 1.0e0_dp
      twolp1 = 2.0e0_dp*l + 1.0e0_dp

      if (kappa<0) then

        ! J = L + 1/2
        do m = 1, nmuetab(k)

          mue = mue + 1.0e0_dp
          ikm = ikm + 1
          cgc(ikm, 1) = sqrt((l-mue+0.5e0_dp)/twolp1)
          cgc(ikm, 2) = sqrt((l+mue+0.5e0_dp)/twolp1)
        end do
      else
        ! J = L - 1/2
        do m = 1, nmuetab(k)

          mue = mue + 1.0e0_dp
          ikm = ikm + 1
          cgc(ikm, 1) = sqrt((l+mue+0.5e0_dp)/twolp1)
          cgc(ikm, 2) = -sqrt((l-mue+0.5e0_dp)/twolp1)

        end do
      end if

    end do

  end subroutine calccgc

end module mod_calccgc
