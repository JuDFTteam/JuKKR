subroutine calccgc(ltab, kaptab, nmuetab, cgc, nkmax, nmuemax, nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
!   *                                                                  *
!   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
!   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
!   *   IS= 1/2  SPIN DOWN/UP                                          *
!   *                                                                  *
!   ********************************************************************

  implicit none

! Dummy arguments
  integer :: nkmax, nkmpmax, nmuemax
  real *8 :: cgc(nkmpmax, 2)
  integer :: kaptab(nmuemax), ltab(nmuemax), nmuetab(nmuemax)

! Local variables
  integer :: ikm, k, kappa, m
  real *8 :: j, l, mue, twolp1

  ikm = 0
  do k = 1, (nkmax+1)
    l = ltab(k)
    kappa = kaptab(k)
    j = abs(kappa) - 0.5d0
    mue = -j - 1.0d0
    twolp1 = 2.0d0*l + 1.0d0

    if (kappa<0) then

!     J = L + 1/2
      do m = 1, nmuetab(k)

        mue = mue + 1.0d0
        ikm = ikm + 1
        cgc(ikm, 1) = dsqrt((l-mue+0.5d0)/twolp1)
        cgc(ikm, 2) = dsqrt((l+mue+0.5d0)/twolp1)
      end do
    else
!     J = L - 1/2
      do m = 1, nmuetab(k)

        mue = mue + 1.0d0
        ikm = ikm + 1
        cgc(ikm, 1) = dsqrt((l+mue+0.5d0)/twolp1)
        cgc(ikm, 2) = -dsqrt((l-mue+0.5d0)/twolp1)

      end do
    end if


  end do

end subroutine
