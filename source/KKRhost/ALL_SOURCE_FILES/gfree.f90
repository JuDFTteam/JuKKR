module mod_gfree13

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of free electron Green's function
  !> Author: 
  !> Category: KKRhost, special-functions, reference-system, structural-greensfunction
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> CALCULATION OF FREE ELECTRON GREEN'S FUNCTION :  G(M)LL'(E0)
  !>
  !> Also:
  !> ---- derivative of free electron green function matrix elements
  !> returned in array DGMLL=dGMLL / dE
  !>
  !> the analytical formula for the derivative of spherical Hankel
  !> functions is used:
  !>
  !> d                     l+1
  !> --  h (x) = h   (x) - --- h (x)
  !> dx   l       l-1       x   l
  !>
  !> which for x = sqrt(E0)*r leads to
  !>
  !> d                       r           rl
  !> --- ( sqrt(E0) h (x) ) = - h   (x) - -- h (x) )
  !> dE0             l        2  l-1      2x  l
  !>
  !> Ported from KKRnano by Phivos Mavropoulos 10.10.2013
  !-------------------------------------------------------------------------------
  subroutine gfree13(rdiff, e0, gmll, dgmll, cleb, icleb, loflm, iend)
    use :: mod_datatypes, only: dp
    use :: global_variables, only: ncleb, lmgf0d, lm2d, lmaxd
    use :: mod_beshan, only: beshan
    use :: mod_ymy, only: ymy
    use :: mod_constants, only: czero, ci, pi
    implicit none

    real (kind=dp), parameter :: fpi = 4.e0_dp*pi, rfpi = sqrt(fpi)

    complex (kind=dp) :: e0
    integer :: iend

    complex (kind=dp) :: gmll(lmgf0d, lmgf0d)
    complex (kind=dp) :: dgmll(lmgf0d, lmgf0d)
    real (kind=dp) :: cleb(ncleb), rdiff(3)
    integer :: icleb(ncleb, 4), loflm(lm2d)

    real (kind=dp) :: rabs, x, y, z
    integer :: ifac, j, lm1, lm2, lm3

    complex (kind=dp) :: bl(lmaxd*2+1), hl(lmaxd*2+1), hyl((lmaxd*2+1)**2), nl(lmaxd*2+1)
    complex (kind=dp) :: dhl(lmaxd*2+1), dhyl((lmaxd*2+1)**2) 
    real (kind=dp) :: yl((lmaxd*2+1)**2)
    integer :: lf((lmaxd*2+1)**2)


    do lm1 = 1, (lmaxd*2+1)**2
      lf(lm1) = loflm(lm1) + 1
    end do
    x = rdiff(1)
    y = rdiff(2)
    z = rdiff(3)
    call ymy(x, y, z, rabs, yl, lmaxd*2)
    call beshan(hl, bl, nl, sqrt(e0)*rabs, lmaxd*2)

    ! Derivatives of Hankel functions ! LLY
    dhl(1) = 0.5e0_dp*ci*rabs*hl(1)
    do lm1 = 2, lmaxd*2 + 1
      dhl(lm1) = 0.5e0_dp*(rabs*hl(lm1-1)-(lm1-1)*hl(lm1)/sqrt(e0))
    end do
    do lm1 = 1, (lmaxd*2+1)**2
      hyl(lm1) = -fpi*ci*sqrt(e0)*yl(lm1)*hl(lf(lm1))
      dhyl(lm1) = -fpi*ci*yl(lm1)*dhl(lf(lm1)) ! LLY
    end do

    do lm1 = 1, lmgf0d
      gmll(lm1, lm1) = hyl(1)/rfpi
      dgmll(lm1, lm1) = dhyl(1)/rfpi ! LLY
      do lm2 = 1, lm1 - 1
        gmll(lm1, lm2) = czero
        dgmll(lm1, lm2) = czero ! LLY
      end do
    end do

    do j = 1, iend
      lm1 = icleb(j, 1)
      lm2 = icleb(j, 2)
      lm3 = icleb(j, 3)
      gmll(lm1, lm2) = gmll(lm1, lm2) + cleb(j)*hyl(lm3)
      dgmll(lm1, lm2) = dgmll(lm1, lm2) + cleb(j)*dhyl(lm3)  ! LLY
    end do
    do lm1 = 1, lmgf0d
      do lm2 = 1, lm1 - 1
        ifac = (-1)**(loflm(lm1)+loflm(lm2))
        gmll(lm2, lm1) = ifac*gmll(lm1, lm2)
        dgmll(lm2, lm1) = ifac*dgmll(lm1, lm2) ! LLY
      end do
    end do
    return

  end subroutine gfree13

end module mod_gfree13
