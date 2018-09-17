module mod_cinthff
  use :: mod_datatypes, only: dp
  private :: dp

contains

  subroutine cinthff(ag, af, bg, bf, rmehf, nka, nkb, jtop, fx, r, drdi, nrmax)
    ! ********************************************************************
    ! *                                                                  *
    ! *  routine to calculate the radial hyperfine matrixelement         *
    ! *  by extrapolating the lower integration boundary to r -> 0       *
    ! *                                                                  *
    ! ********************************************************************
    use :: mod_datatypes
    use :: mod_cint4pts
    implicit none

    real (kind=dp), parameter :: eps = 1.0e-12_dp

    ! Dummy arguments
    integer :: jtop, nka, nkb, nrmax
    complex (kind=dp) :: af(nrmax, nka), ag(nrmax, nka), bf(nrmax, nkb), bg(nrmax, nkb), fx(jtop), rmehf(2, 2)
    real (kind=dp) :: drdi(nrmax), r(nrmax)

    ! Local variables
    real (kind=dp) :: ai, ar, bi, br, delta, hibar, r1bar, reldif, rlim1, rlim2, s, s_i, sr, sx, sxi, sxr, sxx, yi, yr

    complex (kind=dp) :: hi, hif, z(nrmax)
    integer :: i, imax, imin, ka, kb, n

    do kb = 1, nkb
      do ka = 1, nka

        do i = 1, jtop
          fx(i) = (ag(i,ka)*bf(i,kb)+af(i,ka)*bg(i,kb))*drdi(i)
        end do

        call cint4pts(fx, jtop, z)

        rlim1 = 0.5e-5_dp
        rlim2 = 0.5e-4_dp
        ! RLIM1 = 1D-4
        ! RLIM2 = 5D-4
        if (r(1)>rlim1) then
          rlim1 = r(1)
          rlim2 = r(20)
        end if

        imin = 0
        imax = 0
        do i = 1, jtop
          if (r(i)<=rlim1) imin = i
          if (r(i)>=rlim2) then
            imax = i
            go to 100
          end if
        end do
100     continue
        if (imin==0) stop '<CINTHFF> IMIN = 0'
        if (imin>jtop) stop '<CINTHFF> IMIN > JTOP'
        if (imax==0) stop '<CINTHFF> IMAX = 0'
        if (imax>jtop) stop '<CINTHFF> IMAX > JTOP'

        n = 0
        s = 0.0e0_dp
        sx = 0.0e0_dp
        sxx = 0.0e0_dp
        sr = 0.0e0_dp
        s_i = 0.0e0_dp
        sxr = 0.0e0_dp
        sxi = 0.0e0_dp

        do i = imin, imax
          yr = real(z(jtop)-z(i))
          yi = aimag(z(jtop)-z(i))
          n = n + 1
          s = s + 1.0e0_dp
          sx = sx + r(i)
          sxx = sxx + r(i)**2
          sr = sr + yr
          s_i = s_i + yi
          sxr = sxr + yr*r(i)
          sxi = sxi + yi*r(i)
        end do

        delta = s*sxx - sx*sx
        ar = (sxx*sr-sx*sxr)/delta
        ai = (sxx*s_i-sx*sxi)/delta
        br = (s*sxr-sx*sr)/delta
        bi = (s*sxi-sx*s_i)/delta

        hibar = 0.0e0_dp
        r1bar = 0.0e0_dp
        do i = imin, imax
          hibar = hibar + real(z(jtop)-z(i))
          r1bar = r1bar + r(i)
        end do
        hibar = hibar/real(imax-imin+1, kind=dp)
        r1bar = r1bar/real(imax-imin+1, kind=dp)

        sxx = 0.0e0_dp
        sxr = 0.0e0_dp
        sxi = 0.0e0_dp
        reldif = 0.0e0_dp
        do i = imin, imax
          hi = z(jtop) - z(i)
          hif = cmplx(ar, ai, kind=dp) + cmplx(br, bi, kind=dp)*r(i)
          if (abs(hif)>eps) reldif = max(reldif, abs(1.0e0_dp-hi/hif))
        end do

        rmehf(ka, kb) = cmplx(ar, ai, kind=dp)

      end do
    end do

  end subroutine cinthff

end module mod_cinthff
