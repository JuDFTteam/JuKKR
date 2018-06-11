    Subroutine cinthff(ag, af, bg, bf, rmehf, nka, nkb, jtop, fx, r, drdi, &
      nrmax)
!   ********************************************************************
!   *                                                                  *
!   *  routine to calculate the radial hyperfine matrixelement         *
!   *  by extrapolating the lower integration boundary to r -> 0       *
!   *                                                                  *
!   ********************************************************************
      Use mod_datatypes
      Implicit None

!Dummy arguments
      Integer :: jtop, nka, nkb, nrmax
      Complex (Kind=dp) :: af(nrmax, nka), ag(nrmax, nka), bf(nrmax, nkb), &
        bg(nrmax, nkb), fx(jtop), rmehf(2, 2)
      Real (Kind=dp) :: drdi(nrmax), r(nrmax)

! Local variables
      Real (Kind=dp) :: ai, ar, bi, br, delta, hibar, r1bar, reldif, rlim1, &
        rlim2, s, s_i, sr, sx, sxi, sxr, sxx, yi, yr

      Complex (Kind=dp) :: hi, hif, z(nrmax)
      Integer :: i, imax, imin, ka, kb, n

      Do kb = 1, nkb
        Do ka = 1, nka

          Do i = 1, jtop
            fx(i) = (ag(i,ka)*bf(i,kb)+af(i,ka)*bg(i,kb))*drdi(i)
          End Do

          Call cint4pts(fx, jtop, z)

          rlim1 = 0.5E-5_dp
          rlim2 = 0.5E-4_dp
!            RLIM1 = 1D-4
!            RLIM2 = 5D-4
          If (r(1)>rlim1) Then
            rlim1 = r(1)
            rlim2 = r(20)
          End If

          imin = 0
          imax = 0
          Do i = 1, jtop
            If (r(i)<=rlim1) imin = i
            If (r(i)>=rlim2) Then
              imax = i
              Go To 100
            End If
          End Do
100       Continue
          If (imin==0) Stop '<CINTHFF> IMIN = 0'
          If (imin>jtop) Stop '<CINTHFF> IMIN > JTOP'
          If (imax==0) Stop '<CINTHFF> IMAX = 0'
          If (imax>jtop) Stop '<CINTHFF> IMAX > JTOP'

          n = 0
          s = 0.0E0_dp
          sx = 0.0E0_dp
          sxx = 0.0E0_dp
          sr = 0.0E0_dp
          s_i = 0.0E0_dp
          sxr = 0.0E0_dp
          sxi = 0.0E0_dp

          Do i = imin, imax
            yr = real(z(jtop)-z(i))
            yi = aimag(z(jtop)-z(i))
            n = n + 1
            s = s + 1.0E0_dp
            sx = sx + r(i)
            sxx = sxx + r(i)**2
            sr = sr + yr
            s_i = s_i + yi
            sxr = sxr + yr*r(i)
            sxi = sxi + yi*r(i)
          End Do

          delta = s*sxx - sx*sx
          ar = (sxx*sr-sx*sxr)/delta
          ai = (sxx*s_i-sx*sxi)/delta
          br = (s*sxr-sx*sr)/delta
          bi = (s*sxi-sx*s_i)/delta

          hibar = 0.0E0_dp
          r1bar = 0.0E0_dp
          Do i = imin, imax
            hibar = hibar + real(z(jtop)-z(i))
            r1bar = r1bar + r(i)
          End Do
          hibar = hibar/real(imax-imin+1, kind=dp)
          r1bar = r1bar/real(imax-imin+1, kind=dp)

          sxx = 0.0E0_dp
          sxr = 0.0E0_dp
          sxi = 0.0E0_dp
          reldif = 0.0E0_dp
          Do i = imin, imax
            hi = z(jtop) - z(i)
            hif = cmplx(ar, ai, kind=dp) + cmplx(br, bi, kind=dp)*r(i)
            If (abs(hif)/=0.0E0_dp) reldif = max(reldif, abs(1.0E0_dp-hi/hif))
          End Do

          rmehf(ka, kb) = cmplx(ar, ai, kind=dp)

        End Do
      End Do

    End Subroutine
