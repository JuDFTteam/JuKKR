    Subroutine cradwf(eryd, ek, nsra, alpha, ipan, ircut, cvlight, rs, sl, pz, &
      fz, qz, sz, tmat, vm2z, drdi, rmesh, zat, lirrsol, idoldau, lopt, &
      wldauav, cutoff)
      Use mod_datatypes, Only: dp
!-----------------------------------------------------------------------
!  subroutine for radial wave functions of spherical potentials

!             the generalized phase shifts are calculated by
!             a wronski relation :

!                 alpha(z,l) =-sqrt(z)*wronski{hl(r;z),rl(r;z)}; r->0

!             where hl is the free hankel function and rl the regular
!             solution . Using the analytical behaviour of rl at the
!             origin (rl = alphal * r**(l+1)  ; r->0),
!             the generalized phase shifts can be calculated
!             directly with the renormalization alphal .
!                                           b.drittler nov.1987

!   LDA+U added, March 2003 - Dec 2004, Munich/Juelich
!-----------------------------------------------------------------------
      Implicit None
!     .. Parameters ..
      Include 'inc.p'
      Integer :: lmaxp1
      Parameter (lmaxp1=lmaxd+1)
      Complex (Kind=dp) :: ci, czero
      Parameter (ci=(0.E0_dp,1.E0_dp), czero=(0.0E0_dp,0.0E0_dp))
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: eryd, ek
      Real (Kind=dp) :: cvlight, zat
      Real (Kind=dp) :: wldauav
      Integer :: ipan, nsra, idoldau, lopt
      Logical :: lirrsol
!..
!.. Local Arrays ..
      Complex (Kind=dp) :: alpha(0:lmaxd), fz(irmd, 0:lmaxd), &
        pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd), tmat(0:lmaxd)
      Real (Kind=dp) :: drdi(irmd), rmesh(irmd), rs(irmd, 0:lmaxd), &
        sl(0:lmaxd), vm2z(irmd), cutoff(irmd)
      Integer :: ircut(0:ipand)
!..
!.. External Subroutines ..
      Complex (Kind=dp) :: alphal, arg, bl, eklfac, hl, pn, qf, slope, tlsqez, &
        value
      Real (Kind=dp) :: rirc, rirc1, rsirc, s1
      Integer :: i, ir, irc1, l1
!..
!.. Intrinsic Functions ..
      Complex (Kind=dp) :: bessjw(0:lmaxp1), bessyw(0:lmaxp1), &
        dlogdp(0:lmaxd), hamf(irmd, 0:lmaxd), hankws(0:lmaxp1), mass(irmd)
      Real (Kind=dp) :: dror(irmd)


      External :: beshan, irwsol, regsol
!---> calculate regular wavefunctions

      Intrinsic :: real

      irc1 = ircut(ipan)
      Do ir = 2, irc1
        dror(ir) = drdi(ir)/rmesh(ir)
      End Do
      rirc = rmesh(irc1)
      rirc1 = 1E0_dp/rirc
      arg = rirc*ek
      Call beshan(hankws, bessjw, bessyw, arg, lmaxp1)
! ======================================================================

!---> determine t - matrix
      Call regsol(cvlight, eryd, nsra, dlogdp, fz, hamf, mass, pz, dror, &
        rmesh, sl, vm2z, zat, ipan, ircut, idoldau, lopt, wldauav, cutoff, &
        irmd, ipand, lmaxd)

      eklfac = ek

      Do l1 = 0, lmaxd
!---> determine the renormalization


        qf = real(l1, kind=dp)*rirc1
        hl = hankws(l1)*dlogdp(l1)
        bl = bessjw(l1)*dlogdp(l1)
        hl = qf*hankws(l1) - ek*hankws(l1+1) - hl
        bl = bl - qf*bessjw(l1) + ek*bessjw(l1+1)
        hl = hl*ek
        tmat(l1) = ci*bl/hl
!---> determine the alpha matrix


        tlsqez = tmat(l1)*ek
        s1 = sl(l1)
        rsirc = rs(irc1, l1)
        eklfac = eklfac/ek*real(2*l1+1, kind=dp)
        pn = pz(irc1, l1)*rsirc
        alphal = (bessjw(l1)-ci*hankws(l1)*tlsqez)*rirc/pn

! ======================================================================

        alpha(l1) = alphal*eklfac
! -> calculate irregular wavefunctions
        Do i = 2, irc1
          pz(i, l1) = pz(i, l1)*alphal
          fz(i, l1) = fz(i, l1)*alphal
        End Do

        value = -ci*hankws(l1)*rirc*rsirc
        slope = real(l1+1, kind=dp)*hankws(l1) - rirc*ek*hankws(l1+1)
        slope = (-ci*slope*rsirc+s1/rirc*value)
        qz(irc1, l1) = value
        sz(irc1, l1) = (slope*rirc-(s1+1.0E0_dp)*value)/mass(irc1)*dror(irc1)
      End Do
! ======================================================================
! ======================================================================
!-----------------------------------------------------------------------
!  subroutine for radial wave functions of spherical potentials
      If (lirrsol) Call irwsol(ek, fz, hamf, mass, pz, qz, sz, dror, sl, ipan, &
        ircut, irmd, ipand, lmaxd)

      Do l1 = 0, lmaxd
        If (nsra==2) Then
          Do i = 2, irc1
            pz(i, l1) = pz(i, l1)*rs(i, l1)
            qz(i, l1) = qz(i, l1)/rs(i, l1)
            fz(i, l1) = fz(i, l1)*rs(i, l1)/cvlight
            sz(i, l1) = sz(i, l1)/rs(i, l1)/cvlight
          End Do
        Else
          Do i = 2, irc1
            pz(i, l1) = pz(i, l1)*rs(i, l1)
            qz(i, l1) = qz(i, l1)/rs(i, l1)
            fz(i, l1) = czero
            sz(i, l1) = czero
          End Do
        End If
      End Do
!             the generalized phase shifts are calculated by
    End Subroutine
