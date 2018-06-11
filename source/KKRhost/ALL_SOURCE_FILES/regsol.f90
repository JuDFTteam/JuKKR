    Subroutine regsol(cvlight, e, nsra, dlogdp, fz, hamf, mass, pz, dror, r, &
      s, vm2z, z, ipan, ircut, idoldau, lopt, wldauav, cutoff, irmd, ipand, &
      lmaxd)
!-----------------------------------------------------------------------
!  calculates the regular solution of the schroedinger equation or
!    in semi relativistic approximation for a spherically averaged
!    potential and given energy . to archieve greater presion the
!    leading power r**s ( in schroedinger case s = l , in case of sra
!    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
!    from the wavefunction .

!  the t - matrix has to be determined at the mt radius in case of
!    a mt calculation or at the ws radius in case of a ws calcu-
!    lation . therefore the logarithmic derivative is calculated
!    at that point (=ircut(ipan) )

!  the differential equation is solved with a 5 point adams - bashforth
!    and adams - moulton predictor corrector method integrating
!    outwards and extended for potentials with kinks

!                                               b.drittler   nov 1989

!  LDA+U included  March 2003 - Dec 2004, Munich/Juelich
!                                               ph. mavropoulos

!-----------------------------------------------------------------------
      Use mod_datatypes
      Implicit None
!.. Scalar Arguments ..
      Complex (Kind=dp) :: e
      Real (Kind=dp) :: cvlight, z, wldauav
      Integer :: ipan, ipand, irmd, lmaxd, nsra, idoldau, lopt
!..
!.. Array Arguments ..
      Complex (Kind=dp) :: dlogdp(0:lmaxd), fz(irmd, 0:lmaxd), &
        hamf(irmd, 0:lmaxd), mass(irmd), pz(irmd, 0:lmaxd)
      Real (Kind=dp) :: dror(irmd), r(irmd), s(0:lmaxd), vm2z(irmd)
      Real (Kind=dp) :: cutoff(irmd)
      Integer :: ircut(0:ipand)
!..
!.. Local Scalars ..
      Complex (Kind=dp) :: dfd0, dpd0, fip0, fip1, hamf1, k1f, k1p, k2f, k2p, &
        k3f, k3p, k4f, k4p, mass1, pip0, pip1, vme, vmefac, vmetr1
      Real (Kind=dp) :: dror1, drsm1, drsp1, s1, sm1, sp1, srafac
      Integer :: ip, ir, irc, ire, irs, irsp1, j, k, l
!..
!.. Local Arrays ..
      Complex (Kind=dp) :: a(-1:4), b(0:4), dfdi(-4:0), dpdi(-4:0)
      Complex (Kind=dp) :: hamfldau(irmd)
!..
!.. Intrinsic Functions ..
      Intrinsic :: cmplx, real
!     ..
      If (nsra==2) Then

!---> in case of sra  srafac = 1/c - otherwise srafac = 0

        srafac = 1.0E0_dp/cvlight
      Else
        srafac = 0.0E0_dp
      End If

      irc = ircut(ipan)

      Do ir = 2, irc
        vmetr1 = (vm2z(ir)-e)*r(ir) - 2.0E0_dp*z
        hamf(ir, 0) = vmetr1*dror(ir)
        mass(ir) = r(ir) - srafac*srafac*vmetr1
      End Do

      Do l = 1, lmaxd
        Do ir = 7, irc
          hamf(ir, l) = real(l*l+l, kind=dp)/mass(ir)*dror(ir) + hamf(ir, 0)
        End Do
      End Do



! ======================================================================
! LDA+U

!  Account for potential shift in case of LDA+U (averaged over m)
!  by adding the average WLDAUAV to the spherical part of the
!  potential.

      If (idoldau==1 .And. lopt>=0) Then
        s1 = real(lopt*lopt+lopt, kind=dp)
        Do ir = 2, irc
          vmetr1 = (vm2z(ir)-e+wldauav*cutoff(ir))*r(ir)
          hamfldau(ir) = (vmetr1-2.0E0_dp*z)*dror(ir)
        End Do

        Do ir = 7, irc
          hamf(ir, lopt) = s1/mass(ir)*dror(ir) + hamfldau(ir)
        End Do
      End If

! LDA+U
! ======================================================================
      Do ir = 2, irc
        mass(ir) = mass(ir)*dror(ir)
      End Do

      Do l = 0, lmaxd

        s1 = s(l)
        sm1 = s1 - 1.0E0_dp
        sp1 = s1 + 1.0E0_dp

!---> loop over the number of kinks

        Do ip = 1, ipan

          If (ip==1) Then
            irs = 2
            ire = ircut(1)


!---> initial values

            vme = vm2z(2) - e
            vmefac = 1.0E0_dp - vme*srafac*srafac
            If (nsra==2 .And. z>0.0E0_dp) Then
              a(-1) = 0.0E0_dp
              a(0) = 1.0E0_dp
              b(0) = cmplx(sm1*cvlight*cvlight/(2*z), 0.0E0_dp, kind=dp)
              Do j = 1, 3
                a(j) = (0.0E0_dp, 0.E0_dp)
                b(j) = (0.0E0_dp, 0.E0_dp)
              End Do

            Else


              a(0) = 0.0E0_dp
              b(0) = real(l, kind=dp)/vmefac
              a(1) = 1.0E0_dp
              Do j = 2, 4
                a(j) = (vme*vmefac*a(j-2)-2.0E0_dp*z*a(j-1))/ &
                  real((j-1)*(j+2*l), kind=dp)
                b(j-1) = real(l+j-1, kind=dp)*a(j)/vmefac
              End Do

            End If

            k = -4

!---> power series near origin

            Do ir = 2, 6
              pip0 = a(3)
              dpd0 = 3.0E0_dp*a(3)
              fip0 = b(3)
              dfd0 = 3.0E0_dp*b(3)
              Do j = 2, 0, -1
                pip0 = a(j) + pip0*r(ir)
                dpd0 = real(j, kind=dp)*a(j) + dpd0*r(ir)
                fip0 = b(j) + fip0*r(ir)
                dfd0 = real(j, kind=dp)*b(j) + dfd0*r(ir)
              End Do

              pz(ir, l) = pip0
              fz(ir, l) = fip0
              dpdi(k) = dpd0*dror(ir)
              dfdi(k) = dfd0*dror(ir)

              k = k + 1
            End Do

          Else

!---> runge kutta step to restart algorithm

            irs = ircut(ip-1) + 1
            ire = ircut(ip)
            irsp1 = irs + 1
            pip0 = pz(irs, l)
            fip0 = fz(irs, l)
            drsp1 = dror(irs)*sp1
            drsm1 = dror(irs)*sm1
            dpdi(-4) = mass(irs)*fip0 - drsm1*pip0
            dfdi(-4) = hamf(irs, l)*pip0 - drsp1*fip0

!---> first step - 4 point runge kutta with interpolation

            k1p = dpdi(-4)
            k1f = dfdi(-4)

            dror1 = (3.0E0_dp*dror(irs+3)-15.0E0_dp*dror(irs+2)+ &
              45.0E0_dp*dror(irsp1)+15.0E0_dp*dror(irs))/48.0E0_dp
            drsp1 = dror1*sp1
            drsm1 = dror1*sm1
            mass1 = (3.0E0_dp*mass(irs+3)-15.0E0_dp*mass(irs+2)+ &
              45.0E0_dp*mass(irsp1)+15.0E0_dp*mass(irs))/48.0E0_dp
            hamf1 = (3.0E0_dp*hamf(irs+3,l)-15.0E0_dp*hamf(irs+2,l)+ &
              45.0E0_dp*hamf(irsp1,l)+15.0E0_dp*hamf(irs,l))/48.0E0_dp
            k2p = mass1*(fip0+0.5E0_dp*k1f) - drsm1*(pip0+0.5E0_dp*k1p)
            k2f = hamf1*(pip0+0.5E0_dp*k1p) - drsp1*(fip0+0.5E0_dp*k1f)
            k3p = mass1*(fip0+0.5E0_dp*k2f) - drsm1*(pip0+0.5E0_dp*k2p)
            k3f = hamf1*(pip0+0.5E0_dp*k2p) - drsp1*(fip0+0.5E0_dp*k2f)

            drsp1 = dror(irsp1)*sp1
            drsm1 = dror(irsp1)*sm1
            k4p = mass(irsp1)*(fip0+k3f) - drsm1*(pip0+k3p)
            k4f = hamf(irsp1, l)*(pip0+k3p) - drsp1*(fip0+k3f)
            pip0 = pip0 + (k1p+2.0E0_dp*(k2p+k3p)+k4p)/6.0E0_dp
            fip0 = fip0 + (k1f+2.0E0_dp*(k2f+k3f)+k4f)/6.0E0_dp

            pz(irsp1, l) = pip0
            fz(irsp1, l) = fip0
            dpdi(-3) = mass(irsp1)*fip0 - drsm1*pip0
            dfdi(-3) = hamf(irsp1, l)*pip0 - drsp1*fip0

            k = -2

!---> 4 point runge kutta with h = i+2 - i

            Do ir = irs + 2, irs + 4
              pip0 = pz(ir-2, l)
              fip0 = fz(ir-2, l)
              k1p = dpdi(k-2)
              k1f = dfdi(k-2)
              k2p = mass(ir-1)*(fip0+k1f) - drsm1*(pip0+k1p)
              k2f = hamf(ir-1, l)*(pip0+k1p) - drsp1*(fip0+k1f)
              k3p = mass(ir-1)*(fip0+k2f) - drsm1*(pip0+k2p)
              k3f = hamf(ir-1, l)*(pip0+k2p) - drsp1*(fip0+k2f)

              drsp1 = dror(ir)*sp1
              drsm1 = dror(ir)*sm1

              k4p = mass(ir)*(fip0+2.0E0_dp*k3f) - drsm1*(pip0+2.0E0_dp*k3p)
              k4f = hamf(ir, l)*(pip0+2.0E0_dp*k3p) - &
                drsp1*(fip0+2.0E0_dp*k3f)
              pip0 = pip0 + (k1p+2.0E0_dp*(k2p+k3p)+k4p)/3.0E0_dp
              fip0 = fip0 + (k1f+2.0E0_dp*(k2f+k3f)+k4f)/3.0E0_dp

              pz(ir, l) = pip0
              fz(ir, l) = fip0
              dpdi(k) = mass(ir)*fip0 - drsm1*pip0
              dfdi(k) = hamf(ir, l)*pip0 - drsp1*fip0
              k = k + 1
            End Do
          End If

          Do ir = irs + 5, ire
            drsp1 = dror(ir)*sp1
            drsm1 = dror(ir)*sm1

!---> predictor : 5 point adams - bashforth

            pip1 = pip0 + (1901.0E0_dp*dpdi(0)-2774.0E0_dp*dpdi(-1)+ &
              2616.0E0_dp*dpdi(-2)-1274.0E0_dp*dpdi(-3)+251.0E0_dp*dpdi(-4))/ &
              720.0E0_dp
            fip1 = fip0 + (1901.0E0_dp*dfdi(0)-2774.0E0_dp*dfdi(-1)+ &
              2616.0E0_dp*dfdi(-2)-1274.0E0_dp*dfdi(-3)+251.0E0_dp*dfdi(-4))/ &
              720.0E0_dp

            dpdi(-4) = dpdi(-3)
            dpdi(-3) = dpdi(-2)
            dpdi(-2) = dpdi(-1)
            dpdi(-1) = dpdi(0)
            dfdi(-4) = dfdi(-3)
            dfdi(-3) = dfdi(-2)
            dfdi(-2) = dfdi(-1)
            dfdi(-1) = dfdi(0)

            dpdi(0) = mass(ir)*fip1 - drsm1*pip1
            dfdi(0) = hamf(ir, l)*pip1 - drsp1*fip1

!---> corrector : 5 point adams - moulton

            pip0 = pip0 + (251.0E0_dp*dpdi(0)+646.0E0_dp*dpdi(-1)-264.0E0_dp* &
              dpdi(-2)+106.0E0_dp*dpdi(-3)-19.0E0_dp*dpdi(-4))/720.0E0_dp
            fip0 = fip0 + (251.0E0_dp*dfdi(0)+646.0E0_dp*dfdi(-1)-264.0E0_dp* &
              dfdi(-2)+106.0E0_dp*dfdi(-3)-19.0E0_dp*dfdi(-4))/720.0E0_dp

            pz(ir, l) = pip0
            fz(ir, l) = fip0
            dpdi(0) = mass(ir)*fip0 - drsm1*pip0
            dfdi(0) = hamf(ir, l)*pip0 - drsp1*fip0
          End Do

!---> remember that the r - mesh contains the kinks two times
!     store the values of pz and fz to restart the algorithm

          If (ip/=ipan) Then
            pz(ire+1, l) = pip0
            fz(ire+1, l) = fip0
          End If

        End Do


!---> logarithmic derivate of real wavefunction ( r**s *pz / r)

        dlogdp(l) = (dpdi(0)/(pip0*dror(irc))+sm1)/r(irc)
      End Do


    End Subroutine
