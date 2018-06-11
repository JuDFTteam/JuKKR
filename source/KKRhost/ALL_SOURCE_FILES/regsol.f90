subroutine regsol(cvlight, e, nsra, dlogdp, fz, hamf, mass, pz, dror, r, s, &
  vm2z, z, ipan, ircut, idoldau, lopt, wldauav, cutoff, irmd, ipand, lmaxd)
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
  use mod_DataTypes
  implicit none
!.. Scalar Arguments ..
  double complex :: e
  double precision :: cvlight, z, wldauav
  integer :: ipan, ipand, irmd, lmaxd, nsra, idoldau, lopt
!..
!.. Array Arguments ..
  double complex :: dlogdp(0:lmaxd), fz(irmd, 0:lmaxd), hamf(irmd, 0:lmaxd), &
    mass(irmd), pz(irmd, 0:lmaxd)
  double precision :: dror(irmd), r(irmd), s(0:lmaxd), vm2z(irmd)
  double precision :: cutoff(irmd)
  integer :: ircut(0:ipand)
!..
!.. Local Scalars ..
  double complex :: dfd0, dpd0, fip0, fip1, hamf1, k1f, k1p, k2f, k2p, k3f, &
    k3p, k4f, k4p, mass1, pip0, pip1, vme, vmefac, vmetr1
  double precision :: dror1, drsm1, drsp1, s1, sm1, sp1, srafac
  integer :: ip, ir, irc, ire, irs, irsp1, j, k, l
!..
!.. Local Arrays ..
  double complex :: a(-1:4), b(0:4), dfdi(-4:0), dpdi(-4:0)
  double complex :: hamfldau(irmd)
!..
!.. Intrinsic Functions ..
  intrinsic :: cmplx, dble
!     ..
  if (nsra==2) then

!---> in case of sra  srafac = 1/c - otherwise srafac = 0

    srafac = 1.0d0/cvlight
  else
    srafac = 0.0d0
  end if

  irc = ircut(ipan)

  do ir = 2, irc
    vmetr1 = (vm2z(ir)-e)*r(ir) - 2.0d0*z
    hamf(ir, 0) = vmetr1*dror(ir)
    mass(ir) = r(ir) - srafac*srafac*vmetr1
  end do

  do l = 1, lmaxd
    do ir = 7, irc
      hamf(ir, l) = dble(l*l+l)/mass(ir)*dror(ir) + hamf(ir, 0)
    end do
  end do



! ======================================================================
! LDA+U

!  Account for potential shift in case of LDA+U (averaged over m)
!  by adding the average WLDAUAV to the spherical part of the
!  potential.

  if (idoldau==1 .and. lopt>=0) then
    s1 = dble(lopt*lopt+lopt)
    do ir = 2, irc
      vmetr1 = (vm2z(ir)-e+wldauav*cutoff(ir))*r(ir)
      hamfldau(ir) = (vmetr1-2.0d0*z)*dror(ir)
    end do

    do ir = 7, irc
      hamf(ir, lopt) = s1/mass(ir)*dror(ir) + hamfldau(ir)
    end do
  end if

! LDA+U
! ======================================================================
  do ir = 2, irc
    mass(ir) = mass(ir)*dror(ir)
  end do

  do l = 0, lmaxd

    s1 = s(l)
    sm1 = s1 - 1.0d0
    sp1 = s1 + 1.0d0

!---> loop over the number of kinks

    do ip = 1, ipan

      if (ip==1) then
        irs = 2
        ire = ircut(1)


!---> initial values

        vme = vm2z(2) - e
        vmefac = 1.0d0 - vme*srafac*srafac
        if (nsra==2 .and. z>0.0d0) then
          a(-1) = 0.0d0
          a(0) = 1.0d0
          b(0) = cmplx(sm1*cvlight*cvlight/(2*z), 0.0d0, kind=dp)
          do j = 1, 3
            a(j) = (0.0d0, 0.d0)
            b(j) = (0.0d0, 0.d0)
          end do

        else


          a(0) = 0.0d0
          b(0) = dble(l)/vmefac
          a(1) = 1.0d0
          do j = 2, 4
            a(j) = (vme*vmefac*a(j-2)-2.0d0*z*a(j-1))/dble((j-1)*(j+2*l))
            b(j-1) = dble(l+j-1)*a(j)/vmefac
          end do

        end if

        k = -4

!---> power series near origin

        do ir = 2, 6
          pip0 = a(3)
          dpd0 = 3.0d0*a(3)
          fip0 = b(3)
          dfd0 = 3.0d0*b(3)
          do j = 2, 0, -1
            pip0 = a(j) + pip0*r(ir)
            dpd0 = dble(j)*a(j) + dpd0*r(ir)
            fip0 = b(j) + fip0*r(ir)
            dfd0 = dble(j)*b(j) + dfd0*r(ir)
          end do

          pz(ir, l) = pip0
          fz(ir, l) = fip0
          dpdi(k) = dpd0*dror(ir)
          dfdi(k) = dfd0*dror(ir)

          k = k + 1
        end do

      else

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

        dror1 = (3.0d0*dror(irs+3)-15.0d0*dror(irs+2)+45.0d0*dror(irsp1)+ &
          15.0d0*dror(irs))/48.0d0
        drsp1 = dror1*sp1
        drsm1 = dror1*sm1
        mass1 = (3.0d0*mass(irs+3)-15.0d0*mass(irs+2)+45.0d0*mass(irsp1)+ &
          15.0d0*mass(irs))/48.0d0
        hamf1 = (3.0d0*hamf(irs+3,l)-15.0d0*hamf(irs+2,l)+ &
          45.0d0*hamf(irsp1,l)+15.0d0*hamf(irs,l))/48.0d0
        k2p = mass1*(fip0+0.5d0*k1f) - drsm1*(pip0+0.5d0*k1p)
        k2f = hamf1*(pip0+0.5d0*k1p) - drsp1*(fip0+0.5d0*k1f)
        k3p = mass1*(fip0+0.5d0*k2f) - drsm1*(pip0+0.5d0*k2p)
        k3f = hamf1*(pip0+0.5d0*k2p) - drsp1*(fip0+0.5d0*k2f)

        drsp1 = dror(irsp1)*sp1
        drsm1 = dror(irsp1)*sm1
        k4p = mass(irsp1)*(fip0+k3f) - drsm1*(pip0+k3p)
        k4f = hamf(irsp1, l)*(pip0+k3p) - drsp1*(fip0+k3f)
        pip0 = pip0 + (k1p+2.0d0*(k2p+k3p)+k4p)/6.0d0
        fip0 = fip0 + (k1f+2.0d0*(k2f+k3f)+k4f)/6.0d0

        pz(irsp1, l) = pip0
        fz(irsp1, l) = fip0
        dpdi(-3) = mass(irsp1)*fip0 - drsm1*pip0
        dfdi(-3) = hamf(irsp1, l)*pip0 - drsp1*fip0

        k = -2

!---> 4 point runge kutta with h = i+2 - i

        do ir = irs + 2, irs + 4
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

          k4p = mass(ir)*(fip0+2.0d0*k3f) - drsm1*(pip0+2.0d0*k3p)
          k4f = hamf(ir, l)*(pip0+2.0d0*k3p) - drsp1*(fip0+2.0d0*k3f)
          pip0 = pip0 + (k1p+2.0d0*(k2p+k3p)+k4p)/3.0d0
          fip0 = fip0 + (k1f+2.0d0*(k2f+k3f)+k4f)/3.0d0

          pz(ir, l) = pip0
          fz(ir, l) = fip0
          dpdi(k) = mass(ir)*fip0 - drsm1*pip0
          dfdi(k) = hamf(ir, l)*pip0 - drsp1*fip0
          k = k + 1
        end do
      end if

      do ir = irs + 5, ire
        drsp1 = dror(ir)*sp1
        drsm1 = dror(ir)*sm1

!---> predictor : 5 point adams - bashforth

        pip1 = pip0 + (1901.0d0*dpdi(0)-2774.0d0*dpdi(-1)+2616.0d0*dpdi(-2)- &
          1274.0d0*dpdi(-3)+251.0d0*dpdi(-4))/720.0d0
        fip1 = fip0 + (1901.0d0*dfdi(0)-2774.0d0*dfdi(-1)+2616.0d0*dfdi(-2)- &
          1274.0d0*dfdi(-3)+251.0d0*dfdi(-4))/720.0d0

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

        pip0 = pip0 + (251.0d0*dpdi(0)+646.0d0*dpdi(-1)-264.0d0*dpdi(-2)+ &
          106.0d0*dpdi(-3)-19.0d0*dpdi(-4))/720.0d0
        fip0 = fip0 + (251.0d0*dfdi(0)+646.0d0*dfdi(-1)-264.0d0*dfdi(-2)+ &
          106.0d0*dfdi(-3)-19.0d0*dfdi(-4))/720.0d0

        pz(ir, l) = pip0
        fz(ir, l) = fip0
        dpdi(0) = mass(ir)*fip0 - drsm1*pip0
        dfdi(0) = hamf(ir, l)*pip0 - drsp1*fip0
      end do

!---> remember that the r - mesh contains the kinks two times
!     store the values of pz and fz to restart the algorithm

      if (ip/=ipan) then
        pz(ire+1, l) = pip0
        fz(ire+1, l) = fip0
      end if

    end do


!---> logarithmic derivate of real wavefunction ( r**s *pz / r)

    dlogdp(l) = (dpdi(0)/(pip0*dror(irc))+sm1)/r(irc)
  end do


end subroutine
