subroutine irwsol(ek, fz, hamf, mass, pz, qz, sz, dror, s, ipan, ircut, irmd, &
  ipand, lmaxd)
!-----------------------------------------------------------------------
!  calculates the irregular solution of the schroedinger equation or
!    in semi relativistic approximation for a spherically averaged
!    potential and given energy . to achieve greater precision the
!    leading power r**-s ( in schroedinger case s = l , in case of sra
!    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
!    from the wavefunction .


!  the differential equation is solved with a 5 point adams - bashforth
!    and adams - moulton predictor corrector method integrating
!    inwards and extended for potentials with kinks


!                                               b.drittler   nov.1989
!-----------------------------------------------------------------------
!     ..
!.. Scalar Arguments ..
  double complex :: ek
  integer :: ipan, ipand, irmd, lmaxd
!..
!.. Array Arguments ..
  double complex :: fz(irmd, 0:lmaxd), hamf(irmd, 0:lmaxd), mass(irmd), &
    pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd)
  double precision :: dror(irmd), s(0:lmaxd)
  integer :: ircut(0:ipand)
!..
!.. Local Scalars ..
  double complex :: hamf1, k1f, k1p, k2f, k2p, k3f, k3p, k4f, k4p, mass1, &
    qim0, qim1, sim0, sim1
  double precision :: dror1, drsm1, drsp1, s1, sm1, sp1
  integer :: ip, ir, ire, irs, irwsk, k, l
!..
!.. Local Arrays ..
  double complex :: dqdi(0:4), dsdi(0:4)
!..
!.. Intrinsic Functions ..
  intrinsic :: max
!..
!TIMO      IRWSK = IRCUT(1)/30
  irwsk = 41
  irwsk = ircut(1)/9

  do l = 0, lmaxd
    s1 = s(l)
    sm1 = s1 - 1.0d0
    sp1 = s1 + 1.0d0

!---> loop over kinks

    do ip = ipan, 1, -1
      irs = ircut(ip)
      ire = max(ircut(ip-1), irwsk) + 1
      drsp1 = dror(irs)*sp1
      drsm1 = dror(irs)*sm1
      qim0 = qz(irs, l)
      sim0 = sz(irs, l)
      dqdi(4) = mass(irs)*sim0 + drsp1*qim0
      dsdi(4) = hamf(irs, l)*qim0 + drsm1*sim0

!---> start algorithm - 4 point runge kutta with interpolation

      k1p = dqdi(4)
      k1f = dsdi(4)

      dror1 = (3.0d0*dror(irs-3)-15.0d0*dror(irs-2)+45.0d0*dror(irs-1)+ &
        15.0d0*dror(irs))/48.0d0
      drsp1 = dror1*sp1
      drsm1 = dror1*sm1
      mass1 = (3.0d0*mass(irs-3)-15.0d0*mass(irs-2)+45.0d0*mass(irs-1)+ &
        15.0d0*mass(irs))/48.0d0
      hamf1 = (3.0d0*hamf(irs-3,l)-15.0d0*hamf(irs-2,l)+45.0d0*hamf(irs-1,l)+ &
        15.0d0*hamf(irs,l))/48.0d0
      k2p = mass1*(sim0-0.5d0*k1f) + drsp1*(qim0-0.5d0*k1p)
      k2f = hamf1*(qim0-0.5d0*k1p) + drsm1*(sim0-0.5d0*k1f)
      k3p = mass1*(sim0-0.5d0*k2f) + drsp1*(qim0-0.5d0*k2p)
      k3f = hamf1*(qim0-0.5d0*k2p) + drsm1*(sim0-0.5d0*k2f)
      drsp1 = dror(irs-1)*sp1
      drsm1 = dror(irs-1)*sm1
      k4p = mass(irs-1)*(sim0-k3f) + drsp1*(qim0-k3p)
      k4f = hamf(irs-1, l)*(qim0-k3p) + drsm1*(sim0-k3f)
      qim0 = qim0 - (k1p+2.0d0*(k2p+k3p)+k4p)/6.0d0
      sim0 = sim0 - (k1f+2.0d0*(k2f+k3f)+k4f)/6.0d0
      qz(irs-1, l) = qim0
      sz(irs-1, l) = sim0
      dqdi(3) = mass(irs-1)*sim0 + drsp1*qim0
      dsdi(3) = hamf(irs-1, l)*qim0 + drsm1*sim0

      k = 2

!---> 4 point runge kutta with h = i+2 - 1

      do ir = irs - 2, irs - 4, -1
        qim0 = qz(ir+2, l)
        sim0 = sz(ir+2, l)
        k1p = dqdi(k+2)
        k1f = dsdi(k+2)
        k2p = mass(ir+1)*(sim0-k1f) + drsp1*(qim0-k1p)
        k2f = hamf(ir+1, l)*(qim0-k1p) + drsm1*(sim0-k1f)
        k3p = mass(ir+1)*(sim0-k2f) + drsp1*(qim0-k2p)
        k3f = hamf(ir+1, l)*(qim0-k2p) + drsm1*(sim0-k2f)

        drsp1 = dror(ir)*sp1
        drsm1 = dror(ir)*sm1

        k4p = mass(ir)*(sim0-2.0d0*k3f) + drsp1*(qim0-2.0d0*k3p)
        k4f = hamf(ir, l)*(qim0-2.0d0*k3p) + drsm1*(sim0-2.0d0*k3f)
        qim0 = qim0 - (k1p+2.0d0*(k2p+k3p)+k4p)/3.0d0
        sim0 = sim0 - (k1f+2.0d0*(k2f+k3f)+k4f)/3.0d0
        qz(ir, l) = qim0
        sz(ir, l) = sim0
        dqdi(k) = mass(ir)*sim0 + drsp1*qim0
        dsdi(k) = hamf(ir, l)*qim0 + drsm1*sim0
        k = k - 1
      end do

      do ir = irs - 5, ire, -1

!---> predictor : 5 point adams - bashforth

        qim1 = qim0 - (1901.0d0*dqdi(0)-2774.0d0*dqdi(1)+2616.0d0*dqdi(2)- &
          1274.0d0*dqdi(3)+251.0d0*dqdi(4))/720.0d0
        sim1 = sim0 - (1901.0d0*dsdi(0)-2774.0d0*dsdi(1)+2616.0d0*dsdi(2)- &
          1274.0d0*dsdi(3)+251.0d0*dsdi(4))/720.0d0

        dqdi(4) = dqdi(3)
        dqdi(3) = dqdi(2)
        dqdi(2) = dqdi(1)
        dqdi(1) = dqdi(0)
        dsdi(4) = dsdi(3)
        dsdi(3) = dsdi(2)
        dsdi(2) = dsdi(1)
        dsdi(1) = dsdi(0)

        drsp1 = dror(ir)*sp1
        drsm1 = dror(ir)*sm1

        dqdi(0) = mass(ir)*sim1 + drsp1*qim1
        dsdi(0) = hamf(ir, l)*qim1 + drsm1*sim1

!---> corrector : 5 point adams - moulton

        qim0 = qim0 - (251.0d0*dqdi(0)+646.0d0*dqdi(1)-264.0d0*dqdi(2)+106.0d0 &
          *dqdi(3)-19.0d0*dqdi(4))/720.0d0
        sim0 = sim0 - (251.0d0*dsdi(0)+646.0d0*dsdi(1)-264.0d0*dsdi(2)+106.0d0 &
          *dsdi(3)-19.0d0*dsdi(4))/720.0d0

        dqdi(0) = mass(ir)*sim0 + drsp1*qim0
        dsdi(0) = hamf(ir, l)*qim0 + drsm1*sim0

        qz(ir, l) = qim0
        sz(ir, l) = sim0
      end do

      if (ip/=1) then
        qz(ire-1, l) = qim0
        sz(ire-1, l) = sim0
      end if

    end do


!---> use wronski relation near origin

    do ir = irwsk, 2, -1

!---> 2 point corrector - predictor

      qim1 = qim0 - 1.5d0*dqdi(0) + 0.5d0*dqdi(1)

      dqdi(1) = dqdi(0)
      drsp1 = dror(ir)*sp1

      dqdi(0) = mass(ir)*(1.0d0/ek+qim1*fz(ir,l))/pz(ir, l) + drsp1*qim1

      qim0 = qim0 - 0.5d0*dqdi(0) - 0.5d0*dqdi(1)

      dqdi(0) = mass(ir)*(1.0d0/ek+qim0*fz(ir,l))/pz(ir, l) + drsp1*qim0

      qz(ir, l) = qim0
    end do

    do ir = irwsk, 2, -1
      sz(ir, l) = (1.0d0/ek+qz(ir,l)*fz(ir,l))/pz(ir, l)
    end do
  end do

end subroutine
