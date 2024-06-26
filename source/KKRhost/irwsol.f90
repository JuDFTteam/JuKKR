!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_irwsol
  
  private
  public :: irwsol

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculates the irregular solution of the Schroedinger (or SRA) equation
  !> Author: B. Drittler
  !> Date: Nov. 1989
  !> Category: KKRhost, single-site
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Calculates the irregular solution of the Schroedinger equation or
  !> in semi relativistic approximation for a spherically averaged
  !> potential and given energy. To achieve greater precision the
  !> leading power r**-s ( in Schroedinger case s = l , in case of SRA
  !> s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
  !> from the wavefunction.
  !>
  !> The differential equation is solved with a 5 point Adams-Bashforth
  !> and Adams-Moulton predictor corrector method integrating
  !> inwards and extended for potentials with kinks.
  !-------------------------------------------------------------------------------
  subroutine irwsol(ek, fz, hamf, mass, pz, qz, sz, dror, s, ipan, ircut, irmd, ipand, lmaxd)

    use :: mod_datatypes, only: dp
    implicit none
    ! ..
    ! .. Scalar Arguments ..
    complex (kind=dp) :: ek
    integer :: ipan, ipand, irmd, lmaxd
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp) :: fz(irmd, 0:lmaxd), hamf(irmd, 0:lmaxd), mass(irmd), pz(irmd, 0:lmaxd), qz(irmd, 0:lmaxd), sz(irmd, 0:lmaxd)
    real (kind=dp) :: dror(irmd), s(0:lmaxd)
    integer :: ircut(0:ipand)
    ! ..
    ! .. Local Scalars ..
    complex (kind=dp) :: hamf1, k1f, k1p, k2f, k2p, k3f, k3p, k4f, k4p, mass1, qim0, qim1, sim0, sim1
    real (kind=dp) :: dror1, drsm1, drsp1, s1, sm1, sp1
    integer :: ip, ir, ire, irs, irwsk, k, l
    ! ..
    ! .. Local Arrays ..
    complex (kind=dp) :: dqdi(0:4), dsdi(0:4)
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: max


    ! TIMO      IRWSK = IRCUT(1)/30
    irwsk = 41
    irwsk = ircut(1)/9

    do l = 0, lmaxd
      s1 = s(l)
      sm1 = s1 - 1.0e0_dp
      sp1 = s1 + 1.0e0_dp

      ! ---> loop over kinks

      do ip = ipan, 1, -1
        irs = ircut(ip)
        ire = max(ircut(ip-1), irwsk) + 1
        drsp1 = dror(irs)*sp1
        drsm1 = dror(irs)*sm1
        qim0 = qz(irs, l)
        sim0 = sz(irs, l)
        dqdi(4) = mass(irs)*sim0 + drsp1*qim0
        dsdi(4) = hamf(irs, l)*qim0 + drsm1*sim0

        ! ---> start algorithm - 4 point runge kutta with interpolation

        k1p = dqdi(4)
        k1f = dsdi(4)

        dror1 = (3.0e0_dp*dror(irs-3)-15.0e0_dp*dror(irs-2)+45.0e0_dp*dror(irs-1)+15.0e0_dp*dror(irs))/48.0e0_dp
        drsp1 = dror1*sp1
        drsm1 = dror1*sm1
        mass1 = (3.0e0_dp*mass(irs-3)-15.0e0_dp*mass(irs-2)+45.0e0_dp*mass(irs-1)+15.0e0_dp*mass(irs))/48.0e0_dp
        hamf1 = (3.0e0_dp*hamf(irs-3,l)-15.0e0_dp*hamf(irs-2,l)+45.0e0_dp*hamf(irs-1,l)+15.0e0_dp*hamf(irs,l))/48.0e0_dp
        k2p = mass1*(sim0-0.5e0_dp*k1f) + drsp1*(qim0-0.5e0_dp*k1p)
        k2f = hamf1*(qim0-0.5e0_dp*k1p) + drsm1*(sim0-0.5e0_dp*k1f)
        k3p = mass1*(sim0-0.5e0_dp*k2f) + drsp1*(qim0-0.5e0_dp*k2p)
        k3f = hamf1*(qim0-0.5e0_dp*k2p) + drsm1*(sim0-0.5e0_dp*k2f)
        drsp1 = dror(irs-1)*sp1
        drsm1 = dror(irs-1)*sm1
        k4p = mass(irs-1)*(sim0-k3f) + drsp1*(qim0-k3p)
        k4f = hamf(irs-1, l)*(qim0-k3p) + drsm1*(sim0-k3f)
        qim0 = qim0 - (k1p+2.0e0_dp*(k2p+k3p)+k4p)/6.0e0_dp
        sim0 = sim0 - (k1f+2.0e0_dp*(k2f+k3f)+k4f)/6.0e0_dp
        qz(irs-1, l) = qim0
        sz(irs-1, l) = sim0
        dqdi(3) = mass(irs-1)*sim0 + drsp1*qim0
        dsdi(3) = hamf(irs-1, l)*qim0 + drsm1*sim0

        k = 2

        ! ---> 4 point runge kutta with h = i+2 - 1

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

          k4p = mass(ir)*(sim0-2.0e0_dp*k3f) + drsp1*(qim0-2.0e0_dp*k3p)
          k4f = hamf(ir, l)*(qim0-2.0e0_dp*k3p) + drsm1*(sim0-2.0e0_dp*k3f)
          qim0 = qim0 - (k1p+2.0e0_dp*(k2p+k3p)+k4p)/3.0e0_dp
          sim0 = sim0 - (k1f+2.0e0_dp*(k2f+k3f)+k4f)/3.0e0_dp
          qz(ir, l) = qim0
          sz(ir, l) = sim0
          dqdi(k) = mass(ir)*sim0 + drsp1*qim0
          dsdi(k) = hamf(ir, l)*qim0 + drsm1*sim0
          k = k - 1
        end do

        do ir = irs - 5, ire, -1

          ! ---> predictor : 5 point adams - bashforth

          qim1 = qim0 - (1901.0e0_dp*dqdi(0)-2774.0e0_dp*dqdi(1)+2616.0e0_dp*dqdi(2)-1274.0e0_dp*dqdi(3)+251.0e0_dp*dqdi(4))/720.0e0_dp
          sim1 = sim0 - (1901.0e0_dp*dsdi(0)-2774.0e0_dp*dsdi(1)+2616.0e0_dp*dsdi(2)-1274.0e0_dp*dsdi(3)+251.0e0_dp*dsdi(4))/720.0e0_dp

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

          ! ---> corrector : 5 point adams - moulton

          qim0 = qim0 - (251.0e0_dp*dqdi(0)+646.0e0_dp*dqdi(1)-264.0e0_dp*dqdi(2)+106.0e0_dp*dqdi(3)-19.0e0_dp*dqdi(4))/720.0e0_dp
          sim0 = sim0 - (251.0e0_dp*dsdi(0)+646.0e0_dp*dsdi(1)-264.0e0_dp*dsdi(2)+106.0e0_dp*dsdi(3)-19.0e0_dp*dsdi(4))/720.0e0_dp

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


      ! ---> use wronski relation near origin

      do ir = irwsk, 2, -1

        ! ---> 2 point corrector - predictor

        qim1 = qim0 - 1.5e0_dp*dqdi(0) + 0.5e0_dp*dqdi(1)

        dqdi(1) = dqdi(0)
        drsp1 = dror(ir)*sp1

        dqdi(0) = mass(ir)*(1.0e0_dp/ek+qim1*fz(ir,l))/pz(ir, l) + drsp1*qim1

        qim0 = qim0 - 0.5e0_dp*dqdi(0) - 0.5e0_dp*dqdi(1)

        dqdi(0) = mass(ir)*(1.0e0_dp/ek+qim0*fz(ir,l))/pz(ir, l) + drsp1*qim0

        qz(ir, l) = qim0
      end do

      do ir = irwsk, 2, -1
        sz(ir, l) = (1.0e0_dp/ek+qz(ir,l)*fz(ir,l))/pz(ir, l)
      end do
    end do

  end subroutine irwsol

end module mod_irwsol
