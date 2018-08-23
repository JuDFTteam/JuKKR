module mod_setfactl

contains

subroutine setfactl(factl, lmax, krel, lmmaxd)
  use :: mod_datatypes, only: dp
  implicit none
  ! ..
  ! .. Parameters
  complex (kind=dp) :: ci
  parameter (ci=(0e0_dp,1e0_dp))
  ! ..
  ! .. Arguments
  integer :: krel, lmax, lmmaxd
  complex (kind=dp) :: factl(lmmaxd, lmmaxd)
  ! ..
  ! .. Locals
  integer :: ii1, ii2, l1, l2, lm1, lm2, mm1, mm2, imu1, imu2
  integer :: kap1(2), kap2(2), nsol1, nsol2
  real (kind=dp) :: mu1, mu2, mu1m05, mu2m05
  ! ..
  ! .. Externals
  external :: cinit
  ! ..
  call cinit(lmmaxd*lmmaxd, factl)

  ! ----------------------------------------------------------------------
  if (krel==0) then
    lm1 = 0
    do l1 = 0, lmax
      do mm1 = 1, 2*l1 + 1
        lm1 = lm1 + 1
        lm2 = 0
        do l2 = 0, lmax
          do mm2 = 1, 2*l2 + 1
            lm2 = lm2 + 1
            factl(lm1, lm2) = (-1)**(l1+l2)
          end do
        end do
      end do
    end do
    ! ----------------------------------------------------------------------
  else                             ! KREL.EQ.1
    ! ----------------------------------------------------------------------
    lm1 = 0
    ! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
    do l1 = 0, lmax
      kap1(1) = l1
      kap1(2) = -l1 - 1
      nsol1 = 2
      if (l1==0) then
        kap1(1) = kap1(2)
        nsol1 = 1
      end if
      do ii1 = 1, nsol1
        mu1m05 = abs(kap1(ii1)) - 0.5e0_dp
        ! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
        ! DO MU1 = -MU1M05,MU1M05,1D0
        do imu1 = 1, 2*nint(mu1m05) + 1
          mu1 = -mu1m05 + real(imu1-1, kind=dp)
          lm1 = lm1 + 1
          lm2 = 0
          ! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
          do l2 = 0, lmax
            kap2(1) = l2
            kap2(2) = -l2 - 1
            nsol2 = 2
            if (l2==0) then
              kap2(1) = kap2(2)
              nsol2 = 1
            end if
            do ii2 = 1, nsol2
              mu2m05 = abs(kap2(ii2)) - 0.5e0_dp
              ! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
              ! DO MU2 = -MU2M05,MU2M05,1D0
              do imu2 = 1, 2*nint(mu2m05) + 1
                mu2 = -mu2m05 + real(imu2-1, kind=dp)
                lm2 = lm2 + 1
                mm1 = int(mu2-mu1)

                factl(lm1, lm2) = (-1)**(l1+l2)*ci**mm1

              end do
              ! m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2m2
            end do
          end do
          ! l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2l2
        end do
        ! m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1m1
      end do
    end do
    ! l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1l1
  end if
  ! ----------------------------------------------------------------------
end subroutine setfactl

end module mod_setfactl
