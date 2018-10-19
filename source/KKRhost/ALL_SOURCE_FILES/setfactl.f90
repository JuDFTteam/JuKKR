!------------------------------------------------------------------------------------
!> Summary: Determination of the `factl` needed to match the Greens functions in the decimation method
!> Author: People who wrote it
!> Determination of the `factl` needed to match the Greens functions in the decimation method
!------------------------------------------------------------------------------------
module mod_setfactl
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Determination of the `factl` needed to match the Greens functions in the decimation method
  !> Author:
  !> Category: reference-system, structural-greensfunction, KKRhost
  !> Deprecated: False 
  !> Determination of the `factl` needed to match the Greens functions in the decimation 
  !> method.
  !-------------------------------------------------------------------------------
  subroutine setfactl(factl, lmax, krel, lmmaxd)
    use :: mod_cinit
    use :: constants, only: ci
    implicit none
    ! ..
    ! .. Arguments
    integer :: krel, lmax, lmmaxd
    complex (kind=dp) :: factl(lmmaxd, lmmaxd)
    ! ..
    ! .. Locals
    integer :: ii1, ii2, l1, l2, lm1, lm2, mm1, mm2
    integer :: kap1(2), kap2(2), nsol1, nsol2
    real (kind=dp) :: mu1, mu2, mu1m05, mu2m05
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
    else                           ! KREL.EQ.1
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
          do mu1 = -mu1m05, mu1m05, 1.0_dp
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
                do mu2 = -mu2m05, mu2m05, 1.0_dp
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
