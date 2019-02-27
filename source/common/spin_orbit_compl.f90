!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real spherical harmonics
!> Author: 
!> In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real 
!> spherical harmonics
!------------------------------------------------------------------------------------
module mod_spin_orbit_compl
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real spherical harmonics
  !> Author: 
  !> Category: spin-orbit-coupling, special-functions, KKRhost
  !> Deprecated: False 
  !> In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real spherical harmonics
  !-------------------------------------------------------------------------------
  subroutine spin_orbit_compl(lmax, lmmaxd, l_s)

    use :: mod_spin_orbit, only: spin_orbit_one_l
    use :: mod_cinit, only: cinit
    implicit none

    integer, intent (in) :: lmax, lmmaxd
    complex (kind=dp), dimension(lmmaxd*2, lmmaxd*2), intent (out) :: l_s

    ! local variables
    integer :: rl, lm1, lm2
    complex (kind=dp), dimension(:,:), allocatable :: ls_l


    call cinit((2*lmmaxd)**2, l_s)

    do rl = 0, lmax
      allocate (ls_l((2*rl+1)*2,(2*rl+1)*2))
      call cinit(((2*rl+1)*2)**2, ls_l)
      call spin_orbit_one_l(rl, ls_l)
      do lm1 = 1, (2*rl+1)*2
        if (lm1<=2*rl+1) then
          do lm2 = 1, (2*rl+1)
            l_s(rl**2+lm1, rl**2+lm2) = 0.5e0_dp*ls_l(lm1, lm2)
          end do
          do lm2 = (2*rl+1) + 1, (2*rl+1)*2
            l_s(rl**2+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) = 0.5e0_dp*ls_l(lm1, lm2)
          end do
        else
          do lm2 = 1, (2*rl+1)
            l_s(lmmaxd+rl**2-(2*rl+1)+lm1, rl**2+lm2) = 0.5e0_dp*ls_l(lm1, lm2)
          end do
          do lm2 = (2*rl+1) + 1, (2*rl+1)*2
            l_s(lmmaxd+rl**2-(2*rl+1)+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) = 0.5e0_dp*ls_l(lm1, lm2)
          end do
        end if
      end do                       ! lm1
      deallocate (ls_l)
    end do                         ! rl=0,lmax


  end subroutine spin_orbit_compl

end module mod_spin_orbit_compl
