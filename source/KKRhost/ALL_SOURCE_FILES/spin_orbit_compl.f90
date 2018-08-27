module mod_spin_orbit_compl

contains

! ************************************************************************
subroutine spin_orbit_compl(lmax, lmmaxd, l_s)
  ! ************************************************************************
  ! in this subroutine the matrix L*S is calculated for the basis of
  ! real spherical harmonics

  use mod_datatypes, only: dp
  use mod_spin_orbit, only: spin_orbit_one_l
   use mod_cinit
  implicit none

  integer, intent (in) :: lmax, lmmaxd
  complex (kind=dp), intent (out) :: l_s(lmmaxd*2, lmmaxd*2)

  ! local variables
  integer :: rl, lm1, lm2
  complex (kind=dp) :: icompl
  complex (kind=dp), allocatable :: ls_l(:, :)

  icompl = (0e0_dp, 1e0_dp)


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
          l_s(lmmaxd+rl**2-(2*rl+1)+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) &
            = 0.5e0_dp*ls_l(lm1, lm2)
        end do
      end if

    end do                         ! lm1

    deallocate (ls_l)


  end do                           ! rl=0,lmax


end subroutine spin_orbit_compl

end module mod_spin_orbit_compl
