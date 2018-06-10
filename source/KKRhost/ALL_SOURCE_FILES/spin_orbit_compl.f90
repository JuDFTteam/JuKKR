! ************************************************************************
subroutine spin_orbit_compl(lmax, lmmaxd, l_s)
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

  implicit none

  integer, intent (in) :: lmax, lmmaxd
  double complex, intent (out) :: l_s(lmmaxd*2, lmmaxd*2)

! local variables 
  integer :: rl, lm1, lm2
  double complex :: icompl
  double complex, allocatable :: ls_l(:, :)

  icompl = (0d0, 1d0)


  call cinit((2*lmmaxd)**2, l_s)

  do rl = 0, lmax

    allocate (ls_l((2*rl+1)*2,(2*rl+1)*2))
    call cinit(((2*rl+1)*2)**2, ls_l)


    call spin_orbit_one_l(rl, ls_l)

    do lm1 = 1, (2*rl+1)*2

      if (lm1<=2*rl+1) then
        do lm2 = 1, (2*rl+1)
          l_s(rl**2+lm1, rl**2+lm2) = 0.5d0*ls_l(lm1, lm2)
        end do
        do lm2 = (2*rl+1) + 1, (2*rl+1)*2
          l_s(rl**2+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) = 0.5d0*ls_l(lm1, lm2)
        end do
      else
        do lm2 = 1, (2*rl+1)
          l_s(lmmaxd+rl**2-(2*rl+1)+lm1, rl**2+lm2) = 0.5d0*ls_l(lm1, lm2)
        end do
        do lm2 = (2*rl+1) + 1, (2*rl+1)*2
          l_s(lmmaxd+rl**2-(2*rl+1)+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) = 0.5d0* &
            ls_l(lm1, lm2)
        end do
      end if

    end do !lm1

    deallocate (ls_l)


  end do !rl=0,lmax


end subroutine
