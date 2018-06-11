! ************************************************************************
    Subroutine spin_orbit_compl(lmax, lmmaxd, l_s)
      Use mod_datatypes, Only: dp
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

      Implicit None

      Integer, Intent (In) :: lmax, lmmaxd
      Complex (Kind=dp), Intent (Out) :: l_s(lmmaxd*2, lmmaxd*2)

! local variables 
      Integer :: rl, lm1, lm2
      Complex (Kind=dp) :: icompl
      Complex (Kind=dp), Allocatable :: ls_l(:, :)

      icompl = (0E0_dp, 1E0_dp)


      Call cinit((2*lmmaxd)**2, l_s)

      Do rl = 0, lmax

        Allocate (ls_l((2*rl+1)*2,(2*rl+1)*2))
        Call cinit(((2*rl+1)*2)**2, ls_l)


        Call spin_orbit_one_l(rl, ls_l)

        Do lm1 = 1, (2*rl+1)*2

          If (lm1<=2*rl+1) Then
            Do lm2 = 1, (2*rl+1)
              l_s(rl**2+lm1, rl**2+lm2) = 0.5E0_dp*ls_l(lm1, lm2)
            End Do
            Do lm2 = (2*rl+1) + 1, (2*rl+1)*2
              l_s(rl**2+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) = 0.5E0_dp* &
                ls_l(lm1, lm2)
            End Do
          Else
            Do lm2 = 1, (2*rl+1)
              l_s(lmmaxd+rl**2-(2*rl+1)+lm1, rl**2+lm2) = 0.5E0_dp* &
                ls_l(lm1, lm2)
            End Do
            Do lm2 = (2*rl+1) + 1, (2*rl+1)*2
              l_s(lmmaxd+rl**2-(2*rl+1)+lm1, lmmaxd+rl**2-(2*rl+1)+lm2) &
                = 0.5E0_dp*ls_l(lm1, lm2)
            End Do
          End If

        End Do !lm1

        Deallocate (ls_l)


      End Do !rl=0,lmax


    End Subroutine
