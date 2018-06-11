! ************************************************************************
    Subroutine spin_orbit_one_l(lmax, l_s)
      Use mod_datatypes, Only: dp
! ************************************************************************
!      in this subroutine the matrix L*S is calculated for the basis of
!      real spherical harmonics

!      schematically it has the form
!      (  -L_z    L_+  )
!      (  L_-     L_z  )
      Implicit None

      Integer, Intent (In) :: lmax
      Complex (Kind=dp), Intent (Out) :: l_s((2*lmax+1)*2, (2*lmax+1)*2)

!  local variables 
      Integer :: i1, i2, i1l
      Complex (Kind=dp) :: icompl
      Complex (Kind=dp), Allocatable :: l_min(:, :)
      Complex (Kind=dp), Allocatable :: l_up(:, :)
      Real (Kind=dp) :: lfac


      icompl = (0E0_dp, 1E0_dp)


      Allocate (l_min(-lmax:lmax,-lmax:lmax))
      Allocate (l_up(-lmax:lmax,-lmax:lmax))

!  initialize the matrix

      Do i1 = 1, (2*lmax+1)*2
        Do i2 = 1, (2*lmax+1)*2
          l_s(i2, i1) = 0E0_dp
        End Do
      End Do

      Do i1 = -lmax, lmax
        Do i2 = -lmax, lmax
          l_min(i2, i1) = 0E0_dp
          l_up(i2, i1) = 0E0_dp
        End Do
      End Do

!  fill the second and the forth quadrant with L_z
! (-L_z,respectively)


      Do i1 = 1, 2*lmax + 1
        i1l = i1 - lmax - 1 ! the value of m (varies from -l to +l)
        i2 = 2*lmax + 1 - (i1-1)

!         L_S(i2,i1)=icompl*i1l
        l_s(i2, i1) = -icompl*i1l

      End Do

      Do i1 = 2*lmax + 2, (2*lmax+1)*2
        i1l = i1 - lmax - 1 - (2*lmax+1) ! the value of m (varies from -l to +l)
        i2 = (2*lmax+1)*2 - (i1-(2*lmax+2))

!         L_S(i2,i1)=-icompl*i1l
        l_s(i2, i1) = icompl*i1l

      End Do


!  implement now L_- in the third quadrant

      If (lmax>0) Then

        lfac = sqrt(lmax*(lmax+1E0_dp))/sqrt(2E0_dp)
        l_min(0, -1) = -icompl*lfac
!         l_min(0,-1)=icompl*lfac
        l_min(0, 1) = lfac
        l_min(-1, 0) = icompl*lfac
        l_min(1, 0) = -lfac

        If (lmax>1) Then

          Do i1 = 2, lmax

            lfac = 0.5E0_dp*sqrt(lmax*(lmax+1E0_dp)-i1*(i1-1E0_dp))
            l_min(-i1, -i1+1) = -lfac
            l_min(-i1, i1-1) = icompl*lfac
            l_min(i1, -i1+1) = -icompl*lfac
            l_min(i1, i1-1) = -lfac

            lfac = 0.5E0_dp*sqrt(lmax*(lmax+1E0_dp)-(i1-1)*(i1))
            l_min(-i1+1, -i1) = lfac
            l_min(-i1+1, i1) = icompl*lfac
            l_min(i1-1, -i1) = -icompl*lfac
            l_min(i1-1, i1) = lfac

          End Do

        End If
      End If


      Do i1 = -lmax, lmax
        Do i2 = -lmax, lmax
          l_s(i2+3*lmax+2, i1+lmax+1) = l_min(i1, i2)
        End Do
      End Do


!  implement now L_+ in the   quadrant

      If (lmax>0) Then

        lfac = sqrt(lmax*(lmax+1E0_dp))/sqrt(2E0_dp)
        l_up(0, -1) = -icompl*lfac
        l_up(0, 1) = -lfac
        l_up(-1, 0) = icompl*lfac
        l_up(1, 0) = lfac

        If (lmax>1) Then

          Do i1 = 2, lmax

            lfac = 0.5E0_dp*sqrt(lmax*(lmax+1E0_dp)-i1*(i1-1E0_dp))
            l_up(-i1, -i1+1) = lfac
            l_up(-i1, i1-1) = icompl*lfac
            l_up(i1, -i1+1) = -icompl*lfac
            l_up(i1, i1-1) = lfac

            lfac = 0.5E0_dp*sqrt(lmax*(lmax+1E0_dp)-(i1-1)*(i1))
            l_up(-i1+1, -i1) = -lfac
            l_up(-i1+1, i1) = icompl*lfac
            l_up(i1-1, -i1) = -icompl*lfac
            l_up(i1-1, i1) = -lfac

          End Do

        End If
      End If


      Do i1 = -lmax, lmax
        Do i2 = -lmax, lmax
          l_s(i2+lmax+1, i1+3*lmax+2) = l_up(i1, i2)
        End Do
      End Do



      Deallocate (l_min)
      Deallocate (l_up)


    End Subroutine
