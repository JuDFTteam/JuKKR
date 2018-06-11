!-------------------------------------------------------------------------------
!> @brief Generate an angular mesh and spherical harmonics at those
!> mesh points. For an angular integration the weights are generated .
!
!> @author R. Zeller
!> @date Feb. 1996
!> @note - Jonathan Chico: Rewrote to Fortran90
!-------------------------------------------------------------------------------
    Subroutine sphere_nogga(lmax, yr, wtyr, rij, ijd)

      Use constants
      Use mod_datatypes, Only: dp
! ..
! .. Scalar Arguments
      Integer, Intent (In) :: ijd
      Integer, Intent (In) :: lmax !< Maximum l component in wave function expansion
! .. Output variables
      Real (Kind=dp), Dimension (ijd, *), Intent (Out) :: yr
      Real (Kind=dp), Dimension (ijd, 3), Intent (Out) :: rij
      Real (Kind=dp), Dimension (ijd, *), Intent (Out) :: wtyr
! .. Local variables
      Integer :: ij, lm1
      Real (Kind=dp) :: wght
      Real (Kind=dp) :: r, r1, r2, r3
      Real (Kind=dp), Dimension (1000) :: y
! .. External Subroutines
      External :: ymy
!
      Write (1337, *) ' SPHERE : read LEBEDEV mesh'
      If (ijd>1000) Stop ' SPHERE '
!
      Do ij = 1, ijd
        Call lebedev(ij, r1, r2, r3, wght)
        rij(ij, 1) = r1
        rij(ij, 2) = r2
        rij(ij, 3) = r3
        Call ymy(r1, r2, r3, r, y, lmax)
        Do lm1 = 1, (lmax+1)**2
          yr(ij, lm1) = y(lm1)
        End Do ! LM1
!-------------------------------------------------------------------------
! Multiply the spherical harmonics with the weights
!-------------------------------------------------------------------------
        Do lm1 = 1, (lmax+1)**2
          wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.E0_dp
        End Do ! LM1
      End Do ! IJ

    End Subroutine
