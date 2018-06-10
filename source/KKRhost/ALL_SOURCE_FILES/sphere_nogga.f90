!-------------------------------------------------------------------------------
!> @brief Generate an angular mesh and spherical harmonics at those
!> mesh points. For an angular integration the weights are generated .
!
!> @author R. Zeller
!> @date Feb. 1996
!> @note - Jonathan Chico: Rewrote to Fortran90
!-------------------------------------------------------------------------------
subroutine sphere_nogga(lmax, yr, wtyr, rij, ijd)

  use :: constants
! ..
! .. Scalar Arguments
  integer, intent (in) :: ijd
  integer, intent (in) :: lmax !< Maximum l component in wave function expansion
! .. Output variables
  double precision, dimension (ijd, *), intent (out) :: yr
  double precision, dimension (ijd, 3), intent (out) :: rij
  double precision, dimension (ijd, *), intent (out) :: wtyr
! .. Local variables
  integer :: ij, lm1
  double precision :: wght
  double precision :: r, r1, r2, r3
  double precision, dimension (1000) :: y
! .. External Subroutines
  external :: ymy
!
  write (1337, *) ' SPHERE : read LEBEDEV mesh'
  if (ijd>1000) stop ' SPHERE '
!
  do ij = 1, ijd
    call lebedev(ij, r1, r2, r3, wght)
    rij(ij, 1) = r1
    rij(ij, 2) = r2
    rij(ij, 3) = r3
    call ymy(r1, r2, r3, r, y, lmax)
    do lm1 = 1, (lmax+1)**2
      yr(ij, lm1) = y(lm1)
    end do ! LM1
!-------------------------------------------------------------------------
! Multiply the spherical harmonics with the weights
!-------------------------------------------------------------------------
    do lm1 = 1, (lmax+1)**2
      wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.d0
    end do ! LM1
  end do ! IJ

end subroutine
