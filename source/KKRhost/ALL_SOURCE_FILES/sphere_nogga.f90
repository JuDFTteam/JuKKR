module mod_sphere_nogga

contains

! -------------------------------------------------------------------------------
! > @brief Generate an angular mesh and spherical harmonics at those
! > mesh points. For an angular integration the weights are generated .

! > @author R. Zeller
! > @date Feb. 1996
! > @note - Jonathan Chico: Rewrote to Fortran90
! -------------------------------------------------------------------------------
subroutine sphere_nogga(lmax, yr, wtyr, rij, ijd)

  use :: constants
  use :: mod_datatypes, only: dp
   use mod_lebedev
   use mod_ymy
  ! ..
  ! .. Scalar Arguments
  integer, intent (in) :: ijd
  integer, intent (in) :: lmax     ! < Maximum l component in wave function
                                   ! expansion
  ! .. Output variables
  real (kind=dp), dimension (ijd, *), intent (out) :: yr
  real (kind=dp), dimension (ijd, 3), intent (out) :: rij
  real (kind=dp), dimension (ijd, *), intent (out) :: wtyr
  ! .. Local variables
  integer :: ij, lm1
  real (kind=dp) :: wght
  real (kind=dp) :: r, r1, r2, r3
  real (kind=dp), dimension (1000) :: y
  ! .. External Subroutines
  external :: ymy

  write (1337, *) ' SPHERE : read LEBEDEV mesh'
  if (ijd>1000) stop ' SPHERE '

  do ij = 1, ijd
    call lebedev(ij, r1, r2, r3, wght)
    rij(ij, 1) = r1
    rij(ij, 2) = r2
    rij(ij, 3) = r3
    call ymy(r1, r2, r3, r, y, lmax)
    do lm1 = 1, (lmax+1)**2
      yr(ij, lm1) = y(lm1)
    end do                         ! LM1
    ! -------------------------------------------------------------------------
    ! Multiply the spherical harmonics with the weights
    ! -------------------------------------------------------------------------
    do lm1 = 1, (lmax+1)**2
      wtyr(ij, lm1) = yr(ij, lm1)*wght*pi*4.e0_dp
    end do                         ! LM1
  end do                           ! IJ

end subroutine sphere_nogga

end module mod_sphere_nogga
