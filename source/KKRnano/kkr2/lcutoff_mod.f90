module lcutoff_mod
  implicit none

  contains

  !> Calculate distance of two points taking into account the
  !> periodic boundary conditions
  !> Assume that points are in same unit cell
  !> Is this a valid assumption???
  double precision function distance_pbc(point1, point2, bravais)
    implicit none
    double precision, intent(in) :: point1(3)
    double precision, intent(in) :: point2(3)
    double precision, intent(in) :: bravais(3,3)
    !----------------
    double precision :: vec(3), frac(3)
    integer :: ii

    vec = point2 - point1

    do ii = 1,3
      ! calculate fractional coordinates
      frac(ii) = vec(1) * bravais(1,ii) + &
                 vec(2) * bravais(2,ii) + &
                 vec(3) * bravais(3,ii)

      frac(ii) = frac(ii) / sqrt(bravais(1,ii)**2 + &
                                 bravais(2,ii)**2 + &
                                 bravais(3,ii)**2)

      ! apply periodic boundary conditions
      if (frac(ii) < -0.5d0) frac(ii) = frac(ii) + 1.0d0
      if (frac(ii) >= 0.5d0) frac(ii) = frac(ii) - 1.0d0
    end do

    ! back to cartesian coordinates
    vec = frac(1) * bravais(:,1) + &
          frac(2) * bravais(:,2) + &
          frac(3) * bravais(:,3)

    ! calc. distance
    distance_pbc = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)

  end function

  subroutine getLMarray(lmarray, rbasis, center, bravais, dist_cut, lm_high, lm_low)
    implicit none
    integer, dimension(:), intent(out) :: lmarray
    double precision, dimension(:,:), intent(in) :: rbasis
    double precision, dimension(3), intent(in) :: center
    double precision, dimension(3,3), intent(in) :: bravais
    double precision, intent(in) :: dist_cut
    integer, intent(in):: lm_high
    integer, intent(in) :: lm_low
    !-----------------------------

    integer :: ii
    double precision :: dist

    do ii = 1, size(rbasis, 2)
      dist = distance_pbc(rbasis(:,ii), center, bravais)
      if (dist > dist_cut) then
        lmarray(ii) = lm_low
      else
        lmarray(ii) = lm_high
      end if
    end do

  end subroutine
end module lcutoff_mod
