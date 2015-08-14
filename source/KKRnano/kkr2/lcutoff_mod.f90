module lcutoff_mod
  implicit none
  private
  public :: distance_pbc, getLMarray, calcCutoffarray

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
    double precision :: vec(3), vectrans(3), dist_sq
    integer :: nx, ny, nz

    vec = point2 - point1

    dist_sq = vec(1) ** 2 + vec(2) ** 2 + vec(3) ** 2

    ! brute force distance checking
    do nx = -1, 1
      do ny = -1, 1
        do nz = -1, 1
          vectrans = vec + nx*bravais(:, 1) + ny*bravais(:, 2) + nz*bravais(:, 3)
          dist_sq = min(dist_sq, vectrans(1) ** 2 + vectrans(2) ** 2 + vectrans(3) ** 2)
        end do
      end do
    end do

    distance_pbc = sqrt(dist_sq)

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

  !----------------------------------------------------------------------------
  !> Modifies the array 'cutoffarray' on positions that correspond to sites that
  !> are further away from 'center' than 'dist_cut'. The value 'lm_low' is written
  !> at the modified positions.
  !>
  !> If merging = .true.: change the lm value only when it is larger than the
  !> original value -> this merges the truncation zones
  subroutine calcCutoffarray(cutoffarray, rbasis, center, bravais, dist_cut, lm_low, merging)
    implicit none
    integer, dimension(:), intent(inout) :: cutoffarray
    double precision, dimension(:,:), intent(in) :: rbasis
    double precision, dimension(3), intent(in) :: center
    double precision, dimension(3,3), intent(in) :: bravais
    double precision, intent(in) :: dist_cut
    integer, intent(in) :: lm_low
    logical, intent(in), optional :: merging
    !-----------------------------

    integer :: ii
    double precision :: dist
    logical :: do_merge

    do_merge = .false.
    if (present(merging)) then
      do_merge = merging
    end if

    do ii = 1, size(rbasis, 2)
      dist = distance_pbc(rbasis(:,ii), center, bravais)
      if (dist > dist_cut) then
        if (.not. do_merge) then
          cutoffarray(ii) = lm_low
        else
          cutoffarray(ii) = max(cutoffarray(ii), lm_low)
        end if
      end if
    end do

  end subroutine

end module lcutoff_mod
