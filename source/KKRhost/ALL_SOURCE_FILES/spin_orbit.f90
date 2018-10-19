!------------------------------------------------------------------------------------
!> Summary: In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real spherical harmonics
!> Author: 
!> In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real 
!> spherical harmonics. Schematically it has the form
!> \begin{equation}
!> \begin{bmatrix} -L_z & \cdots & L_{+} \\ L_{-} & \cdots & L_z \end{bmatrix}
!> \end{equation}
!------------------------------------------------------------------------------------
module mod_spin_orbit
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real spherical harmonics
  !> Author: 
  !> Category: spin-orbit-coupling, special-functions, KKRhost
  !> Deprecated: False 
  !> In this subroutine the matrix \(L\cdot S\) is calculated for the basis of real 
  !> spherical harmonics. Schematically it has the form
  !> \begin{equation}
  !> \begin{bmatrix} -L_z & \cdots & L_{+} \\ L_{-} & \cdots & L_z \end{bmatrix}
  !> \end{equation}
  !-------------------------------------------------------------------------------
  subroutine spin_orbit_one_l(lmax, l_s)

    use :: mod_constants, only: ci
    implicit none

    integer, intent (in) :: lmax
    complex (kind=dp), intent (out) :: l_s((2*lmax+1)*2, (2*lmax+1)*2)

    ! local variables
    integer :: i1, i2, i1l
    complex (kind=dp), allocatable :: l_min(:, :)
    complex (kind=dp), allocatable :: l_up(:, :)
    real (kind=dp) :: lfac

    allocate (l_min(-lmax:lmax,-lmax:lmax))
    allocate (l_up(-lmax:lmax,-lmax:lmax))

    ! initialize the matrix

    do i1 = 1, (2*lmax+1)*2
      do i2 = 1, (2*lmax+1)*2
        l_s(i2, i1) = 0e0_dp
      end do
    end do

    do i1 = -lmax, lmax
      do i2 = -lmax, lmax
        l_min(i2, i1) = 0e0_dp
        l_up(i2, i1) = 0e0_dp
      end do
    end do

    ! fill the second and the forth quadrant with L_z
    ! (-L_z,respectively)
    do i1 = 1, 2*lmax + 1
      i1l = i1 - lmax - 1          ! the value of m (varies from -l to +l)
      i2 = 2*lmax + 1 - (i1-1)

      ! L_S(i2,i1)=ci*i1l
      l_s(i2, i1) = -ci*i1l
    end do

    do i1 = 2*lmax + 2, (2*lmax+1)*2
      i1l = i1 - lmax - 1 - (2*lmax+1) ! the value of m (varies from -l to +l)
      i2 = (2*lmax+1)*2 - (i1-(2*lmax+2))

      ! L_S(i2,i1)=-ci*i1l
      l_s(i2, i1) = ci*i1l

    end do
    ! implement now L_- in the third quadrant

    if (lmax>0) then

      lfac = sqrt(lmax*(lmax+1e0_dp))/sqrt(2e0_dp)
      l_min(0, -1) = -ci*lfac
      ! l_min(0,-1)=ci*lfac
      l_min(0, 1) = lfac
      l_min(-1, 0) = ci*lfac
      l_min(1, 0) = -lfac

      if (lmax>1) then

        do i1 = 2, lmax

          lfac = 0.5e0_dp*sqrt(lmax*(lmax+1e0_dp)-i1*(i1-1e0_dp))
          l_min(-i1, -i1+1) = -lfac
          l_min(-i1, i1-1) = ci*lfac
          l_min(i1, -i1+1) = -ci*lfac
          l_min(i1, i1-1) = -lfac

          lfac = 0.5e0_dp*sqrt(lmax*(lmax+1e0_dp)-(i1-1)*(i1))
          l_min(-i1+1, -i1) = lfac
          l_min(-i1+1, i1) = ci*lfac
          l_min(i1-1, -i1) = -ci*lfac
          l_min(i1-1, i1) = lfac
        end do
      end if
    end if

    do i1 = -lmax, lmax
      do i2 = -lmax, lmax
        l_s(i2+3*lmax+2, i1+lmax+1) = l_min(i1, i2)
      end do
    end do

    ! implement now L_+ in the   quadrant

    if (lmax>0) then

      lfac = sqrt(lmax*(lmax+1e0_dp))/sqrt(2e0_dp)
      l_up(0, -1) = -ci*lfac
      l_up(0, 1) = -lfac
      l_up(-1, 0) = ci*lfac
      l_up(1, 0) = lfac

      if (lmax>1) then
        do i1 = 2, lmax
          lfac = 0.5e0_dp*sqrt(lmax*(lmax+1e0_dp)-i1*(i1-1e0_dp))
          l_up(-i1, -i1+1) = lfac
          l_up(-i1, i1-1) = ci*lfac
          l_up(i1, -i1+1) = -ci*lfac
          l_up(i1, i1-1) = lfac

          lfac = 0.5e0_dp*sqrt(lmax*(lmax+1e0_dp)-(i1-1)*(i1))
          l_up(-i1+1, -i1) = -lfac
          l_up(-i1+1, i1) = ci*lfac
          l_up(i1-1, -i1) = -ci*lfac
          l_up(i1-1, i1) = -lfac

        end do
      end if
    end if

    do i1 = -lmax, lmax
      do i2 = -lmax, lmax
        l_s(i2+lmax+1, i1+3*lmax+2) = l_up(i1, i2)
      end do
    end do

    deallocate (l_min)
    deallocate (l_up)

  end subroutine spin_orbit_one_l

end module mod_spin_orbit
