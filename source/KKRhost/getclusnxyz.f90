!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_getclusnxyz
  
  private
  public :: getclusnxyz

contains

  !-------------------------------------------------------------------------------
  !> Summary: Find integer boudaries (multiplied with lattice vectors) of spherical cluster
  !> Author: 
  !> Category: KKRhost, geometry, reference-system
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Given a spherical cluster of radius CLURAD it determines the three
  !> integers N1,N2,N3 such that any vector
  !>
  !>   R_i = r_i + SUM_j  N_j * a_j
  !>
  !> with i = 1,NAEZ and a_j the primitive Bravais vectors, is inside 
  !> the cluster. Subroutine also returns the CLURAD**2 value
  !-------------------------------------------------------------------------------
  subroutine getclusnxyz(clurad, bravais, ndim, cluradsq, nbr)

    use :: mod_datatypes, only: dp
    implicit none

    real (kind=dp) :: clurad, cluradsq
    real (kind=dp) :: bravais(3, 3)
    integer :: ndim, nbr(3)

    real (kind=dp) :: dr(3)
    integer :: i, j


    do i = 1, ndim
      dr(i) = 0e0_dp
      do j = 1, ndim
        dr(i) = dr(i) + bravais(j, i)*bravais(j, i)
      end do
      dr(i) = sqrt(dr(i))
    end do

    if (abs(clurad)<1e-6_dp) then
      do i = 1, ndim
        nbr(i) = 0
      end do
      cluradsq = 1e10_dp
    else
      do i = 1, ndim
        nbr(i) = int(clurad/dr(i)) + 2
      end do
      cluradsq = clurad*clurad
    end if

  end subroutine getclusnxyz

end module mod_getclusnxyz
