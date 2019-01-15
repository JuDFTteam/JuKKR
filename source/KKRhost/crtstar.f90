!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_crtstar
  
  private
  public :: crtstar

contains

  !-------------------------------------------------------------------------------
  !> Summary: Apply space-group symmetries to real-space vector
  !> Author: 
  !> Date: 20.07.96
  !> Category: KKRhost, k-points
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> THE SYMMETRY OPERATIONS OF THE SYMMETRY GROUP ARE APPLIED TO THE
  !> INPUT VECTOR RATOM
  !-------------------------------------------------------------------------------
  subroutine crtstar(ratom, nshell, nd, irot, isymindex, rrot)

    use :: mod_datatypes, only: dp
    implicit none

    integer :: irot, nshell
    real (kind=dp) :: nd(64, 3, *), ratom(3, *), rrot(48, 3, *)
    integer :: isymindex(*)

    integer :: i, id, ns, k, j, isym


    do ns = 1, nshell
      do id = 1, irot
        isym = isymindex(id)
        do i = 1, 3
          rrot(id, i, ns) = nd(isym, i, 1)*ratom(1, ns) + nd(isym, i, 2)*ratom(2, ns) + nd(isym, i, 3)*ratom(3, ns)
        end do
      end do
    end do

    return

  end subroutine crtstar

end module mod_crtstar
