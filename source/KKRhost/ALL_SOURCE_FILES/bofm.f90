!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_bofm

  private
  public :: bofm

contains

  !-------------------------------------------------------------------------------
  !> Summary: Take block out of bigger matrix
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> Assumes double complex matrices
  !-------------------------------------------------------------------------------
  subroutine bofm(pl1, pl2, block, nsize, gin, almd)
    use :: mod_datatypes, only: dp
    implicit none
    ! .. Scalar Arguments ..
    integer, intent(in) :: almd, nsize, pl1, pl2
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp), dimension(almd,almd), intent(in) :: gin
    complex (kind=dp), dimension(nsize,nsize), intent(out) :: block
    ! ..
    ! .. Local Scalars ..
    integer :: i1, i1s, i2, i2s
    ! ..
    i1s = (pl1-1)*nsize
    i2s = (pl2-1)*nsize
    do i1 = 1, nsize
      do i2 = 1, nsize
        block(i1, i2) = gin(i1s+i1, i2s+i2)
      end do
    end do

    return

  end subroutine bofm

end module mod_bofm
