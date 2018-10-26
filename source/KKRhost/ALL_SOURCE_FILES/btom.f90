!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_btom

  private
  public :: btom

contains

  !-------------------------------------------------------------------------------
  !> Summary: Copy or substract block from matrix
  !> Author: 
  !> Category: KKRhost, 
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> This subroutine copies or subtracts a block to a matrix
  !-------------------------------------------------------------------------------
  subroutine btom(pl1, pl2, block_mat, nsize, gin, almd, lsub)
    use :: mod_datatypes, only: dp
    implicit none
    ! .. Scalar Arguments ..
    integer :: almd, nsize, pl1, pl2
    logical :: lsub
    ! ..
    ! .. Array Arguments ..
    complex (kind=dp), dimension(almd, almd) :: gin
    complex (kind=dp), dimension(nsize,nsize) :: block_mat
    ! ..
    ! .. Local Scalars ..
    integer :: i1, i1s, i2, i2s
    ! ..
    i1s = (pl1-1)*nsize
    i2s = (pl2-1)*nsize
    if (lsub) then
      do i1 = 1, nsize
        do i2 = 1, nsize
          gin(i1s+i1, i2s+i2) = gin(i1s+i1, i2s+i2) - block_mat(i1, i2)
        end do
      end do
    else
      do i1 = 1, nsize
        do i2 = 1, nsize
          gin(i1s+i1, i2s+i2) = block_mat(i1, i2)
        end do
      end do
    end if

    return

  end subroutine btom

end module mod_btom
