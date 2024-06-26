!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

module mod_dsort
  
  private
  public :: dsort

contains

  !-------------------------------------------------------------------------------
  !> Summary: Sort double precision array returning sorted index array
  !> Author: P. Zahn
  !> Date: April 96
  !> Category: KKRhost, numerical-tools
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> W   is the original array returned unchanged
  !> IND is an array that holds the new positions
  !> max number of ellements to be sorted
  !> pos the position where the first element is found
  !-------------------------------------------------------------------------------
  subroutine dsort(w, ind, max, pos)

    use :: mod_datatypes, only: dp
    implicit none

    integer :: max, pos
    real (kind=dp) :: w(*)
    integer :: ind(*)

    integer :: i, ii, j, jj, k
    real (kind=dp) :: bound, diff
    data bound/1.0e-12_dp/
    ! ------------------------------------------------------------------------
    do i = 1, max
      ind(i) = i
    end do

    j = max
    j = 1
    do while (j<max/3)
      j = 3*j + 1
    end do

    do while (j>1)
      j = j/3
      jj = 1
      do while (jj==1)
        jj = 0
        do k = 1, max - j
          diff = abs(w(ind(k))-w(ind(k+j)))
          if (w(ind(k))>w(ind(k+j)) .and. diff>bound) then
            ii = ind(k)
            ind(k) = ind(k+j)
            ind(k+j) = ii
            jj = 1
          end if
        end do                     ! K=1,MAX-J
      end do
    end do                         ! WHILE (JJ.EQ.1)

    do i = 1, max
      if (ind(i)==1) pos = i
    end do

    return
  end subroutine dsort

end module mod_dsort
