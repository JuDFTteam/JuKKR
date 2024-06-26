!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary:
!> Author: 
!------------------------------------------------------------------------------------
module mod_sname
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: 
  !> Author: 
  !> Category: deprecated, KKRhost
  !> Deprecated: True
  !> 
  !-------------------------------------------------------------------------------
  subroutine sname(name, new, band)

    use :: mod_length

    implicit none
    ! .. scalar arguments
    integer, intent(in) :: band
    character (len=40), intent(inout) :: new
    character (len=40), intent(in) :: name

    ! .. locals
    integer :: i, l, lo
    character (len=1) :: poi
    character (len=1), dimension(50) :: ch
    character (len=10) :: s
    ! ------------------------------------------------------------------------
    poi = '.'
    if (band<0) then
      lo = nint(log(real(-band,kind=dp))/log(10.0e0_dp)) + 1
    else if (band==0) then
      lo = 0
    else
      lo = nint(log(real(band,kind=dp))/log(10.0e0_dp))
    end if

    ! write(6,*) 'LO ',lo

    read (name, fmt='(255a1)')(ch(i), i=1, 40)
    l = length(ch, 40)
    ! write(6,*) 'L  ',l

    ! write(6,*) 'CH ',(CH(I),I=1,25)

    write (s, fmt='(I10)') band
    ! write(6,*) 'S  ',s

    read (s, fmt='(255A1)')(ch(i), i=l+1, l+10)
    ! write(6,*) 'CH ',(CH(I),I=L+1,L+10)

    write (new, fmt='(255A1)')(ch(i), i=1, l), poi, (ch(i), i=l+10-lo, l+10)

    return
  end subroutine sname

end module mod_sname
