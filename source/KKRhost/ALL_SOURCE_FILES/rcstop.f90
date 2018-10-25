!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Subroutine to print where the program stops in case of an error
!> Author:
!> Subroutine to print where the program stops in case of an error
!------------------------------------------------------------------------------------
module mod_rcstop

contains

  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to print where the program stops in case of an error
  !> Author: 
  !> Category: sanity-check, input-output, KKRhost 
  !> Deprecated: False 
  !> Subroutine to print where the program stops in case of an error
  !-------------------------------------------------------------------------------
  subroutine rcstop(c)

    implicit none

    character (len=8), intent (in) :: c !! Input string

    ! ..
    print *, 'ERROR: STOP AT POSITION ', c
    stop

  end subroutine rcstop

end module mod_rcstop
