!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Performs the rotation of the matrix `T1` using the rotation-matrix `ROT`, set up by `CALCROTMAT()`
!> Author: 
!> Performs the rotation of the matrix `T1` using the rotation-matrix `ROT`, set 
!> up by `CALCROTMAT()`
!> \(T2 = ROT T1 ROT_+\) `IF  MODE = 'L->G'`
!>
!> \(T2 = ROT_+ T1 ROT\) `IF  MODE = 'G->L'`
!> see: E.M. ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM
!------------------------------------------------------------------------------------
module mod_rotate
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: Performs the rotation of the matrix `T1` using the rotation-matrix `ROT`, set up by `CALCROTMAT()`
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> Performs the rotation of the matrix `T1` using the rotation-matrix `ROT`, set 
  !> up by `CALCROTMAT()`
  !> \(T2 = ROT T1 ROT_+\) `IF  MODE = 'L->G'`
  !>
  !> \(T2 = ROT_+ T1 ROT\) `IF  MODE = 'G->L'`
  !> see: E.M. ROSE ELEMENTARY THEORY OF ANGULAR MOMENTUM
  !-------------------------------------------------------------------------------
  subroutine rotate(t1, mode, t2, n, rot, nkmmax)
    use :: mod_constants, only: cone, czero

    implicit none

    ! Dummy arguments
    character (len=4) :: mode
    integer :: n, nkmmax
    complex (kind=dp), dimension(nkmmax, nkmmax) :: rot, t1, t2

    ! Local variables
    character (len=1) :: fl1, fl2
    complex (kind=dp), dimension(nkmmax, nkmmax) :: w1


    if (mode=='L->G') then
      fl1 = 'N'
      fl2 = 'C'
    else if (mode=='G->L') then
      fl1 = 'C'
      fl2 = 'N'
    else
      write (*, *) ' MODE = ', mode
      stop 'in <ROTATE>  MODE not allowed'
    end if

    call zgemm(fl1, 'N', n, n, n, cone, rot, nkmmax, t1, nkmmax, czero, w1, nkmmax)

    call zgemm('N', fl2, n, n, n, cone, w1, nkmmax, rot, nkmmax, czero, t2, nkmmax)

  end subroutine rotate

end module mod_rotate
