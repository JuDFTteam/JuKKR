!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: This subroutine defines the rotation matrices for all the 32 point groups
!> Author: 
!> This subroutine defines the rotation matrices for all the 32 point groups and names them after
!> J.F. Cornwell (Group Theory??) second edition Appendix D, p 324-325
!------------------------------------------------------------------------------------
module mod_pointgrp
  use :: mod_datatypes, only: dp
  private :: dp

contains

  !-------------------------------------------------------------------------------
  !> Summary: This subroutine defines the rotation matrices for all the 32 point groups
  !> Author: 
  !> Category: numerical-tools, KKRhost
  !> Deprecated: False 
  !> This subroutine defines the rotation matrices for all the 32 point groups and names them after
  !> J.F. Cornwell (Group Theory??) second edition Appendix D, p 324-325
  !-------------------------------------------------------------------------------
  subroutine pointgrp(rotmat, rotname, writesymfile)

    implicit none

    ! optional input, can trigger writeout of symmetry matrices to sym.out file
    logical, optional :: writesymfile
    integer, parameter :: iou=3514 ! file unit for sym.out file

    ! .. Output variables
    real (kind=dp), dimension(64, 3, 3), intent(out) :: rotmat !! Rotation matrices
    character (len=10), dimension(64), intent(out) :: rotname  !! Name of the space group
    ! .. Local variables
    integer :: i, j, i1, is
    real (kind=dp) :: rthree, half

    rthree = sqrt(3.e0_dp)/2.e0_dp
    half = 0.5e0_dp
    ! set to zero
    do i1 = 1, 64
      do i = 1, 3
        do j = 1, 3
          rotmat(i1, i, j) = 0.e0_dp
        end do
      end do
    end do

    rotmat(1, 1, 1) = 1.e0_dp
    rotmat(1, 2, 2) = 1.e0_dp
    rotmat(1, 3, 3) = 1.e0_dp
    rotname(1) = 'E'

    rotmat(2, 1, 2) = 1.e0_dp
    rotmat(2, 2, 3) = -1.e0_dp
    rotmat(2, 3, 1) = -1.e0_dp
    rotname(2) = 'C3alfa'

    rotmat(3, 1, 2) = -1.e0_dp
    rotmat(3, 2, 3) = -1.e0_dp
    rotmat(3, 3, 1) = 1.e0_dp
    rotname(3) = 'C3beta '

    rotmat(4, 1, 2) = -1.e0_dp
    rotmat(4, 2, 3) = 1.e0_dp
    rotmat(4, 3, 1) = -1.e0_dp
    rotname(4) = 'C3gamma'

    rotmat(5, 1, 2) = 1.e0_dp
    rotmat(5, 2, 3) = 1.e0_dp
    rotmat(5, 3, 1) = 1.e0_dp
    rotname(5) = 'C3delta '

    rotmat(6, 1, 3) = -1.e0_dp
    rotmat(6, 2, 1) = 1.e0_dp
    rotmat(6, 3, 2) = -1.e0_dp
    rotname(6) = 'C3alfa-1'

    rotmat(7, 1, 3) = 1.e0_dp
    rotmat(7, 2, 1) = -1.e0_dp
    rotmat(7, 3, 2) = -1.e0_dp
    rotname(7) = 'C3beta-1 '

    rotmat(8, 1, 3) = -1.e0_dp
    rotmat(8, 2, 1) = -1.e0_dp
    rotmat(8, 3, 2) = 1.e0_dp
    rotname(8) = 'C3gamma-1'

    rotmat(9, 1, 3) = 1.e0_dp
    rotmat(9, 2, 1) = 1.e0_dp
    rotmat(9, 3, 2) = 1.e0_dp
    rotname(9) = 'C3delta-1'

    rotmat(10, 1, 1) = 1.e0_dp
    rotmat(10, 2, 2) = -1.e0_dp
    rotmat(10, 3, 3) = -1.e0_dp
    rotname(10) = 'C2x'

    rotmat(11, 1, 1) = -1.e0_dp
    rotmat(11, 2, 2) = 1.e0_dp
    rotmat(11, 3, 3) = -1.e0_dp
    rotname(11) = 'C2y'

    rotmat(12, 1, 1) = -1.e0_dp
    rotmat(12, 2, 2) = -1.e0_dp
    rotmat(12, 3, 3) = 1.e0_dp
    rotname(12) = 'C2z'

    rotmat(13, 1, 1) = 1.e0_dp
    rotmat(13, 2, 3) = 1.e0_dp
    rotmat(13, 3, 2) = -1.e0_dp
    rotname(13) = 'C4x'

    rotmat(14, 1, 3) = -1.e0_dp
    rotmat(14, 2, 2) = 1.e0_dp
    rotmat(14, 3, 1) = 1.e0_dp
    rotname(14) = 'C4y '

    rotmat(15, 1, 2) = 1.e0_dp
    rotmat(15, 2, 1) = -1.e0_dp
    rotmat(15, 3, 3) = 1.e0_dp
    rotname(15) = 'C4z'

    rotmat(16, 1, 1) = 1.e0_dp
    rotmat(16, 2, 3) = -1.e0_dp
    rotmat(16, 3, 2) = 1.e0_dp
    rotname(16) = 'C4x-1 '

    rotmat(17, 1, 3) = 1.e0_dp
    rotmat(17, 2, 2) = 1.e0_dp
    rotmat(17, 3, 1) = -1.e0_dp
    rotname(17) = 'C4y-1'

    rotmat(18, 1, 2) = -1.e0_dp
    rotmat(18, 2, 1) = 1.e0_dp
    rotmat(18, 3, 3) = 1.e0_dp
    rotname(18) = 'C4z-1'

    rotmat(19, 1, 2) = 1.e0_dp
    rotmat(19, 2, 1) = 1.e0_dp
    rotmat(19, 3, 3) = -1.e0_dp
    rotname(19) = 'C2a'

    rotmat(20, 1, 2) = -1.e0_dp
    rotmat(20, 2, 1) = -1.e0_dp
    rotmat(20, 3, 3) = -1.e0_dp
    rotname(20) = 'C2b'

    rotmat(21, 1, 3) = 1.e0_dp
    rotmat(21, 2, 2) = -1.e0_dp
    rotmat(21, 3, 1) = 1.e0_dp
    rotname(21) = 'C2c'

    rotmat(22, 1, 3) = -1.e0_dp
    rotmat(22, 2, 2) = -1.e0_dp
    rotmat(22, 3, 1) = -1.e0_dp
    rotname(22) = 'C2d'

    rotmat(23, 1, 1) = -1.e0_dp
    rotmat(23, 2, 3) = 1.e0_dp
    rotmat(23, 3, 2) = 1.e0_dp
    rotname(23) = 'C2e'

    rotmat(24, 1, 1) = -1.e0_dp
    rotmat(24, 2, 3) = -1.e0_dp
    rotmat(24, 3, 2) = -1.e0_dp
    rotname(24) = 'C2f'
    do i1 = 1, 24
      do i = 1, 3
        do j = 1, 3
          rotmat(i1+24, i, j) = -rotmat(i1, i, j)
        end do
      end do
      rotname(i1+24) = 'I' // rotname(i1)
    end do


    ! *********************************************
    ! Trigonal and hexagonal groups
    ! *********************************************

    rotmat(49, 1, 1) = -half
    rotmat(49, 1, 2) = rthree
    rotmat(49, 2, 1) = -rthree
    rotmat(49, 2, 2) = -half
    rotmat(49, 3, 3) = 1.e0_dp
    rotname(49) = 'C3z'

    rotmat(50, 1, 1) = -half
    rotmat(50, 1, 2) = -rthree
    rotmat(50, 2, 1) = rthree
    rotmat(50, 2, 2) = -half
    rotmat(50, 3, 3) = 1.e0_dp
    rotname(50) = 'C3z-1'

    rotmat(51, 1, 1) = half
    rotmat(51, 1, 2) = rthree
    rotmat(51, 2, 1) = -rthree
    rotmat(51, 2, 2) = half
    rotmat(51, 3, 3) = 1.e0_dp
    rotname(51) = 'C6z'

    rotmat(52, 1, 1) = half
    rotmat(52, 1, 2) = -rthree
    rotmat(52, 2, 1) = rthree
    rotmat(52, 2, 2) = half
    rotmat(52, 3, 3) = 1.e0_dp
    rotname(52) = 'C6z-1'

    rotmat(53, 1, 1) = -half
    rotmat(53, 1, 2) = rthree
    rotmat(53, 2, 1) = rthree
    rotmat(53, 2, 2) = half
    rotmat(53, 3, 3) = -1.e0_dp
    rotname(53) = 'C2A'

    rotmat(54, 1, 1) = -half
    rotmat(54, 1, 2) = -rthree
    rotmat(54, 2, 1) = -rthree
    rotmat(54, 2, 2) = half
    rotmat(54, 3, 3) = -1.e0_dp
    rotname(54) = 'C2B'

    rotmat(55, 1, 1) = half
    rotmat(55, 1, 2) = -rthree
    rotmat(55, 2, 1) = -rthree
    rotmat(55, 2, 2) = -half
    rotmat(55, 3, 3) = -1.e0_dp
    rotname(55) = 'C2C'

    rotmat(56, 1, 1) = half
    rotmat(56, 1, 2) = rthree
    rotmat(56, 2, 1) = rthree
    rotmat(56, 2, 2) = -half
    rotmat(56, 3, 3) = -1.e0_dp
    rotname(56) = 'C2D'
    do is = 1, 8
      do i = 1, 3
        do j = 1, 3
          rotmat(56+is, i, j) = -rotmat(48+is, i, j)
        end do
      end do
      rotname(56+is) = 'I' // rotname(48+is)
    end do

    ! for FS code we should be able to write out the symmetries to a file
    if (present(writesymfile)) then
      if(writesymfile)then
        open(unit=iou,file='sym.out',form='formatted',action='write')
        write(iou,'(I0)') 64
        do is=1,64
          write(iou,'(A)') ROTNAME(is)
          write(iou,'(3ES25.16)') ROTMAT(is,:,:)
        end do
        close(iou)
        writesymfile = .false.
      end if ! writesymfile
    end if ! present(writesymfile)

    ! ccccccccccccccccccccccccccccccccccccccccccccccccc
  end subroutine pointgrp

end module mod_pointgrp
