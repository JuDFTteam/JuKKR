      subroutine pointgrp(rotmat, rotname)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after 
! J.F. Cornwell (Group Theory??) second edition 
! Appendix D, p 324-325
! 
! *********************************************    
      implicit none
      double precision, intent(out) :: rotmat(64,3,3)
      character(len=*), intent(out) :: rotname(64)
      integer :: is
      double precision, parameter :: RTHREE = SQRT(3.d0)/2.d0, HALF=0.5d0, ONE=1.d0

      rotmat(:,:,:) =  0.d0 ! set all matrices to zero
!
      rotmat(1,1,1) =  ONE
      rotmat(1,2,2) =  ONE
      rotmat(1,3,3) =  ONE
      rotname(1) = 'E'
!
      rotmat(2,1,2) =  ONE
      rotmat(2,2,3) = -ONE
      rotmat(2,3,1) = -ONE
      rotname(2) = 'C3alfa'          
!
      rotmat(3,1,2) = -ONE
      rotmat(3,2,3) = -ONE
      rotmat(3,3,1) =  ONE
      rotname(3) = 'C3beta '
!
      rotmat(4,1,2) = -ONE
      rotmat(4,2,3) =  ONE
      rotmat(4,3,1) = -ONE
      rotname(4) = 'C3gamma'
!
      rotmat(5,1,2) = ONE
      rotmat(5,2,3) = ONE
      rotmat(5,3,1) = ONE
      rotname(5) = 'C3delta '
!
      rotmat(6,1,3) = -ONE
      rotmat(6,2,1) =  ONE
      rotmat(6,3,2) = -ONE
      rotname(6) = 'C3alfa-1'
!
      rotmat(7,1,3) =  ONE
      rotmat(7,2,1) = -ONE
      rotmat(7,3,2) = -ONE
      rotname(7) = 'C3beta-1 '
!
      rotmat(8,1,3) = -ONE
      rotmat(8,2,1) = -ONE
      rotmat(8,3,2) =  ONE
      rotname(8) = 'C3gamma-1'
!
      rotmat(9,1,3) =  ONE
      rotmat(9,2,1) =  ONE
      rotmat(9,3,2) =  ONE
      rotname(9) = 'C3delta-1'
!
      rotmat(10,1,1) =  ONE
      rotmat(10,2,2) = -ONE
      rotmat(10,3,3) = -ONE
      rotname(10) = 'C2x'
!           
      rotmat(11,1,1) = -ONE
      rotmat(11,2,2) =  ONE
      rotmat(11,3,3) = -ONE
      rotname(11) = 'C2y'
!
      rotmat(12,1,1) = -ONE
      rotmat(12,2,2) = -ONE
      rotmat(12,3,3) =  ONE
      rotname(12) = 'C2z'
!
      rotmat(13,1,1) =  ONE
      rotmat(13,2,3) =  ONE
      rotmat(13,3,2) = -ONE
      rotname(13) = 'C4x'
!           
      rotmat(14,1,3) = -ONE
      rotmat(14,2,2) =  ONE
      rotmat(14,3,1) =  ONE
      rotname(14) = 'C4y '
!
      rotmat(15,1,2) =  ONE
      rotmat(15,2,1) = -ONE
      rotmat(15,3,3) =  ONE
      rotname(15) = 'C4z'
!           
      rotmat(16,1,1) =  ONE
      rotmat(16,2,3) = -ONE
      rotmat(16,3,2) =  ONE
      rotname(16) = 'C4x-1 '
!
      rotmat(17,1,3) =  ONE
      rotmat(17,2,2) =  ONE
      rotmat(17,3,1) = -ONE
      rotname(17) = 'C4y-1'
!
      rotmat(18,1,2) = -ONE
      rotmat(18,2,1) =  ONE
      rotmat(18,3,3) =  ONE
      rotname(18) = 'C4z-1'
!           
      rotmat(19,1,2) =  ONE
      rotmat(19,2,1) =  ONE
      rotmat(19,3,3) = -ONE
      rotname(19) = 'C2a'
!
      rotmat(20,1,2) = -ONE
      rotmat(20,2,1) = -ONE
      rotmat(20,3,3) = -ONE
      rotname(20) = 'C2b'
!
      rotmat(21,1,3) =  ONE
      rotmat(21,2,2) = -ONE
      rotmat(21,3,1) =  ONE
      rotname(21) = 'C2c'
!
      rotmat(22,1,3) = -ONE
      rotmat(22,2,2) = -ONE
      rotmat(22,3,1) = -ONE
      rotname(22) = 'C2d'
!
      rotmat(23,1,1) = -ONE
      rotmat(23,2,3) =  ONE
      rotmat(23,3,2) =  ONE
      rotname(23) = 'C2e'
!
      rotmat(24,1,1) = -ONE
      rotmat(24,2,3) = -ONE
      rotmat(24,3,2) = -ONE
      rotname(24) = 'C2f'
      
      do is = 1, 24
        rotmat(is+24,1:3,1:3) = -rotmat(is,1:3,1:3)
        rotname(is+24) = 'I'//rotname(is)
      enddo ! i1      
!
!      
!*********************************************
! Trigonal and hexagonal groups
!*********************************************
!
      rotmat(49,1,1) = -HALF
      rotmat(49,1,2) =  RTHREE
      rotmat(49,2,1) = -RTHREE
      rotmat(49,2,2) = -HALF
      rotmat(49,3,3) =  ONE
      rotname(49) = 'C3z'  
!
      rotmat(50,1,1) = -HALF
      rotmat(50,1,2) = -RTHREE
      rotmat(50,2,1) =  RTHREE
      rotmat(50,2,2) = -HALF
      rotmat(50,3,3) =  ONE
      rotname(50) = 'C3z-1'
!
      rotmat(51,1,1) =  HALF
      rotmat(51,1,2) =  RTHREE
      rotmat(51,2,1) = -RTHREE
      rotmat(51,2,2) =  HALF
      rotmat(51,3,3) =  ONE
      rotname(51) = 'C6z'
!
      rotmat(52,1,1) =  HALF
      rotmat(52,1,2) = -RTHREE
      rotmat(52,2,1) =  RTHREE
      rotmat(52,2,2) =  HALF
      rotmat(52,3,3) =  ONE
      rotname(52) = 'C6z-1'
!
      rotmat(53,1,1) = -HALF
      rotmat(53,1,2) =  RTHREE
      rotmat(53,2,1) =  RTHREE
      rotmat(53,2,2) =  HALF
      rotmat(53,3,3) = -ONE
      rotname(53) = 'C2A'    
!
      rotmat(54,1,1) = -HALF
      rotmat(54,1,2) = -RTHREE
      rotmat(54,2,1) = -RTHREE
      rotmat(54,2,2) =  HALF
      rotmat(54,3,3) = -ONE
      rotname(54) = 'C2B'
!
      rotmat(55,1,1) =  HALF
      rotmat(55,1,2) = -RTHREE
      rotmat(55,2,1) = -RTHREE
      rotmat(55,2,2) = -HALF
      rotmat(55,3,3) = -ONE
      rotname(55) = 'C2C'
!
      rotmat(56,1,1) =  HALF
      rotmat(56,1,2) =  RTHREE
      rotmat(56,2,1) =  RTHREE
      rotmat(56,2,2) = -HALF
      rotmat(56,3,3) = -ONE
      rotname(56) = 'C2D'
      
      do is = 1, 8
        rotmat(56+is,1:3,1:3) = -rotmat(48+is,1:3,1:3)
        rotname(56+is) = 'I'//rotname(48+is) 
      enddo ! is
      
!ccccccccccccccccccccccccccccccccccccccccccccccccc
      endsubroutine pointgrp
