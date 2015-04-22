      subroutine pointgrp(rotmat,rotname)
c **********************************************
c This subroutine defines the rotation matrices for
c all the 32 point groups and names them after 
c J.F. Cornwell (Group Theory??) second edition 
c Appendix D, p 324-325
c 
c *********************************************    
      implicit none
      integer i,j,i1,is
      double precision ROTMAT(64,3,3)
      double precision RTHREE,HALF 
      character*10 ROTNAME(64)

      RTHREE = SQRT(3.d0)/2.d0
      HALF = 0.5d0
c set to zero
      do i1=1,64  
          do i=1,3
             do j=1,3
                ROTMAT(i1,i,j) = 0.d0
             end do
          end do
      end do
c
      ROTMAT(1,1,1) =  1.d0
      ROTMAT(1,2,2) =  1.d0
      ROTMAT(1,3,3) =  1.d0
      ROTNAME(1) = 'E'
c
      ROTMAT(2,1,2) =  1.d0
      ROTMAT(2,2,3) = -1.d0
      ROTMAT(2,3,1) = -1.d0
      ROTNAME(2) = 'C3alfa'          
c
      ROTMAT(3,1,2) = -1.d0
      ROTMAT(3,2,3) = -1.d0
      ROTMAT(3,3,1) =  1.d0
      ROTNAME(3) = 'C3beta '
c
      ROTMAT(4,1,2) = -1.d0
      ROTMAT(4,2,3) =  1.d0
      ROTMAT(4,3,1) = -1.d0
      ROTNAME(4) = 'C3gamma'
c
      ROTMAT(5,1,2) = 1.d0
      ROTMAT(5,2,3) = 1.d0
      ROTMAT(5,3,1) = 1.d0
      ROTNAME(5) = 'C3delta '
c
      ROTMAT(6,1,3) = -1.d0
      ROTMAT(6,2,1) =  1.d0
      ROTMAT(6,3,2) = -1.d0
      ROTNAME(6) = 'C3alfa-1'
c
      ROTMAT(7,1,3) =  1.d0
      ROTMAT(7,2,1) = -1.d0
      ROTMAT(7,3,2) = -1.d0
      ROTNAME(7) = 'C3beta-1 '
c
      ROTMAT(8,1,3) = -1.d0
      ROTMAT(8,2,1) = -1.d0
      ROTMAT(8,3,2) =  1.d0
      ROTNAME(8) = 'C3gamma-1'
c
      ROTMAT(9,1,3) =  1.d0
      ROTMAT(9,2,1) =  1.d0
      ROTMAT(9,3,2) =  1.d0
      ROTNAME(9) = 'C3delta-1'
c
      ROTMAT(10,1,1) =  1.d0
      ROTMAT(10,2,2) = -1.d0
      ROTMAT(10,3,3) = -1.d0
      ROTNAME(10) = 'C2x'
c           
      ROTMAT(11,1,1) = -1.d0
      ROTMAT(11,2,2) =  1.d0
      ROTMAT(11,3,3) = -1.d0
      ROTNAME(11) = 'C2y'
c
      ROTMAT(12,1,1) = -1.d0
      ROTMAT(12,2,2) = -1.d0
      ROTMAT(12,3,3) =  1.d0
      ROTNAME(12) = 'C2z'
c
      ROTMAT(13,1,1) =  1.d0
      ROTMAT(13,2,3) =  1.d0
      ROTMAT(13,3,2) = -1.d0
      ROTNAME(13) = 'C4x'
c           
      ROTMAT(14,1,3) = -1.d0
      ROTMAT(14,2,2) =  1.d0
      ROTMAT(14,3,1) =  1.d0
      ROTNAME(14) = 'C4y '
c
      ROTMAT(15,1,2) =  1.d0
      ROTMAT(15,2,1) = -1.d0
      ROTMAT(15,3,3) =  1.d0
      ROTNAME(15) = 'C4z'
c           
      ROTMAT(16,1,1) =  1.d0
      ROTMAT(16,2,3) = -1.d0
      ROTMAT(16,3,2) =  1.d0
      ROTNAME(16) = 'C4x-1 '
c
      ROTMAT(17,1,3) =  1.d0
      ROTMAT(17,2,2) =  1.d0
      ROTMAT(17,3,1) = -1.d0
      ROTNAME(17) = 'C4y-1'
c
      ROTMAT(18,1,2) = -1.d0
      ROTMAT(18,2,1) =  1.d0
      ROTMAT(18,3,3) =  1.d0
      ROTNAME(18) = 'C4z-1'
c           
      ROTMAT(19,1,2) =  1.d0
      ROTMAT(19,2,1) =  1.d0
      ROTMAT(19,3,3) = -1.d0
      ROTNAME(19) = 'C2a'
c
      ROTMAT(20,1,2) = -1.d0
      ROTMAT(20,2,1) = -1.d0
      ROTMAT(20,3,3) = -1.d0
      ROTNAME(20) = 'C2b'
c
      ROTMAT(21,1,3) =  1.d0
      ROTMAT(21,2,2) = -1.d0
      ROTMAT(21,3,1) =  1.d0
      ROTNAME(21) = 'C2c'
c
      ROTMAT(22,1,3) = -1.d0
      ROTMAT(22,2,2) = -1.d0
      ROTMAT(22,3,1) = -1.d0
      ROTNAME(22) = 'C2d'
c
      ROTMAT(23,1,1) = -1.d0
      ROTMAT(23,2,3) =  1.d0
      ROTMAT(23,3,2) =  1.d0
      ROTNAME(23) = 'C2e'
c
      ROTMAT(24,1,1) = -1.d0
      ROTMAT(24,2,3) = -1.d0
      ROTMAT(24,3,2) = -1.d0
      ROTNAME(24) = 'C2f'
      do i1=1,24
         do i=1,3
            do j=1,3 
              ROTMAT(i1+24,i,j) = -ROTMAT(i1,i,j)
            end do
         end do
      ROTNAME(i1+24) = 'I'//ROTNAME(i1)
      end do
c
c      
c*********************************************
c Trigonal and hexagonal groups
c*********************************************
c
      ROTMAT(49,1,1) = -HALF
      ROTMAT(49,1,2) =  RTHREE
      ROTMAT(49,2,1) = -RTHREE
      ROTMAT(49,2,2) = -HALF
      ROTMAT(49,3,3) =  1.d0
      ROTNAME(49) = 'C3z'  
c
      ROTMAT(50,1,1) = -HALF
      ROTMAT(50,1,2) = -RTHREE
      ROTMAT(50,2,1) =  RTHREE
      ROTMAT(50,2,2) = -HALF
      ROTMAT(50,3,3) =  1.d0
      ROTNAME(50) = 'C3z-1'
c
      ROTMAT(51,1,1) =  HALF
      ROTMAT(51,1,2) =  RTHREE
      ROTMAT(51,2,1) = -RTHREE
      ROTMAT(51,2,2) =  HALF
      ROTMAT(51,3,3) =  1.d0
      ROTNAME(51) = 'C6z'
c
      ROTMAT(52,1,1) =  HALF
      ROTMAT(52,1,2) = -RTHREE
      ROTMAT(52,2,1) =  RTHREE
      ROTMAT(52,2,2) =  HALF
      ROTMAT(52,3,3) =  1.d0
      ROTNAME(52) = 'C6z-1'
c
      ROTMAT(53,1,1) = -HALF
      ROTMAT(53,1,2) =  RTHREE
      ROTMAT(53,2,1) =  RTHREE
      ROTMAT(53,2,2) =  HALF
      ROTMAT(53,3,3) = -1.d0
      ROTNAME(53) = 'C2A'    
c
      ROTMAT(54,1,1) = -HALF
      ROTMAT(54,1,2) = -RTHREE
      ROTMAT(54,2,1) = -RTHREE
      ROTMAT(54,2,2) =  HALF
      ROTMAT(54,3,3) = -1.d0
      ROTNAME(54) = 'C2B'
c
      ROTMAT(55,1,1) =  HALF
      ROTMAT(55,1,2) = -RTHREE
      ROTMAT(55,2,1) = -RTHREE
      ROTMAT(55,2,2) = -HALF
      ROTMAT(55,3,3) = -1.d0
      ROTNAME(55) = 'C2C'
c
      ROTMAT(56,1,1) =  HALF
      ROTMAT(56,1,2) =  RTHREE
      ROTMAT(56,2,1) =  RTHREE
      ROTMAT(56,2,2) = -HALF
      ROTMAT(56,3,3) = -1.d0
      ROTNAME(56) = 'C2D'
      do is=1,8
          do i=1,3
             do j=1,3     
                ROTMAT(56+is,i,j) = -ROTMAT(48+is,i,j)
             end do
          end do
          ROTNAME(56+is) = 'I'//ROTNAME(48+is) 
      end do
      
cccccccccccccccccccccccccccccccccccccccccccccccccc
      END
