program reduce_sym
! This module provides routines regardinf the symmetry of the lattice and the irreducible brillouin zone.

  implicit none

  call find_subgroup

contains

  subroutine find_subgroup

      implicit none

      double precision  :: rotmat(64,3,3)
      character(len=10) :: rotname(64)

      integer :: nsymall, nsymminus, nsymsub
      integer, allocatable :: isymall(:), isymminus(:), isymsub(:), isymsub_clean(:)

      integer :: ii, isy1, isy2, isy3, i1, i2, i3, ierr, iofile
      double precision :: dist, grpelm(3,3), trace

      call pointgrp(rotmat,rotname)

      nsymall=16
      nsymminus=1

      allocate(isymall(nsymall), isymminus(nsymminus), isymsub(nsymall), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating isymall etc.'

!     isymall = (/ (ii, ii=1,48) /)
      isymall = (/ 1, 10, 11, 12, 15, 18, 19, 20, 25, 34, 35, 36, 39, 42, 43, 44 /)

      isymminus = (/ 44 /)

      isymsub = 0!0=unset, 1=in subgroup, -1=excluded from subgroup

      do isy1=1,nsymall
        if(isymsub(isy1)/=0) cycle

        isymsub(isy1) = 1

        !multipliziere dieses Element der Gruppe mit den verbotenen Elementen
        do isy2=1,nsymminus

          !bilde multiplikation
          grpelm=0d0
          do i1=1,3
           do i2=1,3
            do i3=1,3
              grpelm(i2,i1) = grpelm(i2,i1) + rotmat(isymall(isy1),i3,i1)*rotmat(isymminus(isy2),i2,i3)
            end do!i3
           end do!i2
          end do!i1

          !finde correspondierende symmetrieoperation
          inner: do isy3=1,nsymall
            dist = sum(abs( grpelm - rotmat(isy3,:,:) ))
            if(dist<1d-8 .and. isy3/=1)then
              isymsub(isy3) = -1
              exit inner
            elseif(dist<1d-8 .and. isy3==1)then
               isymsub(isy1)= -1
            end if
          end do inner!isy3

        end do!isy2
      end do!isy1

      write(*,'(8I8)') isymsub

      nsymsub = count(isymsub>0)
      allocate(isymsub_clean(nsymsub), STAT=ierr)
      if(ierr/=0) stop 'Problem allocating isymsub_clean etc.'

      isy1=0
      do isy2=1,nsymall
        if(isymsub(isy2)>0)then
          isy1=isy1+1
          isymsub_clean(isy1) = isymall(isy2)
        end if
      end do

      write(*,*) 'nsymsub=', nsymsub
      write(*,1050) (rotname(isymsub_clean(isy1)), isy1=1,nsymsub)

!     write(*,*) 'Calculate traces:'
!     do isy1=1,nsymall
!       trace=0d0
!       do i1=1,3
!         trace = trace+rotmat(isymall(isy1),i1,i1)
!       end do
!       write(*,'(A10,F8.2,L1)') rotname(isymall(isy1)), trace, isymsub(isy1)==1
!     end do

      iofile=135
      write(iofile,*) 'Print the matrices'
      do isy1=1,nsymall
        write(iofile,'(A,I0,A,A10)') 'isym=', isy1, '  name= ', rotname(isy1)
        do i1=1,3
          write(iofile,'(4X,3F4.0)') rotmat(isy1,:,i1)
        end do
        write(iofile,*) ''
      end do


 1050 FORMAT(5(A10,2X))
  end subroutine find_subgroup


  subroutine pointgrp(rotmat,rotname)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after 
! J.F. Cornwell (Group Theory??) second edition 
! Appendix D, p 324-325
! 
! *********************************************    
      implicit none

      double precision,  intent(out) :: ROTMAT(64,3,3)
      character(len=10), intent(out) :: ROTNAME(64)

      !Locals
      integer i,j,i1,is
      double precision RTHREE,HALF 

      RTHREE = sqrt(3.d0)/2.d0
      HALF = 0.5d0
! set to zero
      do i1=1,64  
          do i=1,3
             do j=1,3
                ROTMAT(i1,i,j) = 0.d0
             end do
          end do
      end do

      ROTMAT(1,1,1) =  1.d0
      ROTMAT(1,2,2) =  1.d0
      ROTMAT(1,3,3) =  1.d0
      ROTNAME(1) = 'E'

      ROTMAT(2,1,2) =  1.d0
      ROTMAT(2,2,3) = -1.d0
      ROTMAT(2,3,1) = -1.d0
      ROTNAME(2) = 'C3alfa'          

      ROTMAT(3,1,2) = -1.d0
      ROTMAT(3,2,3) = -1.d0
      ROTMAT(3,3,1) =  1.d0
      ROTNAME(3) = 'C3beta '

      ROTMAT(4,1,2) = -1.d0
      ROTMAT(4,2,3) =  1.d0
      ROTMAT(4,3,1) = -1.d0
      ROTNAME(4) = 'C3gamma'

      ROTMAT(5,1,2) = 1.d0
      ROTMAT(5,2,3) = 1.d0
      ROTMAT(5,3,1) = 1.d0
      ROTNAME(5) = 'C3delta '

      ROTMAT(6,1,3) = -1.d0
      ROTMAT(6,2,1) =  1.d0
      ROTMAT(6,3,2) = -1.d0
      ROTNAME(6) = 'C3alfa-1'

      ROTMAT(7,1,3) =  1.d0
      ROTMAT(7,2,1) = -1.d0
      ROTMAT(7,3,2) = -1.d0
      ROTNAME(7) = 'C3beta-1 '

      ROTMAT(8,1,3) = -1.d0
      ROTMAT(8,2,1) = -1.d0
      ROTMAT(8,3,2) =  1.d0
      ROTNAME(8) = 'C3gamma-1'

      ROTMAT(9,1,3) =  1.d0
      ROTMAT(9,2,1) =  1.d0
      ROTMAT(9,3,2) =  1.d0
      ROTNAME(9) = 'C3delta-1'

      ROTMAT(10,1,1) =  1.d0
      ROTMAT(10,2,2) = -1.d0
      ROTMAT(10,3,3) = -1.d0
      ROTNAME(10) = 'C2x'
           
      ROTMAT(11,1,1) = -1.d0
      ROTMAT(11,2,2) =  1.d0
      ROTMAT(11,3,3) = -1.d0
      ROTNAME(11) = 'C2y'

      ROTMAT(12,1,1) = -1.d0
      ROTMAT(12,2,2) = -1.d0
      ROTMAT(12,3,3) =  1.d0
      ROTNAME(12) = 'C2z'

      ROTMAT(13,1,1) =  1.d0
      ROTMAT(13,2,3) =  1.d0
      ROTMAT(13,3,2) = -1.d0
      ROTNAME(13) = 'C4x'
           
      ROTMAT(14,1,3) = -1.d0
      ROTMAT(14,2,2) =  1.d0
      ROTMAT(14,3,1) =  1.d0
      ROTNAME(14) = 'C4y '

      ROTMAT(15,1,2) =  1.d0
      ROTMAT(15,2,1) = -1.d0
      ROTMAT(15,3,3) =  1.d0
      ROTNAME(15) = 'C4z'
           
      ROTMAT(16,1,1) =  1.d0
      ROTMAT(16,2,3) = -1.d0
      ROTMAT(16,3,2) =  1.d0
      ROTNAME(16) = 'C4x-1 '

      ROTMAT(17,1,3) =  1.d0
      ROTMAT(17,2,2) =  1.d0
      ROTMAT(17,3,1) = -1.d0
      ROTNAME(17) = 'C4y-1'

      ROTMAT(18,1,2) = -1.d0
      ROTMAT(18,2,1) =  1.d0
      ROTMAT(18,3,3) =  1.d0
      ROTNAME(18) = 'C4z-1'
           
      ROTMAT(19,1,2) =  1.d0
      ROTMAT(19,2,1) =  1.d0
      ROTMAT(19,3,3) = -1.d0
      ROTNAME(19) = 'C2a'

      ROTMAT(20,1,2) = -1.d0
      ROTMAT(20,2,1) = -1.d0
      ROTMAT(20,3,3) = -1.d0
      ROTNAME(20) = 'C2b'

      ROTMAT(21,1,3) =  1.d0
      ROTMAT(21,2,2) = -1.d0
      ROTMAT(21,3,1) =  1.d0
      ROTNAME(21) = 'C2c'

      ROTMAT(22,1,3) = -1.d0
      ROTMAT(22,2,2) = -1.d0
      ROTMAT(22,3,1) = -1.d0
      ROTNAME(22) = 'C2d'

      ROTMAT(23,1,1) = -1.d0
      ROTMAT(23,2,3) =  1.d0
      ROTMAT(23,3,2) =  1.d0
      ROTNAME(23) = 'C2e'

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

!*********************************************
! Trigonal and hexagonal groups
!*********************************************

      ROTMAT(49,1,1) = -HALF
      ROTMAT(49,1,2) =  RTHREE
      ROTMAT(49,2,1) = -RTHREE
      ROTMAT(49,2,2) = -HALF
      ROTMAT(49,3,3) =  1.d0
      ROTNAME(49) = 'C3z'  

      ROTMAT(50,1,1) = -HALF
      ROTMAT(50,1,2) = -RTHREE
      ROTMAT(50,2,1) =  RTHREE
      ROTMAT(50,2,2) = -HALF
      ROTMAT(50,3,3) =  1.d0
      ROTNAME(50) = 'C3z-1'

      ROTMAT(51,1,1) =  HALF
      ROTMAT(51,1,2) =  RTHREE
      ROTMAT(51,2,1) = -RTHREE
      ROTMAT(51,2,2) =  HALF
      ROTMAT(51,3,3) =  1.d0
      ROTNAME(51) = 'C6z'

      ROTMAT(52,1,1) =  HALF
      ROTMAT(52,1,2) = -RTHREE
      ROTMAT(52,2,1) =  RTHREE
      ROTMAT(52,2,2) =  HALF
      ROTMAT(52,3,3) =  1.d0
      ROTNAME(52) = 'C6z-1'

      ROTMAT(53,1,1) = -HALF
      ROTMAT(53,1,2) =  RTHREE
      ROTMAT(53,2,1) =  RTHREE
      ROTMAT(53,2,2) =  HALF
      ROTMAT(53,3,3) = -1.d0
      ROTNAME(53) = 'C2A'    

      ROTMAT(54,1,1) = -HALF
      ROTMAT(54,1,2) = -RTHREE
      ROTMAT(54,2,1) = -RTHREE
      ROTMAT(54,2,2) =  HALF
      ROTMAT(54,3,3) = -1.d0
      ROTNAME(54) = 'C2B'

      ROTMAT(55,1,1) =  HALF
      ROTMAT(55,1,2) = -RTHREE
      ROTMAT(55,2,1) = -RTHREE
      ROTMAT(55,2,2) = -HALF
      ROTMAT(55,3,3) = -1.d0
      ROTNAME(55) = 'C2C'

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
      
  end subroutine pointgrp

end program
