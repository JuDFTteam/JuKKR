module Symmetry_mod
  implicit none
  private
  public :: pointgrp, findgroup

  contains
  
      subroutine findgroup(bravais,recbv,rbasis,nbasis, rsymat,rotname,isymindex,nsymat, NAEZD)
! **********************************************************
! This subroutine finds the rotation matrices that leave the
! real lattice unchanged. 
! input:  bravais(i,j)    true bravais lattice vectors
!                         i = x,y,z ; j = A, B, C (a.u.)
!         recbv(i,j)      reciprocal basis vectors 
!         rbasis          coordinates of basis atoms
!         nbasis          number of basis atoms
!         rsymat          all 64 rotation matrices.
!         rotname         names for the rotation matrices
! output: nsymat          number of rotations that restore the lattice.
!         ISYMINDEX       index for the symmeties found
!
! This sub makes all 64 rotations in the basis vectors and bravais
! vectors and checks if the new rotated vectror belongs in the 
! lattice. The proper rotation must bring all vectors to a lattice
! vector. Information about the rotations found is printed in the end.         
! The array ISYMINDEX holds the numbers of the symmetry operations
! that are stored in array RSYMAT
! **********************************************************
      implicit none

      integer NSYMAXD
      parameter (NSYMAXD=48)

      integer nbasis,nsymat
      integer isymindex(NSYMAXD)
      double precision BRAVAIS(3,3),RBASIS(3,NAEZD)
      double precision RSYMAT(64,3,3),recbv(3,3)
      integer NAEZD
!
! Local variables
!
      double precision r(3,4),rotrbas(3,naezd)
      double precision bravais1(3,3)
      integer i,j,isym,nsym,i0,ia
      double precision MDOTMP,MVECQ(3,NAEZD),MVECQP(3,NAEZD)
      double precision MROTR(3,3),SYMDET,SUMMDOTMP
      double precision STET
      double precision DDOT,DDET33,PI
      Character*10 ROTNAME(64)
      CHARACTER*10 CHAR(64)
      logical llatbas,latvec,LBULK,TEST
!     -------------------------------------------------------------
!
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99000)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
!
      NSYM = 0
      PI = 4.D0*ATAN(1.D0)
      do i=1,3
         do j=1,3
            bravais1(j,i) = bravais(j,i) 
         end do
      end do
!     Check for surface mode. If so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found. Not checked, be careful of z--> -z!
      LBULK=.TRUE.
!     Now check the bravais vectors if they have a z component 
      if ((bravais(1,3).eq.0.d0).and.(bravais(2,3).eq.0.d0).and.(bravais(3,3).eq.0.d0)) THEN
         LBULK=.FALSE.
      END IF
!     
      do isym=1,64
!
!--------------------------------- store rotation matrix
!
         DO I=1,3    
            DO J=1,3 
               MROTR(I,J) = RSYMAT(ISYM,I,J)
            END DO
         END DO

         SUMMDOTMP = 0D0

         SYMDET = DDET33(MROTR)
!
!     rotate bravais lattice vectors
!     
!     In the case of slab/interface geometry look only for 
!     symmetry opperations that preserve the z axis..
!     
         IF (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) ) THEN 
!     do rotation only in case bulk or if slab and z axis is restored..  
            
            
            do i=1,3            ! Loop on bravais vectors
               do j=1,3         ! Loop on coordinates
                  r(j,i) = rsymat(isym,j,1)*bravais1(1,i) + rsymat(isym,j,2)*bravais1(2,i) + rsymat(isym,j,3)*bravais1(3,i)
               enddo
            enddo
!     
!     rotate the basis atoms p and take RSYMAT.p - p then
!     find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
!     lattice. This is done by function latvec by checking
!     if R.q = integer (q reciprocal lattice vector)
!     
            llatbas = .true.
            do ia=1,nbasis      ! Loop on basis atoms
               do j=1,3         ! Loop on coordinates
                  rotrbas(j,ia) = rsymat(isym,j,1)*rbasis(1,ia) + rsymat(isym,j,2)*rbasis(2,ia) + rsymat(isym,j,3)*rbasis(3,ia)
!     
                  rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
                  r(j,4) = rotrbas(j,ia) 
               enddo
!     do j=1,4
!     do i=1,3
!     write(6,*) 'rrr',i,j,r(i,j)
!     end do
!     end do
!     write(6,*) 'latvec',latvec(4,recbv,r)
               if (.not.latvec(4,recbv,r)) llatbas=.false.
      IF(TEST('Oh-symm ').AND.ISYM.LE.48)  LLATBAS=.TRUE.
      IF(TEST('Td-symm ').AND.ISYM.LE.12)  LLATBAS=.TRUE.
      IF(TEST('Td-symm ').AND.ISYM.GE.37 .AND. ISYM.LE.48)  LLATBAS=.TRUE.

            enddo               ! ia=1,nbasis

!     
!     if llatbas=.true. the rotation does not change the lattice 
!     
!     write(6,*) 'llatbas',llatbas
            if (llatbas) then
               NSYM = NSYM + 1
               ISYMINDEX(NSYM) = ISYM
            end if
         END IF                 ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
      end do                    ! isym=1,nmatd
!     nsym symmetries were found
!     the ISYMINDEX array has the numbers of the symmetries found
!     
!     
      NSYMAT = NSYM 
!
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE(6,'(8X,60(1H-))')
      IF ( LBULK ) THEN
         WRITE(6,99001) 
      ELSE
         WRITE(6,99002) 
      END IF
      WRITE(6,99003) NSYMAT
      DO I=1,NSYMAT
         I0 = ISYMINDEX(I)
         CHAR(I) =  ROTNAME(I0) 
      END DO
      NSYM = NSYMAT/5
      DO I=1,NSYM + 1
         ISYM = MIN(5,NSYMAT-(I-1)*5)
         WRITE(6,99004) (CHAR(J),J=(I-1)*5+1,(I-1)*5+ISYM)
!        WRITE(6,99004) (CHAR(J),J=(I-1)*NSYM+1,(I-1)*NSYM+ISYM)
      END DO
      WRITE(6,99005)
99000 FORMAT (5X,'< FINDGROUP > : Finding symmetry operations',/)
99001 FORMAT (8X,'3D symmetries',$)
99002 FORMAT (8X,'surface symmetries',$)
99003 FORMAT(' found for this lattice: ',I2,/,8X,60(1H-))
99004 FORMAT(8X,5(A10,2X))
99005 FORMAT(8X,60(1H-),/)
      endsubroutine
  
  
  
  subroutine pointgrp(rotmat, rotname)
! **********************************************
! This subroutine defines the rotation matrices for
! all the 32 point groups and names them after 
! J.F. Cornwell (Group Theory??) second edition 
! Appendix D, p 324-325
! 
! *********************************************    
    double precision, intent(out) :: rotmat(64,3,3)
    character(len=*), intent(out) :: rotname(64)
    integer :: is
    double precision, parameter :: RTHREE = SQRT(3.d0)/2.d0, HALF=0.5d0, ONE=1.d0

    rotmat(:,:,:) =  0.d0 ! set all matrices to zero

    rotmat(1,1,1) =  ONE
    rotmat(1,2,2) =  ONE
    rotmat(1,3,3) =  ONE
    rotname(1) = 'E'

    rotmat(2,1,2) =  ONE
    rotmat(2,2,3) = -ONE
    rotmat(2,3,1) = -ONE
    rotname(2) = 'C3alfa'          

    rotmat(3,1,2) = -ONE
    rotmat(3,2,3) = -ONE
    rotmat(3,3,1) =  ONE
    rotname(3) = 'C3beta '

    rotmat(4,1,2) = -ONE
    rotmat(4,2,3) =  ONE
    rotmat(4,3,1) = -ONE
    rotname(4) = 'C3gamma'

    rotmat(5,1,2) = ONE
    rotmat(5,2,3) = ONE
    rotmat(5,3,1) = ONE
    rotname(5) = 'C3delta '

    rotmat(6,1,3) = -ONE
    rotmat(6,2,1) =  ONE
    rotmat(6,3,2) = -ONE
    rotname(6) = 'C3alfa-1'

    rotmat(7,1,3) =  ONE
    rotmat(7,2,1) = -ONE
    rotmat(7,3,2) = -ONE
    rotname(7) = 'C3beta-1 '

    rotmat(8,1,3) = -ONE
    rotmat(8,2,1) = -ONE
    rotmat(8,3,2) =  ONE
    rotname(8) = 'C3gamma-1'

    rotmat(9,1,3) =  ONE
    rotmat(9,2,1) =  ONE
    rotmat(9,3,2) =  ONE
    rotname(9) = 'C3delta-1'

    rotmat(10,1,1) =  ONE
    rotmat(10,2,2) = -ONE
    rotmat(10,3,3) = -ONE
    rotname(10) = 'C2x'
          
    rotmat(11,1,1) = -ONE
    rotmat(11,2,2) =  ONE
    rotmat(11,3,3) = -ONE
    rotname(11) = 'C2y'

    rotmat(12,1,1) = -ONE
    rotmat(12,2,2) = -ONE
    rotmat(12,3,3) =  ONE
    rotname(12) = 'C2z'

    rotmat(13,1,1) =  ONE
    rotmat(13,2,3) =  ONE
    rotmat(13,3,2) = -ONE
    rotname(13) = 'C4x'
          
    rotmat(14,1,3) = -ONE
    rotmat(14,2,2) =  ONE
    rotmat(14,3,1) =  ONE
    rotname(14) = 'C4y '

    rotmat(15,1,2) =  ONE
    rotmat(15,2,1) = -ONE
    rotmat(15,3,3) =  ONE
    rotname(15) = 'C4z'
          
    rotmat(16,1,1) =  ONE
    rotmat(16,2,3) = -ONE
    rotmat(16,3,2) =  ONE
    rotname(16) = 'C4x-1 '

    rotmat(17,1,3) =  ONE
    rotmat(17,2,2) =  ONE
    rotmat(17,3,1) = -ONE
    rotname(17) = 'C4y-1'

    rotmat(18,1,2) = -ONE
    rotmat(18,2,1) =  ONE
    rotmat(18,3,3) =  ONE
    rotname(18) = 'C4z-1'
          
    rotmat(19,1,2) =  ONE
    rotmat(19,2,1) =  ONE
    rotmat(19,3,3) = -ONE
    rotname(19) = 'C2a'

    rotmat(20,1,2) = -ONE
    rotmat(20,2,1) = -ONE
    rotmat(20,3,3) = -ONE
    rotname(20) = 'C2b'

    rotmat(21,1,3) =  ONE
    rotmat(21,2,2) = -ONE
    rotmat(21,3,1) =  ONE
    rotname(21) = 'C2c'

    rotmat(22,1,3) = -ONE
    rotmat(22,2,2) = -ONE
    rotmat(22,3,1) = -ONE
    rotname(22) = 'C2d'

    rotmat(23,1,1) = -ONE
    rotmat(23,2,3) =  ONE
    rotmat(23,3,2) =  ONE
    rotname(23) = 'C2e'

    rotmat(24,1,1) = -ONE
    rotmat(24,2,3) = -ONE
    rotmat(24,3,2) = -ONE
    rotname(24) = 'C2f'
    
    do is = 1, 24
      rotmat(is+24,1:3,1:3) = -rotmat(is,1:3,1:3)
      rotname(is+24) = 'I'//rotname(is)
    enddo ! i1      

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

    rotmat(50,1,1) = -HALF
    rotmat(50,1,2) = -RTHREE
    rotmat(50,2,1) =  RTHREE
    rotmat(50,2,2) = -HALF
    rotmat(50,3,3) =  ONE
    rotname(50) = 'C3z-1'

    rotmat(51,1,1) =  HALF
    rotmat(51,1,2) =  RTHREE
    rotmat(51,2,1) = -RTHREE
    rotmat(51,2,2) =  HALF
    rotmat(51,3,3) =  ONE
    rotname(51) = 'C6z'

    rotmat(52,1,1) =  HALF
    rotmat(52,1,2) = -RTHREE
    rotmat(52,2,1) =  RTHREE
    rotmat(52,2,2) =  HALF
    rotmat(52,3,3) =  ONE
    rotname(52) = 'C6z-1'

    rotmat(53,1,1) = -HALF
    rotmat(53,1,2) =  RTHREE
    rotmat(53,2,1) =  RTHREE
    rotmat(53,2,2) =  HALF
    rotmat(53,3,3) = -ONE
    rotname(53) = 'C2A'    

    rotmat(54,1,1) = -HALF
    rotmat(54,1,2) = -RTHREE
    rotmat(54,2,1) = -RTHREE
    rotmat(54,2,2) =  HALF
    rotmat(54,3,3) = -ONE
    rotname(54) = 'C2B'

    rotmat(55,1,1) =  HALF
    rotmat(55,1,2) = -RTHREE
    rotmat(55,2,1) = -RTHREE
    rotmat(55,2,2) = -HALF
    rotmat(55,3,3) = -ONE
    rotname(55) = 'C2C'

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
      
  endsubroutine pointgrp

endmodule Symmetry_mod      