      subroutine findgroup(bravais,recbv,rbasis,nbasis,
     &                     rsymat,rotname,isymindex,nsymat)
c **********************************************************
c This subroutine finds the rotation matrices that leave the
c real lattice unchanged. 
c input:  bravais(i,j)    true bravais lattice vectors
c                         i = x,y,z ; j = A, B, C (a.u.)
c         recbv(i,j)      reciprocal basis vectors 
c         rbasis          coordinates of basis atoms
c         nbasis          number of basis atoms
c         rsymat          all 64 rotation matrices.
c         rotname         names for the rotation matrices
c output: nsymat          number of rotations that restore the lattice.
c         ISYMINDEX       index for the symmeties found
c
c This sub makes all 64 rotations in the basis vectors and bravais
c vectors and checks if the new rotated vectror belongs in the 
c lattice. The proper rotation must bring all vectors to a lattice
c vector. Information about the rotations found is printed in the end.         
c The array ISYMINDEX holds the numbers of the symmetry operations
c that are stored in array RSYMAT
c **********************************************************
      implicit none
      include 'inc.p'
      integer NSYMAXD
      parameter (NSYMAXD=48)
      integer nbasis,nsymat
      integer isymindex(NSYMAXD)
      double precision BRAVAIS(3,3),RBASIS(3,NAEZD)
      double precision RSYMAT(64,3,3),recbv(3,3)
c
c Local variables
c
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
c     -------------------------------------------------------------
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,99000)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      NSYM = 0
      PI = 4.D0*ATAN(1.D0)
      do i=1,3
         do j=1,3
            bravais1(j,i) = bravais(j,i) 
         end do
      end do
c     Check for surface mode. If so, set bravais1(3,3) very large, so
c     that only the in-plane symmetries are found. Not checked, be careful of z--> -z!
      LBULK=.TRUE.
c     Now check the bravais vectors if they have a z component 
      if ((bravais(1,3).eq.0.d0).and.(bravais(2,3).eq.0.d0).and.
     &     (bravais(3,3).eq.0.d0)) THEN
         LBULK=.FALSE.
      END IF
c     
      do isym=1,64
c
c--------------------------------- store rotation matrix
c
         DO I=1,3    
            DO J=1,3 
               MROTR(I,J) = RSYMAT(ISYM,I,J)
            END DO
         END DO

         SUMMDOTMP = 0D0

         SYMDET = DDET33(MROTR)
c
c     rotate bravais lattice vectors
c     
c     In the case of slab/interface geometry look only for 
c     symmetry opperations that preserve the z axis..
c     
         IF (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) ) THEN 
c     do rotation only in case bulk or if slab and z axis is restored..  
            
            
            do i=1,3            ! Loop on bravais vectors
               do j=1,3         ! Loop on coordinates
                  r(j,i) = rsymat(isym,j,1)*bravais1(1,i) +
     &                 rsymat(isym,j,2)*bravais1(2,i) +
     &                 rsymat(isym,j,3)*bravais1(3,i)
               enddo
            enddo
c     
c     rotate the basis atoms p and take RSYMAT.p - p then
c     find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
c     lattice. This is done by function latvec by checking
c     if R.q = integer (q reciprocal lattice vector)
c     
            llatbas = .true.
            do ia=1,nbasis      ! Loop on basis atoms
               do j=1,3         ! Loop on coordinates
                  rotrbas(j,ia) = rsymat(isym,j,1)*rbasis(1,ia) +
     &                 rsymat(isym,j,2)*rbasis(2,ia) +
     &                 rsymat(isym,j,3)*rbasis(3,ia)
c     
                  rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
                  r(j,4) = rotrbas(j,ia) 
               enddo
c     do j=1,4
c     do i=1,3
c     write(6,*) 'rrr',i,j,r(i,j)
c     end do
c     end do
c     write(6,*) 'latvec',latvec(4,recbv,r)
               if (.not.latvec(4,recbv,r)) llatbas=.false.
      IF(TEST('Oh-symm ').AND.ISYM.LE.48)  LLATBAS=.TRUE.
      IF(TEST('Td-symm ').AND.ISYM.LE.12)  LLATBAS=.TRUE.
      IF(TEST('Td-symm ').AND.ISYM.GE.37
     +                   .AND.ISYM.LE.48)  LLATBAS=.TRUE.

            enddo               ! ia=1,nbasis

c     
c     if llatbas=.true. the rotation does not change the lattice 
c     
c     write(6,*) 'llatbas',llatbas
            if (llatbas) then
               NSYM = NSYM + 1
               ISYMINDEX(NSYM) = ISYM
            end if
         END IF                 ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
      end do                    ! isym=1,nmatd
c     nsym symmetries were found
c     the ISYMINDEX array has the numbers of the symmetries found
c     
c     
      NSYMAT = NSYM 
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
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
c        WRITE(6,99004) (CHAR(J),J=(I-1)*NSYM+1,(I-1)*NSYM+ISYM)
      END DO
      WRITE(6,99005)
99000 FORMAT (5X,'< FINDGROUP > : Finding symmetry operations',/)
99001 FORMAT (8X,'3D symmetries',$)
99002 FORMAT (8X,'surface symmetries',$)
99003 FORMAT(' found for this lattice: ',I2,/,8X,60(1H-))
99004 FORMAT(8X,5(A10,2X))
99005 FORMAT(8X,60(1H-),/)
      END
