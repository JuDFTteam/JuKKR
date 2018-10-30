      SUBROUTINE FINDGROUP(BRAVAIS,RECBV,RBASIS,NBASIS,
     &                     RSYMAT,ROTNAME,ISYMINDEX,NSYMAT,
     &                     PARA,QMTET,QMPHI,SYMUNITARY,
     &                     KREL,NAEZD,NEMBD,NSYMAXD)
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
c----------------------------------------------------------------
c in case of relativistic calculation: take account of 
c direction of the magnetic moment specified by (QMTET,QMPHI)
c if the PARA(magnetic) flag is set to .FALSE.
c
c **********************************************************
      IMPLICIT NONE
      INTEGER KREL
      INTEGER NAEZD,NEMBD,NSYMAXD
C     ..
      INTEGER NBASIS,NSYMAT
      INTEGER ISYMINDEX(NSYMAXD)
      DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,NAEZD+NEMBD)
      DOUBLE PRECISION RSYMAT(64,3,3),RECBV(3,3)
      DOUBLE PRECISION QMTET(NAEZD), QMPHI(NAEZD)
      LOGICAL SYMUNITARY(NSYMAXD),PARA
C     ..
C     .. Local variables
      DOUBLE PRECISION R(3,4),ROTRBAS(3,NAEZD+NEMBD)
      DOUBLE PRECISION BRAVAIS1(3,3)
      INTEGER I,J,ISYM,NSYM,I0,IA
      DOUBLE PRECISION MDOTMP,MVECQ(3,NAEZD),MVECQP(3,NAEZD)
      DOUBLE PRECISION MROTR(3,3),SYMDET,SUMMDOTMP
      DOUBLE PRECISION STET,DDOT,DDET33,PI
      CHARACTER*10 ROTNAME(64)
      CHARACTER*10 CHAR(64)
      LOGICAL LLATBAS,LATVEC,LBULK
C     ..................................................................
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,99000)
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      NSYM = 0
      PI = 4.D0*ATAN(1.D0)
      DO ISYM = 1,NSYMAXD
          SYMUNITARY(ISYM) = .TRUE.
      END DO
c     - ---------------------------------
      do i=1,3
         do j=1,3
            bravais1(j,i) = bravais(j,i) 
         end do
      end do
C     Check for surface mode. If so, set bravais1(3,3) very large, so
C     that only the in-plane symmetries are found. 
C     not checked, be careful of z--> -z!
C
      LBULK=.TRUE.
C     Now check the bravais vectors if they have a z component 
      IF ((BRAVAIS(1,3).EQ.0.D0).AND.(BRAVAIS(2,3).EQ.0.D0).AND.
     &     (BRAVAIS(3,3).EQ.0.D0)) THEN
         LBULK=.FALSE.
      END IF
C     
      DO ISYM=1,64
C
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

               IF (.NOT.LATVEC(4,RECBV,R)) LLATBAS=.FALSE.

               IF( (KREL.EQ.1) .AND. (.NOT.PARA) ) THEN
                  STET = SIN(QMTET(IA)*PI/180D0)
                  MVECQ(1,IA) = STET*COS(QMPHI(IA)*PI/180D0)
                  MVECQ(2,IA) = STET*SIN(QMPHI(IA)*PI/180D0)
                  MVECQ(3,IA) = COS(QMTET(IA)*PI/180D0)
C
                  CALL DGEMV('N',3,3,1D0,MROTR,3,MVECQ(1,IA),1,0D0,
     &                       MVECQP(1,IA),1)
C
                  CALL DSCAL(3,DBLE(SYMDET),MVECQP(1,IA),1)

                  MDOTMP = DDOT(3,MVECQ(1,IA),1,MVECQP(1,IA),1)
                  SUMMDOTMP = SUMMDOTMP + MDOTMP
                  
               END IF

            enddo               ! ia=1,nbasis

            IF( (KREL.EQ.1) .AND. (.NOT.PARA) ) THEN
                IF( ABS(ABS(SUMMDOTMP)-NBASIS) .GT. 0.00001D0 ) THEN
                   LLATBAS=.FALSE.   
                ELSE 
                   IF( SUMMDOTMP .GT. 0.00001D0 ) THEN
                      SYMUNITARY(NSYM + 1) = .TRUE.
                   ELSE          
                      SYMUNITARY(NSYM + 1) = .FALSE.
                   END IF   
                END IF
            END IF
c     
c     if llatbas=.true. the rotation does not change the lattice 
c     
            IF (LLATBAS) THEN
               NSYM = NSYM + 1
               ISYMINDEX(NSYM) = ISYM
            END IF
         END IF                 ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
      end do                    ! isym=1,nmatd
c     nsym symmetries were found
c     the ISYMINDEX array has the numbers of the symmetries found
c     
c     
      NSYMAT = NSYM 
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE(1337,'(8X,60(1H-))')
      IF ( LBULK ) THEN
         WRITE(1337,99001) 
      ELSE
         WRITE(1337,99002) 
      END IF
      WRITE(1337,99003) NSYMAT
      DO I=1,NSYMAT
         I0 = ISYMINDEX(I)
         CHAR(I) =  ROTNAME(I0) 
      END DO
      NSYM = NSYMAT/5
      DO I=1,NSYM + 1
         I0 = (I-1)*5
         ISYM = MIN(5,NSYMAT-I0)
         WRITE(1337,99004) (CHAR(J),J=I0+1,I0+ISYM)
      END DO
      WRITE(1337,99005)
99000 FORMAT (5X,'< FINDGROUP > : Finding symmetry operations',/)
99001 FORMAT (8X,'3D symmetries:')
99002 FORMAT (8X,'surface symmetries:')
99003 FORMAT(' found for this lattice: ',I2,/,8X,60(1H-))
99004 FORMAT(8X,5(A10,2X))
99005 FORMAT(8X,60(1H-),/)
      END
