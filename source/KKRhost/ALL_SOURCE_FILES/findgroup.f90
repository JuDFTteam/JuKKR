SUBROUTINE findgroup(bravais,recbv,rbasis,nbasis,  &
        rsymat,rotname,isymindex,nsymat,  &
        para,qmtet,qmphi,symunitary,  &
        krel,naezd,nembd,nsymaxd)
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

! This sub makes all 64 rotations in the basis vectors and bravais
! vectors and checks if the new rotated vectror belongs in the
! lattice. The proper rotation must bring all vectors to a lattice
! vector. Information about the rotations found is printed in the end.
! The array ISYMINDEX holds the numbers of the symmetry operations
! that are stored in array RSYMAT
!----------------------------------------------------------------
! in case of relativistic calculation: take account of
! direction of the magnetic moment specified by (QMTET,QMPHI)
! if the PARA(magnetic) flag is set to .FALSE.

! **********************************************************
IMPLICIT NONE
INTEGER KREL
INTEGER NAEZD,NEMBD,NSYMAXD
!..
INTEGER NBASIS,NSYMAT
INTEGER ISYMINDEX(NSYMAXD)
DOUBLE PRECISION BRAVAIS(3,3),RBASIS(3,NAEZD+NEMBD)
DOUBLE PRECISION RSYMAT(64,3,3),RECBV(3,3)
DOUBLE PRECISION QMTET(NAEZD), QMPHI(NAEZD)
LOGICAL SYMUNITARY(NSYMAXD),PARA
!..
!.. Local variables
DOUBLE PRECISION R(3,4),ROTRBAS(3,NAEZD+NEMBD)
DOUBLE PRECISION BRAVAIS1(3,3)
INTEGER I,J,ISYM,NSYM,I0,IA
DOUBLE PRECISION MDOTMP,MVECQ(3,NAEZD),MVECQP(3,NAEZD)
DOUBLE PRECISION MROTR(3,3),SYMDET,SUMMDOTMP
DOUBLE PRECISION STET,DDOT,DDET33,PI
CHARACTER*10 ROTNAME(64)
CHARACTER*10 CHAR(64)
LOGICAL LLATBAS,LATVEC,LBULK
!..................................................................

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE (1337,99000)
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT

nsym = 0
pi = 4.d0*ATAN(1.d0)
DO isym = 1,nsymaxd
  symunitary(isym) = .true.
END DO
!     - ---------------------------------
DO i=1,3
  DO j=1,3
    bravais1(j,i) = bravais(j,i)
  END DO
END DO
!     Check for surface mode. If so, set bravais1(3,3) very large, so
!     that only the in-plane symmetries are found.
!     not checked, be careful of z--> -z!

lbulk=.true.
!     Now check the bravais vectors if they have a z component
IF ((bravais(1,3) == 0.d0).AND.(bravais(2,3) == 0.d0).AND.  &
      (bravais(3,3) == 0.d0)) THEN
  lbulk=.false.
END IF

DO isym=1,64
  
!--------------------------------- store rotation matrix
  
  DO i=1,3
    DO j=1,3
      mrotr(i,j) = rsymat(isym,i,j)
    END DO
  END DO
  
  summdotmp = 0D0
  
  symdet = ddet33(mrotr)
  
!     rotate bravais lattice vectors
  
!     In the case of slab/interface geometry look only for
!     symmetry opperations that preserve the z axis..
  
  IF (lbulk .OR. (rsymat(isym,3,3) == 1) ) THEN
!     do rotation only in case bulk or if slab and z axis is restored..
    
    
    DO i=1,3            ! Loop on bravais vectors
      DO j=1,3         ! Loop on coordinates
        r(j,i) = rsymat(isym,j,1)*bravais1(1,i) +  &
            rsymat(isym,j,2)*bravais1(2,i) + rsymat(isym,j,3)*bravais1(3,i)
      END DO
    END DO
    
!     rotate the basis atoms p and take RSYMAT.p - p then
!     find if R = (RSYMAT.bravais + RSYMAT.p - p) belongs to the
!     lattice. This is done by function latvec by checking
!     if R.q = integer (q reciprocal lattice vector)
    
    llatbas = .true.
    DO ia=1,nbasis      ! Loop on basis atoms
      DO j=1,3         ! Loop on coordinates
        rotrbas(j,ia) = rsymat(isym,j,1)*rbasis(1,ia) +  &
            rsymat(isym,j,2)*rbasis(2,ia) + rsymat(isym,j,3)*rbasis(3,ia)
        
        rotrbas(j,ia) = rotrbas(j,ia) - rbasis(j,ia)
        r(j,4) = rotrbas(j,ia)
      END DO
      
      IF (.NOT.latvec(4,recbv,r)) llatbas=.false.
      
      IF( (krel == 1) .AND. (.NOT.para) ) THEN
        stet = SIN(qmtet(ia)*pi/180D0)
        mvecq(1,ia) = stet*COS(qmphi(ia)*pi/180D0)
        mvecq(2,ia) = stet*SIN(qmphi(ia)*pi/180D0)
        mvecq(3,ia) = COS(qmtet(ia)*pi/180D0)
        
        CALL dgemv('N',3,3,1D0,mrotr,3,mvecq(1,ia),1,0D0, mvecqp(1,ia),1)
        
        CALL dscal(3,DBLE(symdet),mvecqp(1,ia),1)
        
        mdotmp = ddot(3,mvecq(1,ia),1,mvecqp(1,ia),1)
        summdotmp = summdotmp + mdotmp
        
      END IF
      
    END DO               ! ia=1,nbasis
    
    IF( (krel == 1) .AND. (.NOT.para) ) THEN
      IF( ABS(ABS(summdotmp)-nbasis) > 0.00001D0 ) THEN
        llatbas=.false.
      ELSE
        IF( summdotmp > 0.00001D0 ) THEN
          symunitary(nsym + 1) = .true.
        ELSE
          symunitary(nsym + 1) = .false.
        END IF
      END IF
    END IF
    
!     if llatbas=.true. the rotation does not change the lattice
    
    IF (llatbas) THEN
      nsym = nsym + 1
      isymindex(nsym) = isym
    END IF
  END IF                 ! (LBULK .OR. (RSYMAT(ISYM,3,3).EQ.1) )
END DO                    ! isym=1,nmatd
!     nsym symmetries were found
!     the ISYMINDEX array has the numbers of the symmetries found


nsymat = nsym

! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
WRITE(1337,'(8X,60(1H-))')
IF ( lbulk ) THEN
  WRITE(1337,99001)
ELSE
  WRITE(1337,99002)
END IF
WRITE(1337,99003) nsymat
DO i=1,nsymat
  i0 = isymindex(i)
  CHAR(i) =  rotname(i0)
END DO
nsym = nsymat/5
DO i=1,nsym + 1
  i0 = (i-1)*5
  isym = MIN(5,nsymat-i0)
  WRITE(1337,99004) (CHAR(j),j=i0+1,i0+isym)
END DO
WRITE(1337,99005)
99000 FORMAT (5X,'< FINDGROUP > : Finding symmetry operations',/)
99001 FORMAT (8X,'3D symmetries:')
99002 FORMAT (8X,'surface symmetries:')
99003 FORMAT(' found for this lattice: ',i2,/,8X,60(1H-))
99004 FORMAT(8X,5(a10,2X))
99005 FORMAT(8X,60(1H-),/)
END SUBROUTINE findgroup
