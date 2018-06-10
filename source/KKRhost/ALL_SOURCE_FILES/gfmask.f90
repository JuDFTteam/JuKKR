SUBROUTINE gfmask(linterface,icheck,icc,invmod,nsh1,nsh2,naez,  &
        nshell,naezd,nprincd)
! **********************************************************************
! *                                                                    *
! * This subroutine prepares the ICHECK matrix that is used for        *
! * calculating the proper off-diagonal GF matrix elements ( e.g.      *
! * impurity) in case of no full inversion algorithm                   *
! *                                                                    *
! * ICHECK(I,J) points to the block (I,J) of the GF matrix having the  *
! * size NPRINCD                                                       *
! *                                                                    *
! *                                            29.02.2000              *
! *                                                                    *
! **********************************************************************
      IMPLICIT NONE
!     ..
!     .. Scalar arguments
      INTEGER NAEZD,NPRINCD
      INTEGER ICC,INVMOD,NLAYER,NAEZ,NSHELL
      LOGICAL LINTERFACE
!     ..
!     .. Array arguments
      INTEGER ICHECK(NAEZD/NPRINCD,NAEZD/NPRINCD)    
      INTEGER NSH1(*),NSH2(*)
!     .. Local variables
      INTEGER ICOUPLE(NAEZD,NAEZD)
      INTEGER I,J,K,II,ISTEP1,ILT1,ISTEP2,ILT2,IL2,IL1,LFCHK
      CHARACTER*80 FMTCHK
      CHARACTER*35 INVALG(0:2)
!     ..
!     .. External functions
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
!     ..
!     .. Data statements
DATA invalg /'FULL MATRIX                        ',  &
    'BANDED MATRIX (slab)               ', 'BANDED + CORNERS MATRIX (supercell)' /

WRITE (1337,99000)

! --> set default inversion to SUPERCELL mode = banded matrix + corners

invmod = 2

! --> LINTERFACE = use band diagonal mode

IF (linterface) invmod = 1

! --> full inversion is performed ONLY BY EXPLICIT request

IF ( opt('full inv') )  invmod = 0

IF ( ( invmod /= 0 ).AND. ( MOD(naez,nprincd) /= 0 ) ) THEN
  WRITE(6,99001) naez,nprincd
  STOP
END IF

WRITE (1337,99002) invalg(invmod)

nlayer=naez/nprincd
! ----------------------------------------------------------- INVMOD = 1
!                                                   band-diagonal matrix
IF (invmod == 1) THEN
  DO i=1,nlayer
    DO j=1,nlayer
      IF (i == j) THEN
        icheck(i,j)=1
      ELSE
        icheck(i,j)=0
      END IF
    END DO
  END DO
END IF
! ----------------------------------------------------------- INVMOD = 2
!              band-diagonal matrix with corners (slab periodic along z)
IF (invmod == 2) THEN
  DO i=1,nlayer
    DO j=1,nlayer
      IF ( ( i == j ).OR. ( j == nlayer ) .OR. ( i == nlayer ) ) THEN
        icheck(i,j)=1
      ELSE
        icheck(i,j)=0
      END IF
    END DO
  END DO
END IF
! ================================================= INVMOD = 1, ICC <> 0
!                   band-diagonal matrix, off-diagonal G elements needed

! --> prepare the matrix ICOUPLE which has 1 in all nn' blocks
!     (atomic sites) that are needed

! ======================================================================
IF ( ( icc /= 0 ) .AND. ( invmod == 1 ) ) THEN
  DO i=1,naez
    DO j=1,naez
      icouple(i,j) = 0
      
      DO ii=1,nshell
        IF ( ( ( nsh1(ii) == i ) .AND. ( nsh2(ii) == j ) )  &
            .OR. ( ( nsh1(ii) == j ) .AND. ( nsh2(ii) == i ) ) ) icouple(i,j) = 1
      END DO
    END DO
  END DO
!cccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!cccC                                               conductivity calculation
!ccc         IF (OPT('CONDUCT ')) THEN
!ccc            DO I=1,NLAYER
!ccc               DO J=1,NLAYER
!ccc                  ICHECK(I,J)=0
!ccc               ENDDO
!ccc            ENDDO
!ccc            DO I=1,NAEZ
!ccc               DO J=1,NAEZ
!ccc                  ICOUPLE(I,J) = 0
!ccc               END DO
!ccc            END DO
!ccc            DO II=1,NCONDPAIR
!ccc               I = IATCONDL(II)
!ccc               J = IATCONDR(II)
!ccc               ICOUPLE(I,J) = 1
!ccc               ICOUPLE(J,I) = 1
!ccc            ENDDO
!ccc         END IF                 ! Conductivity calculation
!cccC ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
! ----------------------------------------------------------------------
! Now given the matrix ICOUPLE prepare the matrix ICHECK which has 1 in
! all principal-layer blocks that we need -- this will be used in the
! matrix inversion
! ----------------------------------------------------------------------
  istep1=0
  ilt1=1
! ----------------------------------------------------------------------
  DO il1=1,naez
    istep1=istep1+1
    
    IF ( istep1 > nprincd ) THEN
      ilt1=ilt1+1
      istep1=1
    END IF
    
    ilt2=1
    istep2=0
! ......................................................................
    DO il2=1,naez
      istep2=istep2+1
      
      IF ( istep2 > nprincd ) THEN
        ilt2=ilt2+1
        istep2=1
      END IF
      
      IF ( icouple(il1,il2) == 1 ) icheck(ilt1,ilt2)=1
    END DO
! ......................................................................
  END DO
! ----------------------------------------------------------------------
! in the case of calculation of single blocks it has to put the correct
! value to ICHECK in order to calculate all the elements also necessary
! to calculate that single block          ?????
! ----------------------------------------------------------------------
  DO j=1,nlayer
    
! --> loop over the element ICHECK(I,J) with fixed J and I < J
    
    IF (j /= 1) THEN
      DO i=1,j-1
        IF (icheck(i,j) == 1) THEN
          DO k=i+1,j
            icheck(k,j)=1
          END DO
          DO k=j,nlayer
            icheck(k,k)=1
          END DO
        END IF
      END DO
    END IF
    
    IF ( .NOT.opt('CONDUCT ') ) THEN
      
! --> loop over the element ICHECK(I,J) with fixed J and I > J
      
      IF (j /= nlayer) THEN
        DO i=nlayer,j+1,-1
          IF (icheck(i,j) == 1) THEN
            DO k=i-1,j,-1
              icheck(k,j)=1
            END DO
          END IF
        END DO
      END IF
    END IF
  END DO
! ----------------------------------------------------------------------
END IF
! ======================================================================

IF ( test('ICHECK  ') ) THEN
  
  fmtchk=' '
  lfchk = 1
  DO i = 1,MIN(35,nlayer)
    fmtchk=fmtchk(1:lfchk)//'--'
    lfchk = lfchk+2
  END DO
  
  WRITE (1337,'(8X,A,/,8X,A)') 'ICHECK matrix :',fmtchk(1:lfchk)
  DO i=1,nlayer
    WRITE (1337,'(9X,35I2)') (icheck(i,j),j=1,MIN(35,nlayer))
  END DO
  WRITE (1337,'(8X,A,/)') fmtchk(1:lfchk)
END IF

99000 FORMAT (5X,'< GFMASK > : set KKR matrix inversion algorithm',/)
99001 FORMAT(6X,"ERROR: Number of sites (NAEZ) =",i3,  &
    " not an integer multiplier",/,6X,  &
    "of principal layers (NPRINCD) =",i3,/,6X,  &
    "Use ONLY  full inversion in this case")
99002 FORMAT (8X,'INVERSION algorithm used : ',a35,/)
END SUBROUTINE gfmask
