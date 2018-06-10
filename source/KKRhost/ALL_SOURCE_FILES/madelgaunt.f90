SUBROUTINE madelgaunt(lpot,yrg,wg,cleb,icleb,iend,lassld,nclebd)
      IMPLICIT NONE
!..
!.. Scalar arguments
      INTEGER LPOT,IEND
      INTEGER LASSLD,NCLEBD
!..
!.. Array arguments
!.. Attention: Dimension NCLEBD appears sometimes as NCLEB1
!..            an empirical factor - it has to be optimized
      DOUBLE PRECISION YRG(LASSLD,0:LASSLD,0:LASSLD),WG(LASSLD)
      DOUBLE PRECISION CLEB(NCLEBD)
      INTEGER ICLEB(NCLEBD,3)
!..
!.. Local scalars
      DOUBLE PRECISION CLECG,FACTOR,S
      INTEGER I,J,L1,L2,L3,M1,M1A,M1S,M2,M2A,M2S,M3,M3A,M3S
!..
!.. Intrinsic functions
      INTRINSIC ABS,ATAN,DBLE,SIGN

! --> set up of the gaunt coefficients with an index field
!     recognize that they are needed here only for l3=l1+l2

IF ( 2*lpot > lassld ) THEN
  WRITE (6,*) 'Dim ERROR in MADELGAUNT -- 2*LPOT > LASSLD', 2*lpot,lassld
  STOP
END IF

i = 1
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
DO l1 = 0,lpot
  DO l2 = 0,lpot
    l3 = l1 + l2
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
    DO m1 = -l1,l1
      DO m2 = -l2,l2
        DO m3 = -l3,l3
          m1s = SIGN(1,m1)
          m2s = SIGN(1,m2)
          m3s = SIGN(1,m3)
! **********************************************************************
          IF ( m1s*m2s*m3s >= 0 ) THEN
            m1a = ABS(m1)
            m2a = ABS(m2)
            m3a = ABS(m3)
            
            factor = 0.0D0
            IF ( m1a+m2a == m3a ) factor = factor +  &
                DBLE(3*m3s+SIGN(1,-m3))/8.0D0
            IF ( m1a-m2a == m3a ) factor = factor + DBLE(m1s)/4.0D0
            IF ( m2a-m1a == m3a ) factor = factor + DBLE(m2s)/4.0D0
! ======================================================================
            IF ( factor /= 0.0D0 ) THEN
              IF ( m1s*m2s /= 1 .OR. m2s*m3s /= 1 .OR.  &
                  m1s*m3s /= 1 ) factor = -factor
              
              s = 0.0D0
              DO j = 1,lassld
                s = s + wg(j)*yrg(j,l1,m1a)*yrg(j,l2,m2a) *yrg(j,l3,m3a)
              END DO
              
              clecg = s*factor
! ----------------------------------------------------------------------
              IF ( ABS(clecg) > 1.d-10 ) THEN
                cleb(i) = clecg
                icleb(i,1) = l1*(l1+1) + m1 + 1
                icleb(i,2) = l2*(l2+1) + m2 + 1
                icleb(i,3) = l3*(l3+1) + m3 + 1
                i = i + 1
                IF ( i > nclebd ) THEN
                  WRITE (6,FMT='(2I10)') i,nclebd
                  STOP ' Dim stop in MADELGAUNT '
                END IF
              END IF
! ----------------------------------------------------------------------
            END IF
! ======================================================================
          END IF
! **********************************************************************
        END DO
      END DO
    END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  END DO
END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
iend = i - 1
END SUBROUTINE madelgaunt
