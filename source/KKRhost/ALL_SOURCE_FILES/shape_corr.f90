SUBROUTINE shape_corr(lpot,natyp,gsh,ilm_map,imaxsh,lmsp,ntcell,  &
        w,yr,lassld,lmpotd,natypd,ngshd)
! **********************************************************************
! *  Prepares shape corrections using gaussian quadrature as given by  *
! *  m. abramowitz and i.a. stegun, handbook of mathematical functions *
! *  nbs applied mathematics series 55 (1968), pages 887 and 916       *
! *                                                                    *
! *  the parameter LASSLD has to be chosen such that                   *
! *                        l1+l2+l3 .le. 2*LASSLD                      *
! *                                                                    *
! **********************************************************************

      IMPLICIT NONE
!..
!.. Scalar Arguments ..
      INTEGER LASSLD,LMPOTD,NATYPD,NGSHD
      INTEGER LPOT,NATYP
!..
!.. Array Arguments ..
      DOUBLE PRECISION GSH(*),W(LASSLD),YR(LASSLD,0:LASSLD,0:LASSLD)
      INTEGER ILM_MAP(NGSHD,3),IMAXSH(0:LMPOTD),LMSP(NATYPD,*),NTCELL(*)
!..
!.. Local Scalars ..
      DOUBLE PRECISION FACTOR,GAUNT,S
      INTEGER I,IAT,ICELL,ISUM,J,L1,L2,L3,LM1,LM2,LM3,M1,M1A,M1S,M2,M2A, &
              M2S,M3,M3A,M3S
      LOGICAL TRIANGLE
!..
!.. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,SIGN
!..
!.. External Subroutines ..
      EXTERNAL RCSTOP,TRIANGLE
!..

! -> set up of the gaunt coefficients with an index field
!    so that  c(lm,lm',lm'') is mapped to c(i)
i = 1
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
DO l1 = 0,lpot
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
  DO m1 = -l1,l1
    
    lm1 = l1*l1 + l1 + m1 + 1
    imaxsh(lm1-1) = i - 1
! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
    DO l3 = 0,lpot*2
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
      DO m3 = -l3,l3
        
        lm3 = l3*l3 + l3 + m3 + 1
        isum = 0
        
        DO iat = 1,natyp
          icell = ntcell(iat)
!                      write(*,*) 'test icell=ntcell(iat) in shape_corr.f',
!      +                icell,iat
          isum = isum + lmsp(icell,lm3)
        END DO
        
! ======================================================================
        IF ( isum > 0 ) THEN
          DO l2 = 0,lpot
! ----------------------------------------------------------------------
            IF ( triangle(l1,l2,l3) ) THEN
              DO m2 = -l2,l2
                
                lm2 = l2*l2 + l2 + m2 + 1
                
! -> use the m-conditions for the gaunt coefficients not to be 0
                
                m1s = SIGN(1,m1)
                m2s = SIGN(1,m2)
                m3s = SIGN(1,m3)
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
                IF ( m1s*m2s*m3s >= 0 ) THEN
                  m1a = ABS(m1)
                  m2a = ABS(m2)
                  m3a = ABS(m3)
                  factor = 0.0D0
                  
                  IF (m1a+m2a == m3a) factor = factor +  &
                      DBLE(3*m3s+SIGN(1,-m3))/8.0D0
                  
                  IF (m1a-m2a == m3a) factor = factor + DBLE(m1s)/4.0D0
                  
                  IF (m2a-m1a == m3a) factor = factor + DBLE(m2s)/4.0D0
! ......................................................................
                  IF (factor /= 0.0D0) THEN
                    
                    IF ( m1s*m2s /= 1 .OR. m2s*m3s /= 1  &
                        .OR.m1s*m3s /= 1 ) factor = -factor
                    
                    s = 0.0D0
                    DO j = 1,lassld
                      s = s + w(j) * yr(j,l1,m1a)  &
                          * yr(j,l2,m2a) * yr(j,l3,m3a)
                    END DO
                    
                    gaunt = s*factor
                    IF ( ABS(gaunt) > 1D-10 ) THEN
                      gsh(i) = gaunt
                      ilm_map(i,1) = lm1
                      ilm_map(i,2) = lm2
                      ilm_map(i,3) = lm3
                      i = i + 1
                    END IF
                  END IF
! ......................................................................
                END IF
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
              END DO
            END IF
! ----------------------------------------------------------------------
          END DO
        END IF
! ======================================================================
      END DO
! mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
    END DO
! llllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllllll
  END DO
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
END DO
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL

imaxsh(lm1) = i - 1
WRITE (1337,FMT=9000) imaxsh(lm1),ngshd
IF ( imaxsh(lm1) > ngshd ) CALL rcstop('SHAPE   ')

9000 FORMAT(' >>> SHAPE : IMAXSH(',i4,'),NGSHD :',2I6)

END SUBROUTINE shape_corr

FUNCTION triangle(l1,l2,l3)
      IMPLICIT NONE
      INTEGER L1,L2,L3
      LOGICAL TRIANGLE
      INTRINSIC MOD
!     ..
triangle = (l1 >= ABS(l3-l2)) .AND. (l1 <= (l3+l2))  &
    .AND. (MOD((l1+l2+l3),2) == 0)
END FUNCTION triangle
