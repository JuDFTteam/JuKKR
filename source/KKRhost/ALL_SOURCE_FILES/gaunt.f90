! ************************************************************************
SUBROUTINE gaunt(lmax,lpot,w,yr,cleb,loflm,icleb,iend,jend,  &
    ncleb,lmaxd,lmgf0d,lmpotd)
! ************************************************************************

!   - fills the array cleb with the gaunt coeffients ,i.e.
!      the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)
!      but only for lm2.le.lm1 and lm3>1
!   - calculate the pointer array jend  to project the indices
!      array cleb with the same lm3,l1,l2 values - because of
!      the special ordering of array cleb only the last index
!      has to be determined .
!     (the parameter n has to be chosen that l1+l2+l3 .lt. 2*n)
!     using gaussian quadrature as given by
!     m. abramowitz and i.a. stegun, handbook of mathematical functions,
!     nbs applied mathematics series 55 (1968), pages 887 and 916
!     m. weinert and e. wimmer
!     northwestern university march 1980

!     an index array -icleb- is used to save storage place .
!     fills the array loflm which is used to determine the
!     l-value of a given lm-value .
!     this subroutine has to be called only once !

!                               b.drittler   november 1987

!     modified gaunt coefficients are als calculated defined by
!     the integral of y(l1,m1)*y(l2,m2)*y(l3,m3)*i**(l2-l1+l3)
!-----------------------------------------------------------------------

!---> attention : ncleb is an empirical factor - it has to be optimized

implicit none
!..
DOUBLE COMPLEX CI
PARAMETER (CI= (0.0D0,1.0D0))
!..
INTEGER LMPOTD,LMGF0D,LMAXD,NCLEB
!..
!.. Scalar Arguments ..
INTEGER IEND,LMAX,LPOT
!..
!.. Array Arguments ..
DOUBLE PRECISION CLEB(NCLEB,2),W(*), &
                  YR(4*LMAXD,0:4*LMAXD,0:4*LMAXD)
INTEGER ICLEB(NCLEB,4),JEND(LMPOTD,0:LMAXD,0:LMAXD),LOFLM(*)
!..
!.. Local Scalars ..
DOUBLE PRECISION CLECG,FACTOR,FCI,S
INTEGER I,J,L,L1,L1P,L2,L2P,L3,LM1,LM2,LM3,LM3P,LMPOT,M,M1,M1A, &
        M1S,M2,M2A,M2S,M3,M3A,M3S
!..
!.. Intrinsic Functions ..
 INTRINSIC ABS,DBLE,MOD,REAL,SIGN
!..
!.. External Subroutines ..
 EXTERNAL RCSTOP
!..

i = 1
DO  l = 0,2*lmax
  DO  m = -l,l
    loflm(i) = l
    i = i + 1
  END DO
END DO

icleb=0
cleb=0D0
IF (lpot == 0) THEN
  iend = 1
  icleb(1,1) = (lmax+1)**2
  icleb(1,3) = 1
END IF

IF (lpot /= 0) THEN
  
!---> set up of the gaunt coefficients with an index field
  
  i = 1
  DO  l3 = 1,lpot
    DO  m3 = -l3,l3
      
      DO  l1 = 0,lmax
        DO  l2 = 0,l1
          
          IF (MOD((l1+l2+l3),2) /= 1 .AND. (l1+l2-l3) >= 0 .AND.  &
                (l1-l2+l3) >= 0 .AND. (l2-l1+l3) >= 0) THEN
            
            fci = DBLE(ci** (l2-l1+l3))
            DO  m1 = -l1,l1
              DO  m2 = -l2,l2
                
!---> store only gaunt coeffients for lm2.le.lm1
                
                lm1 = l1* (l1+1) + m1 + 1
                lm2 = l2* (l2+1) + m2 + 1
                IF (lm2 <= lm1) THEN
                  
                  m1s = SIGN(1,m1)
                  m2s = SIGN(1,m2)
                  m3s = SIGN(1,m3)
                  
                  IF (m1s*m2s*m3s >= 0) THEN
                    
                    m1a = ABS(m1)
                    m2a = ABS(m2)
                    m3a = ABS(m3)
                    
                    factor = 0.0
                    
                    IF (m1a+m2a == m3a) factor = factor +  &
                        REAL(3*M3S+SIGN(1,-M3))/8.0d0
                    IF (m1a-m2a == m3a) factor = factor + REAL(m1s)/4.0D0
                    IF (m2a-m1a == m3a) factor = factor + REAL(m2s)/4.0D0
                    
                    IF (factor /= 0.0) THEN
                      
                      IF (m1s*m2s /= 1 .OR. m2s*m3s /= 1 .OR.  &
                          m1s*m3s /= 1) factor = -factor
                      
                      s = 0.0
                      DO  j = 1,4*lmaxd
                        s = s + w(j)*yr(j,l1,m1a)*yr(j,l2,m2a)* yr(j,l3,m3a)
                      END DO
                      clecg = s*factor
                      IF (ABS(clecg) > 1.d-10) THEN
                        cleb(i,1) = clecg
                        cleb(i,2) = fci*clecg
                        icleb(i,1) = lm1
                        icleb(i,2) = lm2
                        icleb(i,3) = l3* (l3+1) + m3 + 1
                        icleb(i,4) = lm2*lmgf0d - (lm2*lm2-lm2)/2 + lm1 -  &
                            lmgf0d
                        i = i + 1
                      END IF
                      
                    END IF
                    
                  END IF
                  
                END IF
                
              END DO
            END DO
          END IF
          
        END DO
      END DO
    END DO
  END DO
  iend = i - 1
  IF (ncleb < iend) THEN
    WRITE (6,FMT=9000) ncleb,iend
    CALL rcstop('33      ')
    
  ELSE
    
!---> set up of the pointer array jend,use explicitly
!     the ordering of the gaunt coeffients
    
    lmpot = (lpot+1)* (lpot+1)
    DO  l1 = 0,lmax
      DO  l2 = 0,l1
        DO  lm3 = 2,lmpot
          jend(lm3,l1,l2) = 0
        END DO
      END DO
    END DO
    
    lm3 = icleb(1,3)
    l1 = loflm(icleb(1,1))
    l2 = loflm(icleb(1,2))
    
    DO  j = 2,iend
      lm3p = icleb(j,3)
      l1p = loflm(icleb(j,1))
      l2p = loflm(icleb(j,2))
      
      IF (lm3 /= lm3p .OR. l1 /= l1p .OR. l2 /= l2p) THEN
        jend(lm3,l1,l2) = j - 1
        lm3 = lm3p
        l1 = l1p
        l2 = l2p
      END IF
      
    END DO
    jend(lm3,l1,l2) = iend
    
    
  END IF
  
END IF



9000 FORMAT (13X,'error stop in gaunt : dimension of NCLEB = ',i10,  &
    ' too small ',/, 13X,'change NCLEB to ',i6)
END SUBROUTINE gaunt
