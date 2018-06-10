SUBROUTINE calcrotmat( nk,irel, alfdeg,betdeg,gamdeg,  &
        rot, fact, nkmmax )
!   ********************************************************************
!   *                                                                  *
!   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
!   *           ( ALFDEG, BETDEG, GAMDEG )                             *
!   *                                                                  *
!   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
!   *            EQS. (4.8), (4.12) AND (4.13)                         *
!   *                                                                  *
!   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
!   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
!   *                                                                  *
!   *   12/11/96  HE  deal with beta = 0                               *
!   ********************************************************************

IMPLICIT NONE

DOUBLE COMPLEX CI, C0
PARAMETER ( CI = (0.0D0,1.0D0), C0 = (0.0D0,0.0D0) )
DOUBLE PRECISION PI
PARAMETER ( PI = 3.141592653589793238462643D0 )   

DOUBLE PRECISION     NUM,MSB05,MSB05SQ,MSB05PW,J,M1,M2,RFAC,DOM, & 
           X, CB05, CB05SQ, ALFDEG, BETDEG, GAMDEG, SUM, CB05PW
DOUBLE PRECISION     FACT(0:100)

INTEGER    S, SLOW, SHIGH, OFF, NKMMAX, IREL, NK, I1, I2, K, L, &
           IM1, IM2, NMUE
DOUBLE COMPLEX EMIM2A, EMIM1G, ROT(NKMMAX,NKMMAX) 

! INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
rfac(x) = fact( nint(x) )

IF( irel == 2 ) CALL errortrap('calcrotmat',12,1)
IF( irel == 3 .AND. MOD(nk,2) == 0) CALL errortrap('CALCROTMAT',13,1)

DO  i2=1,nkmmax
  DO  i1=1,nkmmax
    rot(i1,i2) = c0
  END DO
END DO

cb05   =   DCOS( betdeg*0.5D0*pi/180.0D0 )
cb05sq =   cb05 *  cb05
msb05   = - DSIN( betdeg*0.5D0*pi/180.0D0 )
msb05sq =  msb05 * msb05

off = 0
DO  k=1,nk
  IF( irel < 2 ) THEN
    l = k - 1
    j = l
  ELSE
    l = k/2
    IF( l*2 == k ) THEN
      j = l - 0.5D0
    ELSE
      j = l + 0.5D0
    END IF
  END IF
  
  nmue = nint( 2*j + 1 )
  
  DO  im2 = 1, nmue
    m2 = - j + (im2-1.0D0)
    emim2a = CDEXP( -ci*m2*alfdeg*pi/180.0D0 )
    
    DO  im1 = 1, nmue
      m1 = - j + (im1-1.0D0)
      emim1g = CDEXP( -ci*m1*gamdeg*pi/180.0D0 )
      
      IF( DABS(betdeg) < 1D-8 ) THEN
        IF( im1 == im2 ) THEN
          sum = 1.0D0
        ELSE
          sum = 0.0D0
        END IF
      ELSE
        slow   = MAX(          0, nint(m1-m2) )
        shigh  = MIN( nint(j-m2), nint( j+m1) )
        cb05pw =  cb05**nint(2*j+m1-m2-2*slow    +2)
        msb05pw = msb05**nint(    m2-m1+2*slow    -2)
        dom = (-1.0D0)**(slow-1) *  &
            DSQRT( rfac(j+m1)*rfac(j-m1)*rfac(j+m2)*rfac(j-m2) )
        sum = 0.0D0
        
        DO s=slow,shigh
          dom = -dom
          num =    fact(s) * rfac(j-m2-s) * rfac(j+m1-s) * rfac(m2-m1+s)
          cb05pw =  cb05pw /  cb05sq
          msb05pw = msb05pw * msb05sq
          sum = sum + (dom/num) * cb05pw * msb05pw
        END DO
      END IF
      
      rot(off+im2,off+im1) = emim1g * sum * emim2a
    END DO
    
  END DO
  
  off = off + nmue
END DO

RETURN
END SUBROUTINE calcrotmat
