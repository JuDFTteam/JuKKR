SUBROUTINE calccgc(ltab,kaptab,nmuetab,cgc,nkmax,nmuemax,nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *   CLEBSCH-GORDON-COEFFICIENTS     CGC(IKM,IS)                    *
!   *                                                                  *
!   *   IKM NUMBERS  CGC  FOR INCREASING  K  AND  MUE                  *
!   *   IKM  = L*2*(J+1/2) + J + MUE + 1                               *
!   *   IS= 1/2  SPIN DOWN/UP                                          *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! Dummy arguments
INTEGER NKMAX,NKMPMAX,NMUEMAX
REAL*8 CGC(NKMPMAX,2)
INTEGER KAPTAB(NMUEMAX),LTAB(NMUEMAX),NMUETAB(NMUEMAX)

! Local variables
INTEGER IKM,K,KAPPA,M
REAL*8 J,L,MUE,TWOLP1

ikm = 0
DO k = 1,(nkmax+1)
  l = ltab(k)
  kappa = kaptab(k)
  j = ABS(kappa) - 0.5D0
  mue = -j - 1.0D0
  twolp1 = 2.0D0*l + 1.0D0
  
  IF ( kappa < 0 ) THEN
    
!     J = L + 1/2
    DO m = 1,nmuetab(k)
      
      mue = mue + 1.0D0
      ikm = ikm + 1
      cgc(ikm,1) = DSQRT((l-mue+0.5D0)/twolp1)
      cgc(ikm,2) = DSQRT((l+mue+0.5D0)/twolp1)
    END DO
  ELSE
!     J = L - 1/2
    DO m = 1,nmuetab(k)
      
      mue = mue + 1.0D0
      ikm = ikm + 1
      cgc(ikm,1) = DSQRT((l+mue+0.5D0)/twolp1)
      cgc(ikm,2) = -DSQRT((l-mue+0.5D0)/twolp1)
      
    END DO
  END IF
  
  
END DO

END SUBROUTINE calccgc
