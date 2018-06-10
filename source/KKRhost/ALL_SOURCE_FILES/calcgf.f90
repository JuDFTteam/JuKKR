SUBROUTINE calcgf(nk,cgc,gdia,gmdia,goff,gmoff,fdia,fmdia,foff,  &
        fmoff,ltab,lbtab,kaptab,nmuetab,nmuemax,nkmmax,  &
        nkmpmax)
!   ********************************************************************
!   *                                                                  *
!   *   G- AND F-COEFFICIENTS                                          *
!   *                                                                  *
!   *   G(K,K',MUE) =                                                  *
!   *     CGC(K,MUE,2)*CGC(K',MUE,2) - CGC(K,MUE,1)*CGC(K',MUE,1)      *
!   *                                                                  *
!   *   GM = G(-K,-K',MUE)                                             *
!   *                                                                  *
!   *   F(K,K',MUE) =    (MUE-1/2) * CGC(K,MUE,2)*CGC(K',MUE,2)        *
!   *                  + (MUE+1/2) * CGC(K,MUE,1)*CGC(K',MUE,1)        *
!   *                                                                  *
!   *   FM = F(-K,-K',MUE)                                             *
!   *                                                                  *
!   *   ..DIA/..OFF ARE THE ELEMENTS FOR  K=K'/K=-K'-1                 *
!   *   IG NUMBERS THE G'S   COLUMN-WISE STARTING WITH COLUMN  1       *
!   *                                                                  *
!   ********************************************************************

IMPLICIT NONE

! Dummy arguments
INTEGER NK,NKMMAX,NKMPMAX,NMUEMAX
REAL*8 CGC(NKMPMAX,2),FDIA(NKMMAX),FMDIA(NKMMAX),FMOFF(NKMMAX), &
       FOFF(NKMMAX),GDIA(NKMMAX),GMDIA(NKMMAX),GMOFF(NKMMAX), &
       GOFF(NKMMAX)
INTEGER KAPTAB(NMUEMAX),LBTAB(NMUEMAX),LTAB(NMUEMAX), &
        NMUETAB(NMUEMAX)

! Local variables
INTEGER IG,IKM1,IKM2,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1, &
        L1BAR,L2,L2BAR,M1,M2,MLDN,MLUP,MUE1M05,MUE2M05


!     JP05 = J +0.5     MUEP05 = MUE + 0.5
!     IKM  = L*2*(J+1/2) + J + MUE + 1

ikm2 = 0
DO k2 = 1,nk
  l2 = ltab(k2)
  kap2 = kaptab(k2)
  l2bar = lbtab(k2)
  j2p05 = ABS(kap2)
  mue2m05 = -j2p05 - 1
  
  DO m2 = 1,nmuetab(k2)
    mue2m05 = mue2m05 + 1
    ikm2 = ikm2 + 1
    ig = ikm2
    goff(ig) = 0.0D0
    gmoff(ig) = 0.0D0
    foff(ig) = 0.0D0
    fmoff(ig) = 0.0D0
    
    imkm2 = l2bar*2*j2p05 + j2p05 + mue2m05 + 1
    
    ikm1 = 0
    DO k1 = 1,nk
      l1 = ltab(k1)
      kap1 = kaptab(k1)
      l1bar = lbtab(k1)
      j1p05 = ABS(kap1)
      mue1m05 = -j1p05 - 1
      
      DO m1 = 1,nmuetab(k1)
        mue1m05 = mue1m05 + 1
        ikm1 = ikm1 + 1
        mlup = mue1m05
        mldn = mlup + 1
        
        IF ( l1 == l2 ) THEN
          IF ( mue1m05 == mue2m05 ) THEN
            
            imkm1 = l1bar*2*j1p05 + j1p05 + mue1m05 + 1
            
            IF ( kap1 == kap2 ) THEN
              
              gdia(ig) = cgc(ikm1,2)*cgc(ikm2,2) - cgc(ikm1,1)*cgc(ikm2,1)
              gmdia(ig) = cgc(imkm1,2)*cgc(imkm2,2)  &
                  - cgc(imkm1,1)*cgc(imkm2,1)
              fdia(ig) = mlup*cgc(ikm1,2)*cgc(ikm2,2)  &
                  + mldn*cgc(ikm1,1)*cgc(ikm2,1)
              fmdia(ig) = mlup*cgc(imkm1,2)*cgc(imkm2,2)  &
                  + mldn*cgc(imkm1,1)*cgc(imkm2,1)
              
            ELSE
              
              goff(ig) = cgc(ikm1,2)*cgc(ikm2,2) - cgc(ikm1,1)*cgc(ikm2,1)
              gmoff(ig) = cgc(imkm1,2)*cgc(imkm2,2)  &
                  - cgc(imkm1,1)*cgc(imkm2,1)
              foff(ig) = mlup*cgc(ikm1,2)*cgc(ikm2,2)  &
                  + mldn*cgc(ikm1,1)*cgc(ikm2,1)
              fmoff(ig) = mlup*cgc(imkm1,2)*cgc(imkm2,2)  &
                  + mldn*cgc(imkm1,1)*cgc(imkm2,1)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO
END DO


END SUBROUTINE calcgf
