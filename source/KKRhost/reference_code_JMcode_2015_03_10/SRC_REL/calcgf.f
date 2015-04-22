      SUBROUTINE CALCGF(NK,CGC,GDIA,GMDIA,GOFF,GMOFF,FDIA,FMDIA,FOFF,
     &                  FMOFF,LTAB,LBTAB,KAPTAB,NMUETAB,NMUEMAX,NKMMAX,
     &                  NKMPMAX)
C   ********************************************************************
C   *                                                                  *
C   *   G- AND F-COEFFICIENTS                                          *
C   *                                                                  *
C   *   G(K,K',MUE) =                                                  *
C   *     CGC(K,MUE,2)*CGC(K',MUE,2) - CGC(K,MUE,1)*CGC(K',MUE,1)      *
C   *                                                                  *
C   *   GM = G(-K,-K',MUE)                                             *
C   *                                                                  *
C   *   F(K,K',MUE) =    (MUE-1/2) * CGC(K,MUE,2)*CGC(K',MUE,2)        *
C   *                  + (MUE+1/2) * CGC(K,MUE,1)*CGC(K',MUE,1)        *
C   *                                                                  *
C   *   FM = F(-K,-K',MUE)                                             *
C   *                                                                  *
C   *   ..DIA/..OFF ARE THE ELEMENTS FOR  K=K'/K=-K'-1                 *
C   *   IG NUMBERS THE G'S   COLUMN-WISE STARTING WITH COLUMN  1       *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER NK,NKMMAX,NKMPMAX,NMUEMAX
      REAL*8 CGC(NKMPMAX,2),FDIA(NKMMAX),FMDIA(NKMMAX),FMOFF(NKMMAX),
     &       FOFF(NKMMAX),GDIA(NKMMAX),GMDIA(NKMMAX),GMOFF(NKMMAX),
     &       GOFF(NKMMAX)
      INTEGER KAPTAB(NMUEMAX),LBTAB(NMUEMAX),LTAB(NMUEMAX),
     &        NMUETAB(NMUEMAX)
C
C Local variables
C
      INTEGER IG,IKM1,IKM2,IMKM1,IMKM2,J1P05,J2P05,K1,K2,KAP1,KAP2,L1,
     &        L1BAR,L2,L2BAR,M1,M2,MLDN,MLUP,MUE1M05,MUE2M05
C
C
C     JP05 = J +0.5     MUEP05 = MUE + 0.5
C     IKM  = L*2*(J+1/2) + J + MUE + 1
C
      IKM2 = 0
      DO K2 = 1,NK
         L2 = LTAB(K2)
         KAP2 = KAPTAB(K2)
         L2BAR = LBTAB(K2)
         J2P05 = ABS(KAP2)
         MUE2M05 = -J2P05 - 1
C
         DO M2 = 1,NMUETAB(K2)
            MUE2M05 = MUE2M05 + 1
            IKM2 = IKM2 + 1
            IG = IKM2
            GOFF(IG) = 0.0D0
            GMOFF(IG) = 0.0D0
            FOFF(IG) = 0.0D0
            FMOFF(IG) = 0.0D0
C
            IMKM2 = L2BAR*2*J2P05 + J2P05 + MUE2M05 + 1
C
            IKM1 = 0
            DO K1 = 1,NK
               L1 = LTAB(K1)
               KAP1 = KAPTAB(K1)
               L1BAR = LBTAB(K1)
               J1P05 = ABS(KAP1)
               MUE1M05 = -J1P05 - 1
C
               DO M1 = 1,NMUETAB(K1)
                  MUE1M05 = MUE1M05 + 1
                  IKM1 = IKM1 + 1
                  MLUP = MUE1M05
                  MLDN = MLUP + 1
C
                  IF ( L1.EQ.L2 ) THEN
                     IF ( MUE1M05.EQ.MUE2M05 ) THEN
C
                        IMKM1 = L1BAR*2*J1P05 + J1P05 + MUE1M05 + 1
C
                        IF ( KAP1.EQ.KAP2 ) THEN
C
                           GDIA(IG) = CGC(IKM1,2)*CGC(IKM2,2)
     &                                - CGC(IKM1,1)*CGC(IKM2,1)
                           GMDIA(IG) = CGC(IMKM1,2)*CGC(IMKM2,2)
     &                                 - CGC(IMKM1,1)*CGC(IMKM2,1)
                           FDIA(IG) = MLUP*CGC(IKM1,2)*CGC(IKM2,2)
     &                                + MLDN*CGC(IKM1,1)*CGC(IKM2,1)
                           FMDIA(IG) = MLUP*CGC(IMKM1,2)*CGC(IMKM2,2)
     &                                 + MLDN*CGC(IMKM1,1)*CGC(IMKM2,1)
C
                        ELSE
C
                           GOFF(IG) = CGC(IKM1,2)*CGC(IKM2,2)
     &                                - CGC(IKM1,1)*CGC(IKM2,1)
                           GMOFF(IG) = CGC(IMKM1,2)*CGC(IMKM2,2)
     &                                 - CGC(IMKM1,1)*CGC(IMKM2,1)
                           FOFF(IG) = MLUP*CGC(IKM1,2)*CGC(IKM2,2)
     &                                + MLDN*CGC(IKM1,1)*CGC(IKM2,1)
                           FMOFF(IG) = MLUP*CGC(IMKM1,2)*CGC(IMKM2,2)
     &                                 + MLDN*CGC(IMKM1,1)*CGC(IMKM2,1)
                        END IF
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO
C
C
      END
