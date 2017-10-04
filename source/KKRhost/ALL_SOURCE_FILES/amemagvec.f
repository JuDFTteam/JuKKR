C*==amemagvec.f    processed by SPAG 6.05Rc at 18:24 on  6 Apr 2003
      SUBROUTINE AMEMAGVEC(IREL,IPRINT,NKM,AMEMVEC,IKMLLIM1,IKMLLIM2,
     &                     IMKMTAB,CGC,NLMAX,NKMMAX,NKMPMAX,NMVECMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the angular matrix elements connected with           *
C   *                                                                  *
C   *   1:       < LAM | sigma(ipol) | LAM' >    spin moment           *
C   *   2:       < LAM |     l(ipol) | LAM' >    orbital moment        *
C   *   3:       < LAM |     T(ipol) | LAM' >    spin dipole moment    *
C   *   4:       < LAM |  B_hf(ipol) | LAM' >    hyperfine field       *
C   *                                                                  *
C   *   ipol= 1,2,3  ==  (+),(-),(z):                                  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
C
C Dummy arguments
C
      INTEGER IPRINT,IREL,NKM,NKMMAX,NKMPMAX,NLMAX,NMVECMAX
      DOUBLE PRECISION AMEMVEC(NKMMAX,NKMMAX,3,NMVECMAX),CGC(NKMPMAX,2)
      INTEGER IKMLLIM1(NKMMAX),IKMLLIM2(NKMMAX),IMKMTAB(NKMMAX)
C
C Local variables
C
      CHARACTER*1 CHPOL(3)
      DOUBLE PRECISION DBLE
      INTEGER I,IKM1,IKM2,IMKM,IMV,IMVEC,IPOL,J1P05,J2P05,K,K1,K2,KAP1,
     &        KAP2,L,L1,L2,LB1,LB2,M2,MSM05,MUE1M05,MUE2M05,NK,NMVEC
      INTEGER IABS,NINT
      CHARACTER*20 STR20
      DOUBLE PRECISION SUM,XJ,XJM,XJP,XM,XYNORM
      CHARACTER*4 TXTMVEC(4)
C
      DATA CHPOL/'+','-','z'/
      DATA TXTMVEC/'spin','orb ','T_z ','B_hf'/
C
      NK = 2*NLMAX - 1
C     XYNORM = DSQRT(2.0D0)
      XYNORM = 2.0D0
C
      CALL RINIT(NKMMAX*NKMMAX*3*NMVECMAX,AMEMVEC)
C
      IF ( IREL.LE.1 ) RETURN
C
C ----------------------------------------------------------------------
C     find the bounding indices  IKMLLIM1  and  IKMLLIM2  for IKM-loops
C     assuming that the matrix elements are diagonal with respect to l
C     this does not hold for B_hf for which there are l-(l+/-2)-terms
C ----------------------------------------------------------------------
C
      I = 0
      DO K = 1,NK
         L = K/2
         XJM = L - 0.5D0
         XJP = L + 0.5D0
         IF ( MOD(K,2).EQ.1 ) THEN
            XJ = L + 0.5D0
         ELSE
            XJ = L - 0.5D0
         END IF
         DO XM = -XJ, + XJ
            I = I + 1
            IKMLLIM1(I) = NINT(L*2*(XJM+0.5D0)+1)
            IKMLLIM2(I) = NINT(L*2*(XJP+0.5D0)+2*XJP+1)
         END DO
      END DO
C
C ----------------------------------------------------------------------
C
      IKM1 = 0
      DO K1 = 1,NK
         L1 = K1/2
         IF ( MOD(K1,2).EQ.0 ) THEN
            KAP1 = L1
            LB1 = L1 - 1
         ELSE
            KAP1 = -L1 - 1
            LB1 = L1 + 1
         END IF
         J1P05 = IABS(KAP1)
C
         DO MUE1M05 = -J1P05,J1P05 - 1
            IKM1 = IKM1 + 1
            IMKM = LB1*2*J1P05 + J1P05 + MUE1M05 + 1
            IMKMTAB(IKM1) = IMKM
C
            IKM2 = 0
            DO K2 = 1,NK
               L2 = K2/2
               IF ( MOD(K2,2).EQ.0 ) THEN
                  KAP2 = L2
                  LB2 = L2 - 1
               ELSE
                  KAP2 = -L2 - 1
                  LB2 = L2 + 1
               END IF
               J2P05 = IABS(KAP2)
C
               DO MUE2M05 = -J2P05,J2P05 - 1
                  IKM2 = IKM2 + 1
C ----------------------------------------------------------------------
                  IF ( (MUE1M05-MUE2M05).EQ.+1 ) THEN
                     AMEMVEC(IKM1,IKM2,1,1) = XYNORM*CGC(IKM1,2)
     &                  *CGC(IKM2,1)
C
                     SUM = 0D0
                     DO MSM05 = -1,0
                        M2 = MUE2M05 - MSM05
                        IF (ABS(M2).LE.L2) SUM = SUM + 
     &                       CGC(IKM1,MSM05+2) * CGC(IKM2,MSM05+2)
     &                       * SQRT(DBLE((L2-M2)*(L2+M2+1)))
                     END DO
                     AMEMVEC(IKM1,IKM2,1,2) = SUM
                  END IF
C
                  IF ( (MUE1M05-MUE2M05).EQ.-1 ) THEN
                     AMEMVEC(IKM1,IKM2,2,1) = XYNORM*CGC(IKM1,1)
     &                  *CGC(IKM2,2)
C
                     SUM = 0D0
                     DO MSM05 = -1,0
                        M2 = MUE2M05 - MSM05
                        IF (ABS(M2).LE.L2) SUM = SUM + 
     &                       CGC(IKM1,MSM05+2) * CGC(IKM2,MSM05+2)
     &                       * SQRT(DBLE((L2+M2)*(L2-M2+1)))
                        
                     END DO
                     AMEMVEC(IKM1,IKM2,2,2) = SUM
                  END IF
C
                  IF ( (MUE1M05-MUE2M05).EQ.0 ) THEN
                     AMEMVEC(IKM1,IKM2,3,1) = CGC(IKM1,2)*CGC(IKM2,2)
     &                  - CGC(IKM1,1)*CGC(IKM2,1)
C
                     SUM = 0D0
                     DO MSM05 = -1,0
                        M2 = MUE2M05 - MSM05
                        SUM = SUM + CGC(IKM1,MSM05+2)*CGC(IKM2,MSM05+2)
     &                        *M2
                     END DO
                     AMEMVEC(IKM1,IKM2,3,2) = SUM
                  END IF
C
C ----------------------------------------------------------------------
               END DO
            END DO
         END DO
      END DO
C
      NMVEC = 2
C ----------------------------------------------------------------------
      IF ( IPRINT.LE.90 ) RETURN
C
      DO IMVEC = 1,NMVEC
         IMV = MIN(IMVEC,2)
         DO IPOL = 1,3
            STR20 = 'A  '//TXTMVEC(IMV)//'  ('//CHPOL(IPOL)//')'
            CALL RMATSTR(STR20,12,AMEMVEC(1,1,IPOL,IMVEC),NKM,NKMMAX,3,
     &                   3,1D-8,6)
         END DO
      END DO
      END
