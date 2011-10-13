      SUBROUTINE LDAUOVRLP(RESULT,PHI,PZ,QZ,PQNS,ACR,DR,LIRREG,
     &                     IPAN,IRCUT,DRDI,IRMIN,LPHI,
     &                     IPAND,LMAXD,LMMAXD,MMAXD,LMPOTD,IRMIND,IRMD)
C **********************************************************************
C *                                                                    *
C * Calculates the overlap integral of test function PHI with regular  *
C * or irregular wavefunction (PZ, QZ: spherical part of the large     *
C * component, PQNS: nonspherical part of the large component) for     *
C * LDA+U.                                                             *
C *                                                                    *
C * THETAS are the shape functions (needed only in case of cell        *
C *        integration)                                                *
C * LPHI is the angular momentum of PHI.                               *
C * LIRREG is true if the irrregular wavefunction is to be used.       *
C * The overlap is given out in array RESULT                           *
C *                                                                    *
C * In the inner region the non-spherical wavefunctions are            *
C * approximated as:                                                   *
C *                                                                    *
C *  * the regular one (ir < irmin = irws-irns) :                      *
C *                                                                    *
C *          pns(ir,lm1,lm2) = pz(ir,l1) * ar(lm1,lm2)                 *
C *                                                                    *
C *   where pz is the regular wavefct of the spherically symmetric     *
C *   part of the potential and ar the alpha matrix (see sub. regns)   *
C *                                                                    *
C *  * the irregular one (ir < irmin) :                                *
C *                                                                    *
C *          qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)                 *
C *                                + qz(ir,l1) * dr(lm1,lm2)           *
C *                                                                    *
C *   where pz is the regular and qz is the irregular                  *
C *   wavefunction of the spherically symmetric part of the            *
C *   potential and cr, dr the matrices calculated at the point irmin  *
C *   (see sub. irwns)                                                 *
C *                                                                    *
C *  Attention: last index of pns,qns is not the spin index but        *
C *             sra index.                                             *
C *                                                                    *
C *  The integrand is convoluted with the shape functions:             *
C *  (commented out now)                                               *
C *    Result = Int conjg(phi_l)                                       *
C *                 sum_{L''L'''} R_L''L' Theta_L''' GAUNT_LL''L'''    *
C *                                                                    *
C *  Finally, the result should then be transformed to complex         *
C *  spherical harmonics basis.                                        *
C *                                                                    *
C *  Querry: Here only large component is used, but was PHI normalised *
C *  according to large + small component? (qldau)                     *
C *                                                                    *
C *                             ph. mavropoulos, juelich, 2002         *
C *                                                                    *
C **********************************************************************
      IMPLICIT NONE
C     ..
C     .. Scalar Arguments ..
      INTEGER  IPAND,LMAXD,LMMAXD,MMAXD,LMPOTD,IRMIND,IRMD
      INTEGER IRMIN,IPAN,LPHI   ! l-value for LDA+U
      LOGICAL LIRREG
C     ..
C     .. Array arguments ..
      DOUBLE COMPLEX PHI(IRMD),PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD)
      DOUBLE COMPLEX PQNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
      DOUBLE COMPLEX ACR(LMMAXD,LMMAXD),DR(LMMAXD,LMMAXD)
      DOUBLE COMPLEX RESULT(MMAXD,MMAXD)
      DOUBLE PRECISION DRDI(IRMD)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Locals ..
      INTEGER LPHISQ,MMAX,IRS1,IRC1
      INTEGER LM1,LM2,LM3,MM1,MM2,MM3,IR
      DOUBLE PRECISION GAUNT(LMMAXD,LMMAXD,LMPOTD)
      DOUBLE COMPLEX WINT(IRMD)
C     ..
      MMAX = 2*LPHI + 1
      LPHISQ = LPHI*LPHI
      IRS1 = IRCUT(1)
      IRC1 = IRCUT(IPAN)
C      
      CALL RINIT(LMMAXD*LMMAXD*LMPOTD,GAUNT)
C
      DO LM1 = 1,LMMAXD
         GAUNT(LM1,LM1,1) = 1.D0
      ENDDO
      CALL CINIT(MMAXD*MMAXD,RESULT)
C
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
C Loop over mm1,mm2:
      DO MM1 = 1,MMAX
         LM1 = LPHISQ + MM1
         DO MM2 = 1,MMAX
            LM2 = LPHISQ + MM2
C Set up integrand
            CALL CINIT(IRMD,WINT)
C ----------------------------------------------------------------------
            DO MM3 = 1,MMAX
               LM3 = LPHISQ + MM3
C              
C -> First inner part (up to irmin, no thetas)
               DO IR = 2,IRMIN
                  WINT(IR) = WINT(IR) + 
     &                     PHI(IR) * PZ(IR,LPHI) * ACR(LM3,LM2) 
     &                     * GAUNT(LM1,LM3,1)
               ENDDO
C
               IF (LIRREG) THEN
                  DO IR = 2,IRMIN
                     WINT(IR) = WINT(IR) + 
     &                        PHI(IR) * QZ(IR,LPHI) * DR(LM3,LM2) 
     &                        * GAUNT(LM1,LM3,1)      
                  ENDDO
               ENDIF
C
C -> Next middle part (from irmin+1 to irc1, no thetas)
               DO IR = IRMIN+1,IRS1
                  WINT(IR) = WINT(IR) 
     &                     + PHI(IR) * PQNS(LM3,LM2,IR,1)
     &                     * GAUNT(LM1,LM3,1)       
               ENDDO
C
C -> Finally last part - from irc1+1 to irs1 - with THETAS and proper
C    GAUNTS if we integrate in cell:
Cccc               DO  LM4 = 1,LMPOTD
Cccc                  IF (LMSP(LM4).GT.0) THEN
Cccc                     DO IR = IRS1+1,IRC1
Cccc                        WINT(IR) = WINT(IR) 
Cccc     &                           + PHI(IR) * PQNS(LM3,LM2,IR,1)
Cccc     &                           * GAUNT(LM1,LM3,LM4) 
Cccc     &                           * THETAS(IR-IRS1,LM4)
Cccc                     ENDDO
Cccc                  END IF
Cccc               ENDDO
C
C    or still without THETAS if we integrate in sphere
C
               DO IR = IRS1+1,IRC1
                  WINT(IR) = WINT(IR) 
     &                     + PHI(IR) * PQNS(LM3,LM2,IR,1)
     &                     * GAUNT(LM1,LM3,1)       
               ENDDO
            ENDDO
C ----------------------------------------------------------------------
C
            CALL CSIMPK(WINT,RESULT(MM1,MM2),IPAN,IRCUT,DRDI(1))
         ENDDO                 
      ENDDO                    
C AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      END
