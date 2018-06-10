SUBROUTINE overlap(result,phi,pz,qz,pqns,acr,dr,lirreg,  &
        ipan,ircut,drdi,irmin,lphi,  &
        ipand,lmaxd,lmmaxd,mmaxd,lmpotd,irmind,irmd)
! **********************************************************************
! *                                                                    *
! * Calculates the overlap integral of test function PHI with regular  *
! * or irregular wavefunction (PZ, QZ: spherical part of the large     *
! * component, PQNS: nonspherical part of the large component) for     *
! * LDA+U.                                                             *
! *                                                                    *
! * THETAS are the shape functions (needed only in case of cell        *
! *        integration)                                                *
! * LPHI is the angular momentum of PHI.                               *
! * LIRREG is true if the irrregular wavefunction is to be used.       *
! * The overlap is given out in array RESULT                           *
! *                                                                    *
! * In the inner region the non-spherical wavefunctions are            *
! * approximated as:                                                   *
! *                                                                    *
! *  * the regular one (ir < irmin = irws-irns) :                      *
! *                                                                    *
! *          pns(ir,lm1,lm2) = pz(ir,l1) * ar(lm1,lm2)                 *
! *                                                                    *
! *   where pz is the regular wavefct of the spherically symmetric     *
! *   part of the potential and ar the alpha matrix (see sub. regns)   *
! *                                                                    *
! *  * the irregular one (ir < irmin) :                                *
! *                                                                    *
! *          qns(ir,lm1,lm2) = pz(ir,l1) * cr(lm1,lm2)                 *
! *                                + qz(ir,l1) * dr(lm1,lm2)           *
! *                                                                    *
! *   where pz is the regular and qz is the irregular                  *
! *   wavefunction of the spherically symmetric part of the            *
! *   potential and cr, dr the matrices calculated at the point irmin  *
! *   (see sub. irwns)                                                 *
! *                                                                    *
! *  Attention: last index of pns,qns is not the spin index but        *
! *             sra index.                                             *
! *                                                                    *
! *  The integrand is convoluted with the shape functions:             *
! *  (commented out now)                                               *
! *    Result = Int conjg(phi_l)                                       *
! *                 sum_{L''L'''} R_L''L' Theta_L''' GAUNT_LL''L'''    *
! *                                                                    *
! *  Finally, the result should then be transformed to complex         *
! *  spherical harmonics basis.                                        *
! *                                                                    *
! *  Querry: Here only large component is used, but was PHI normalised *
! *  according to large + small component? (qldau)                     *
! *                                                                    *
! *                             ph. mavropoulos, juelich, 2002         *
! *                                                                    *
! **********************************************************************
IMPLICIT NONE
!..
!.. Scalar Arguments ..
INTEGER  IPAND,LMAXD,LMMAXD,MMAXD,LMPOTD,IRMIND,IRMD
INTEGER IRMIN,IPAN,LPHI   ! l-value for LDA+U
LOGICAL LIRREG
!..
!.. Array arguments ..
DOUBLE COMPLEX PHI(IRMD),PZ(IRMD,0:LMAXD),QZ(IRMD,0:LMAXD)
DOUBLE COMPLEX PQNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
DOUBLE COMPLEX ACR(LMMAXD,LMMAXD),DR(LMMAXD,LMMAXD)
DOUBLE COMPLEX RESULT(MMAXD,MMAXD)
DOUBLE PRECISION DRDI(IRMD)
INTEGER IRCUT(0:IPAND)
!..
!.. Locals ..
INTEGER LPHISQ,MMAX,IRS1,IRC1
INTEGER LM1,LM2,LM3,MM1,MM2,MM3,IR
DOUBLE PRECISION GAUNT(LMMAXD,LMMAXD,LMPOTD)
DOUBLE COMPLEX WINT(IRMD)
!..
mmax = 2*lphi + 1
lphisq = lphi*lphi
irs1 = ircut(1)
irc1 = ircut(ipan)

CALL rinit(lmmaxd*lmmaxd*lmpotd,gaunt)

DO lm1 = 1,lmmaxd
  gaunt(lm1,lm1,1) = 1.d0
END DO
CALL cinit(mmaxd*mmaxd,result)

! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! Loop over mm1,mm2:
DO mm1 = 1,mmax
  lm1 = lphisq + mm1
  DO mm2 = 1,mmax
    lm2 = lphisq + mm2
! Set up integrand
    CALL cinit(irmd,wint)
! ----------------------------------------------------------------------
    DO mm3 = 1,mmax
      lm3 = lphisq + mm3
      
! -> First inner part (up to irmin, no thetas)
      DO ir = 2,irmin
        wint(ir) = wint(ir) + phi(ir) * pz(ir,lphi) * acr(lm3,lm2)  &
            * gaunt(lm1,lm3,1)
      END DO
      
      IF (lirreg) THEN
        DO ir = 2,irmin
          wint(ir) = wint(ir) + phi(ir) * qz(ir,lphi) * dr(lm3,lm2)  &
              * gaunt(lm1,lm3,1)
        END DO
      END IF
      
! -> Next middle part (from irmin+1 to irc1, no thetas)
      DO ir = irmin+1,irs1
        wint(ir) = wint(ir) + phi(ir) * pqns(lm3,lm2,ir,1)  &
            * gaunt(lm1,lm3,1)
      END DO
      
! -> Finally last part - from irc1+1 to irs1 - with THETAS and proper
!    GAUNTS if we integrate in cell:
!ccc               DO  LM4 = 1,LMPOTD
!ccc                  IF (LMSP(LM4).GT.0) THEN
!ccc                     DO IR = IRS1+1,IRC1
!ccc                        WINT(IR) = WINT(IR)
!ccc     &                           + PHI(IR) * PQNS(LM3,LM2,IR,1)
!ccc     &                           * GAUNT(LM1,LM3,LM4)
!ccc     &                           * THETAS(IR-IRS1,LM4)
!ccc                     ENDDO
!ccc                  END IF
!ccc               ENDDO
      
!    or still without THETAS if we integrate in sphere
      
      DO ir = irs1+1,irc1
        wint(ir) = wint(ir) + phi(ir) * pqns(lm3,lm2,ir,1)  &
            * gaunt(lm1,lm3,1)
      END DO
    END DO
! ----------------------------------------------------------------------
    
    CALL csimpk(wint,result(mm1,mm2),ipan,ircut,drdi(1))
  END DO
END DO
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
END SUBROUTINE overlap
