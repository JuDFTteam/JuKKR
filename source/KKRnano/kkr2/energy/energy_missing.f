C>     calculates some missing energy terms.
C>
C>     @author E. Rabel
C>
C>    The correction is:
C>    \sum_M (Z * \int_{r_{MT}}^{r_{BS}} rho_L \Theta_L r dr
C>            -\delta_{L,(0,0)} Z**2 / R)
C>
      SUBROUTINE energy_missing(e_miss,LPOT,RHO2NS,Z,R,DRDI,
     +                          IRCUT,IPAN,IFUNM,
     +                          THETAS,LMSP,
     &                          irmd, irid, nfund, ipand)

      IMPLICIT NONE

      double precision, intent(out) :: e_miss(0:LPOT)

      INTEGER irmd
      INTEGER irid
      INTEGER nfund
      INTEGER ipand

      INTEGER LPOT

      DOUBLE PRECISION DRDI(IRMD),
     &                 R(IRMD),
     &                 RHO2NS(IRMD,(LPOT + 1)**2,2),
     &                 THETAS(IRID,NFUND)

      DOUBLE PRECISION Z

      INTEGER IFUNM(*),IPAN,
     &        IRCUT(0:IPAND),LMSP(*)

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION RFPI,RHOSP
      INTEGER IFUN,IPAN1,IR,IRC1,IRH,IRS1,L,LM,M
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION ER(IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL SIMPK
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SQRT
C
      RFPI = SQRT(16.D0*ATAN(1.0D0))
c
      IPAN1 = IPAN
      IRS1 = IRCUT(1)
      IRC1 = IRCUT(IPAN1)

      e_miss = 0.0d0

C term Z * \int_{r_{MT}}^{r_{BS}} rho_L \Theta_L r dr
C      (summed over M but not L)

      DO 80 L = 0,LPOT

          ER = 0.0d0

          DO 60 M = -L,L
            LM = L*L + L + M + 1

              IF (LMSP(LM).GT.0) THEN
                IFUN = IFUNM(LM)

                  DO 40 IR = IRS1 + 1,IRC1
                    IRH = IR - IRS1
                    RHOSP = RHO2NS(IR,LM,1)
                    ER(IR) = ER(IR) + Z * RHOSP*THETAS(IRH,IFUN) / R(IR)
   40             CONTINUE
              END IF   

   60     CONTINUE

c--->     now integrate
c
        CALL SIMPK(ER,e_miss(L),IPAN1,IRCUT,DRDI)

   80   CONTINUE                    ! L = 0,LPOT

C     Additional term -\delta_{L,(0,0)} Z**2 / R  where R is the
C     reference radius chosen to calculate the gen. Madelung potential
C     Here it is equal to the muffin-tin radius!!!

        e_miss(0) = e_miss(0) - Z**2 / R(IRS1)

      END SUBROUTINE
