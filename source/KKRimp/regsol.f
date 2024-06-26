      !------------------------------------------------------------------------------------
      !> Summary: Calculates the regular solution of the schroedinger equation or in semi relativistic approximation for a spherically averaged potential and given energy
      !> Author: B. Drittler
      !> Calculates the regular solution of the schroedinger equation or in semi relativistic 
      !> approximation for a spherically averaged potential and given energy.
      !> To archieve greater presion the leading power \(r^s\) (in schroedinger case s = l,
      !> in case of sra \(s = \sqrt{ (l^2+l-1) - \frac{4z^2}{c^2} } )\) is analytically separated
      !> from the wavefunction.
      !> The t-matrix has to be determined at the mt radius in case of a mt calculation 
      !> or at the ws radius in case of a ws calculation. Therefore the logarithmic 
      !> derivative is calculated at that point (`=ircut(ipan)` )
      !------------------------------------------------------------------------------------
      !> @note Ph. Mavropoulos March 2003 - Dec 2004, Munich/Juelich: LDA+U included
      !> @endnote
      !------------------------------------------------------------------------------------
      MODULE mod_regsol
      CONTAINS

      !-------------------------------------------------------------------------------
      !> Summary: Calculates the regular solution of the schroedinger equation or in semi relativistic approximation for a spherically averaged potential and given energy
      !> Author: B. Drittler
      !> Category: single-site, old-mesh, lda+u, KKRimp, KKRhost
      !> Deprecated: False 
      !> Calculates the regular solution of the schroedinger equation or in semi relativistic 
      !> approximation for a spherically averaged potential and given energy.
      !> To archieve greater presion the leading power \(r^s\) (in schroedinger case s = l,
      !> in case of sra \(s = \sqrt{ (l^2+l-1) - \frac{4z^2}{c^2} } )\) is analytically separated
      !> from the wavefunction.
      !> The t-matrix has to be determined at the mt radius in case of a mt calculation 
      !> or at the ws radius in case of a ws calculation. Therefore the logarithmic 
      !> derivative is calculated at that point (`=ircut(ipan)` )
      !-------------------------------------------------------------------------------
      !> @note Ph. Mavropoulos March 2003 - Dec 2004, Munich/Juelich: LDA+U included
      !> @endnote
      !-------------------------------------------------------------------------------
      SUBROUTINE REGSOL(CVLIGHT,E,NSRA,DLOGDP,FZ,HAMF,MASS,PZ,DROR,R,
     +                  S,VM2Z,Z,IPAN,IRCUT,IDOLDAU,LOPT,WLDAUAV,CUTOFF,
     +                  IRMD,IPAND,LMAXATOM)
      use :: nrtype, only: dp
      IMPLICIT NONE

C     .. Scalar Arguments ..
      complex(kind=dp), intent(in) :: E
      real(kind=dp), intent(in) :: CVLIGHT !!Speed of light
      real(kind=dp) Z
      real(kind=dp) WLDAUAV
      INTEGER IPAN,IPAND,IRMD,LMAXATOM,IDOLDAU,LOPT
      integer, intent(in) :: NSRA
C     ..
C     .. Array Arguments ..
      complex(kind=dp), intent(out) :: DLOGDP(0:LMAXATOM) !!logarithmic derivate of real wavefunction
      complex(kind=dp), intent(out) :: FZ(IRMD,0:LMAXATOM) 
      complex(kind=dp), intent(out) :: HAMF(IRMD,0:LMAXATOM)
      complex(kind=dp), intent(out) :: MASS(IRMD)
      complex(kind=dp), intent(out) :: PZ(IRMD,0:LMAXATOM)
      real(kind=dp), intent(in) :: DROR(IRMD)
      real(kind=dp), intent(in) :: R(IRMD)
      real(kind=dp), intent(in) :: S(0:LMAXATOM)
      real(kind=dp), intent(in) :: VM2Z(IRMD)
      real(kind=dp), intent(in) :: CUTOFF(IRMD)
      integer, intent(in) :: IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      complex(kind=dp) :: DF_dp,DP_dp,FIP0,FIP1,HAMF1,K1F,K1P,K2F,K2P,
     +   K3F, K3P, K4F,K4P,MASS1,PIP0,PIP1,VME,VMEFAC,VMETR1
      real(kind=dp) DROR1,DRSM1,DRSP1,S1,SM1,SP1,SRAFAC
      INTEGER IP,IR,IRC,IRE,IRS,IRSP1,J,K,L
C     ..
C     .. Local Arrays ..
      complex(kind=dp) A(-1:4),B(0:4),DFDI(-4:0),DPDI(-4:0)
      complex(kind=dp) HAMFLDAU(IRMD)
C     ..

      IF (NSRA.EQ.2) THEN
c
c---> in case of sra  srafac = 1/c - otherwise srafac = 0
c
        SRAFAC = 1.0_dp/CVLIGHT
      ELSE
        SRAFAC = 0.0_dp
      END IF
c
      IRC = IRCUT(IPAN)
c

      DO 10 IR = 2,IRC
        VMETR1 = (VM2Z(IR)-E)*R(IR) - 2.0_dp*Z
        HAMF(IR,0) = VMETR1*DROR(IR)
        MASS(IR) = R(IR) - SRAFAC*SRAFAC*VMETR1
   10 CONTINUE
c
      DO 30 L = 1,LMAXATOM
        DO 20 IR = 7,IRC
          HAMF(IR,L) = real(L*L+L, kind=dp)/MASS(IR)*DROR(IR)+HAMF(IR,0)
   20   CONTINUE
   30 CONTINUE
c


C ======================================================================
C LDA+U
C
C  Account for potential shift in case of LDA+U (averaged over m)
C  by adding the average WLDAUAV to the spherical part of the 
C  potential.
C
      IF ( IDOLDAU.EQ.1.AND.LOPT.GE.0 ) THEN
         S1 = real(LOPT*LOPT+LOPT, kind=dp)
         DO IR = 2,IRC
            VMETR1 = ( VM2Z(IR) - E + WLDAUAV*CUTOFF(IR) )*R(IR) 
            HAMFLDAU(IR) = (VMETR1-2.0_dp*Z)*DROR(IR)
         END DO
C
         DO IR = 7,IRC
            HAMF(IR,LOPT) = S1/MASS(IR)*DROR(IR) + HAMFLDAU(IR)
         END DO
      END IF
C
C LDA+U
C ======================================================================

      DO 40 IR = 2,IRC
        MASS(IR) = MASS(IR)*DROR(IR)
   40 CONTINUE
c
      DO 120 L = 0,LMAXATOM
c
        S1 = S(L)
        SM1 = S1 - 1.0_dp
        SP1 = S1 + 1.0_dp
c
c---> loop over the number of kinks
c
        DO 110 IP = 1,IPAN
c
          IF (IP.EQ.1) THEN
            IRS = 2
            IRE = IRCUT(1)

c
c---> initial values
c
            VME = VM2Z(2) - E
            VMEFAC = 1.0_dp - VME*SRAFAC*SRAFAC
            IF (NSRA.EQ.2 .AND. Z.GT.0.0_dp) THEN
              A(-1) = 0.0_dp
              A(0) = 1.0_dp
              B(0) = CMPLX(SM1*CVLIGHT*CVLIGHT/ (2*Z),0.0_dp, kind=dp)
              DO 50 J = 1,3
                A(J) = (0.0_dp,0._dp)
                B(J) = (0.0_dp,0._dp)
   50         CONTINUE

            ELSE


              A(0) = 0.0_dp
              B(0) = real(L, kind=dp)/VMEFAC
              A(1) = 1.0_dp
              DO 60 J = 2,4
                A(J) = (VME*VMEFAC*A(J-2)-2.0_dp*Z*A(J-1))/
     +                 real((J-1)* (J+2*L), kind=dp)
                B(J-1) = real(L+J-1, kind=dp)*A(J)/VMEFAC
   60         CONTINUE

            END IF
c
            K = -4
c
c---> power series near origin
c
            DO 80 IR = 2,6
              PIP0 = A(3)
              DP_dp = 3.0_dp*A(3)
              FIP0 = B(3)
              DF_dp = 3.0_dp*B(3)
              DO 70 J = 2,0,-1
                PIP0 = A(J) + PIP0*R(IR)
                DP_dp = real(J, kind=dp)*A(J) + DP_dp*R(IR)
                FIP0 = B(J) + FIP0*R(IR)
                DF_dp = real(J, kind=dp)*B(J) + DF_dp*R(IR)
   70         CONTINUE
c
              PZ(IR,L) = PIP0
              FZ(IR,L) = FIP0
              DPDI(K) = DP_dp*DROR(IR)
              DFDI(K) = DF_dp*DROR(IR)
c
              K = K + 1
   80       CONTINUE

          ELSE
c
c---> runge kutta step to restart algorithm
c
            IRS = IRCUT(IP-1) + 1
            IRE = IRCUT(IP)
            IRSP1 = IRS + 1
            PIP0 = PZ(IRS,L)
            FIP0 = FZ(IRS,L)
            DRSP1 = DROR(IRS)*SP1
            DRSM1 = DROR(IRS)*SM1
            DPDI(-4) = MASS(IRS)*FIP0 - DRSM1*PIP0
            DFDI(-4) = HAMF(IRS,L)*PIP0 - DRSP1*FIP0
c
c---> first step - 4 point runge kutta with interpolation
c
            K1P = DPDI(-4)
            K1F = DFDI(-4)
c
            DROR1 = (3.0_dp*DROR(IRS+3)-15.0_dp*DROR(IRS+2)+
     +              45.0_dp*DROR(IRSP1)+15.0_dp*DROR(IRS))/48.0_dp
            DRSP1 = DROR1*SP1
            DRSM1 = DROR1*SM1
            MASS1 = (3.0_dp*MASS(IRS+3)-15.0_dp*MASS(IRS+2)+
     +              45.0_dp*MASS(IRSP1)+15.0_dp*MASS(IRS))/48.0_dp
            HAMF1 = (3.0_dp*HAMF(IRS+3,L)-15.0_dp*HAMF(IRS+2,L)+
     +              45.0_dp*HAMF(IRSP1,L)+15.0_dp*HAMF(IRS,L))/48.0_dp
            K2P = MASS1* (FIP0+0.5_dp*K1F) - DRSM1* (PIP0+0.5_dp*K1P)
            K2F = HAMF1* (PIP0+0.5_dp*K1P) - DRSP1* (FIP0+0.5_dp*K1F)
            K3P = MASS1* (FIP0+0.5_dp*K2F) - DRSM1* (PIP0+0.5_dp*K2P)
            K3F = HAMF1* (PIP0+0.5_dp*K2P) - DRSP1* (FIP0+0.5_dp*K2F)
c
            DRSP1 = DROR(IRSP1)*SP1
            DRSM1 = DROR(IRSP1)*SM1
            K4P = MASS(IRSP1)* (FIP0+K3F) - DRSM1* (PIP0+K3P)
            K4F = HAMF(IRSP1,L)* (PIP0+K3P) - DRSP1* (FIP0+K3F)
            PIP0 = PIP0 + (K1P+2.0_dp* (K2P+K3P)+K4P)/6.0_dp
            FIP0 = FIP0 + (K1F+2.0_dp* (K2F+K3F)+K4F)/6.0_dp
c
            PZ(IRSP1,L) = PIP0
            FZ(IRSP1,L) = FIP0
            DPDI(-3) = MASS(IRSP1)*FIP0 - DRSM1*PIP0
            DFDI(-3) = HAMF(IRSP1,L)*PIP0 - DRSP1*FIP0
c
            K = -2
c
c---> 4 point runge kutta with h = i+2 - i
c
            DO 90 IR = IRS + 2,IRS + 4
              PIP0 = PZ(IR-2,L)
              FIP0 = FZ(IR-2,L)
              K1P = DPDI(K-2)
              K1F = DFDI(K-2)
              K2P = MASS(IR-1)* (FIP0+K1F) - DRSM1* (PIP0+K1P)
              K2F = HAMF(IR-1,L)* (PIP0+K1P) - DRSP1* (FIP0+K1F)
              K3P = MASS(IR-1)* (FIP0+K2F) - DRSM1* (PIP0+K2P)
              K3F = HAMF(IR-1,L)* (PIP0+K2P) - DRSP1* (FIP0+K2F)
c
              DRSP1 = DROR(IR)*SP1
              DRSM1 = DROR(IR)*SM1
c
              K4P = MASS(IR)*(FIP0+2.0_dp*K3F)-DRSM1*(PIP0+2.0_dp*K3P)
              K4F = HAMF(IR,L)* (PIP0+2.0_dp*K3P) -
     +              DRSP1* (FIP0+2.0_dp*K3F)
              PIP0 = PIP0 + (K1P+2.0_dp* (K2P+K3P)+K4P)/3.0_dp
              FIP0 = FIP0 + (K1F+2.0_dp* (K2F+K3F)+K4F)/3.0_dp
c
              PZ(IR,L) = PIP0
              FZ(IR,L) = FIP0
              DPDI(K) = MASS(IR)*FIP0 - DRSM1*PIP0
              DFDI(K) = HAMF(IR,L)*PIP0 - DRSP1*FIP0
              K = K + 1
   90       CONTINUE
          END IF
c
          DO 100 IR = IRS + 5,IRE
            DRSP1 = DROR(IR)*SP1
            DRSM1 = DROR(IR)*SM1
c
c---> predictor : 5 point adams - bashforth
c
            PIP1 = PIP0 + (1901.0_dp*DPDI(0)-2774.0_dp*DPDI(-1)+
     +             2616.0_dp*DPDI(-2)-1274.0_dp*DPDI(-3)+
     +             251.0_dp*DPDI(-4))/720.0_dp
            FIP1 = FIP0 + (1901.0_dp*DFDI(0)-2774.0_dp*DFDI(-1)+
     +             2616.0_dp*DFDI(-2)-1274.0_dp*DFDI(-3)+
     +             251.0_dp*DFDI(-4))/720.0_dp
c
            DPDI(-4) = DPDI(-3)
            DPDI(-3) = DPDI(-2)
            DPDI(-2) = DPDI(-1)
            DPDI(-1) = DPDI(0)
            DFDI(-4) = DFDI(-3)
            DFDI(-3) = DFDI(-2)
            DFDI(-2) = DFDI(-1)
            DFDI(-1) = DFDI(0)
c
            DPDI(0) = MASS(IR)*FIP1 - DRSM1*PIP1
            DFDI(0) = HAMF(IR,L)*PIP1 - DRSP1*FIP1
c
c---> corrector : 5 point adams - moulton
c
            PIP0 = PIP0 + (251.0_dp*DPDI(0)+646.0_dp*DPDI(-1)-
     +           264.0_dp*DPDI(-2)+106.0_dp*DPDI(-3)-19.0_dp*DPDI(-4))/
     +           720.0_dp
            FIP0 = FIP0 + (251.0_dp*DFDI(0)+646.0_dp*DFDI(-1)-
     +           264.0_dp*DFDI(-2)+106.0_dp*DFDI(-3)-19.0_dp*DFDI(-4))/
     +           720.0_dp
c
            PZ(IR,L) = PIP0
            FZ(IR,L) = FIP0
            DPDI(0) = MASS(IR)*FIP0 - DRSM1*PIP0
            DFDI(0) = HAMF(IR,L)*PIP0 - DRSP1*FIP0
  100     CONTINUE
c
c---> remember that the r - mesh contains the kinks two times
c     store the values of pz and fz to restart the algorithm
c
          IF (IP.NE.IPAN) THEN
            PZ(IRE+1,L) = PIP0
            FZ(IRE+1,L) = FIP0
          END IF

  110   CONTINUE

c
c---> logarithmic derivate of real wavefunction ( r**s *pz / r)
c
        DLOGDP(L) = (DPDI(0)/ (PIP0*DROR(IRC))+SM1)/R(IRC)
  120 CONTINUE


      END SUBROUTINE
      END MODULE mod_regsol
