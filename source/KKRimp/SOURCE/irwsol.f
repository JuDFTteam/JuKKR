!------------------------------------------------------------------------------------
!> Summary:  Calculates the irregular solution of the schroedinger equation 
!> Author: B. Drittler Nov. 1989
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!>
!> One can write Latex comments like this \(i\hbar\frac{\partial \psi}{\partial t}=-\mathcal{H}\psi\)
!> or add labeled equations using the standard latex way
!> \begin{equation}
!> \mathbf{A} \mathbf{v}= \eta\mathbf{v}
!> \end{equation}
!> **FORd** also accepts markdown style so you can _write with style_ 
!> 
!> **IMPORTANT**
!> The JM-KKR follows the coding conventions noted in this example, one that is
!> not obvious is that **each level of indentation consists of two spaces**. Please keep this:
!> _do it for the children_.
!> So please keep the conventions.
! These are special boxes for ford, notice that this comment does not appear in the html file.
! These boxes contain important information and should be added when necessary. ALWAYS remember to close the box
! BEFORE oppening a new one or they will be nested.
!------------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!------------------------------------------------------------------------------------



      MODULE mod_IRWSOL

      CONTAINS
!-------------------------------------------------------------------------------
!> Summary:   calculates the irregular solution of the schroedinger equation or
!>    in semi relativistic approximation for a spherically averaged
!>    potential and given energy . to achieve greater precision the
!>    leading power r**-s ( in schroedinger case s = l , in case of sra
!>    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
!>    from the wavefunction .
!>
!>
!>  the differential equation is solved with a 5 point adams - bashforth
!>    and adams - moulton predictor corrector method integrating
!>   inwards and extended for potentials with kinks
!> Author: B. Drittler   Nov. 1989
!> Category: Wavefunction, physical-observables, kkrimp
!> Deprecated: False ! This needs to be set to True for deprecated subroutines
!> A More detailed explanation with the math, concepts, etc necessary to understand the routine
!-------------------------------------------------------------------------------
!> @note Notes on the code
!> @endnote
!> @todo things that must be checked
!> @endtodo
!> @warning Important precautions
!> @endwarning
!> @bug If nasty things are found
!> @endbug
!-------------------------------------------------------------------------------



      SUBROUTINE IRWSOL(EK,FZ,HAMF,MASS,PZ,QZ,SZ,DROR,S,IPAN,IRCUT,
     +                    IRMD,IPAND,LMAXATOM)
      IMPLICIT NONE
c-----------------------------------------------------------------------
c  calculates the irregular solution of the schroedinger equation or
c    in semi relativistic approximation for a spherically averaged
c    potential and given energy . to achieve greater precision the
c    leading power r**-s ( in schroedinger case s = l , in case of sra
c    s = sqrt( (l*l+l-1) - 4*z*z/c/c ) ) is analytically separated
c    from the wavefunction .
c
c
c  the differential equation is solved with a 5 point adams - bashforth
c    and adams - moulton predictor corrector method integrating
c    inwards and extended for potentials with kinks
c
c
c                                               b.drittler   nov.1989
c-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE COMPLEX EK
      INTEGER IPAN,IPAND,IRMD,LMAXATOM
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX FZ(IRMD,0:LMAXATOM),HAMF(IRMD,0:LMAXATOM),
     +               MASS(IRMD),
     +               PZ(IRMD,0:LMAXATOM),QZ(IRMD,0:LMAXATOM),
     +               SZ(IRMD,0:LMAXATOM)
      DOUBLE PRECISION DROR(IRMD),S(0:LMAXATOM)
      INTEGER IRCUT(0:IPAND)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX HAMF1,K1F,K1P,K2F,K2P,K3F,K3P,K4F,K4P,MASS1,QIM0,
     +               QIM1,SIM0,SIM1
      DOUBLE PRECISION DROR1,DRSM1,DRSP1,S1,SM1,SP1
      INTEGER IP,IR,IRE,IRS,IRWSK,K,L
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX DQDI(0:4),DSDI(0:4)
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
CTIMO      IRWSK = IRCUT(1)/30
      IRWSK = 41
      IRWSK = IRCUT(1)/9
c
      DO 60 L = 0,LMAXATOM
        S1 = S(L)
        SM1 = S1 - 1.0D0
        SP1 = S1 + 1.0D0
c
c---> loop over kinks
c
        DO 30 IP = IPAN,1,-1
          IRS = IRCUT(IP)
          IRE = MAX(IRCUT(IP-1),IRWSK) + 1
          DRSP1 = DROR(IRS)*SP1
          DRSM1 = DROR(IRS)*SM1
          QIM0 = QZ(IRS,L)
          SIM0 = SZ(IRS,L)
          DQDI(4) = MASS(IRS)*SIM0 + DRSP1*QIM0
          DSDI(4) = HAMF(IRS,L)*QIM0 + DRSM1*SIM0
c
c---> start algorithm - 4 point runge kutta with interpolation
c
          K1P = DQDI(4)
          K1F = DSDI(4)
c
          DROR1 = (3.0D0*DROR(IRS-3)-15.0D0*DROR(IRS-2)+
     +            45.0D0*DROR(IRS-1)+15.0D0*DROR(IRS))/48.0D0
          DRSP1 = DROR1*SP1
          DRSM1 = DROR1*SM1
          MASS1 = (3.0D0*MASS(IRS-3)-15.0D0*MASS(IRS-2)+
     +            45.0D0*MASS(IRS-1)+15.0D0*MASS(IRS))/48.0D0
          HAMF1 = (3.0D0*HAMF(IRS-3,L)-15.0D0*HAMF(IRS-2,L)+
     +            45.0D0*HAMF(IRS-1,L)+15.0D0*HAMF(IRS,L))/48.0D0
          K2P = MASS1* (SIM0-0.5D0*K1F) + DRSP1* (QIM0-0.5D0*K1P)
          K2F = HAMF1* (QIM0-0.5D0*K1P) + DRSM1* (SIM0-0.5D0*K1F)
          K3P = MASS1* (SIM0-0.5D0*K2F) + DRSP1* (QIM0-0.5D0*K2P)
          K3F = HAMF1* (QIM0-0.5D0*K2P) + DRSM1* (SIM0-0.5D0*K2F)
          DRSP1 = DROR(IRS-1)*SP1
          DRSM1 = DROR(IRS-1)*SM1
          K4P = MASS(IRS-1)* (SIM0-K3F) + DRSP1* (QIM0-K3P)
          K4F = HAMF(IRS-1,L)* (QIM0-K3P) + DRSM1* (SIM0-K3F)
          QIM0 = QIM0 - (K1P+2.0D0* (K2P+K3P)+K4P)/6.0D0
          SIM0 = SIM0 - (K1F+2.0D0* (K2F+K3F)+K4F)/6.0D0
          QZ(IRS-1,L) = QIM0
          SZ(IRS-1,L) = SIM0
          DQDI(3) = MASS(IRS-1)*SIM0 + DRSP1*QIM0
          DSDI(3) = HAMF(IRS-1,L)*QIM0 + DRSM1*SIM0
c
          K = 2
c
c---> 4 point runge kutta with h = i+2 - 1
c
          DO 10 IR = IRS - 2,IRS - 4,-1
            QIM0 = QZ(IR+2,L)
            SIM0 = SZ(IR+2,L)
            K1P = DQDI(K+2)
            K1F = DSDI(K+2)
            K2P = MASS(IR+1)* (SIM0-K1F) + DRSP1* (QIM0-K1P)
            K2F = HAMF(IR+1,L)* (QIM0-K1P) + DRSM1* (SIM0-K1F)
            K3P = MASS(IR+1)* (SIM0-K2F) + DRSP1* (QIM0-K2P)
            K3F = HAMF(IR+1,L)* (QIM0-K2P) + DRSM1* (SIM0-K2F)
c
            DRSP1 = DROR(IR)*SP1
            DRSM1 = DROR(IR)*SM1
c
            K4P = MASS(IR)* (SIM0-2.0D0*K3F) + DRSP1* (QIM0-2.0D0*K3P)
            K4F = HAMF(IR,L)* (QIM0-2.0D0*K3P) + DRSM1* (SIM0-2.0D0*K3F)
            QIM0 = QIM0 - (K1P+2.0D0* (K2P+K3P)+K4P)/3.0D0
            SIM0 = SIM0 - (K1F+2.0D0* (K2F+K3F)+K4F)/3.0D0
            QZ(IR,L) = QIM0
            SZ(IR,L) = SIM0
            DQDI(K) = MASS(IR)*SIM0 + DRSP1*QIM0
            DSDI(K) = HAMF(IR,L)*QIM0 + DRSM1*SIM0
            K = K - 1
   10     CONTINUE
c
          DO 20 IR = IRS - 5,IRE,-1
c
c---> predictor : 5 point adams - bashforth
c
            QIM1 = QIM0 - (1901.0D0*DQDI(0)-2774.0D0*DQDI(1)+
     +             2616.0D0*DQDI(2)-1274.0D0*DQDI(3)+251.0D0*DQDI(4))/
     +             720.0D0
            SIM1 = SIM0 - (1901.0D0*DSDI(0)-2774.0D0*DSDI(1)+
     +             2616.0D0*DSDI(2)-1274.0D0*DSDI(3)+251.0D0*DSDI(4))/
     +             720.0D0
c
            DQDI(4) = DQDI(3)
            DQDI(3) = DQDI(2)
            DQDI(2) = DQDI(1)
            DQDI(1) = DQDI(0)
            DSDI(4) = DSDI(3)
            DSDI(3) = DSDI(2)
            DSDI(2) = DSDI(1)
            DSDI(1) = DSDI(0)
c
            DRSP1 = DROR(IR)*SP1
            DRSM1 = DROR(IR)*SM1
c
            DQDI(0) = MASS(IR)*SIM1 + DRSP1*QIM1
            DSDI(0) = HAMF(IR,L)*QIM1 + DRSM1*SIM1
c
c---> corrector : 5 point adams - moulton
c
            QIM0 = QIM0 - (251.0D0*DQDI(0)+646.0D0*DQDI(1)-
     +             264.0D0*DQDI(2)+106.0D0*DQDI(3)-19.0D0*DQDI(4))/
     +             720.0D0
            SIM0 = SIM0 - (251.0D0*DSDI(0)+646.0D0*DSDI(1)-
     +             264.0D0*DSDI(2)+106.0D0*DSDI(3)-19.0D0*DSDI(4))/
     +             720.0D0

            DQDI(0) = MASS(IR)*SIM0 + DRSP1*QIM0
            DSDI(0) = HAMF(IR,L)*QIM0 + DRSM1*SIM0

            QZ(IR,L) = QIM0
            SZ(IR,L) = SIM0
   20     CONTINUE
c
          IF (IP.NE.1) THEN
            QZ(IRE-1,L) = QIM0
            SZ(IRE-1,L) = SIM0
          END IF

   30   CONTINUE

c
c---> use wronski relation near origin
c
        DO 40 IR = IRWSK,2,-1
c
c---> 2 point corrector - predictor
c
          QIM1 = QIM0 - 1.5D0*DQDI(0) + 0.5D0*DQDI(1)
c
          DQDI(1) = DQDI(0)
          DRSP1 = DROR(IR)*SP1
c
          DQDI(0) = MASS(IR)* (1.0D0/EK+QIM1*FZ(IR,L))/PZ(IR,L) +
     +              DRSP1*QIM1

          QIM0 = QIM0 - 0.5D0*DQDI(0) - 0.5D0*DQDI(1)

          DQDI(0) = MASS(IR)* (1.0D0/EK+QIM0*FZ(IR,L))/PZ(IR,L) +
     +              DRSP1*QIM0

          QZ(IR,L) = QIM0
   40   CONTINUE
c
        DO 50 IR = IRWSK,2,-1
          SZ(IR,L) = (1.0D0/EK+QZ(IR,L)*FZ(IR,L))/PZ(IR,L)
   50   CONTINUE
   60 CONTINUE

      END SUBROUTINE
      END MODULE mod_IRWSOL
