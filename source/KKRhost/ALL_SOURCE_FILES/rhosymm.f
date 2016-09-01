c 30.08.96 ***************************************************************
      SUBROUTINE RHOSYMM(LMPOT,NSPIN,NSTART,NEND,RHO2NS,IXIPOL,
     +                 IRWS,IRCUT,IPAN,KSHAPE)
c ************************************************************************
c     symmetrize the charge densities and magnetic moments of
c     atoms which are magnetic 'antisymmetric' 
c     (dependencies in IXIPOL(*))
c
c     p. zahn, aug. 96
c-----------------------------------------------------------------------
      IMPLICIT NONE
C     .. Parameters ..
      include 'inc.p'
      INTEGER LMPOTD, NSPIN
      PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER KSHAPE,LMPOT,NEND,NSTART
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION  RHO2NS(IRMD,LMPOTD,NATYPD,*)
      INTEGER IXIPOL(*),IRCUT(0:IPAND,*),IPAN(*),IRWS(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,IATYP,IATYP1,IRC,IRC1,LM
      DOUBLE PRECISION FAC
C     ..
C     .. External Subroutines ..
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS
C     ..
c ------------------------------------------------------------------------
c
      DO 110 IATYP = NSTART,NEND

        IATYP1 = ABS(IXIPOL(IATYP))

        FAC = 1.D0
        IF (IXIPOL(IATYP).LT.0) FAC = -1.d0

        IF (IATYP1.GE.IATYP) THEN

          write(1337,*) 'Symmetrize atom ',IATYP,' with ',IATYP1,'.'
          IF (KSHAPE.NE.0) THEN
            IRC  = IRCUT(IPAN(IATYP),IATYP)
            IRC1 = IRCUT(IPAN(IATYP1),IATYP1)
          ELSE
            IRC  = IRWS(IATYP)
            IRC1 = IRWS(IATYP1)
          END IF

          IF (IRC.NE.IRC1) THEN
            write(6,*) 'Error in RHOSYMM : ***********************'
            write(6,*) 'Radial mesh of atoms ',iatyp,
     +           ' and ',iatyp1,' are not equal.'
          END IF
          
          DO 10 LM = 1,LMPOT
            DO 20 I = 1,IRC1
              RHO2NS(I,LM,IATYP,1) = ( RHO2NS(I,LM,IATYP,1) + 
     +             RHO2NS(I,LM,IATYP1,1) )/2.d0
              RHO2NS(I,LM,IATYP1,1) = RHO2NS(I,LM,IATYP,1)
              IF (NSPIN.GT.1) THEN
                RHO2NS(I,LM,IATYP,2) = ( RHO2NS(I,LM,IATYP,2) +
     +               FAC*RHO2NS(I,LM,IATYP1,2) )/2.d0
                RHO2NS(I,LM,IATYP1,2) = FAC*RHO2NS(I,LM,IATYP,2)
              END IF
 20         END DO                  ! I = 1,IRC1
 10       END DO                    ! LM = 1,LMPOT

        END IF                      ! (IATYP1.GT.IATYP)
 110  CONTINUE                      ! IATYP = NSTART,NEND

      RETURN

      END                           ! RHOSYMM
