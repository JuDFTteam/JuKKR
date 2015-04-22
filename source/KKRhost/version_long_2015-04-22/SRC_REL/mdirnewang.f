      SUBROUTINE MDIRNEWANG(IT,NMVEC,MVEVI,MVPHI,MVTET,MVGAM,
     &                      NATYPD,LMAXD,NMVECMAX)
C   ********************************************************************
C   *                                                                  *
C   *  this routine has been build up from the last part of the        *
C   *  original Munich CALCMVEC routine.                               *
C   *  After correcting MVEVI with the Fermi energy value MVEVIEF      *
C   *  (outside this routine) it calculates the new angles of the      *
C   *  LOCAL FRAME quantisation axis with respect to the GLOBAL FRAME  *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C ..
C ..  Parameter definitions
      INTEGER LMAXDLOC
      PARAMETER (LMAXDLOC=8)
C ..
C ..  Scalar Arguments
      INTEGER IT,NMVEC,NATYPD,LMAXD,NMVECMAX
C ..
C ..  Array Arguments
      DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX)
      DOUBLE PRECISION MVPHI(NATYPD,NMVECMAX),MVTET(NATYPD,NMVECMAX),
     &                 MVGAM(NATYPD,NMVECMAX)
C ..
C ..  Local Scalars
      DOUBLE PRECISION MV,MVX,MVXY,MVY,MVZ,PI
      INTEGER I,IMV,ICALL
C ..
C ..  Local Arrays
      DOUBLE PRECISION MVGLO(3,NMVECMAX)
C ..
C ..  Intrinsic Functions
      INTRINSIC ABS,ATAN
C ..
C ..  Data Statements
      DATA ICALL /0/
C ..
C ..  Save Statements
      SAVE ICALL,PI
C ..
      ICALL = ICALL + 1
C=======================================================================
      IF (ICALL.EQ.1) THEN
C
         IF (LMAXD.GT.LMAXDLOC) THEN
            WRITE (6,*) 
            WRITE (6,*) ' Please increase parameter LMAXDLOC to ',LMAXD
            WRITE (6,*) ' in the < MVECGLOBAL > routine.'
            STOP ' < TBKKR2 > '
         END IF
C
         PI = 4.D0 * ATAN(1.D0)
C
      END IF
C=======================================================================
C
      DO IMV = 1,NMVEC
C
         DO I = 1,3
            MVGLO(I,IMV) = DIMAG(MVEVI(IT,I,IMV))
         END DO
C
         MVX = MVGLO(1,IMV)
         MVY = MVGLO(2,IMV)
         MVZ = MVGLO(3,IMV)
C
         MV = SQRT(MVX**2+MVY**2+MVZ**2)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ( MV.LT.1D-8 ) THEN
            MVPHI(IT,IMV) = 0D0
            MVTET(IT,IMV) = 0D0
            MVGAM(IT,IMV) = 0D0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ELSE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            MVXY = SQRT(MVX**2+MVY**2)
C ----------------------------------------------------------------------
            IF ( ABS(MVXY).LT.1D-8 ) THEN
               MVPHI(IT,IMV) = 0D0
C ----------------------------------------------------------------------
            ELSE
C ----------------------------------------------------------------------
               IF ( MVY.GE.0D0 ) THEN
                  MVPHI(IT,IMV) = ACOS(MVX/MVXY)
               ELSE IF ( MVX.LT.0D0 ) THEN
                  MVPHI(IT,IMV) = PI + ACOS(-MVX/MVXY)
               ELSE
                  MVPHI(IT,IMV) = 2*PI - ACOS(MVX/MVXY)
               END IF
               MVPHI(IT,IMV) = MVPHI(IT,IMV)*180D0/PI
               IF ( ABS(MVPHI(IT,IMV)-360.0D0).LT.1D-8 )
     &              MVPHI(IT,IMV) = 0D0
            END IF
C ----------------------------------------------------------------------
            IF (MVPHI(IT,IMV).GE.345.D0) 
     &           MVPHI(IT,IMV) = 360.D0 - MVPHI(IT,IMV)
            MVTET(IT,IMV) = ACOS(MVZ/MV)*180D0/PI
            MVGAM(IT,IMV) = 0D0
         END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END DO
C=======================================================================
      END
