      SUBROUTINE MVECGLOBAL(IT,IQ,NATYP,QMPHI,QMTET,
     &                      MVEVI,MVEVIL,MVEVIEF,
     &                      NATYPD,LMAXD,NMVECMAX)
C   ********************************************************************
C   *                                                                  *
C   *  this routine has been build up from the last part of the        *
C   *  original Munich CALCMVEC routine.                               *
C   *  on exit, MVEVI,MVEVIL,MVEVIEF are in the CARTESIAN GLOBAL       *
C   *                                           frame of reference     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C ..
C ..  Parameter definitions
      INTEGER LMAXDLOC
      PARAMETER (LMAXDLOC=8)
      COMPLEX*16 CI,CZERO
      PARAMETER (CI = (0.0D0,1.0D0), CZERO = (0.0D0,0.0D0))
C ..
C ..  Scalar Arguments
      INTEGER IT,IQ,NATYP,NATYPD,LMAXD,NMVECMAX
      DOUBLE PRECISION QMPHI,QMTET
C ..
C ..  Array Arguments
      DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX),
     &           MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX) 
      DOUBLE COMPLEX MVEVIEF(NATYPD,3,NMVECMAX)
C ..
C ..  Local Scalars
      INTEGER ICALL,I,J,K,L,IMV,NMVEC
      DOUBLE COMPLEX CS
      DOUBLE COMPLEX AMIN,APLS
      DOUBLE PRECISION PI,MV,MVX,MVXY,MVY,MVZ,WSQ2
C ..
C ..  Local Arrays
      DOUBLE COMPLEX USC(3,3),DROT4(4,4),W3X3(3,3)
      DOUBLE COMPLEX MVG(3,NMVECMAX),MVGEF(3,NMVECMAX)
      DOUBLE COMPLEX MVGL(0:LMAXD,3,NMVECMAX)
      DOUBLE PRECISION MROT(3,3),FACT(0:100)
      DOUBLE PRECISION MVGLO(3,NMVECMAX),MVGLOL(0:LMAXD,3,NMVECMAX)
      DOUBLE PRECISION MVPHI(NMVECMAX),MVTET(NMVECMAX)
      CHARACTER*1  TXTL(0:LMAXDLOC)
C ..
C ..  Intrinsic Functions
      INTRINSIC ABS,ACOS,ATAN,DREAL,DIMAG,DCONJG
C ..
C ..  External subroutines
      EXTERNAL CALCROTMAT      
C ..
C ..  Data Statements
      DATA ICALL /0/
C ..
C ..  Save Statements
      SAVE ICALL,NMVEC,USC,FACT,TXTL,PI
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
         TXTL(0) = 's'
         TXTL(1) = 'p'
         TXTL(2) = 'd'
         IF( LMAXD .GE. 3 ) THEN
            DO L=3,LMAXD
               TXTL(L) = CHAR( ICHAR('f') + L - 3 )
            ENDDO
         ENDIF
C     
         WRITE (6,'(78(1H#))')
         WRITE (6,99001)
         WRITE (6,'(78(1H#))')
         WRITE (6,*)
         WRITE (6,99002)
C
         NMVEC = 2
         PI = 4.D0 * ATAN(1.D0)
C
         FACT(0) = 1.0D0
         DO I=1,100
            FACT(I) = FACT(I-1) * DBLE(I)
         END DO
C-----------------------------------------------------------------------
C  create transformation matrix   U  cartesian/sperical ccordinates
C-----------------------------------------------------------------------
C  RC,RCP  vectors in cartesian coordinates
C  RS,RSP  vectors in spherical coordinates
C         RS  = USC * RC                                 (4.40)
C         RSP = MS  * RS                                 (4.37)
C     MS(i,j) = D(j,i)                                   (4.42)
C     D  rotation matrix for complex spherical harmonics
C
C ordering of: m=-1,0,+1 >>> row 1 and 3 interchanged compared to (4.44)
C
         WSQ2 = 1.0D0/SQRT(2.0D0)
C
         USC(1,1) = WSQ2
         USC(1,2) = -CI*WSQ2
         USC(1,3) = 0.0D0
         USC(2,1) = 0.0D0
         USC(2,2) = 0.0D0
         USC(2,3) = 1.0D0
         USC(3,1) = -WSQ2
         USC(3,2) = -CI*WSQ2
         USC(3,3) = 0.0D0
C-----------------------------------------------------------------------
      END IF
C=======================================================================
C
C-----------------------------------------------------------------------
C   create the rotation matrices  DROT4 for complex spherical harmonics
C-----------------------------------------------------------------------
C
      CALL CALCROTMAT(2,1,QMPHI,QMTET,0.0D0,DROT4,FACT,4)
C
C-----------------------------------------------------------------------
C create the rotation matrix  MROT for vectors in cartesian coordinates
C NOTE:  U^+ D^T U gives the inverse of the real matrix  M
C        for that reason  the transposed matrix is stored as  MROT(J,I)
C-----------------------------------------------------------------------
C
      DO I = 1,3
         DO J = 1,3
            CS = 0.0D0
            DO K = 1,3
               CS = CS + DROT4(K+1,I+1)*USC(K,J)
            END DO
            W3X3(I,J) = CS
         END DO
      END DO
C
      DO I = 1,3
         DO J = 1,3
            CS = 0.0D0
            DO K = 1,3
               CS = CS + DCONJG(USC(K,I))*W3X3(K,J)
            END DO
            IF ( DIMAG(CS).GT.1D-8 ) WRITE (*,*) ' MROT',I,J,CS,
     &           ' ???????????'
C     see above >> MROT(I,J) = DREAL(CS)
            MROT(J,I) = DREAL(CS)
         END DO
      END DO
C-----------------------------------------------------------------------
C
C **********************************************************************
      DO IMV = 1,NMVEC
C-----------------------------------------------------------------------
C     transform from (+,-,z) to cartesian coordinates  (x,y,z)
C     note the convention
C-----------------------------------------------------------------------
         APLS = MVEVI(IT,1,IMV)
         AMIN = MVEVI(IT,2,IMV)
         MVEVI(IT,1,IMV) = (AMIN+APLS)*0.5D0
         MVEVI(IT,2,IMV) = (AMIN-APLS)*0.5D0*CI
C
         APLS = MVEVIEF(IT,1,IMV)
         AMIN = MVEVIEF(IT,2,IMV)
         MVEVIEF(IT,1,IMV) = (AMIN+APLS)*0.5D0
         MVEVIEF(IT,2,IMV) = (AMIN-APLS)*0.5D0*CI
C     
         DO L = 0,LMAXD
            APLS = MVEVIL(L,IT,1,IMV)
            AMIN = MVEVIL(L,IT,2,IMV)
            MVEVIL(L,IT,1,IMV) = (AMIN+APLS)*0.5D0
            MVEVIL(L,IT,2,IMV) = (AMIN-APLS)*0.5D0*CI
         END DO
C-----------------------------------------------------------------------
C     transform from LOCAL cartesian coordinates (x,y,z)
C               to  GLOBAL cartesian coordinates 
C-----------------------------------------------------------------------
         DO I = 1,3
            MVG(I,IMV) = CZERO
            MVGEF(I,IMV) = CZERO
            DO J = 1,3
               MVG  (I,IMV) = MVG  (I,IMV) + MROT(I,J)*MVEVI(IT,J,IMV)
               MVGEF(I,IMV) = MVGEF(I,IMV) + MROT(I,J)*MVEVIEF(IT,J,IMV)
            END DO
            MVGLO(I,IMV) = DIMAG(MVG(I,IMV))
C
            DO L = 0,LMAXD
               MVGL(L,I,IMV) = CZERO
               DO J = 1,3
                  MVGL(L,I,IMV) = MVGL(L,I,IMV) 
     &                 + MROT(I,J)*MVEVIL(L,IT,J,IMV)
               END DO
               MVGLOL(L,I,IMV) = DIMAG(MVGL(L,I,IMV))
            END DO
C
         END DO
C ......................................................................
         DO I = 1,3
            MVEVI(IT,I,IMV)   = MVG  (I,IMV) 
            MVEVIEF(IT,I,IMV) = MVGEF(I,IMV) 
C
            DO L = 0,LMAXD
               MVEVIL(L,IT,I,IMV)   = MVGL  (L,I,IMV)
            END DO
         END DO
C-----------------------------------------------------------------------
C        calculate the angles
C-----------------------------------------------------------------------
         MVX = MVGLO(1,IMV)
         MVY = MVGLO(2,IMV)
         MVZ = MVGLO(3,IMV)
C
         MV = SQRT(MVX**2+MVY**2+MVZ**2)
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF ( MV.LT.1D-8 ) THEN
            MVPHI(IMV) = 0D0
            MVTET(IMV) = 0D0
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         ELSE
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            MVXY = SQRT(MVX**2+MVY**2)
C ======================================================================
            IF ( ABS(MVXY).LT.1D-8 ) THEN
               MVPHI(IMV) = 0D0
C ======================================================================
            ELSE
C ======================================================================
               IF ( MVY.GE.0D0 ) THEN
                  MVPHI(IMV) = ACOS(MVX/MVXY)
               ELSE IF ( MVX.LT.0D0 ) THEN
                  MVPHI(IMV) = PI + ACOS(-MVX/MVXY)
               ELSE
                  MVPHI(IMV) = 2*PI - ACOS(MVX/MVXY)
               END IF
               MVPHI(IMV) = MVPHI(IMV)*180D0/PI
               IF ( ABS(MVPHI(IMV)-360.0D0).LT.1D-8 )
     &              MVPHI(IMV) = 0D0
            END IF
C ======================================================================
            IF (MVPHI(IMV).GE.345.D0) 
     &           MVPHI(IMV) = 360.D0 - MVPHI(IMV)
            MVTET(IMV) = ACOS(MVZ/MV)*180D0/PI
         END IF
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------------------------------------------------------
      END DO
C **********************************************************************
C output vector components, in and out angles
C ----------------------------------------------------------------------
      L = 0 
      WRITE (6,99003) IT,IQ,TXTL(L),((MVGLOL(L,I,IMV),I=1,3),IMV=1,2)
      WRITE (6,99004) (TXTL(L),((MVGLOL(L,I,IMV),I=1,3),IMV=1,2),
     &     L=1,LMAXD)
      WRITE (6,99005) ((MVGLO(I,IMV),I=1,3),IMV=1,2)
      WRITE (6,99006) QMPHI,QMTET,(MVPHI(IMV),MVTET(IMV),IMV=1,2)
C ----------------------------------------------------------------------
      IF ( IT.LT.NATYP ) THEN
         WRITE (6,'(3X,75(1H=))')
      ELSE
         WRITE (6,*)
         WRITE (6,'(78(1H#))')
      END IF
C ----------------------------------------------------------------------
C
99001 FORMAT (15X,'vectorial magnetic properties given with respect',/,
     &        15X,'   to the GLOBAL (crystal) frame of reference')
99002 FORMAT (29X,'m_spin',27X,'m_orb',/,
     &     3X,'ATOM/SITE     ',
     &     '    x         y         z',8X,'    x         y         z',/,
     &     3X,75(1H=))
99003 FORMAT (3X,I3,'/',I3,2X,A1,' =',3F10.5,3X,3F10.5)
99004 FORMAT (12X,A1,' =',3F10.5,3X,3F10.5)
99005 FORMAT (12X,66(1H-),/,12X,'sum',3F10.5,3X,3F10.5,/)
99006 FORMAT (3X,'angles (IN)   TET =',F9.4,' PHI =',F9.4,/,
     &        3X,'angles (calc) TET =',F9.4,' PHI =',F9.4,
     &           '   TET =',F9.4,' PHI =',F9.4)
      END
