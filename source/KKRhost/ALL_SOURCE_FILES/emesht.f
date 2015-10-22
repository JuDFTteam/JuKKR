C 20.09.95 *************************************************************
      SUBROUTINE EMESHT(EZ,DF,NPNT,EBOT,EMU,EFERMI,TK,
     &                  NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
C **********************************************************************
C *                                                                    *
C * This subroutine provides the energy mesh in array EZ and the       *
C * appropriate integration weights in array DF.                       *
C *                                                                    *
C * Poles of the Fermi function C (Matsubara frequencies) and          *
C * a contour in the complex energy are used as described in (????).   *
C *                                                                    *
C * The contour consists of three straight lines with                  *
C * NPNT1, NPNT2, and NPNT3 integration points and is determined by    *
C * the input arguments: EBOT, EMU, TK, and NPOL.                      *
C *                                                                    *
C *            TK   = temperature in K                                 *
C *            EMU  = chemical potential in Ry                         *
C *            EBOT = bottom of contour in Ry                          *
C *            NPOL = number of Matsubara frequencies                  *
C *                                                                    *
C * The three lines are defined by:                                    *
C *                                                                    *
C *  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK                    *
C *              with NPNT1 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  2. the line from EBOT+2*NPOL*pi*i*k*TK to                         *
C *                   EMU+(2*NPOL*pi*i-30)*k*TK                        *
C *              with NPNT2 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity            *
C *              with NPNT3 integration points (Gauss-Fermi-Dirac rule)*
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT1+NPNT2+NPNT3+NPOL                           *
C *                                                                    *
C *  The integration points and weights on three lines are chosen      *
C *  according to Gauss integration rules. Only in third interval      *
C *  the Fermi function matters since exp(x) < 10**(-10) for x < -25.  *
C *                                                                    *
C *  There are two special cases determined by NPOL = 0 and NPOL < 0.  *
C *                                                                    *
C *  a) NPOL = 0 leads to density-of-states calculations               *
C *  with constant integration weights and equally distributed points  *
C *  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.                     *
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT2                                            *
C *                                                                    *
C *  b) NPOL < 0 is meant for calculations where the Fermi-Dirac       *
C *  function is replaced by a step function with step at EMU. When    *
C *  this option is used no poles of the Fermi-Dirac function are used *
C *  and the contour consists of the three straight lines:             *
C *                                                                    *
C *  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK                    *
C *              with NPNT1 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK    *
C *              with NPNT2 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU                      *
C *              with NPNT3 integration points (Gauss-Legendre rule)   *
C *                                                                    *
C *  The total number of integration points is given by:               *
C *              NPNT=NPNT1+NPNT2+NPNT3                                *
C *                                                                    *
C **********************************************************************
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EBOT,EMU,TK,EFERMI
      INTEGER NPNT,NPNT1,NPNT2,NPNT3,NPOL,IEMXD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DF(*),EZ(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DE
      DOUBLE PRECISION ER,ETK,KB,PI,RYD
      INTEGER I
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WI(128),XI(128)
C     ..
C     .. External Subroutines ..
      LOGICAL OPT
      EXTERNAL GAUFD,GAULEG,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
C     ..
C     .. Data statements ..
      DATA PI /3.14159265358979312D0/
      DATA KB /0.6333659D-5/
      DATA RYD/13.6058D0/
C     ..
C
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (1337,'(5X,A,F12.6," (Ry)",8X,A,F12.6," (Ry)")')
     &     'E min = ',EBOT,'Fermi energy = ',EFERMI
      WRITE (1337,
     &       '(5X,A,F12.6," (Ry)",8X,A,F12.6," (K )",/,5X,62(1H-))') 
     &     'E max = ',EMU,'Temperature  = ',TK
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      ETK = PI*KB*TK
C ======================================================================
      IF (NPOL.EQ.0) THEN
         DE = (EMU-EBOT)
         IF (NPNT2.GT.1) THEN
            DE = DE/(NPNT2-1)
         ELSE
            DE=DCMPLX(1.0D0,0.0D0)
         END IF
         NPNT = 0
         DO 10 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            ER = EBOT + (I-1)*DE
            EZ(NPNT) = DCMPLX(ER,ETK)
            DF(NPNT) = DE
 10      CONTINUE
         WRITE (1337,FMT=9000) NPNT,ETK,ETK*RYD
C ------------------------------------------------------------- NPOL > 0
      ELSE IF (NPOL.GT.0) THEN
         CALL GAULEG(XI,WI,NPNT1)
         DE = NPOL*DCMPLX(0.0D0,ETK)
         NPNT = 0
         DO 20 I = 1,NPNT1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT
            DF(NPNT) = WI(I)*DE
 20      CONTINUE
         CALL GAULEG(XI,WI,NPNT2)
         DE = (EMU-30*KB*TK-EBOT)*0.5D0
         DO 30 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 30      CONTINUE
         CALL GAUFD(XI,WI,NPNT3)
         DE = 30*KB*TK
         DO 40 I = 1,NPNT3
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 40      CONTINUE
         DO 50 I = NPOL,1,-1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
            DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
 50      CONTINUE
         WRITE(1337,9090) NPNT,NPOL,NPNT1,NPNT2,NPNT3
C ------------------------------------------------------------- NPOL < 0
      ELSE
         IF (NPNT1.GT.0) CALL GAULEG(XI,WI,NPNT1)
         DE = -NPOL*DCMPLX(0.0D0,ETK)
         NPNT = 0
         DO 60 I = 1,NPNT1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            NPNT = NPNT + 1
            EZ(NPNT) = XI(I)*DE + DE + EBOT
            DF(NPNT) = WI(I)*DE
 60      CONTINUE
         CALL GAULEG(XI,WI,NPNT2)
         DE = (EMU-EBOT)*0.5D0
         DO 70 I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
            IF (OPT('GF-EF   ')) EZ(NPNT) = EMU + NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
 70      CONTINUE
         IF (NPNT3.GT.0) CALL GAULEG(XI,WI,NPNT3)
         DE = -NPOL*DCMPLX(0.0D0,ETK)
         DO 80 I = NPNT3,1,-1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EMU
            DF(NPNT) = -WI(I)*DE
 80      CONTINUE
         WRITE(1337,9091) NPNT,-NPOL,NPNT1,NPNT2,NPNT3
       END IF
C ======================================================================
       WRITE(1337,*)
c*****************************************************************
c************Correction Factor for the weight in the integration
c*********according to Phivos Idea******************************
        DO 90 I=1,NPNT
 90       DF(I)=DF(I)
c*(8.5D0/8.48686D0)*(8.75D0/8.74083D0)
c GaCrN*(8.5D0/8.49286D0)
c*(8.5D0/8.48969D0)
c*(8.5D0/8.48823D0)
c*(8.75D0/8.73983D0)
c*(8.5D0/8.48686D0)
c*(8.75D0/8.75659D0)
c*(6.5D0/6.55253D0)*(7.5D0/7.47798D0)
c     *(8.75D0/8.75659D0)
c*(8.5D0/8.54963D0)
c*(6.5D0/6.41299D0)
c*(8.5D0/8.47767D0)
c*(6.5D0/6.45787D0)
c*(4.D0/4.01579D0)
c*(8.8D0/8.80272D0)*(8.8D0/8.78691D0)
c*(4.0D0/4.0419213D0)*(17.5D0/17.508D0)
c     &  *(8.0D0/7.9885D0)*(8.75D0/8.74682D0)*(8.0D0/7.9246D0)
c     &  *(8.25D0/8.24085)
c 90        write(*,*)'DF=',I,DF(I)
c*************************************************************
c**********************************************************       

      RETURN
 9000 FORMAT (5X,'Density-of-States calculation',/,
     &        5X,'Number of energy points :',I4,4X,'broadening =',
     &        3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')
 9090 FORMAT (5X,'GF integration rectangular contour ( ImE > 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
 9091 FORMAT (5X,'GF integration rectangular contour ( ImE < 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
      END
