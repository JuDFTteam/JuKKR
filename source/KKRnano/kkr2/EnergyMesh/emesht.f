C 06.10.09 *************************************************************
      SUBROUTINE EMESHT(EZ,DF,NPNT,EBOT,EMU,EFERMI,TK,
     &                  NPOL,NPNT1,NPNT2,NPNT3,IEMXD)
C **********************************************************************
C *                                                                    *
C * This subroutine provides the energy mesh in array EZ and the       *
C * appropriate integration weights in array DF for contour            *
C * integration in the complex energy plane.                           *
C *                                                                    *
C * The contour consists of two straight lines determined by two       *
C * real valued input arguments: TK and EBOT and by a parameter NPOL   *
C * which is determined in this subroutine (see below).                *
C *                                                                    *
C *         TK   = temperature in K                                    *
C *         EBOT = bottom of contour in Ry                             *
C *         NPOL = number of used poles of the Fermi-Dirac function    *
C *         EMU =  chemical potential in Ry                            *
C *                                                                    *
C * One line begins at EBOT on the real energy axis and ends at        *
C * EBOT + 2*NPOL*i*pi*k*TK. On this line NPNT1 mesh points are used   *
C * which are determined according to a Gauss-Legendre rule.           *
C *                                                                    *
C * The second line begins at EBOT+2*NPOL*i*pi*k*TK and goes parallel  *
C * to the real axis to infinity. This line is divided into            *
C * NLEG + 2 parts (NLEG can be specified below, usually NLEG = 1 is   *
C * enough, but NLEG can be increased for large values of EMU-EBOT.    *
C * On the first NLEG parts Gauss-Legendre integration with NPNT2      *
C * mesh points each are used and the Fermi-Dirac function is replaced *
C * by one since it differs from one by less than 10**(-13).           *
C * On the remaining two parts Gauss-Fermi Dirac integrations rules    *
C * with NPNT3 and NPNT4 mesh points are used. These parts are         *
C * separated by EMU + 2*NPOL*i*pi*k*TK or by the point where the      *
C * weight function (the difference of two Fermi-Dirac functions)      *
C * changes sign. (see value of XZERO in GAUFDT2 and GAUFDT3)          *
C *                                                                    *
C * useful values for NPNT1, NPNT2 and NPNT3 depend on EBOT, TK and on *
C * the option EMESHT1,  EMESHT2 or EMESHT3. For Ni, Cu and Pd useful  *
C * values are NPNT1 = 7 and specifically                              *
C * for EMESHT1                                                        *
C        TK = 3200 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
C        TK = 2400 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
C        TK = 1600 K    EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
C        TK = 1200 K    EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
C        TK = 800 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
C        TK = 400 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
C        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
C * for EMESHT2                                                        *
C        TK = 3200 K    EBOT = -1.36   NPNT2 =  6   NPNT3 = 20         *
C        TK = 2400 K    EBOT = -0.96   NPNT2 =  6   NPNT3 = 16         *
C        TK = 1600 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
C        TK = 1200 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
C        TK = 800 K     EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
C        TK = 400 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
C        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
C * for EMESHT3                                                        *
C        TK = 3200 K    EBOT = -1.36   NPNT2 =  6   NPNT3 = 20         *
C        TK = 2400 K    EBOT = -0.96   NPNT2 =  6   NPNT3 = 20         *
C        TK = 1600 K    EBOT = -0.56   NPNT2 =  8   NPNT3 = 12         *
C        TK = 1200 K    EBOT = -0.56   NPNT2 = 10   NPNT3 =  8         *
C        TK = 800 K     EBOT = -0.56   NPNT2 = 14   NPNT3 =  8         *
C        TK = 400 K     EBOT = -0.56   NPNT2 = 16   NPNT3 =  8         *
C        TK = 200 K     EBOT = -0.56   NPNT2 = 18   NPNT3 =  8         *
C * 
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
C *                                                                    *
C **********************************************************************
C     ..
      IMPLICIT NONE
C     .. Parameters ..
      DOUBLE PRECISION KB,PI,RYD
      PARAMETER (PI = 3.14159265358979312D0)
      PARAMETER (KB = 0.6333659D-5)
      PARAMETER (RYD = 13.6058D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION EBOT,FEBOT,EMU,FEMU,TK,EFERMI
      INTEGER NPNT,NPNT1,NPNT2,NPNT3,NPOL,NPNT4,IEMXD
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX DF(*),EZ(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX DE,DESUB
      DOUBLE PRECISION EA,EB,EFDBOT,ER,ETK,XSCALE
      INTEGER I,NLEG,NSUB
      LOGICAL FINEGRID
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION WI(228),XI(228)
C     ..
C     .. External Subroutines ..
      LOGICAL OPT
      EXTERNAL GAUFD,GAULEG,OPT
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
C     ..
C
      NLEG = 1
      NPNT4 = 3
      IF(TK.LT.2000.D0) NPNT4 = 2
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
      WRITE (6,'(5X,A,F10.6," (Ry)",8X,A,F10.6," (Ry)")')
     &     'E min = ',EBOT,'Fermi energy = ',EFERMI
      WRITE (6,'(5X,A,F10.6," (Ry)",8X,A,F12.6," (K )",/,5X,62(1H-))')
     &     'E max = ',EMU,'Temperature  = ',TK
C OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO OUTPUT
C
      ETK = PI*KB*TK
C ------------------------------------------------------------- NPOL = 0
      IF (NPOL.EQ.0) THEN
C===================================================================
         INQUIRE(FILE='FGRID',EXIST=FINEGRID)
C===================================================================
C===================================================================
         IF (FINEGRID) THEN
C===================================================================
           OPEN(87,FILE='FGRID',FORM='formatted')
             READ(87,*) FEBOT,FEMU
           CLOSE(87)
C
           DE = (FEMU-FEBOT)
           IF (NPNT2.GT.1) THEN
             DE = DE/(NPNT2-3)
           ELSE
             DE=DCMPLX(1.0D0,0.0D0)
           ENDIF
C
           EZ(1) = DCMPLX(EBOT,ETK)
           DF(1) = DE
           NPNT = 1
C
           DO I = 2,NPNT2-1
             NPNT = NPNT + 1
             IF ( (NPNT+1).GT.IEMXD ) THEN
               WRITE(6,'(/,5X,2A,I4)')
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
             ENDIF
             ER = FEBOT + (I-1)*DE
             EZ(NPNT) = DCMPLX(ER,ETK)
             DF(NPNT) = DE
           ENDDO
C
           EZ(NPNT+1) = DCMPLX(EMU,ETK)
           DF(NPNT+1) = DE
           NPNT       = NPNT + 1
C
           WRITE (6,FMT=9000) NPNT,ETK,ETK*RYD
C===================================================================
         ELSE
C===================================================================
           DE = (EMU-EBOT)
           IF (NPNT2.GT.1) THEN
             DE = DE/(NPNT2-1)
           ELSE
             DE=DCMPLX(1.0D0,0.0D0)
           ENDIF
           NPNT = 0
           DO 10 I = 1,NPNT2
             NPNT = NPNT + 1
             IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
             ENDIF
             ER = EBOT + (I-1)*DE
             EZ(NPNT) = DCMPLX(ER,ETK)
             DF(NPNT) = DE
 10        CONTINUE
           WRITE (6,FMT=9000) NPNT,ETK,ETK*RYD
C===================================================================
         ENDIF
C===================================================================
C===================================================================
      ENDIF
C ------------------------------------------------------------- NPOL > 0
      IF (NPOL.GT.0) THEN

C The following statements for NPOL > 0 are kept for backwards compatibility.
C They were used before the year 2009 in older versions of the 
C KKR-GF programs. (remark E.R.: kkrnano - this code is used if NPOL>0)
C *                                                                    *
C * The three lines in the backwards compatibility code are defined by:*
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
         CALL GAULEG(XI,WI,NPNT1)
         DE = NPOL*DCMPLX(0.0D0,ETK)
         NPNT = 0
         DO I = 1,NPNT1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EBOT
            DF(NPNT) = WI(I)*DE
         END DO   
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
         WRITE(6,9090) NPNT,NPOL,NPNT1,NPNT2,NPNT3
      END IF 

C ------------------------------------------------------------- NPOL < 0
      IF (NPOL.LT.0) THEN
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
         WRITE(6,9091) NPNT,-NPOL,NPNT1,NPNT2,NPNT3
      END IF
C ======================================================================
      WRITE(6,*)
      RETURN
 9000 FORMAT (5X,'Density-of-States calculation',/,
     &        5X,'Number of energy points :',I4,4X,'broadening =',
     &        3P,F9.3,' ( mRy )',/,48X,' =',3P,F9.3,' ( meV )')
 9010 FORMAT (/,5X,'GF integration rectangular contour ( ImE < 0 )',/,
     &     5X,'Number of energy points :',I4,/,
     &     13X,'poles: NPOL  =',I3,/,
     &     13X,'points on upwards contour: NPT1 =',I3)
 9020 FORMAT (/,5X,'GF integration rectangular contour ( ImE < 0 )',/,
     &     5X,'Number of energy points :',I4,/,
     &     13X,'poles: NPOL  =',I3,3X,'NPOL2  =',I3,/,
     &     13X,'points on upwards contour: NPT1 =',I3)
 9030 FORMAT (13X,'points used with Gauss-Legendre rule: NPT2 =',I3,/,
     &     13X,'points near Fermi level: NPT3 =',I3,' and NPT4 =',I3)
 9040 FORMAT (13X,'points used with Gauss-Legendre rule: ',
     &                                 I1,'*NPT2 =',I2,/,
     &     13X,'points near Fermi level: NPT3 =',I3,' and NPT4 =',I3)
 9090 FORMAT (5X,'GF integration rectangular contour ( ImE > 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
 9091 FORMAT (5X,'GF integration rectangular contour ( ImE < 0 )',/,
     &     5X,'Number of energy points :',I4,13X,'poles =',I2,/,
     &     23X,'contour: N1 =',I2,', N2 =',I4,', N3 =',I2)
      END