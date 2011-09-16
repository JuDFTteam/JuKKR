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
      IF (OPT('EMESHT1 ') .OR. OPT('EMESHT2 ')
     &                    .OR. OPT('EMESHT3 ')) THEN
C The input value for NPOL is ignored and replaced by a standard one,
C which depends on the temperature TK. 
         NPOL = 2
         IF(TK.LT.3000.D0) NPOL = 3
         IF(TK.LT.2000.D0) NPOL = 4
         IF(TK.LT.1500.D0) NPOL = 6
         IF(TK.LT.1000.D0) NPOL = 8
         IF(TK.LT.500.D0) NPOL = 16
         IF(TK.LT.250.D0) NPOL = 32
C Note that NPOL must be even for option EMESHT2 and divisible by three
C for option EMESHT3
         IF(OPT('EMESHT2 ').AND.TK.LT.3000.D0
     &                     .AND.TK.GT.1500.D0) NPOL = 4
         IF(OPT('EMESHT3 ')) THEN
         NPOL = 3
         IF(TK.LT.3000.D0) NPOL = 3
         IF(TK.LT.2000.D0) NPOL = 6
         IF(TK.LT.1000.D0) NPOL = 9
         IF(TK.LT.500.D0) NPOL = 15
         IF(TK.LT.250.D0) NPOL = 30
         END IF
         IF(TK.LT.0.D0) THEN
            WRITE(6,'(/,5X,A)') 
     &             'Temperature ERROR: temperature is negative'
            STOP '     < EMESHT >'
         END IF
      IF (OPT('EMESHT1 ')) CALL GAUFDT1(XI,WI,EFDBOT,0,-1)
      IF (OPT('EMESHT2 ')) CALL GAUFDT2(XI,WI,EFDBOT,0,-1)
      IF (OPT('EMESHT3 ')) CALL GAUFDT3(XI,WI,EFDBOT,0,-1)
      IF (OPT('EMESHT2 ')) EFDBOT = 2.D0 * EFDBOT
      IF (OPT('EMESHT3 ')) EFDBOT = 3.D0 * EFDBOT
         IF(EFDBOT*KB*TK.GE.(EMU-EBOT)*0.99D0) THEN 
            WRITE(6,'(/,5X,A)') 
     &             'Temperature ERROR: temperature is too high. '
            WRITE(6,'(/,5X,2A)') 
     &             'Fermi-Dirac function should be used everywhere ',
     &             'on the contour, which has not been implemented'
            STOP '     < EMESHT >'
         END IF
         IF(EMU.LT.EBOT) THEN 
            WRITE(6,'(/,5X,2A)') 
     &             'Fermi level is smaller than energy at the bottom',
     &             'of the contour, which has not been implemented'
            STOP '     < EMESHT >'
         END IF
C determine points and weights on upwards contour at EBOT
         CALL GAULEGE(XI,WI,NPNT1)
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
C
C determine Gauss-Legendre points and weights on contour parallel to
C real energy axis.
C note that NLEG is fixed as given above
      DO NSUB = 1,NLEG 
      DESUB = (EMU-EBOT-EFDBOT*KB*TK)/NLEG
      EB = EBOT + NSUB*DESUB
      EA = EB   - DESUB
      DE = (EB-EA)*0.5D0
      WRITE(6,*) 'EA,EB,DE',EA,EB,DE
         CALL GAULEGE(XI,WI,NPNT2)
         DO I = 1,NPNT2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + DE + EA + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
         END DO
      END DO
C
C determine Gauss-Fermi-Dirac points and weights on contour parallel to
C real energy axis below 0, ln(2+sqrt(7) or ln(1/2+sqrt(33)/2)).
      IF (OPT('EMESHT1 '))
     &                CALL GAUFDT1(XI,WI,XSCALE,NPNT3,1)
      IF (OPT('EMESHT2 '))
     &                CALL GAUFDT2(XI,WI,XSCALE,NPNT3,1)
      IF (OPT('EMESHT3 '))
     &                CALL GAUFDT3(XI,WI,XSCALE,NPNT3,1)
      DE = EFDBOT*KB*TK
         DO I = 1,NPNT3
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
         END DO  
C
C determine Gauss-Fermi-Dirac points and weights on contour parallel to
C real energy axis above 0, ln(2+sqrt(7) or ln(1/2+sqrt(33)/2)).
C note that NPNT4 is fixed as given above
      IF (OPT('EMESHT1 '))
     &                CALL GAUFDT1(XI,WI,XSCALE,NPNT4,0)
      IF (OPT('EMESHT2 '))
     &                CALL GAUFDT2(XI,WI,XSCALE,NPNT4,0)
      IF (OPT('EMESHT3 '))
     &                CALL GAUFDT3(XI,WI,XSCALE,NPNT4,0)
      DE = EFDBOT*KB*TK
         DO I = 1,NPNT4
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
            DF(NPNT) = WI(I)*DE
         END DO  
C
C determine Matsubara points and weights for 2*TK
         IF(OPT('EMESHT2 ')) THEN
          DO I = NPOL,2,-2
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN
               WRITE(6,'(/,5X,2A,I4)')
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = EMU + (2*I-2)*DCMPLX(0.0D0,ETK)
            DF(NPNT) = -4*DCMPLX(0.0D0,ETK)
         DF(NPNT) = -1.D0/3.D0*DF(NPNT)
         END DO
         END IF
C
C determine Matsubara points and weights for TK and 3*TK
         DO I = NPOL,1,-1
            NPNT = NPNT + 1
            IF ( NPNT.GT.IEMXD ) THEN 
               WRITE(6,'(/,5X,2A,I4)') 
     &              'Dimension ERROR: Increase IEMXD in inc.p to ',
     &              'at least ',NPNT
               STOP '     < EMESHT >'
            END IF
            EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
            DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
            IF(OPT('EMESHT2 ')) DF(NPNT) = 4.D0/3.D0*DF(NPNT)
            IF(OPT('EMESHT3 ')) THEN
            IF(MOD(2*I-1,3).EQ.0) THEN
            DF(NPNT) = 6.D0/8.D0*DF(NPNT)
            ELSE 
            DF(NPNT) = 9.D0/8.D0*DF(NPNT)
            END IF
            END IF
         END DO
C
c print number of points
         IF(OPT('EMESHT2 ')) THEN
         WRITE(6,9020) NPNT,NPOL,NPOL/2,NPNT1
         ELSE
         WRITE(6,9010) NPNT,NPOL,NPNT1
         END IF
         IF(NLEG.EQ.1) THEN
         WRITE(6,9030) NPNT2,NPNT3,NPNT4
         ELSE
         WRITE(6,9040) NLEG,NLEG*NPNT2,NPNT3,NPNT4
         END IF
         WRITE(6,'(/,5X,A)') 'WARNING  WARNING  WARNING  WARNING'
         WRITE(6,'(5X,2A,I4)') 'in EMESHT ',
     &     'the input value of NPOL is replaced by NPOL = ',NPOL
      ELSE
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
c     ************************************************ 
      SUBROUTINE GAULEGE(XI,WI,N) 
c     ************************************************ 
C     .. Scalar Arguments .. 
      INTEGER N 
C     .. 
C     .. Array Arguments .. 
      DOUBLE PRECISION WI(*),XI(*) 
C     .. 
C     .. Local Scalars .. 
      INTEGER I 
C     .. 
      IF (N.GT.32) N = ((N-1)/4+1)*4 
      IF (N.GT.64) N = ((N-1)/8+1)*8 
      IF(NCASE.EQ.0) THEN 
        IF(NC.EQ.0) THEN 
          IF(N.EQ.1) THEN 
            XI(1) = 0.D0 
            WI(1) = 2.D0 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) = -.57735026918962576451D0 
            WI(1) =  .10000000000000000000D1 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) = -.77459666924148337704D0 
            XI(2) = 0.D0 
            WI(1) =  .55555555555555555556D0 
            WI(2) =  .88888888888888888889D0 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) = -.86113631159405257522D0 
            XI(2) = -.33998104358485626480D0 
            WI(1) =  .34785484513745385737D0 
            WI(2) =  .65214515486254614263D0 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) = -.90617984593866399280D0 
            XI(2) = -.53846931010568309104D0 
            XI(3) = 0.D0 
            WI(1) =  .23692688505618908751D0 
            WI(2) =  .47862867049936646804D0 
            WI(3) =  .56888888888888888889D0 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) = -.93246951420315202781D0 
            XI(2) = -.66120938646626451366D0 
            XI(3) = -.23861918608319690863D0 
            WI(1) =  .17132449237917034504D0 
            WI(2) =  .36076157304813860757D0 
            WI(3) =  .46791393457269104739D0 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) = -.94910791234275852453D0 
            XI(2) = -.74153118559939443986D0 
            XI(3) = -.40584515137739716691D0 
            XI(4) = 0.D0 
            WI(1) =  .12948496616886969327D0 
            WI(2) =  .27970539148927666790D0 
            WI(3) =  .38183005050511894495D0 
            WI(4) =  .41795918367346938776D0 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) = -.96028985649753623168D0 
            XI(2) = -.79666647741362673959D0 
            XI(3) = -.52553240991632898582D0 
            XI(4) = -.18343464249564980494D0 
            WI(1) =  .10122853629037625915D0 
            WI(2) =  .22238103445337447054D0 
            WI(3) =  .31370664587788728734D0 
            WI(4) =  .36268378337836198297D0 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) = -.96816023950762608984D0 
            XI(2) = -.83603110732663579430D0 
            XI(3) = -.61337143270059039731D0 
            XI(4) = -.32425342340380892904D0 
            XI(5) = 0.D0 
            WI(1) =  .81274388361574411972D-1 
            WI(2) =  .18064816069485740406D0 
            WI(3) =  .26061069640293546232D0 
            WI(4) =  .31234707704000284007D0 
            WI(5) =  .33023935500125976316D0 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) = -.97390652851717172008D0 
            XI(2) = -.86506336668898451073D0 
            XI(3) = -.67940956829902440623D0 
            XI(4) = -.43339539412924719080D0 
            XI(5) = -.14887433898163121088D0 
            WI(1) =  .66671344308688137594D-1 
            WI(2) =  .14945134915058059315D0 
            WI(3) =  .21908636251598204400D0 
            WI(4) =  .26926671930999635509D0 
            WI(5) =  .29552422471475287017D0 
            GO TO 10 
          END IF 
          IF(N.EQ.11) THEN 
            XI(1) = -.97822865814605699280D0 
            XI(2) = -.88706259976809529908D0 
            XI(3) = -.73015200557404932409D0 
            XI(4) = -.51909612920681181593D0 
            XI(5) = -.26954315595234497233D0 
            XI(6) = 0.D0 
            WI(1) =  .55668567116173666483D-1 
            WI(2) =  .12558036946490462463D0 
            WI(3) =  .18629021092773425143D0 
            WI(4) =  .23319376459199047992D0 
            WI(5) =  .26280454451024666218D0 
            WI(6) =  .27292508677790063071D0 
            GO TO 10 
          END IF 
          IF(N.EQ.12) THEN 
            XI(1) = -.98156063424671925069D0 
            XI(2) = -.90411725637047485668D0 
            XI(3) = -.76990267419430468704D0 
            XI(4) = -.58731795428661744730D0 
            XI(5) = -.36783149899818019375D0 
            XI(6) = -.12523340851146891547D0 
            WI(1) =  .47175336386511827195D-1 
            WI(2) =  .10693932599531843096D0 
            WI(3) =  .16007832854334622633D0 
            WI(4) =  .20316742672306592175D0 
            WI(5) =  .23349253653835480876D0 
            WI(6) =  .24914704581340278500D0 
            GO TO 10 
          END IF 
          IF(N.EQ.13) THEN 
            XI(1) = -.98418305471858814947D0 
            XI(2) = -.91759839922297796521D0 
            XI(3) = -.80157809073330991279D0 
            XI(4) = -.64234933944034022064D0 
            XI(5) = -.44849275103644685288D0 
            XI(6) = -.23045831595513479407D0 
            XI(7) = 0.D0 
            WI(1) =  .40484004765315879520D-1 
            WI(2) =  .92121499837728447914D-1 
            WI(3) =  .13887351021978723846D0 
            WI(4) =  .17814598076194573828D0 
            WI(5) =  .20781604753688850231D0 
            WI(6) =  .22628318026289723841D0 
            WI(7) =  .23255155323087391019D0 
            GO TO 10 
          END IF 
          IF(N.EQ.14) THEN 
            XI(1) = -.98628380869681233884D0 
            XI(2) = -.92843488366357351734D0 
            XI(3) = -.82720131506976499319D0 
            XI(4) = -.68729290481168547015D0 
            XI(5) = -.51524863635815409197D0 
            XI(6) = -.31911236892788976044D0 
            XI(7) = -.10805494870734366207D0 
            WI(1) =  .35119460331751863032D-1 
            WI(2) =  .80158087159760209806D-1 
            WI(3) =  .12151857068790318469D0 
            WI(4) =  .15720316715819353457D0 
            WI(5) =  .18553839747793781374D0 
            WI(6) =  .20519846372129560397D0 
            WI(7) =  .21526385346315779020D0 
            GO TO 10 
          END IF 
          IF(N.EQ.15) THEN 
            XI(1) = -.98799251802048542849D0 
            XI(2) = -.93727339240070590431D0 
            XI(3) = -.84820658341042721620D0 
            XI(4) = -.72441773136017004742D0 
            XI(5) = -.57097217260853884754D0 
            XI(6) = -.39415134707756336990D0 
            XI(7) = -.20119409399743452230D0 
            XI(8) = 0.D0 
            WI(1) =  .30753241996117268355D-1 
            WI(2) =  .70366047488108124709D-1 
            WI(3) =  .10715922046717193501D0 
            WI(4) =  .13957067792615431445D0 
            WI(5) =  .16626920581699393355D0 
            WI(6) =  .18616100001556221103D0 
            WI(7) =  .19843148532711157646D0 
            WI(8) =  .20257824192556127288D0 
            GO TO 10 
          END IF 
          IF(N.EQ.16) THEN 
            XI(1) = -.98940093499164993260D0 
            XI(2) = -.94457502307323257608D0 
            XI(3) = -.86563120238783174388D0 
            XI(4) = -.75540440835500303390D0 
            XI(5) = -.61787624440264374845D0 
            XI(6) = -.45801677765722738634D0 
            XI(7) = -.28160355077925891323D0 
            XI(8) = -.95012509837637440185D-1 
            WI(1) =  .27152459411754094852D-1 
            WI(2) =  .62253523938647892863D-1 
            WI(3) =  .95158511682492784810D-1 
            WI(4) =  .12462897125553387205D0 
            WI(5) =  .14959598881657673208D0 
            WI(6) =  .16915651939500253819D0 
            WI(7) =  .18260341504492358887D0 
            WI(8) =  .18945061045506849629D0 
            GO TO 10 
          END IF 
          IF(N.EQ.17) THEN 
            XI(1) = -.99057547531441733568D0 
            XI(2) = -.95067552176876776122D0 
            XI(3) = -.88023915372698590212D0 
            XI(4) = -.78151400389680140693D0 
            XI(5) = -.65767115921669076585D0 
            XI(6) = -.51269053708647696789D0 
            XI(7) = -.35123176345387631530D0 
            XI(8) = -.17848418149584785585D0 
            XI(9) = 0.D0 
            WI(1) =  .24148302868547931960D-1 
            WI(2) =  .55459529373987201129D-1 
            WI(3) =  .85036148317179180884D-1 
            WI(4) =  .11188384719340397109D0 
            WI(5) =  .13513636846852547329D0 
            WI(6) =  .15404576107681028808D0 
            WI(7) =  .16800410215645004451D0 
            WI(8) =  .17656270536699264633D0 
            WI(9) =  .17944647035620652546D0 
            GO TO 10 
          END IF 
          IF(N.EQ.18) THEN 
            XI(1) = -.99156516842093094673D0 
            XI(2) = -.95582394957139775518D0 
            XI(3) = -.89260246649755573921D0 
            XI(4) = -.80370495897252311568D0 
            XI(5) = -.69168704306035320787D0 
            XI(6) = -.55977083107394753461D0 
            XI(7) = -.41175116146284264604D0 
            XI(8) = -.25188622569150550959D0 
            XI(9) = -.84775013041735301242D-1 
            WI(1) =  .21616013526483310313D-1 
            WI(2) =  .49714548894969796453D-1 
            WI(3) =  .76425730254889056529D-1 
            WI(4) =  .10094204410628716556D0 
            WI(5) =  .12255520671147846018D0 
            WI(6) =  .14064291467065065120D0 
            WI(7) =  .15468467512626524493D0 
            WI(8) =  .16427648374583272299D0 
            WI(9) =  .16914238296314359184D0 
            GO TO 10 
          END IF 
          IF(N.EQ.19) THEN 
            XI(1) = -.99240684384358440319D0 
            XI(2) = -.96020815213483003085D0 
            XI(3) = -.90315590361481790164D0 
            XI(4) = -.82271465653714282498D0 
            XI(5) = -.72096617733522937862D0 
            XI(6) = -.60054530466168102347D0 
            XI(7) = -.46457074137596094572D0 
            XI(8) = -.31656409996362983199D0 
            XI(9) = -.16035864564022537587D0 
            XI(10) = 0.D0 
            WI(1) =  .19461788229726477036D-1 
            WI(2) =  .44814226765699600333D-1 
            WI(3) =  .69044542737641226581D-1 
            WI(4) =  .91490021622449999464D-1 
            WI(5) =  .11156664554733399472D0 
            WI(6) =  .12875396253933622768D0 
            WI(7) =  .14260670217360661178D0 
            WI(8) =  .15276604206585966678D0 
            WI(9) =  .15896884339395434765D0 
            WI(10) =  .16105444984878369598D0 
            GO TO 10 
          END IF 
          IF(N.EQ.20) THEN 
            XI(1) = -.99312859918509492479D0 
            XI(2) = -.96397192727791379127D0 
            XI(3) = -.91223442825132590587D0 
            XI(4) = -.83911697182221882339D0 
            XI(5) = -.74633190646015079261D0 
            XI(6) = -.63605368072651502545D0 
            XI(7) = -.51086700195082709800D0 
            XI(8) = -.37370608871541956067D0 
            XI(9) = -.22778585114164507808D0 
            XI(10) = -.76526521133497333755D-1 
            WI(1) =  .17614007139152118312D-1 
            WI(2) =  .40601429800386941331D-1 
            WI(3) =  .62672048334109063570D-1 
            WI(4) =  .83276741576704748725D-1 
            WI(5) =  .10193011981724043504D0 
            WI(6) =  .11819453196151841731D0 
            WI(7) =  .13168863844917662690D0 
            WI(8) =  .14209610931838205133D0 
            WI(9) =  .14917298647260374679D0 
            WI(10) =  .15275338713072585070D0 
            GO TO 10 
          END IF 
          IF(N.EQ.21) THEN 
            XI(1) = -.99375217062038950026D0 
            XI(2) = -.96722683856630629432D0 
            XI(3) = -.92009933415040082879D0 
            XI(4) = -.85336336458331728365D0 
            XI(5) = -.76843996347567790862D0 
            XI(6) = -.66713880419741231931D0 
            XI(7) = -.55161883588721980706D0 
            XI(8) = -.42434212020743878357D0 
            XI(9) = -.28802131680240109660D0 
            XI(10) = -.14556185416089509094D0 
            XI(11) = 0.D0 
            WI(1) =  .16017228257774333324D-1 
            WI(2) =  .36953789770852493800D-1 
            WI(3) =  .57134425426857208284D-1 
            WI(4) =  .76100113628379302017D-1 
            WI(5) =  .93444423456033861553D-1 
            WI(6) =  .10879729916714837766D0 
            WI(7) =  .12183141605372853420D0 
            WI(8) =  .13226893863333746178D0 
            WI(9) =  .13988739479107315472D0 
            WI(10) =  .14452440398997005906D0 
            WI(11) =  .14608113364969042719D0 
            GO TO 10 
          END IF 
          IF(N.EQ.22) THEN 
            XI(1) = -.99429458548239929207D0 
            XI(2) = -.97006049783542872712D0 
            XI(3) = -.92695677218717400052D0 
            XI(4) = -.86581257772030013654D0 
            XI(5) = -.78781680597920816200D0 
            XI(6) = -.69448726318668278005D0 
            XI(7) = -.58764040350691159296D0 
            XI(8) = -.46935583798675702641D0 
            XI(9) = -.34193582089208422516D0 
            XI(10) = -.20786042668822128548D0 
            XI(11) = -.69739273319722221214D-1 
            WI(1) =  .14627995298272200685D-1 
            WI(2) =  .33774901584814154793D-1 
            WI(3) =  .52293335152683285940D-1 
            WI(4) =  .69796468424520488095D-1 
            WI(5) =  .85941606217067727414D-1 
            WI(6) =  .10041414444288096493D0 
            WI(7) =  .11293229608053921839D0 
            WI(8) =  .12325237681051242429D0 
            WI(9) =  .13117350478706237073D0 
            WI(10) =  .13654149834601517135D0 
            WI(11) =  .13925187285563199338D0 
            GO TO 10 
          END IF 
          IF(N.EQ.23) THEN 
            XI(1) = -.99476933499755212352D0 
            XI(2) = -.97254247121811523196D0 
            XI(3) = -.93297108682601610235D0 
            XI(4) = -.87675235827044166738D0 
            XI(5) = -.80488840161883989215D0 
            XI(6) = -.71866136313195019446D0 
            XI(7) = -.61960987576364615639D0 
            XI(8) = -.50950147784600754969D0 
            XI(9) = -.39030103803029083142D0 
            XI(10) = -.26413568097034493053D0 
            XI(11) = -.13325682429846611093D0 
            XI(12) = 0.D0 
            WI(1) =  .13411859487141772081D-1 
            WI(2) =  .30988005856979444311D-1 
            WI(3) =  .48037671731084668572D-1 
            WI(4) =  .64232421408525852127D-1 
            WI(5) =  .79281411776718954923D-1 
            WI(6) =  .92915766060035147477D-1 
            WI(7) =  .10489209146454141007D0 
            WI(8) =  .11499664022241136494D0 
            WI(9) =  .12304908430672953047D0 
            WI(10) =  .12890572218808214998D0 
            WI(11) =  .13246203940469661737D0 
            WI(12) =  .13365457218610617535D0 
            GO TO 10 
          END IF 
          IF(N.EQ.24) THEN 
            XI(1) = -.99518721999702136018D0 
            XI(2) = -.97472855597130949820D0 
            XI(3) = -.93827455200273275852D0 
            XI(4) = -.88641552700440103421D0 
            XI(5) = -.82000198597390292195D0 
            XI(6) = -.74012419157855436424D0 
            XI(7) = -.64809365193697556925D0 
            XI(8) = -.54542147138883953566D0 
            XI(9) = -.43379350762604513849D0 
            XI(10) = -.31504267969616337439D0 
            XI(11) = -.19111886747361630916D0 
            XI(12) = -.64056892862605626085D-1 
            WI(1) =  .12341229799987199547D-1 
            WI(2) =  .28531388628933663181D-1 
            WI(3) =  .44277438817419806169D-1 
            WI(4) =  .59298584915436780746D-1 
            WI(5) =  .73346481411080305734D-1 
            WI(6) =  .86190161531953275917D-1 
            WI(7) =  .97618652104113888270D-1 
            WI(8) =  .10744427011596563478D0 
            WI(9) =  .11550566805372560135D0 
            WI(10) =  .12167047292780339120D0 
            WI(11) =  .12583745634682829612D0 
            WI(12) =  .12793819534675215697D0 
            GO TO 10 
          END IF 
          IF(N.EQ.25) THEN 
            XI(1) = -.99555696979049809791D0 
            XI(2) = -.97666392145951751150D0 
            XI(3) = -.94297457122897433941D0 
            XI(4) = -.89499199787827536885D0 
            XI(5) = -.83344262876083400142D0 
            XI(6) = -.75925926303735763058D0 
            XI(7) = -.67356636847346836449D0 
            XI(8) = -.57766293024122296772D0 
            XI(9) = -.47300273144571496052D0 
            XI(10) = -.36117230580938783774D0 
            XI(11) = -.24386688372098843205D0 
            XI(12) = -.12286469261071039639D0 
            XI(13) = 0.D0 
            WI(1) =  .11393798501026287948D-1 
            WI(2) =  .26354986615032137262D-1 
            WI(3) =  .40939156701306312656D-1 
            WI(4) =  .54904695975835191926D-1 
            WI(5) =  .68038333812356917207D-1 
            WI(6) =  .80140700335001018013D-1 
            WI(7) =  .91028261982963649811D-1 
            WI(8) =  .10053594906705064420D0 
            WI(9) =  .10851962447426365312D0 
            WI(10) =  .11485825914571164834D0 
            WI(11) =  .11945576353578477223D0 
            WI(12) =  .12224244299031004169D0 
            WI(13) =  .12317605372671545120D0 
            GO TO 10 
          END IF 
          IF(N.EQ.26) THEN 
            XI(1) = -.99588570114561692900D0 
            XI(2) = -.97838544595647099110D0 
            XI(3) = -.94715906666171425014D0 
            XI(4) = -.90263786198430707422D0 
            XI(5) = -.84544594278849801880D0 
            XI(6) = -.77638594882067885619D0 
            XI(7) = -.69642726041995726486D0 
            XI(8) = -.60669229301761806323D0 
            XI(9) = -.50844071482450571770D0 
            XI(10) = -.40305175512348630648D0 
            XI(11) = -.29200483948595689514D0 
            XI(12) = -.17685882035689018397D0 
            XI(13) = -.59230093429313207094D-1 
            WI(1) =  .10551372617343007156D-1 
            WI(2) =  .24417851092631908790D-1 
            WI(3) =  .37962383294362763950D-1 
            WI(4) =  .50975825297147811998D-1 
            WI(5) =  .63274046329574835539D-1 
            WI(6) =  .74684149765659745887D-1 
            WI(7) =  .85045894313485239210D-1 
            WI(8) =  .94213800355914148464D-1 
            WI(9) =  .10205916109442542324D0 
            WI(10) =  .10847184052857659066D0 
            WI(11) =  .11336181654631966655D0 
            WI(12) =  .11666044348529658204D0 
            WI(13) =  .11832141527926227652D0 
            GO TO 10 
          END IF 
          IF(N.EQ.27) THEN 
            XI(1) = -.99617926288898856694D0 
            XI(2) = -.97992347596150122286D0 
            XI(3) = -.95090055781470500685D0 
            XI(4) = -.90948232067749110430D0 
            XI(5) = -.85620790801829449030D0 
            XI(6) = -.79177163907050822714D0 
            XI(7) = -.71701347373942369929D0 
            XI(8) = -.63290797194649514093D0 
            XI(9) = -.54055156457945689490D0 
            XI(10) = -.44114825175002688059D0 
            XI(11) = -.33599390363850889973D0 
            XI(12) = -.22645936543953685886D0 
            XI(13) = -.11397258560952996693D0 
            XI(14) = 0.D0 
            WI(1) =  .97989960512943602612D-2 
            WI(2) =  .22686231596180623196D-1 
            WI(3) =  .35297053757419711023D-1 
            WI(4) =  .47449412520615062704D-1 
            WI(5) =  .58983536859833599110D-1 
            WI(6) =  .69748823766245592984D-1 
            WI(7) =  .79604867773057771263D-1 
            WI(8) =  .88423158543756950194D-1 
            WI(9) =  .96088727370028507566D-1 
            WI(10) =  .10250163781774579867D0 
            WI(11) =  .10757828578853318721D0 
            WI(12) =  .11125248835684519267D0 
            WI(13) =  .11347634610896514862D0 
            WI(14) =  .11422086737895698905D0 
            GO TO 10 
          END IF 
          IF(N.EQ.28) THEN 
            XI(1) = -.99644249757395444995D0 
            XI(2) = -.98130316537087275369D0 
            XI(3) = -.95425928062893819725D0 
            XI(4) = -.91563302639213207387D0 
            XI(5) = -.86589252257439504894D0 
            XI(6) = -.80564137091717917145D0 
            XI(7) = -.73561087801363177203D0 
            XI(8) = -.65665109403886496122D0 
            XI(9) = -.56972047181140171931D0 
            XI(10) = -.47587422495511826103D0 
            XI(11) = -.37625151608907871022D0 
            XI(12) = -.27206162763517807768D0 
            XI(13) = -.16456928213338077128D0 
            XI(14) = -.55079289884034270427D-1 
            WI(1) =  .91242825930945177388D-2 
            WI(2) =  .21132112592771259752D-1 
            WI(3) =  .32901427782304379978D-1 
            WI(4) =  .44272934759004227840D-1 
            WI(5) =  .55107345675716745431D-1 
            WI(6) =  .65272923966999595793D-1 
            WI(7) =  .74646214234568779024D-1 
            WI(8) =  .83113417228901218390D-1 
            WI(9) =  .90571744393032840942D-1 
            WI(10) =  .96930657997929915850D-1 
            WI(11) =  .10211296757806076981D0 
            WI(12) =  .10605576592284641791D0 
            WI(13) =  .10871119225829413525D0 
            WI(14) =  .11004701301647519628D0 
            GO TO 10 
          END IF 
          IF(N.EQ.29) THEN 
            XI(1) = -.99667944226059658616D0 
            XI(2) = -.98254550526141317487D0 
            XI(3) = -.95728559577808772580D0 
            XI(4) = -.92118023295305878509D0 
            XI(5) = -.87463780492010279042D0 
            XI(6) = -.81818548761525244499D0 
            XI(7) = -.75246285173447713391D0 
            XI(8) = -.67821453760268651516D0 
            XI(9) = -.59628179713822782038D0 
            XI(10) = -.50759295512422764210D0 
            XI(11) = -.41315288817400866389D0 
            XI(12) = -.31403163786763993495D0 
            XI(13) = -.21135228616600107451D0 
            XI(14) = -.10627823013267923017D0 
            XI(15) = 0.D0 
            WI(1) =  .85169038787464096543D-2 
            WI(2) =  .19732085056122705984D-1 
            WI(3) =  .30740492202093622644D-1 
            WI(4) =  .41402062518682836105D-1 
            WI(5) =  .51594826902497923913D-1 
            WI(6) =  .61203090657079138542D-1 
            WI(7) =  .70117933255051278570D-1 
            WI(8) =  .78238327135763783828D-1 
            WI(9) =  .85472257366172527545D-1 
            WI(10) =  .91737757139258763348D-1 
            WI(11) =  .96963834094408606302D-1 
            WI(12) =  .10109127375991496612D0 
            WI(13) =  .10407331007772937391D0 
            WI(14) =  .10587615509732094141D0 
            WI(15) =  .10647938171831424425D0 
            GO TO 10 
          END IF 
          IF(N.EQ.30) THEN 
            XI(1) = -.99689348407464954027D0 
            XI(2) = -.98366812327974720997D0 
            XI(3) = -.96002186496830751222D0 
            XI(4) = -.92620004742927432588D0 
            XI(5) = -.88256053579205268154D0 
            XI(6) = -.82956576238276839744D0 
            XI(7) = -.76777743210482619492D0 
            XI(8) = -.69785049479331579693D0 
            XI(9) = -.62052618298924286114D0 
            XI(10) = -.53662414814201989926D0 
            XI(11) = -.44703376953808917678D0 
            XI(12) = -.35270472553087811347D0 
            XI(13) = -.25463692616788984644D0 
            XI(14) = -.15386991360858354696D0 
            XI(15) = -.51471842555317695833D-1 
            WI(1) =  .79681924961666056155D-2 
            WI(2) =  .18466468311090959142D-1 
            WI(3) =  .28784707883323369350D-1 
            WI(4) =  .38799192569627049597D-1 
            WI(5) =  .48402672830594052903D-1 
            WI(6) =  .57493156217619066482D-1 
            WI(7) =  .65974229882180495128D-1 
            WI(8) =  .73755974737705206268D-1 
            WI(9) =  .80755895229420215355D-1 
            WI(10) =  .86899787201082979802D-1 
            WI(11) =  .92122522237786128718D-1 
            WI(12) =  .96368737174644259639D-1 
            WI(13) =  .99593420586795267063D-1 
            WI(14) =  .10176238974840550460D0 
            WI(15) =  .10285265289355884034D0 
            GO TO 10 
          END IF 
          IF(N.EQ.31) THEN 
            XI(1) = -.99708748181947707406D0 
            XI(2) = -.98468590966515248400D0 
            XI(3) = -.96250392509294966179D0 
            XI(4) = -.93075699789664816496D0 
            XI(5) = -.88976002994827104337D0 
            XI(6) = -.83992032014626734009D0 
            XI(7) = -.78173314841662494041D0 
            XI(8) = -.71577678458685328391D0 
            XI(9) = -.64270672292426034618D0 
            XI(10) = -.56324916140714926272D0 
            XI(11) = -.47819378204490248044D0 
            XI(12) = -.38838590160823294306D0 
            XI(13) = -.29471806998170161662D0 
            XI(14) = -.19812119933557062877D0 
            XI(15) = -.99555312152341520325D-1 
            XI(16) = 0.D0 
            WI(1) =  .74708315792487758587D-2 
            WI(2) =  .17318620790310582463D-1 
            WI(3) =  .27009019184979421801D-1 
            WI(4) =  .36432273912385464024D-1 
            WI(5) =  .45493707527201102902D-1 
            WI(6) =  .54103082424916853712D-1 
            WI(7) =  .62174786561028426910D-1 
            WI(8) =  .69628583235410366168D-1 
            WI(9) =  .76390386598776616426D-1 
            WI(10) =  .82392991761589263904D-1 
            WI(11) =  .87576740608477876126D-1 
            WI(12) =  .91890113893641478215D-1 
            WI(13) =  .95290242912319512807D-1 
            WI(14) =  .97743335386328725093D-1 
            WI(15) =  .99225011226672307875D-1 
            WI(16) =  .99720544793426451428D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.32) THEN 
            XI(1) = -.99726386184948156354D0 
            XI(2) = -.98561151154526833540D0 
            XI(3) = -.96476225558750643077D0 
            XI(4) = -.93490607593773968917D0 
            XI(5) = -.89632115576605212397D0 
            XI(6) = -.84936761373256997013D0 
            XI(7) = -.79448379596794240696D0 
            XI(8) = -.73218211874028968039D0 
            XI(9) = -.66304426693021520098D0 
            XI(10) = -.58771575724076232904D0 
            XI(11) = -.50689990893222939002D0 
            XI(12) = -.42135127613063534536D0 
            XI(13) = -.33186860228212764978D0 
            XI(14) = -.23928736225213707454D0 
            XI(15) = -.14447196158279649349D0 
            XI(16) = -.48307665687738316235D-1 
            WI(1) =  .70186100094700966004D-2 
            WI(2) =  .16274394730905670605D-1 
            WI(3) =  .25392065309262059456D-1 
            WI(4) =  .34273862913021433103D-1 
            WI(5) =  .42835898022226680657D-1 
            WI(6) =  .50998059262376176196D-1 
            WI(7) =  .58684093478535547145D-1 
            WI(8) =  .65822222776361846838D-1 
            WI(9) =  .72345794108848506225D-1 
            WI(10) =  .78193895787070306472D-1 
            WI(11) =  .83311924226946755222D-1 
            WI(12) =  .87652093004403811143D-1 
            WI(13) =  .91173878695763884713D-1 
            WI(14) =  .93844399080804565639D-1 
            WI(15) =  .95638720079274859419D-1 
            WI(16) =  .96540088514727800567D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.36) THEN 
            XI(1) = -.99783046248408583620D0 
            XI(2) = -.98858647890221223807D0 
            XI(3) = -.97202769104969794934D0 
            XI(4) = -.94827298439950754520D0 
            XI(5) = -.91749777451565906608D0 
            XI(6) = -.87992980089039713198D0 
            XI(7) = -.83584716699247530642D0 
            XI(8) = -.78557623013220651283D0 
            XI(9) = -.72948917159355658209D0 
            XI(10) = -.66800123658552106210D0 
            XI(11) = -.60156765813598053508D0 
            XI(12) = -.53068028592624516164D0 
            XI(13) = -.45586394443342026721D0 
            XI(14) = -.37767254711968921632D0 
            XI(15) = -.29668499534402827050D0 
            XI(16) = -.21350089231686557894D0 
            XI(17) = -.12873610380938478865D0 
            XI(18) = -.43018198473708607227D-1 
            WI(1) =  .55657196642450453613D-2 
            WI(2) =  .12915947284065574405D-1 
            WI(3) =  .20181515297735471532D-1 
            WI(4) =  .27298621498568779094D-1 
            WI(5) =  .34213810770307229921D-1 
            WI(6) =  .40875750923644895474D-1 
            WI(7) =  .47235083490265978417D-1 
            WI(8) =  .53244713977759919092D-1 
            WI(9) =  .58860144245324817310D-1 
            WI(10) =  .64039797355015489556D-1 
            WI(11) =  .68745323835736442614D-1 
            WI(12) =  .72941885005653061354D-1 
            WI(13) =  .76598410645870674529D-1 
            WI(14) =  .79687828912071601909D-1 
            WI(15) =  .82187266704339709517D-1 
            WI(16) =  .84078218979661934933D-1 
            WI(17) =  .85346685739338627492D-1 
            WI(18) =  .85983275670394747490D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.40) THEN 
            XI(1) = -.99823770971055920035D0 
            XI(2) = -.99072623869945700645D0 
            XI(3) = -.97725994998377426266D0 
            XI(4) = -.95791681921379165580D0 
            XI(5) = -.93281280827867653336D0 
            XI(6) = -.90209880696887429673D0 
            XI(7) = -.86595950321225950382D0 
            XI(8) = -.82461223083331166320D0 
            XI(9) = -.77830565142651938769D0 
            XI(10) = -.72731825518992710328D0 
            XI(11) = -.67195668461417954838D0 
            XI(12) = -.61255388966798023795D0 
            XI(13) = -.54946712509512820208D0 
            XI(14) = -.48307580168617871291D0 
            XI(15) = -.41377920437160500152D0 
            XI(16) = -.34199409082575847301D0 
            XI(17) = -.26815218500725368114D0 
            XI(18) = -.19269758070137109972D0 
            XI(19) = -.11608407067525520848D0 
            XI(20) = -.38772417506050821933D-1 
            WI(1) =  .45212770985331912585D-2 
            WI(2) =  .10498284531152813615D-1 
            WI(3) =  .16421058381907888713D-1 
            WI(4) =  .22245849194166957262D-1 
            WI(5) =  .27937006980023401098D-1 
            WI(6) =  .33460195282547847393D-1 
            WI(7) =  .38782167974472017640D-1 
            WI(8) =  .43870908185673271992D-1 
            WI(9) =  .48695807635072232061D-1 
            WI(10) =  .53227846983936824355D-1 
            WI(11) =  .57439769099391551367D-1 
            WI(12) =  .61306242492928939167D-1 
            WI(13) =  .64804013456601038075D-1 
            WI(14) =  .67912045815233903826D-1 
            WI(15) =  .70611647391286779695D-1 
            WI(16) =  .72886582395804059061D-1 
            WI(17) =  .74723169057968264200D-1 
            WI(18) =  .76110361900626242372D-1 
            WI(19) =  .77039818164247965588D-1 
            WI(20) =  .77505947978424811264D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.44) THEN 
            XI(1) = -.99854020063677422494D0 
            XI(2) = -.99231639213851580848D0 
            XI(3) = -.98115183307791396666D0 
            XI(4) = -.96509965042249313939D0 
            XI(5) = -.94423950911819409920D0 
            XI(6) = -.91867525998417577432D0 
            XI(7) = -.88853423828604320234D0 
            XI(8) = -.85396659500471037873D0 
            XI(9) = -.81514453964513501049D0 
            XI(10) = -.77226147924875589902D0 
            XI(11) = -.72553105366071700261D0 
            XI(12) = -.67518607066612236533D0 
            XI(13) = -.62147734590357584780D0 
            XI(14) = -.56467245318547076842D0 
            XI(15) = -.50505439138820231798D0 
            XI(16) = -.44292017452541148383D0 
            XI(17) = -.37857935201470713251D0 
            XI(18) = -.31235246650278581224D0 
            XI(19) = -.24456945692820125151D0 
            XI(20) = -.17556801477551678575D0 
            XI(21) = -.10569190170865324712D0 
            XI(22) = -.35289236964135359058D-1 
            WI(1) =  .37454048031127775152D-2 
            WI(2) =  .87004813675248441226D-2 
            WI(3) =  .13619586755579985520D-1 
            WI(4) =  .18471481736814749172D-1 
            WI(5) =  .23231481902019210629D-1 
            WI(6) =  .27875782821281010081D-1 
            WI(7) =  .32381222812069820881D-1 
            WI(8) =  .36725347813808873643D-1 
            WI(9) =  .40886512310346218908D-1 
            WI(10) =  .44843984081970031446D-1 
            WI(11) =  .48578046448352037528D-1 
            WI(12) =  .52070096091704461881D-1 
            WI(13) =  .55302735563728052549D-1 
            WI(14) =  .58259859877595495334D-1 
            WI(15) =  .60926736701561968039D-1 
            WI(16) =  .63290079733203854950D-1 
            WI(17) =  .65338114879181434984D-1 
            WI(18) =  .67060638906293652396D-1 
            WI(19) =  .68449070269366660985D-1 
            WI(20) =  .69496491861572578037D-1 
            WI(21) =  .70197685473558212587D-1 
            WI(22) =  .70549157789354068811D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.48) THEN 
            XI(1) = -.99877100725242611860D0 
            XI(2) = -.99353017226635075755D0 
            XI(3) = -.98412458372282685774D0 
            XI(4) = -.97059159254624725046D0 
            XI(5) = -.95298770316043086072D0 
            XI(6) = -.93138669070655433311D0 
            XI(7) = -.90587913671556967282D0 
            XI(8) = -.87657202027424788591D0 
            XI(9) = -.84358826162439353071D0 
            XI(10) = -.80706620402944262708D0 
            XI(11) = -.76715903251574033925D0 
            XI(12) = -.72403413092381465467D0 
            XI(13) = -.67787237963266390521D0 
            XI(14) = -.62886739677651362400D0 
            XI(15) = -.57722472608397270382D0 
            XI(16) = -.52316097472223303368D0 
            XI(17) = -.46690290475095840454D0 
            XI(18) = -.40868648199071672992D0 
            XI(19) = -.34875588629216073816D0 
            XI(20) = -.28736248735545557674D0 
            XI(21) = -.22476379039468906122D0 
            XI(22) = -.16122235606889171806D0 
            XI(23) = -.97004699209462698930D-1 
            XI(24) = -.32380170962869362033D-1 
            WI(1) =  .31533460523058386327D-2 
            WI(2) =  .73275539012762621024D-2 
            WI(3) =  .11477234579234539490D-1 
            WI(4) =  .15579315722943848728D-1 
            WI(5) =  .19616160457355527814D-1 
            WI(6) =  .23570760839324379141D-1 
            WI(7) =  .27426509708356948200D-1 
            WI(8) =  .31167227832798088902D-1 
            WI(9) =  .34777222564770438893D-1 
            WI(10) =  .38241351065830706317D-1 
            WI(11) =  .41545082943464749214D-1 
            WI(12) =  .44674560856694280419D-1 
            WI(13) =  .47616658492490474826D-1 
            WI(14) =  .50359035553854474958D-1 
            WI(15) =  .52890189485193667096D-1 
            WI(16) =  .55199503699984162868D-1 
            WI(17) =  .57277292100403215705D-1 
            WI(18) =  .59114839698395635746D-1 
            WI(19) =  .60704439165893880053D-1 
            WI(20) =  .62039423159892663904D-1 
            WI(21) =  .63114192286254025657D-1 
            WI(22) =  .63924238584648186624D-1 
            WI(23) =  .64466164435950082207D-1 
            WI(24) =  .64737696812683922503D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.52) THEN 
            XI(1) = -.99895111110395027809D0 
            XI(2) = -.99447759092921602925D0 
            XI(3) = -.98644619565154984065D0 
            XI(4) = -.97488388422174450314D0 
            XI(5) = -.95983182693308655253D0 
            XI(6) = -.94134385364135905684D0 
            XI(7) = -.91948612891642453989D0 
            XI(8) = -.89433689053449532252D0 
            XI(9) = -.86598616284606758524D0 
            XI(10) = -.83453543232673453496D0 
            XI(11) = -.80009728343046832433D0 
            XI(12) = -.76279499519374496028D0 
            XI(13) = -.72276209974998319368D0 
            XI(14) = -.68014190422716770209D0 
            XI(15) = -.63508697769524592430D0 
            XI(16) = -.58775860497957906990D0 
            XI(17) = -.53832620928582743838D0 
            XI(18) = -.48696674569809607778D0 
            XI(19) = -.43386406771876167031D0 
            XI(20) = -.37920826911609366925D0 
            XI(21) = -.32319500343480782550D0 
            XI(22) = -.26602478360500182747D0 
            XI(23) = -.20790226415636605969D0 
            XI(24) = -.14903550860694918049D0 
            XI(25) = -.89635244648900565489D-1 
            XI(26) = -.29914109797338766044D-1 
            WI(1) =  .26913169500471111190D-2 
            WI(2) =  .62555239629732768997D-2 
            WI(3) =  .98026345794627520620D-2 
            WI(4) =  .13315114982340960657D-1 
            WI(5) =  .16780023396300735678D-1 
            WI(6) =  .20184891507980792203D-1 
            WI(7) =  .23517513553984461590D-1 
            WI(8) =  .26765953746504013449D-1 
            WI(9) =  .29918581147143946641D-1 
            WI(10) =  .32964109089718797915D-1 
            WI(11) =  .35891634835097232942D-1 
            WI(12) =  .38690678310423978985D-1 
            WI(13) =  .41351219500560271679D-1 
            WI(14) =  .43863734259000407995D-1 
            WI(15) =  .46219228372784793508D-1 
            WI(16) =  .48409269744074896854D-1 
            WI(17) =  .50426018566342377218D-1 
            WI(18) =  .52262255383906993034D-1 
            WI(19) =  .53911406932757264751D-1 
            WI(20) =  .55367569669302652549D-1 
            WI(21) =  .56625530902368597191D-1 
            WI(22) =  .57680787452526827654D-1 
            WI(23) =  .58529561771813868550D-1 
            WI(24) =  .59168815466042970369D-1 
            WI(25) =  .59596260171248158258D-1 
            WI(26) =  .59810365745291860248D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.56) THEN 
            XI(1) = -.99909434380146558435D0 
            XI(2) = -.99523122608106974722D0 
            XI(3) = -.98829371554016151109D0 
            XI(4) = -.97830170914025638338D0 
            XI(5) = -.96528590190549018363D0 
            XI(6) = -.94928647956196263565D0 
            XI(7) = -.93035288024749630055D0 
            XI(8) = -.90854362042065549085D0 
            XI(9) = -.88392610832782754079D0 
            XI(10) = -.85657643376274863540D0 
            XI(11) = -.82657913214288165167D0 
            XI(12) = -.79402692289386649803D0 
            XI(13) = -.75902042270512890220D0 
            XI(14) = -.72166783445018808352D0 
            XI(15) = -.68208461269447045550D0 
            XI(16) = -.64039310680700689427D0 
            XI(17) = -.59672218277066332010D0 
            XI(18) = -.55120682485553461875D0 
            XI(19) = -.50398771838438171420D0 
            XI(20) = -.45521081487845957895D0 
            XI(21) = -.40502688092709127812D0 
            XI(22) = -.35359103217495452097D0 
            XI(23) = -.30106225386722066905D0 
            XI(24) = -.24760290943433720397D0 
            XI(25) = -.19337823863527525824D0 
            XI(26) = -.13855584681037624201D0 
            XI(27) = -.83305186822435374440D-1 
            XI(28) = -.27797035287275437094D-1 
            WI(1) =  .23238553757732155011D-2 
            WI(2) =  .54025222460153377613D-2 
            WI(3) =  .84690631633078876616D-2 
            WI(4) =  .11509824340383382174D-1 
            WI(5) =  .14515089278021471808D-1 
            WI(6) =  .17475512911400946505D-1 
            WI(7) =  .20381929882402572635D-1 
            WI(8) =  .23225351562565316937D-1 
            WI(9) =  .25996987058391952192D-1 
            WI(10) =  .28688268473822741730D-1 
            WI(11) =  .31290876747310447868D-1 
            WI(12) =  .33796767115611761295D-1 
            WI(13) =  .36198193872315186036D-1 
            WI(14) =  .38487734259247662487D-1 
            WI(15) =  .40658311384744517880D-1 
            WI(16) =  .42703216084667086511D-1 
            WI(17) =  .44616127652692283213D-1 
            WI(18) =  .46391133373001896762D-1 
            WI(19) =  .48022746793600258121D-1 
            WI(20) =  .49505924683047578920D-1 
            WI(21) =  .50836082617798480560D-1 
            WI(22) =  .52009109151741399843D-1 
            WI(23) =  .53021378524010763968D-1 
            WI(24) =  .53869761865714485709D-1 
            WI(25) =  .54551636870889421062D-1 
            WI(26) =  .55064895901762425796D-1 
            WI(27) =  .55407952503245123218D-1 
            WI(28) =  .55579746306514395846D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.60) THEN 
            XI(1) = -.99921012322743602203D0 
            XI(2) = -.99584052511883817388D0 
            XI(3) = -.98978789522222171737D0 
            XI(4) = -.98106720175259818562D0 
            XI(5) = -.96970178876505273372D0 
            XI(6) = -.95572225583999610740D0 
            XI(7) = -.93916627611642324950D0 
            XI(8) = -.92007847617762755286D0 
            XI(9) = -.89851031081004594194D0 
            XI(10) = -.87451992264689831513D0 
            XI(11) = -.84817198478592963249D0 
            XI(12) = -.81953752616214575937D0 
            XI(13) = -.78869373993226405457D0 
            XI(14) = -.75572377530658568687D0 
            XI(15) = -.72071651335573039944D0 
            XI(16) = -.68376632738135543722D0 
            XI(17) = -.64497282848947706781D0 
            XI(18) = -.60444059704851036344D0 
            XI(19) = -.56227890075394453918D0 
            XI(20) = -.51860140005856974742D0 
            XI(21) = -.47352584176170711111D0 
            XI(22) = -.42717374158307838931D0 
            XI(23) = -.37967005657679797715D0 
            XI(24) = -.33114284826844819425D0 
            XI(25) = -.28172293742326169169D0 
            XI(26) = -.23154355137602933801D0 
            XI(27) = -.18073996487342541724D0 
            XI(28) = -.12944913539694500315D0 
            XI(29) = -.77809333949536569419D-1 
            XI(30) = -.25959772301247798589D-1 
            WI(1) =  .20268119688737584964D-2 
            WI(2) =  .47127299269535686409D-2 
            WI(3) =  .73899311633454555315D-2 
            WI(4) =  .10047557182287984358D-1 
            WI(5) =  .12678166476815960131D-1 
            WI(6) =  .15274618596784799307D-1 
            WI(7) =  .17829901014207720260D-1 
            WI(8) =  .20337120729457286775D-1 
            WI(9) =  .22789516943997819864D-1 
            WI(10) =  .25180477621521248380D-1 
            WI(11) =  .27503556749924791635D-1 
            WI(12) =  .29752491500788945241D-1 
            WI(13) =  .31921219019296328949D-1 
            WI(14) =  .34003892724946422835D-1 
            WI(15) =  .35994898051084503067D-1 
            WI(16) =  .37888867569243444031D-1 
            WI(17) =  .39680695452380799470D-1 
            WI(18) =  .41365551235584755613D-1 
            WI(19) =  .42938892835935641954D-1 
            WI(20) =  .44396478795787113328D-1 
            WI(21) =  .45734379716114486647D-1 
            WI(22) =  .46948988848912204847D-1 
            WI(23) =  .48037031819971180964D-1 
            WI(24) =  .48995575455756835389D-1 
            WI(25) =  .49822035690550181011D-1 
            WI(26) =  .50514184532509374598D-1 
            WI(27) =  .51070156069855627405D-1 
            WI(28) =  .51488451500980933995D-1 
            WI(29) =  .51767943174910187544D-1 
            WI(30) =  .51907877631220639733D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.64) THEN 
            XI(1) = -.99930504173577213946D0 
            XI(2) = -.99634011677195527935D0 
            XI(3) = -.99101337147674432074D0 
            XI(4) = -.98333625388462595693D0 
            XI(5) = -.97332682778991096374D0 
            XI(6) = -.96100879965205371892D0 
            XI(7) = -.94641137485840281606D0 
            XI(8) = -.92956917213193957582D0 
            XI(9) = -.91052213707850280576D0 
            XI(10) = -.88931544599511410585D0 
            XI(11) = -.86599939815409281976D0 
            XI(12) = -.84062929625258036275D0 
            XI(13) = -.81326531512279755974D0 
            XI(14) = -.78397235894334140761D0 
            XI(15) = -.75281990726053189661D0 
            XI(16) = -.71988185017161082685D0 
            XI(17) = -.68523631305423324256D0 
            XI(18) = -.64896547125465733986D0 
            XI(19) = -.61115535517239325025D0 
            XI(20) = -.57189564620263403428D0 
            XI(21) = -.53127946401989454566D0 
            XI(22) = -.48940314570705295748D0 
            XI(23) = -.44636601725346408798D0 
            XI(24) = -.40227015796399160370D0 
            XI(25) = -.35722015833766811595D0 
            XI(26) = -.31132287199021095616D0 
            XI(27) = -.26468716220876741637D0 
            XI(28) = -.21742364374000708415D0 
            XI(29) = -.16964442042399281804D0 
            XI(30) = -.12146281929612055447D0 
            XI(31) = -.72993121787799039450D-1 
            XI(32) = -.24350292663424432509D-1 
            WI(1) =  .17832807216964329473D-2 
            WI(2) =  .41470332605624676353D-2 
            WI(3) =  .65044579689783628561D-2 
            WI(4) =  .88467598263639477230D-2 
            WI(5) =  .11168139460131128819D-1 
            WI(6) =  .13463047896718642598D-1 
            WI(7) =  .15726030476024719322D-1 
            WI(8) =  .17951715775697343085D-1 
            WI(9) =  .20134823153530209372D-1 
            WI(10) =  .22270173808383254159D-1 
            WI(11) =  .24352702568710873338D-1 
            WI(12) =  .26377469715054658672D-1 
            WI(13) =  .28339672614259483228D-1 
            WI(14) =  .30234657072402478868D-1 
            WI(15) =  .32057928354851553585D-1 
            WI(16) =  .33805161837141609392D-1 
            WI(17) =  .35472213256882383811D-1 
            WI(18) =  .37055128540240046040D-1 
            WI(19) =  .38550153178615629129D-1 
            WI(20) =  .39953741132720341387D-1 
            WI(21) =  .41262563242623528610D-1 
            WI(22) =  .42473515123653589007D-1 
            WI(23) =  .43583724529323453377D-1 
            WI(24) =  .44590558163756563060D-1 
            WI(25) =  .45491627927418144480D-1 
            WI(26) =  .46284796581314417296D-1 
            WI(27) =  .46968182816210017325D-1 
            WI(28) =  .47540165714830308662D-1 
            WI(29) =  .47999388596458307728D-1 
            WI(30) =  .48344762234802957170D-1 
            WI(31) =  .48575467441503426935D-1 
            WI(32) =  .48690957009139720383D-1 
            GO TO 10 
          END IF 
        END IF 
      END IF 
               
      WRITE (6,FMT=*) 'CASE N=',N,' IS NOT PROVIDED' 
      STOP 'GAULEGE' 
               
 10   CONTINUE 
      DO I = N, (N+3)/2,-1 
        XI(I) = -XI(N+1-I) 
        WI(I) = WI(N+1-I) 
      END DO 
      RETURN   
               
      END 
c     ************************************************ 
      SUBROUTINE GAUFDT1(XI,WI,XSCALE,N,NCASE) 
c     ************************************************ 
C     .. Scalar Arguments .. 
      INTEGER N,NCASE 
      DOUBLE PRECISION XSCALE 
C     .. 
C     .. Array Arguments .. 
      DOUBLE PRECISION WI(*),XI(*) 
C     .. 
      XSCALE = 13.D0 * LOG(10.D0) 
      IF(NCASE.LT.0) RETURN 
      IF(NCASE.EQ.1) THEN 
          IF(N.EQ.1) THEN 
            XI(1) = -.5109128660907D0 
            WI(1) =  .9768438464874D0 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) = -.7923105545187D0 
            XI(2) = -.2261406359238D0 
            WI(1) =  .4913330637194D0 
            WI(2) =  .4855107827680D0 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) = -.8887753610196D0 
            XI(2) = -.5066546172925D0 
            XI(3) = -.1262095386730D0 
            WI(1) =  .2741283742272D0 
            WI(2) =  .4383133292069D0 
            WI(3) =  .2644021430533D0 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) = -.9312627444884D0 
            XI(2) = -.6733132754208D0 
            XI(3) = -.3369401943852D0 
            XI(4) = -.8039183115211D-1 
            WI(1) =  .1721853585682D0 
            WI(2) =  .3227514235415D0 
            WI(3) =  .3222708935931D0 
            WI(4) =  .1596361707845D0 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) = -.9534611319602D0 
            XI(2) = -.7710680603806D0 
            XI(3) = -.5040198559294D0 
            XI(4) = -.2372801245670D0 
            XI(5) = -.5534141311161D-1 
            WI(1) =  .1175252568581D0 
            WI(2) =  .2374010163917D0 
            WI(3) =  .2820759710686D0 
            WI(4) =  .2363146276642D0 
            WI(5) =  .1035269745049D0 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) = -.9664553095016D0 
            XI(2) = -.8317143942144D0 
            XI(3) = -.6218218168207D0 
            XI(4) = -.3848603191969D0 
            XI(5) = -.1755092818482D0 
            XI(6) = -.4010259704260D-1 
            WI(1) =  .8510236119667D-1 
            WI(2) =  .1791940779968D0 
            WI(3) =  .2323849772049D0 
            WI(4) =  .2322331422528D0 
            WI(5) =  .1770532097484D0 
            WI(6) =  .7087607808775D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) = -.9746959183855D0 
            XI(2) = -.8714881943800D0 
            XI(3) = -.7045921763649D0 
            XI(4) = -.5028422288757D0 
            XI(5) = -.3012024790023D0 
            XI(6) = -.1350069020325D0 
            XI(7) = -.3018676730141D-1 
            WI(1) =  .6438111280531D-1 
            WI(2) =  .1390682763472D0 
            WI(3) =  .1898298263499D0 
            WI(4) =  .2077427561299D0 
            WI(5) =  .1895352952046D0 
            WI(6) =  .1355406480747D0 
            WI(7) =  .5074593157573D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) = -.9802418462361D0 
            XI(2) = -.8988303733697D0 
            XI(3) = -.7639306564317D0 
            XI(4) = -.5937365157922D0 
            XI(5) = -.4112573894568D0 
            XI(6) = -.2412561664959D0 
            XI(7) = -.1070797467014D0 
            XI(8) = -.2342347462062D-1 
            WI(1) =  .5036711023742D-1 
            WI(2) =  .1106455124983D0 
            WI(3) =  .1560773181905D0 
            WI(4) =  .1804237450921D0 
            WI(5) =  .1803551283255D0 
            WI(6) =  .1555676018609D0 
            WI(7) =  .1056592165344D0 
            WI(8) =  .3774821374819D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) = -.9841492890000D0 
            XI(2) = -.9183722658097D0 
            XI(3) = -.8075292971039D0 
            XI(4) = -.6636091313827D0 
            XI(5) = -.5022164681903D0 
            XI(6) = -.3408758982388D0 
            XI(7) = -.1972261807832D0 
            XI(8) = -.8698549956060D-1 
            XI(9) = -.1863810149862D-1 
            WI(1) =  .4046057647775D-1 
            WI(2) =  .8993034157261D-1 
            WI(3) =  .1297333942514D0 
            WI(4) =  .1554775737104D0 
            WI(5) =  .1643557878031D0 
            WI(6) =  .1553533406395D0 
            WI(7) =  .1288804690946D0 
            WI(8) =  .8364377243437D-1 
            WI(9) =  .2900859050366D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) = -.9870043750732D0 
            XI(2) = -.9327962980849D0 
            XI(3) = -.8403349479830D0 
            XI(4) = -.7178160324805D0 
            XI(5) = -.5761292194322D0 
            XI(6) = -.4278780351856D0 
            XI(7) = -.2862815567961D0 
            XI(8) = -.1640994505867D0 
            XI(9) = -.7201582312137D-1 
            XI(10) = -.1514721816275D-1 
            WI(1) =  .3320504491102D-1 
            WI(2) =  .7443214903849D-1 
            WI(3) =  .1091105873073D0 
            WI(4) =  .1340959158476D0 
            WI(5) =  .1471584549392D0 
            WI(6) =  .1471219511931D0 
            WI(7) =  .1339026068657D0 
            WI(8) =  .1077678855356D0 
            WI(9) =  .6713481476306D-1 
            WI(10) =  .2291443608624D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.11) THEN 
            XI(1) = -.9891531714106D0 
            XI(2) = -.9437329869853D0 
            XI(3) = -.8655588313823D0 
            XI(4) = -.7604113310305D0 
            XI(5) = -.6360898792672D0 
            XI(6) = -.5018214740063D0 
            XI(7) = -.3675818473679D0 
            XI(8) = -.2433862943135D0 
            XI(9) = -.1386168426609D0 
            XI(10) = -.6054879102978D-1 
            XI(11) = -.1253307389298D-1 
            WI(1) =  .2773494321399D-1 
            WI(2) =  .6256563690393D-1 
            WI(3) =  .9281057734757D-1 
            WI(4) =  .1161747115738D0 
            WI(5) =  .1309189029415D0 
            WI(6) =  .1359431384004D0 
            WI(7) =  .1308539489913D0 
            WI(8) =  .1158772140469D0 
            WI(9) =  .9085581892426D-1 
            WI(10) =  .5458381451816D-1 
            WI(11) =  .1852513962568D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.12) THEN 
            XI(1) = -.9908105293775D0 
            XI(2) = -.9522158575527D0 
            XI(3) = -.8853292573928D0 
            XI(4) = -.7943385756151D0 
            XI(5) = -.6849611680598D0 
            XI(6) = -.5640731036544D0 
            XI(7) = -.4392784621383D0 
            XI(8) = -.3184404393381D0 
            XI(9) = -.2092239166856D0 
            XI(10) = -.1186168071938D0 
            XI(11) = -.5156627870571D-1 
            XI(12) = -.1053047075463D-1 
            WI(1) =  .2351035859183D-1 
            WI(2) =  .5329411729364D-1 
            WI(3) =  .7977546583693D-1 
            WI(4) =  .1012468946998D0 
            WI(5) =  .1163546654199D0 
            WI(6) =  .1241461855420D0 
            WI(7) =  .1241244994860D0 
            WI(8) =  .1162589295619D0 
            WI(9) =  .1007914502439D0 
            WI(10) =  .7713953077551D-1 
            WI(11) =  .4492933759994D-1 
            WI(12) =  .1527241143603D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.13) THEN 
            XI(1) = -.9921154910566D0 
            XI(2) = -.9589241312511D0 
            XI(3) = -.9010902877663D0 
            XI(4) = -.8217188459938D0 
            XI(5) = -.7250883500889D0 
            XI(6) = -.6164101228650D0 
            XI(7) = -.5015474067491D0 
            XI(8) = -.3867022916899D0 
            XI(9) = -.2780931896007D0 
            XI(10) = -.1816557615014D0 
            XI(11) = -.1026354171714D0 
            XI(12) = -.4440059661172D-1 
            XI(13) = -.8965611883246D-2 
            WI(1) =  .2018065660224D-1 
            WI(2) =  .4592095555824D-1 
            WI(3) =  .6922539552576D-1 
            WI(4) =  .8880045981813D-1 
            WI(5) =  .1035872649346D0 
            WI(6) =  .1127867130337D0 
            WI(7) =  .1158995346085D0 
            WI(8) =  .1127481825878D0 
            WI(9) =  .1034499744811D0 
            WI(10) =  .8812206313248D-1 
            WI(11) =  .6590030510630D-1 
            WI(12) =  .3742111078364D-1 
            WI(13) =  .1280123031503D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.14) THEN 
            XI(1) = -.9931612320482D0 
            XI(2) = -.9643183469376D0 
            XI(3) = -.9138445840139D0 
            XI(4) = -.8440886977263D0 
            XI(5) = -.7583117925642D0 
            XI(6) = -.6605258343545D0 
            XI(7) = -.5553057117455D0 
            XI(8) = -.4475764892286D0 
            XI(9) = -.3423871865445D0 
            XI(10) = -.2446891798152D0 
            XI(11) = -.1591306153016D0 
            XI(12) = -.8966080371538D-1 
            XI(13) = -.3859676725500D-1 
            XI(14) = -.7721364952157D-2 
            WI(1) =  .1751023584177D-1 
            WI(2) =  .3996593453129D-1 
            WI(3) =  .6058738530112D-1 
            WI(4) =  .7837825436357D-1 
            WI(5) =  .9250374636668D-1 
            WI(6) =  .1023021160934D0 
            WI(7) =  .1073134681607D0 
            WI(8) =  .1072995393680D0 
            WI(9) =  .1022468197540D0 
            WI(10) =  .9230562260229D-1 
            WI(11) =  .7741127705855D-1 
            WI(12) =  .5661802444588D-1 
            WI(13) =  .3151873287980D-1 
            WI(14) =  .1088268972033D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.15) THEN 
            XI(1) = -.9940120743912D0 
            XI(2) = -.9687193593080D0 
            XI(3) = -.9243035328417D0 
            XI(4) = -.8625729853654D0 
            XI(5) = -.7860542835655D0 
            XI(6) = -.6978807851576D0 
            XI(7) = -.6016638741644D0 
            XI(8) = -.5013455962649D0 
            XI(9) = -.4010388792979D0 
            XI(10) = -.3048643643591D0 
            XI(11) = -.2167975624259D0 
            XI(12) = -.1405115878143D0 
            XI(13) = -.7897947092873D-1 
            XI(14) = -.3383461770351D-1 
            XI(15) = -.6716782786797D-2 
            WI(1) =  .1533611007732D-1 
            WI(2) =  .3509023012556D-1 
            WI(3) =  .5343799234184D-1 
            WI(4) =  .6960022318918D-1 
            WI(5) =  .8291280000308D-1 
            WI(6) =  .9282982829193D-1 
            WI(7) =  .9894430145634D-1 
            WI(8) =  .1010039809610D0 
            WI(9) =  .9891945953798D-1 
            WI(10) =  .9275407885453D-1 
            WI(11) =  .8262667476670D-1 
            WI(12) =  .6828883591397D-1 
            WI(13) =  .4890597711721D-1 
            WI(14) =  .2682855048191D-1 
            WI(15) =  .9364803368826D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.16) THEN 
            XI(1) = -.9947135731351D0 
            XI(2) = -.9723560763431D0 
            XI(3) = -.9329819733666D0 
            XI(4) = -.8780054762164D0 
            XI(5) = -.8094128679959D0 
            XI(6) = -.7296836651445D0 
            XI(7) = -.6417004843807D0 
            XI(8) = -.5486450565983D0 
            XI(9) = -.4538839480282D0 
            XI(10) = -.3608489273991D0 
            XI(11) = -.2729193098985D0 
            XI(12) = -.1933145507386D0 
            XI(13) = -.1249556757741D0 
            XI(14) = -.7007833678774D-1 
            XI(15) = -.2988278726349D-1 
            XI(16) = -.5894647441833D-2 
            WI(1) =  .1354265229821D-1 
            WI(2) =  .3104970368477D-1 
            WI(3) =  .4746124641101D-1 
            WI(4) =  .6215945259804D-1 
            WI(5) =  .7461100529784D-1 
            WI(6) =  .8436522222399D-1 
            WI(7) =  .9106893437111D-1 
            WI(8) =  .9447880932027D-1 
            WI(9) =  .9446933158486D-1 
            WI(10) =  .9103374369405D-1 
            WI(11) =  .8426145123807D-1 
            WI(12) =  .7420389905464D-1 
            WI(13) =  .6046505531216D-1 
            WI(14) =  .4246785913222D-1 
            WI(15) =  .2306157178408D-1 
            WI(16) =  .8143908482067D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.17) THEN 
            XI(1) = -.9952987192274D0 
            XI(2) = -.9753952585126D0 
            XI(3) = -.9402592922878D0 
            XI(4) = -.8910122020707D0 
            XI(5) = -.8292361304776D0 
            XI(6) = -.7569168452321D0 
            XI(7) = -.6763794091393D0 
            XI(8) = -.5902135190917D0 
            XI(9) = -.5011906200898D0 
            XI(10) = -.4121757271684D0 
            XI(11) = -.3260378902085D0 
            XI(12) = -.2455652304542D0 
            XI(13) = -.1733856468680D0 
            XI(14) = -.1118304243829D0 
            XI(15) = -.6258130568763D-1 
            XI(16) = -.2657040737960D-1 
            XI(17) = -.5213681651034D-2 
            WI(1) =  .1204601080660D-1 
            WI(2) =  .2766508112030D-1 
            WI(3) =  .4241873114889D-1 
            WI(4) =  .5581086677154D-1 
            WI(5) =  .6740924182599D-1 
            WI(6) =  .7684056836996D-1 
            WI(7) =  .8380129918484D-1 
            WI(8) =  .8806709293648D-1 
            WI(9) =  .8949966712180D-1 
            WI(10) =  .8805009556148D-1 
            WI(11) =  .8375442769845D-1 
            WI(12) =  .7669712802496D-1 
            WI(13) =  .6684564852296D-1 
            WI(14) =  .5371543769204D-1 
            WI(15) =  .3707079189357D-1 
            WI(16) =  .2000420299119D-1 
            WI(17) =  .7147554816345D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.18) THEN 
            XI(1) = -.9957918772368D0 
            XI(2) = -.9779606638027D0 
            XI(3) = -.9464196855176D0 
            XI(4) = -.9020692127485D0 
            XI(5) = -.8461844562110D0 
            XI(6) = -.7803732474450D0 
            XI(7) = -.7065293063493D0 
            XI(8) = -.6267777679327D0 
            XI(9) = -.5434142258700D0 
            XI(10) = -.4588391292972D0 
            XI(11) = -.3754898320740D0 
            XI(12) = -.2957735609909D0 
            XI(13) = -.2220055918778D0 
            XI(14) = -.1563452772053D0 
            XI(15) = -.1006558779721D0 
            XI(16) = -.5620760099923D-1 
            XI(17) = -.2376897350239D-1 
            XI(18) = -.4643554482141D-2 
            WI(1) =  .1078418937266D-1 
            WI(2) =  .2480245515909D-1 
            WI(3) =  .3812846857929D-1 
            WI(4) =  .5035929136611D-1 
            WI(5) =  .6114147113856D-1 
            WI(6) =  .7016441510723D-1 
            WI(7) =  .7716825753836D-1 
            WI(8) =  .8195112103542D-1 
            WI(9) =  .8437474055153D-1 
            WI(10) =  .8436799897717D-1 
            WI(11) =  .8192718812873D-1 
            WI(12) =  .7710645394706D-1 
            WI(13) =  .6996576886497D-1 
            WI(14) =  .6038709832572D-1 
            WI(15) =  .4786503186210D-1 
            WI(16) =  .3252844112608D-1 
            WI(17) =  .1749747798085D-1 
            WI(18) =  .6323977426486D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.19) THEN 
            XI(1) = -.9962113557307D0 
            XI(2) = -.9801456671909D0 
            XI(3) = -.9516792387203D0 
            XI(4) = -.9115429247262D0 
            XI(5) = -.8607755982614D0 
            XI(6) = -.8006921677757D0 
            XI(7) = -.7328490767180D0 
            XI(8) = -.6590039672325D0 
            XI(9) = -.5810702536883D0 
            XI(10) = -.5010677960056D0 
            XI(11) = -.4210711138840D0 
            XI(12) = -.3431570077916D0 
            XI(13) = -.2693543254524D0 
            XI(14) = -.2015979599826D0 
            XI(15) = -.1416712818570D0 
            XI(16) = -.9106372244552D-1 
            XI(17) = -.5074398756342D-1 
            XI(18) = -.2138034085009D-1 
            XI(19) = -.4161603871680D-2 
            WI(1) =  .9710557612878D-2 
            WI(2) =  .2236025283043D-1 
            WI(3) =  .3444997966060D-1 
            WI(4) =  .4564901975228D-1 
            WI(5) =  .5566590773911D-1 
            WI(6) =  .6424087156394D-1 
            WI(7) =  .7115161019716D-1 
            WI(8) =  .7621888066416D-1 
            WI(9) =  .7931103335952D-1 
            WI(10) =  .8034722253580D-1 
            WI(11) =  .7929887104261D-1 
            WI(12) =  .7618754621042D-1 
            WI(13) =  .7106952614829D-1 
            WI(14) =  .6396833221986D-1 
            WI(15) =  .5469108071671D-1 
            WI(16) =  .4277543273777D-1 
            WI(17) =  .2869008196708D-1 
            WI(18) =  .1542219208802D-1 
            WI(19) =  .5635447440753D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.20) THEN 
            XI(1) = -.9965711269466D0 
            XI(2) = -.9820217717486D0 
            XI(3) = -.9562044934769D0 
            XI(4) = -.9197186155290D0 
            XI(5) = -.8734187374598D0 
            XI(6) = -.8183901066575D0 
            XI(7) = -.7559227696909D0 
            XI(8) = -.6874812945194D0 
            XI(9) = -.6146704931715D0 
            XI(10) = -.5391979344934D0 
            XI(11) = -.4628341938511D0 
            XI(12) = -.3873719877592D0 
            XI(13) = -.3145857694320D0 
            XI(14) = -.2461939456003D0 
            XI(15) = -.1838230095580D0 
            XI(16) = -.1289509262038D0 
            XI(17) = -.8276821690707D-1 
            XI(18) = -.4602596265197D-1 
            XI(19) = -.1932854388469D-1 
            XI(20) = -.3750639642494D-2 
            WI(1) =  .8789501225222D-2 
            WI(2) =  .2026034582014D-1 
            WI(3) =  .3127363748215D-1 
            WI(4) =  .4155532030506D-1 
            WI(5) =  .5086313183303D-1 
            WI(6) =  .5897861601500D-1 
            WI(7) =  .6571139004202D-1 
            WI(8) =  .7090346450999D-1 
            WI(9) =  .7443287737962D-1 
            WI(10) =  .7621645121425D-1 
            WI(11) =  .7621148422456D-1 
            WI(12) =  .7441579293797D-1 
            WI(13) =  .7086322970433D-1 
            WI(14) =  .6560146887913D-1 
            WI(15) =  .5861156308723D-1 
            WI(16) =  .4964530863395D-1 
            WI(17) =  .3833498985328D-1 
            WI(18) =  .2543307964494D-1 
            WI(19) =  .1368822574952D-1 
            WI(20) =  .5053967946010D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.21) THEN 
            XI(1) = -.9968820028221D0 
            XI(2) = -.9836444684178D0 
            XI(3) = -.9601254036393D0 
            XI(4) = -.9268207849972D0 
            XI(5) = -.8844398820500D0 
            XI(6) = -.8338859589619D0 
            XI(7) = -.7762366435858D0 
            XI(8) = -.7127209127856D0 
            XI(9) = -.6446929251024D0 
            XI(10) = -.5736032339105D0 
            XI(11) = -.5009680233466D0 
            XI(12) = -.4283371171620D0 
            XI(13) = -.3572617091415D0 
            XI(14) = -.2892631661744D0 
            XI(15) = -.2258043007821D0 
            XI(16) = -.1682592856596D0 
            XI(17) = -.1178555177545D0 
            XI(18) = -.7554530574536D-1 
            XI(19) = -.4192480685144D-1 
            XI(20) = -.1755408617320D-1 
            XI(21) = -.3397441726851D-2 
            WI(1) =  .7993442734942D-2 
            WI(2) =  .1844187280174D-1 
            WI(3) =  .2851300731084D-1 
            WI(4) =  .3797773847248D-1 
            WI(5) =  .4663319778282D-1 
            WI(6) =  .5429465931274D-1 
            WI(7) =  .6079870452964D-1 
            WI(8) =  .6600658263041D-1 
            WI(9) =  .6980711902319D-1 
            WI(10) =  .7211902784400D-1 
            WI(11) =  .7289252345224D-1 
            WI(12) =  .7211001494838D-1 
            WI(13) =  .6978498550780D-1 
            WI(14) =  .6595481669245D-1 
            WI(15) =  .6065125853299D-1 
            WI(16) =  .5381203267431D-1 
            WI(17) =  .4515836295941D-1 
            WI(18) =  .3445184404889D-1 
            WI(19) =  .2265737189205D-1 
            WI(20) =  .1222685948931D-1 
            WI(21) =  .4558423846755D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.22) THEN 
            XI(1) = -.9971524552065D0 
            XI(2) = -.9850573461686D0 
            XI(3) = -.9635445261874D0 
            XI(4) = -.9330278804235D0 
            XI(5) = -.8941008566956D0 
            XI(6) = -.8475210953112D0 
            XI(7) = -.7941953325884D0 
            XI(8) = -.7351617030131D0 
            XI(9) = -.6715695451241D0 
            XI(10) = -.6046570758411D0 
            XI(11) = -.5357273790158D0 
            XI(12) = -.4661232193119D0 
            XI(13) = -.3972012903267D0 
            XI(14) = -.3303067093150D0 
            XI(15) = -.2667488768297D0 
            XI(16) = -.2077790999062D0 
            XI(17) = -.1545631327269D0 
            XI(18) = -.1081213640962D0 
            XI(19) = -.6921747151514D-1 
            XI(20) = -.3833851460233D-1 
            XI(21) = -.1600987371125D-1 
            XI(22) = -.3091714969025D-2 
            WI(1) =  .7300760925976D-2 
            WI(2) =  .1685687249612D-1 
            WI(3) =  .2609927773617D-1 
            WI(4) =  .3483488792776D-1 
            WI(5) =  .4289263618147D-1 
            WI(6) =  .5011546837092D-1 
            WI(7) =  .5636269718770D-1 
            WI(8) =  .6151263236048D-1 
            WI(9) =  .6546491239441D-1 
            WI(10) =  .6814242328898D-1 
            WI(11) =  .6949273671998D-1 
            WI(12) =  .6948897110296D-1 
            WI(13) =  .6812976741426D-1 
            WI(14) =  .6543700454467D-1 
            WI(15) =  .6144546102307D-1 
            WI(16) =  .5616621622188D-1 
            WI(17) =  .4949729449505D-1 
            WI(18) =  .4115561438595D-1 
            WI(19) =  .3104912925567D-1 
            WI(20) =  .2028120207134D-1 
            WI(21) =  .1098522581996D-1 
            WI(22) =  .4132654562619D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.23) THEN 
            XI(1) = -.9973891981343D0 
            XI(2) = -.9862950240082D0 
            XI(3) = -.9665436402971D0 
            XI(4) = -.9384830721034D0 
            XI(5) = -.9026135715068D0 
            XI(6) = -.8595751885808D0 
            XI(7) = -.8101360267203D0 
            XI(8) = -.7551784865737D0 
            XI(9) = -.6956835232545D0 
            XI(10) = -.6327131686988D0 
            XI(11) = -.5673916341111D0 
            XI(12) = -.5008853509417D0 
            XI(13) = -.4343823621412D0 
            XI(14) = -.3690715748428D0 
            XI(15) = -.3061225813703D0 
            XI(16) = -.2466668705563D0 
            XI(17) = -.1917796374833D0 
            XI(18) = -.1424528980577D0 
            XI(19) = -.9953539649606D-1 
            XI(20) = -.6364265569825D-1 
            XI(21) = -.3518533128276D-1 
            XI(22) = -.1465826893249D-1 
            XI(23) = -.2825349540366D-2 
            WI(1) =  .6694312186004D-2 
            WI(2) =  .1546714798046D-1 
            WI(3) =  .2397716851274D-1 
            WI(4) =  .3206042058030D-1 
            WI(5) =  .3957171150840D-1 
            WI(6) =  .4637680752616D-1 
            WI(7) =  .5235418640894D-1 
            WI(8) =  .5739710930640D-1 
            WI(9) =  .6141549613542D-1 
            WI(10) =  .6433751038823D-1 
            WI(11) =  .6611080572828D-1 
            WI(12) =  .6670338269974D-1 
            WI(13) =  .6610393549530D-1 
            WI(14) =  .6432121997344D-1 
            WI(15) =  .6138044107771D-1 
            WI(16) =  .5730935800617D-1 
            WI(17) =  .5209585632587D-1 
            WI(18) =  .4560554624751D-1 
            WI(19) =  .3757561264108D-1 
            WI(20) =  .2806169658539D-1 
            WI(21) =  .1823770846119D-1 
            WI(22) =  .9922287325162D-2 
            WI(23) =  .3764125387471D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.24) THEN 
            XI(1) = -.9975976078705D0 
            XI(2) = -.9873852749448D0 
            XI(3) = -.9691885813113D0 
            XI(4) = -.9433022290429D0 
            XI(5) = -.9101508141358D0 
            XI(6) = -.8702786322963D0 
            XI(7) = -.8243404429988D0 
            XI(8) = -.7730906703643D0 
            XI(9) = -.7173710142661D0 
            XI(10) = -.6580966466647D0 
            XI(11) = -.5962412189195D0 
            XI(12) = -.5328209368221D0 
            XI(13) = -.4688779922991D0 
            XI(14) = -.4054636917847D0 
            XI(15) = -.3436217238687D0 
            XI(16) = -.2843721688213D0 
            XI(17) = -.2286966826771D0 
            XI(18) = -.1775227977414D0 
            XI(19) = -.1316965379570D0 
            XI(20) = -.9192424120684D-1 
            XI(21) = -.5870608571019D-1 
            XI(22) = -.3239906948671D-1 
            XI(23) = -.1346892498469D-1 
            XI(24) = -.2591889549462D-2 
            WI(1) =  .6160362794625D-2 
            WI(2) =  .1424198243582D-1 
            WI(3) =  .2210189150556D-1 
            WI(4) =  .2959990959432D-1 
            WI(5) =  .3661205195499D-1 
            WI(6) =  .4302300498546D-1 
            WI(7) =  .4872742932934D-1 
            WI(8) =  .5363160201832D-1 
            WI(9) =  .5765493099690D-1 
            WI(10) =  .6073126268513D-1 
            WI(11) =  .6280994696181D-1 
            WI(12) =  .6385662769463D-1 
            WI(13) =  .6385370445484D-1 
            WI(14) =  .6280029255897D-1 
            WI(15) =  .6071097767137D-1 
            WI(16) =  .5761063396944D-1 
            WI(17) =  .5351683323851D-1 
            WI(18) =  .4839379507093D-1 
            WI(19) =  .4208457469555D-1 
            WI(20) =  .3436712880412D-1 
            WI(21) =  .2543384910905D-1 
            WI(22) =  .1647215317859D-1 
            WI(23) =  .9005905611284D-2 
            WI(24) =  .3442995167831D-2 
            GO TO 10 
          END IF 
      END IF 
      IF(NCASE.EQ.0) THEN 
          IF(N.EQ.1) THEN 
            XI(1) =  .3964003208080D-1 
            WI(1) =  .2315615351261D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) =  .2306902985843D-1 
            XI(2) =  .1200117368196D0 
            WI(1) =  .1919793238029D-1 
            WI(2) =  .3958221132328D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) =  .1617090409278D-1 
            XI(2) =  .8173906117058D-1 
            XI(3) =  .2146200566810D0 
            WI(1) =  .1547961006390D-1 
            WI(2) =  .7374634784621D-2 
            WI(3) =  .3019086640932D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) =  .1240330548527D-1 
            XI(2) =  .6265308655909D-1 
            XI(3) =  .1557270131679D0 
            XI(4) =  .3176060243367D0 
            WI(1) =  .1277957667463D-1 
            WI(2) =  .9230469587975D-2 
            WI(3) =  .1130018655620D-2 
            WI(4) =  .1608859438496D-4 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) =  .1003811897495D-1 
            XI(2) =  .5091778322262D-1 
            XI(3) =  .1240647210271D0 
            XI(4) =  .2401968367533D0 
            XI(5) =  .4255933309307D0 
            WI(1) =  .1081755215890D-1 
            WI(2) =  .1003824335015D-1 
            WI(3) =  .2191229326947D-2 
            WI(4) =  .1084216750917D-3 
            WI(5) =  .7070015242700D-6 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) =  .8419214312526D-2 
            XI(2) =  .4290542677693D-1 
            XI(3) =  .1036626910036D0 
            XI(4) =  .1962272714557D0 
            XI(5) =  .3317209836761D0 
            XI(6) =  .5369344820508D0 
            WI(1) =  .9349914474458D-2 
            WI(2) =  .1025049956540D-1 
            WI(3) =  .3234815331047D-2 
            WI(4) =  .3129509769995D-3 
            WI(5) =  .7945723331439D-5 
            WI(6) =  .2744137727977D-7 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) =  .7243475076271D-2 
            XI(2) =  .3706847890336D-1 
            XI(3) =  .8924302923359D-1 
            XI(4) =  .1668639382503D0 
            XI(5) =  .2762542247328D0 
            XI(6) =  .4282226640820D0 
            XI(7) =  .6507157421613D0 
            WI(1) =  .8218765972860D-2 
            WI(2) =  .1014209807392D-1 
            WI(3) =  .4141342055809D-2 
            WI(4) =  .6206212919979D-3 
            WI(5) =  .3283804171470D-4 
            WI(6) =  .4871012408837D-6 
            WI(7) =  .9750685824348D-9 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) =  .6351873322875D-2 
            XI(2) =  .3262146033475D-1 
            XI(3) =  .7844226063563D-1 
            XI(4) =  .1455803314071D0 
            XI(5) =  .2381660990067D0 
            XI(6) =  .3620749282058D0 
            XI(7) =  .5284591487012D0 
            XI(8) =  .7663675856836D0 
            WI(1) =  .7323763572600D-2 
            WI(2) =  .9868175129239D-2 
            WI(3) =  .4875698162453D-2 
            WI(4) =  .1000004185527D-2 
            WI(5) =  .8568783458906D-4 
            WI(6) =  .2798380204007D-5 
            WI(7) =  .2621558921958D-7 
            WI(8) =  .3241276155673D-10 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) =  .5653098803194D-2 
            XI(2) =  .2911927269433D-1 
            XI(3) =  .7001895793699D-1 
            XI(4) =  .1293343849758D0 
            XI(5) =  .2099700239391D0 
            XI(6) =  .3156690870141D0 
            XI(7) =  .4523502299163D0 
            XI(8) =  .6316238080591D0 
            XI(9) =  .8835071510579D0 
            WI(1) =  .6599654768092D-2 
            WI(2) =  .9514888558077D-2 
            WI(3) =  .5443426135563D-2 
            WI(4) =  .1416574205783D-2 
            WI(5) =  .1719618158725D-3 
            WI(6) =  .9442343528659D-5 
            WI(7) =  .2044084591031D-6 
            WI(8) =  .1276215702028D-8 
            WI(9) =  .1022457667154D-11 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) =  .5091053022353D-2 
            XI(2) =  .2628952735379D-1 
            XI(3) =  .6325067567147D-1 
            XI(4) =  .1164749757764D0 
            XI(5) =  .1880896430538D0 
            XI(6) =  .2807327791030D0 
            XI(7) =  .3980424467959D0 
            XI(8) =  .5461811324939D0 
            XI(9) =  .7371573222687D0 
            XI(10) =  .1001863079614D1 
            WI(1) =  .6002680000016D-2 
            WI(2) =  .9129927484784D-2 
            WI(3) =  .5865847563664D-2 
            WI(4) =  .1841565679595D-2 
            WI(5) =  .2918877656779D-3 
            WI(6) =  .2335275955249D-4 
            WI(7) =  .8789670704640D-6 
            WI(8) =  .1323488080632D-7 
            WI(9) =  .5734129620123D-10 
            WI(10) = 0.D0 
            GO TO 10 
          END IF 
      END IF 
               
      WRITE (6,FMT=*) 'CASE N=',N,' IS NOT PROVIDED' 
      STOP 'GAUFDT1' 
               
 10   CONTINUE 
      RETURN   
               
      END 
c     ************************************************ 
      SUBROUTINE GAUFDT2(XI,WI,XSCALE,N,NCASE) 
c     ************************************************ 
C     .. Scalar Arguments .. 
      INTEGER N,NCASE 
      DOUBLE PRECISION XSCALE 
C     .. 
C     .. Array Arguments .. 
      DOUBLE PRECISION WI(*),XI(*) 
C     .. 
C     XZERO = LOG(2.D0+SQRT(7.D0)) 
      XSCALE = 13.D0 * LOG(10.D0) 
      IF(NCASE.LT.0) RETURN 
      IF(NCASE.EQ.1) THEN 
          IF(N.EQ.1) THEN 
            XI(1) = -.4993009779093D0 
            WI(1) =  .1001162208915D1 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) = -.7880759830485D0 
            XI(2) = -.2095990894951D0 
            WI(1) =  .5013831766899D0 
            WI(2) =  .4997790322248D0 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) = -.8866909948658D0 
            XI(2) = -.4973978439804D0 
            XI(3) = -.1094516029170D0 
            WI(1) =  .2792662253525D0 
            WI(2) =  .4465715397876D0 
            WI(3) =  .2753244437746D0 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) = -.9299508686019D0 
            XI(2) = -.6670818296087D0 
            XI(3) = -.3242829383200D0 
            XI(4) = -.6392998689631D-1 
            WI(1) =  .1754711601365D0 
            WI(2) =  .3289029699960D0 
            WI(3) =  .3285419325967D0 
            WI(4) =  .1682461461856D0 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) = -.9524893676049D0 
            XI(2) = -.7662889445836D0 
            XI(3) = -.4936683091676D0 
            XI(4) = -.2212476679926D0 
            XI(5) = -.3844723786276D-1 
            WI(1) =  .1199791281210D0 
            WI(2) =  .2423553679130D0 
            WI(3) =  .2879657963321D0 
            WI(4) =  .2420066934614D0 
            WI(5) =  .1088552230873D0 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) = -.9656743548789D0 
            XI(2) = -.8277968724190D0 
            XI(3) = -.6130190647058D0 
            XI(4) = -.3705323433475D0 
            XI(5) = -.1560556190247D0 
            XI(6) = -.2178272803859D-1 
            WI(1) =  .8708358905675D-1 
            WI(2) =  .1833650910957D0 
            WI(3) =  .2377940002784D0 
            WI(4) =  .2376804454143D0 
            WI(5) =  .1831902224764D0 
            WI(6) =  .7204886059319D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) = -.9740455733798D0 
            XI(2) = -.8681855127110D0 
            XI(3) = -.6970012436153D0 
            XI(4) = -.4900677120834D0 
            XI(5) = -.2832197446702D0 
            XI(6) = -.1125571780738D0 
            XI(7) = -.9603570160278D-2 
            WI(1) =  .6603575765123D-1 
            WI(2) =  .1426419384380D0 
            WI(3) =  .1947068180815D0 
            WI(4) =  .2130826741420D0 
            WI(5) =  .1945354078656D0 
            WI(6) =  .1424083691680D0 
            WI(7) =  .4775124356835D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) = -.9796906792938D0 
            XI(2) = -.8960083873431D0 
            XI(3) = -.7573467744178D0 
            XI(4) = -.5824084231812D0 
            XI(5) = -.3948430615122D0 
            XI(6) = -.2200527346799D0 
            XI(7) = -.8234333530655D-1 
            XI(8) = -.5541407233260D-4 
            WI(1) =  .5177211081855D-1 
            WI(2) =  .1137315144861D0 
            WI(3) =  .1604290568028D0 
            WI(4) =  .1854519572095D0 
            WI(5) =  .1853883297038D0 
            WI(6) =  .1602423918527D0 
            WI(7) =  .1127498893929D0 
            WI(8) =  .3139695864836D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) = -.9836766835521D0 
            XI(2) = -.9159386304066D0 
            XI(3) = -.8017918467368D0 
            XI(4) = -.6535839003164D0 
            XI(5) = -.4873864286809D0 
            XI(6) = -.3212420612649D0 
            XI(7) = -.1732450864015D0 
            XI(8) = -.6062557785559D-1 
            XI(9) =  .7709435283282D-2 
            WI(1) =  .4166693031117D-1 
            WI(2) =  .9261124851787D-1 
            WI(3) =  .1335995953742D0 
            WI(4) =  .1601082187183D0 
            WI(5) =  .1692469484319D0 
            WI(6) =  .1599992684558D0 
            WI(7) =  .1334619773277D0 
            WI(8) =  .8996698872029D-1 
            WI(9) =  .2050103305748D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) = -.9865952274160D0 
            XI(2) = -.9306806361880D0 
            XI(3) = -.8353091560024D0 
            XI(4) = -.7089356658857D0 
            XI(5) = -.5627944670489D0 
            XI(6) = -.4098875383234D0 
            XI(7) = -.2638393766872D0 
            XI(8) = -.1377621067848D0 
            XI(9) = -.4446653675940D-1 
            XI(10) =  .1411773515063D-1 
            WI(1) =  .3425043929881D-1 
            WI(2) =  .7677516245830D-1 
            WI(3) =  .1125441999854D0 
            WI(4) =  .1383134947474D0 
            WI(5) =  .1517825704352D0 
            WI(6) =  .1517407570685D0 
            WI(7) =  .1381765228428D0 
            WI(8) =  .1124988901720D0 
            WI(9) =  .7168876505002D-1 
            WI(10) =  .1339140685642D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.11) THEN 
            XI(1) = -.9887959587466D0 
            XI(2) = -.9418800896061D0 
            XI(3) = -.8611321499978D0 
            XI(4) = -.7525240276533D0 
            XI(5) = -.6241134198016D0 
            XI(6) = -.4854332462528D0 
            XI(7) = -.3467885402025D0 
            XI(8) = -.2185050962872D0 
            XI(9) = -.1103275967907D0 
            XI(10) = -.3200316191782D-1 
            XI(11) =  .1942768197804D-1 
            WI(1) =  .2864831011840D-1 
            WI(2) =  .6462578785042D-1 
            WI(3) =  .9586583453619D-1 
            WI(4) =  .1199973270754D0 
            WI(5) =  .1352232253941D0 
            WI(6) =  .1404069513878D0 
            WI(7) =  .1351506074161D0 
            WI(8) =  .1198582528308D0 
            WI(9) =  .9586023621513D-1 
            WI(10) =  .5669748468077D-1 
            WI(11) =  .8828191409758D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.12) THEN 
            XI(1) = -.9904962857740D0 
            XI(2) = -.9505819160985D0 
            XI(3) = -.8814085843249D0 
            XI(4) = -.7873080269929D0 
            XI(5) = -.6741942215939D0 
            XI(6) = -.5491800500045D0 
            XI(7) = -.4201317128759D0 
            XI(8) = -.2951799480017D0 
            XI(9) = -.1822273690247D0 
            XI(10) = -.8875691167923D-1 
            XI(11) = -.2205009961021D-1 
            XI(12) =  .2383058382770D-1 
            WI(1) =  .2431431012437D-1 
            WI(2) =  .5511634270092D-1 
            WI(3) =  .8250253066600D-1 
            WI(4) =  .1047066043859D0 
            WI(5) =  .1203280363320D0 
            WI(6) =  .1283810504528D0 
            WI(7) =  .1283522461669D0 
            WI(8) =  .1202331802244D0 
            WI(9) =  .1046006368225D0 
            WI(10) =  .8233730737622D-1 
            WI(11) =  .4437449975846D-1 
            WI(12) =  .5915463904325D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.13) THEN 
            XI(1) = -.9918371374714D0 
            XI(2) = -.9574740621204D0 
            XI(3) = -.8975988657113D0 
            XI(4) = -.8154265686616D0 
            XI(5) = -.7153876043374D0 
            XI(6) = -.6028785515621D0 
            XI(7) = -.4839714695449D0 
            XI(8) = -.3650891658639D0 
            XI(9) = -.2526648781734D0 
            XI(10) = -.1528283102516D0 
            XI(11) = -.7154082494724D-1 
            XI(12) = -.1385762476436D-1 
            XI(13) =  .2748483342760D-1 
            WI(1) =  .2089310400765D-1 
            WI(2) =  .4754197055979D-1 
            WI(3) =  .7166859555264D-1 
            WI(4) =  .9193351843756D-1 
            WI(5) =  .1072401034051D0 
            WI(6) =  .1167605358337D0 
            WI(7) =  .1199774568966D0 
            WI(8) =  .1167100490252D0 
            WI(9) =  .1071337884729D0 
            WI(10) =  .9189577819659D-1 
            WI(11) =  .7100597460159D-1 
            WI(12) =  .3435520062962D-1 
            WI(13) =  .4046133295765D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.14) THEN 
            XI(1) = -.9929131132590D0 
            XI(2) = -.9630238255904D0 
            XI(3) = -.9107191261768D0 
            XI(4) = -.8384333840191D0 
            XI(5) = -.7495466571020D0 
            XI(6) = -.6482174281093D0 
            XI(7) = -.5391879663726D0 
            XI(8) = -.4275638684350D0 
            XI(9) = -.3185782346997D0 
            XI(10) = -.2173538479462D0 
            XI(11) = -.1287264887830D0 
            XI(12) = -.5759938404969D-1 
            XI(13) = -.6951890108897D-2 
            XI(14) =  .3052486090959D-1 
            WI(1) =  .1814552250386D-1 
            WI(2) =  .4141581405800D-1 
            WI(3) =  .6278501076101D-1 
            WI(4) =  .8122042153445D-1 
            WI(5) =  .9585669114268D-1 
            WI(6) =  .1060076886579D0 
            WI(7) =  .1111963462135D0 
            WI(8) =  .1111757574733D0 
            WI(9) =  .1059406958805D0 
            WI(10) =  .9575448130199D-1 
            WI(11) =  .8126690304366D-1 
            WI(12) =  .6121467735320D-1 
            WI(13) =  .2635203086025D-1 
            WI(14) =  .2830168130460D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.15) THEN 
            XI(1) = -.9937896289360D0 
            XI(2) = -.9675573554556D0 
            XI(3) = -.9214917685365D0 
            XI(4) = -.8574687462088D0 
            XI(5) = -.7781092006500D0 
            XI(6) = -.6866636125119D0 
            XI(7) = -.5868785733710D0 
            XI(8) = -.4828440441122D0 
            XI(9) = -.3788274311357D0 
            XI(10) = -.2791018991091D0 
            XI(11) = -.1877812258669D0 
            XI(12) = -.1087617448119D0 
            XI(13) = -.4614058844272D-1 
            XI(14) = -.1029795771349D-2 
            XI(15) =  .3306355157604D-1 
            WI(1) =  .1590582663785D-1 
            WI(2) =  .3639369358972D-1 
            WI(3) =  .5542272687265D-1 
            WI(4) =  .7218464873481D-1 
            WI(5) =  .8599046299116D-1 
            WI(6) =  .9627372050064D-1 
            WI(7) =  .1026119153122D0 
            WI(8) =  .1047430467413D0 
            WI(9) =  .1025755723746D0 
            WI(10) =  .9619570631762D-1 
            WI(11) =  .8591278977141D-1 
            WI(12) =  .7228971633521D-1 
            WI(13) =  .5255457996851D-1 
            WI(14) =  .2008319727665D-1 
            WI(15) =  .2024605490416D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.16) THEN 
            XI(1) = -.9945130900241D0 
            XI(2) = -.9713077366701D0 
            XI(3) = -.9304405957302D0 
            XI(4) = -.8733797514408D0 
            XI(5) = -.8021871896065D0 
            XI(6) = -.7194370094081D0 
            XI(7) = -.6281219360288D0 
            XI(8) = -.5315454796568D0 
            XI(9) = -.4332033615514D0 
            XI(10) = -.3366587794290D0 
            XI(11) = -.2454168531566D0 
            XI(12) = -.1628145578776D0 
            XI(13) = -.9207011758903D-1 
            XI(14) = -.3657725720625D-1 
            XI(15) =  .4108279240244D-2 
            XI(16) =  .3519405802483D-1 
            WI(1) =  .1405624191191D-1 
            WI(2) =  .3222715532812D-1 
            WI(3) =  .4926083342262D-1 
            WI(4) =  .6451589605242D-1 
            WI(5) =  .7743864922779D-1 
            WI(6) =  .8756111711458D-1 
            WI(7) =  .9451645760255D-1 
            WI(8) =  .9805174817400D-1 
            WI(9) =  .9803659186707D-1 
            WI(10) =  .9446778897153D-1 
            WI(11) =  .8747937759711D-1 
            WI(12) =  .7740779997805D-1 
            WI(13) =  .6459558061814D-1 
            WI(14) =  .4480699044642D-1 
            WI(15) =  .1526014521574D-1 
            WI(16) =  .1479835386700D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.17) THEN 
            XI(1) = -.9951171544652D0 
            XI(2) = -.9744450418513D0 
            XI(3) = -.9379522652786D0 
            XI(4) = -.8868037023709D0 
            XI(5) = -.8226428902958D0 
            XI(6) = -.7475327104182D0 
            XI(7) = -.6638886181251D0 
            XI(8) = -.5744011647012D0 
            XI(9) = -.4819499835555D0 
            XI(10) = -.3895121339332D0 
            XI(11) = -.3000680628440D0 
            XI(12) = -.2165097502091D0 
            XI(13) = -.1415768112484D0 
            XI(14) = -.7799616142058D-1 
            XI(15) = -.2847375988604D-1 
            XI(16) =  .8600130642497D-2 
            XI(17) =  .3699213597069D-1 
            WI(1) =  .1251122820765D-1 
            WI(2) =  .2873345001516D-1 
            WI(3) =  .4405668660904D-1 
            WI(4) =  .5796559025622D-1 
            WI(5) =  .7001110630450D-1 
            WI(6) =  .7980538058785D-1 
            WI(7) =  .8703293833088D-1 
            WI(8) =  .9146049533144D-1 
            WI(9) =  .9294413245089D-1 
            WI(10) =  .9143354984563D-1 
            WI(11) =  .8697515982533D-1 
            WI(12) =  .7972989735795D-1 
            WI(13) =  .7004510145447D-1 
            WI(14) =  .5787419779385D-1 
            WI(15) =  .3787817527505D-1 
            WI(16) =  .1160156270487D-1 
            WI(17) =  .1103556563995D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.18) THEN 
            XI(1) = -.9956267140399D0 
            XI(2) = -.9770956694323D0 
            XI(3) = -.9443168696304D0 
            XI(4) = -.8982260684201D0 
            XI(5) = -.8401487563854D0 
            XI(6) = -.7717561994724D0 
            XI(7) = -.6950169031250D0 
            XI(8) = -.6121400481913D0 
            XI(9) = -.5255122052594D0 
            XI(10) = -.4376291948023D0 
            XI(11) = -.3510252205277D0 
            XI(12) = -.2682016677771D0 
            XI(13) = -.1915608805832D0 
            XI(14) = -.1233848203415D0 
            XI(15) = -.6603293933369D-1 
            XI(16) = -.2150817645706D-1 
            XI(17) =  .1254661925479D-1 
            XI(18) =  .3851881037551D-1 
            WI(1) =  .1120745211439D-1 
            WI(2) =  .2577586723579D-1 
            WI(3) =  .3962474704417D-1 
            WI(4) =  .5233526458446D-1 
            WI(5) =  .6354000645443D-1 
            WI(6) =  .7291607684686D-1 
            WI(7) =  .8019325980233D-1 
            WI(8) =  .8516154787826D-1 
            WI(9) =  .8767698591063D-1 
            WI(10) =  .8766553785517D-1 
            WI(11) =  .8512520060457D-1 
            WI(12) =  .8013023043555D-1 
            WI(13) =  .7285952080744D-1 
            WI(14) =  .6364433384666D-1 
            WI(15) =  .5187778482122D-1 
            WI(16) =  .3173788009040D-1 
            WI(17) =  .8852283128467D-2 
            WI(18) =  .8382294539525D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.19) THEN 
            XI(1) = -.9960604974419D0 
            XI(2) = -.9793551127996D0 
            XI(3) = -.9497552888093D0 
            XI(4) = -.9080211005648D0 
            XI(5) = -.8552329699114D0 
            XI(6) = -.7927584222058D0 
            XI(7) = -.7222162337756D0 
            XI(8) = -.6454345226949D0 
            XI(9) = -.5644035610435D0 
            XI(10) = -.4812245364812D0 
            XI(11) = -.3980556807238D0 
            XI(12) = -.3170573210558D0 
            XI(13) = -.2403378338681D0 
            XI(14) = -.1699083257251D0 
            XI(15) = -.1077020766154D0 
            XI(16) = -.5578189215882D-1 
            XI(17) = -.1544354097586D-1 
            XI(18) =  .1602599942889D-1 
            XI(19) =  .3982300902589D-1 
            WI(1) =  .1009721601622D-1 
            WI(2) =  .2325056437975D-1 
            WI(3) =  .3582155669824D-1 
            WI(4) =  .4746624740149D-1 
            WI(5) =  .5788149647756D-1 
            WI(6) =  .6679710214611D-1 
            WI(7) =  .7398179894510D-1 
            WI(8) =  .7924905395573D-1 
            WI(9) =  .8246177120539D-1 
            WI(10) =  .8353566952255D-1 
            WI(11) =  .8244127695723D-1 
            WI(12) =  .7920539311613D-1 
            WI(13) =  .7391865089469D-1 
            WI(14) =  .6677370714763D-1 
            WI(15) =  .5804059901685D-1 
            WI(16) =  .4642351046795D-1 
            WI(17) =  .2637329136788D-1 
            WI(18) =  .6795871817994D-2 
            WI(19) =  .6474313802446D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.20) THEN 
            XI(1) = -.9964328144551D0 
            XI(2) = -.9812965863714D0 
            XI(3) = -.9544379789957D0 
            XI(4) = -.9164805870407D0 
            XI(5) = -.8683136288102D0 
            XI(6) = -.8110663341645D0 
            XI(7) = -.7460810685631D0 
            XI(8) = -.6748818574276D0 
            XI(9) = -.5991387605499D0 
            XI(10) = -.5206289172589D0 
            XI(11) = -.4411952232012D0 
            XI(12) = -.3627036992377D0 
            XI(13) = -.2870007193885D0 
            XI(14) = -.2158721837258D0 
            XI(15) = -.1510167919492D0 
            XI(16) = -.9410299151605D-1 
            XI(17) = -.4692615776487D-1 
            XI(18) = -.1010551670312D-1 
            XI(19) =  .1910190010718D-1 
            XI(20) =  .4094392275284D-1 
            WI(1) =  .9144047074604D-2 
            WI(2) =  .2107756615145D-1 
            WI(3) =  .3253500256207D-1 
            WI(4) =  .4323119888914D-1 
            WI(5) =  .5291406608769D-1 
            WI(6) =  .6135626960759D-1 
            WI(7) =  .6835966160645D-1 
            WI(8) =  .7375976586669D-1 
            WI(9) =  .7742954720591D-1 
            WI(10) =  .7928228179922D-1 
            WI(11) =  .7927343970276D-1 
            WI(12) =  .7740173006240D-1 
            WI(13) =  .7371107469567D-1 
            WI(14) =  .6830315580895D-1 
            WI(15) =  .6137950100193D-1 
            WI(16) =  .5308487452804D-1 
            WI(17) =  .4139095002135D-1 
            WI(18) =  .2176200342052D-1 
            WI(19) =  .5258380352508D-2 
            WI(20) =  .5076924697866D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.21) THEN 
            XI(1) = -.9967547503547D0 
            XI(2) = -.9829769731387D0 
            XI(3) = -.9584981059478D0 
            XI(4) = -.9238344452756D0 
            XI(5) = -.8797243258811D0 
            XI(6) = -.8271080418346D0 
            XI(7) = -.7671074254919D0 
            XI(8) = -.7010019103099D0 
            XI(9) = -.6302013202459D0 
            XI(10) = -.5562159420513D0 
            XI(11) = -.4806245425301D0 
            XI(12) = -.4050410660893D0 
            XI(13) = -.3310808080497D0 
            XI(14) = -.2603270186781D0 
            XI(15) = -.1943007444580D0 
            XI(16) = -.1344518772205D0 
            XI(17) = -.8224609514254D-1 
            XI(18) = -.3921242091803D-1 
            XI(19) = -.5365195973929D-2 
            XI(20) =  .2182769936625D-1 
            XI(21) =  .4191299513309D-1 
            WI(1) =  .8319671869631D-2 
            WI(2) =  .1919449988937D-1 
            WI(3) =  .2967657450958D-1 
            WI(4) =  .3952739982595D-1 
            WI(5) =  .4853578456468D-1 
            WI(6) =  .5650941161624D-1 
            WI(7) =  .6327812728926D-1 
            WI(8) =  .6869743205556D-1 
            WI(9) =  .7265149898530D-1 
            WI(10) =  .7505557429251D-1 
            WI(11) =  .7585768395516D-1 
            WI(12) =  .7503964268791D-1 
            WI(13) =  .7261780893192D-1 
            WI(14) =  .6864657993513D-1 
            WI(15) =  .6323662980630D-1 
            WI(16) =  .5658862839663D-1 
            WI(17) =  .4864524299053D-1 
            WI(18) =  .3671350567468D-1 
            WI(19) =  .1786143625561D-1 
            WI(20) =  .4105472410975D-2 
            WI(21) =  .4036029718319D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.22) THEN 
            XI(1) = -.9970349999801D0 
            XI(2) = -.9844410019559D0 
            XI(3) = -.9620408755562D0 
            XI(4) = -.9302656417675D0 
            XI(5) = -.8897333233807D0 
            XI(6) = -.8412329460788D0 
            XI(7) = -.7857088273479D0 
            XI(8) = -.7242421607567D0 
            XI(9) = -.6580300066565D0 
            XI(10) = -.5883620706924D0 
            XI(11) = -.5165957333790D0 
            XI(12) = -.4441298479434D0 
            XI(13) = -.3723778679169D0 
            XI(14) = -.3027409137003D0 
            XI(15) = -.2365817189912D0 
            XI(16) = -.1752036655103D0 
            XI(17) = -.1198593190304D0 
            XI(18) = -.7185417093593D-1 
            XI(19) = -.3243793508519D-1 
            XI(20) = -.1126122195353D-2 
            XI(21) =  .2424896853244D-1 
            XI(22) =  .4275553925603D-1 
            WI(1) =  .7601900940991D-2 
            WI(2) =  .1755216067188D-1 
            WI(3) =  .2717572548862D-1 
            WI(4) =  .3627151149410D-1 
            WI(5) =  .4466136085889D-1 
            WI(6) =  .5218170001074D-1 
            WI(7) =  .5868598742664D-1 
            WI(8) =  .6404744757192D-1 
            WI(9) =  .6816149235119D-1 
            WI(10) =  .7094771204188D-1 
            WI(11) =  .7235137625556D-1 
            WI(12) =  .7234441368985D-1 
            WI(13) =  .7092596819739D-1 
            WI(14) =  .6812340370420D-1 
            WI(15) =  .6399818256485D-1 
            WI(16) =  .5866897077280D-1 
            WI(17) =  .5231828203661D-1 
            WI(18) =  .4460895254589D-1 
            WI(19) =  .3236514083297D-1 
            WI(20) =  .1460944779784D-1 
            WI(21) =  .3236218014421D-2 
            WI(22) =  .3248536455104D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.23) THEN 
            XI(1) = -.9972804605622D0 
            XI(2) = -.9857242312406D0 
            XI(3) = -.9651502630771D0 
            XI(4) = -.9359211288452D0 
            XI(5) = -.8985579949866D0 
            XI(6) = -.8537276788930D0 
            XI(7) = -.8022304210588D0 
            XI(8) = -.7449855649549D0 
            XI(9) = -.6830151714864D0 
            XI(10) = -.6174258317664D0 
            XI(11) = -.5493890069050D0 
            XI(12) = -.4801202637008D0 
            XI(13) = -.4108578083454D0 
            XI(14) = -.3428407487945D0 
            XI(15) = -.2772875841635D0 
            XI(16) = -.2153760817868D0 
            XI(17) = -.1582308155456D0 
            XI(18) = -.1069484525982D0 
            XI(19) = -.6269971646623D-1 
            XI(20) = -.2644066784965D-1 
            XI(21) =  .2685248261883D-2 
            XI(22) =  .2640493787353D-1 
            XI(23) =  .4349202409773D-1 
            WI(1) =  .6973123498381D-2 
            WI(2) =  .1611132245290D-1 
            WI(3) =  .2497571976601D-1 
            WI(4) =  .3339550857915D-1 
            WI(5) =  .4121942082792D-1 
            WI(6) =  .4830759894219D-1 
            WI(7) =  .5453341906360D-1 
            WI(8) =  .5978564575562D-1 
            WI(9) =  .6397038081099D-1 
            WI(10) =  .6701270831759D-1 
            WI(11) =  .6885799000659D-1 
            WI(12) =  .6947278024829D-1 
            WI(13) =  .6884537008658D-1 
            WI(14) =  .6698620589083D-1 
            WI(15) =  .6392963432217D-1 
            WI(16) =  .5974275322097D-1 
            WI(17) =  .5455044756824D-1 
            WI(18) =  .4849120456602D-1 
            WI(19) =  .4088451634935D-1 
            WI(20) =  .2834522118356D-1 
            WI(21) =  .1193071865608D-1 
            WI(22) =  .2576098039481D-2 
            WI(23) =  .2644207622250D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.24) THEN 
            XI(1) = -.9974966598769D0 
            XI(2) = -.9868552130057D0 
            XI(3) = -.9678939341301D0 
            XI(4) = -.9409199524316D0 
            XI(5) = -.9063757709855D0 
            XI(6) = -.8648286471538D0 
            XI(7) = -.8169609730280D0 
            XI(8) = -.7635590291011D0 
            XI(9) = -.7055000842822D0 
            XI(10) = -.6437380256626D0 
            XI(11) = -.5792877542002D0 
            XI(12) = -.5132086129401D0 
            XI(13) = -.4465871392439D0 
            XI(14) = -.3805194539726D0 
            XI(15) = -.3160936238020D0 
            XI(16) = -.2543724668689D0 
            XI(17) = -.1963784499565D0 
            XI(18) = -.1430895265525D0 
            XI(19) = -.9547911741088D-1 
            XI(20) = -.5459438304082D-1 
            XI(21) = -.2109143135335D-1 
            XI(22) =  .6126361712807D-2 
            XI(23) =  .2832948324453D-1 
            XI(24) =  .4413908647076D-1 
            WI(1) =  .6419219207366D-2 
            WI(2) =  .1484041206053D-1 
            WI(3) =  .2303054455719D-1 
            WI(4) =  .3084352516309D-1 
            WI(5) =  .3815013619484D-1 
            WI(6) =  .4483019183161D-1 
            WI(7) =  .5077389420670D-1 
            WI(8) =  .5588354209318D-1 
            WI(9) =  .6007510567005D-1 
            WI(10) =  .6327958462182D-1 
            WI(11) =  .6544411365197D-1 
            WI(12) =  .6653279037177D-1 
            WI(13) =  .6652721475205D-1 
            WI(14) =  .6542680363255D-1 
            WI(15) =  .6324934040312D-1 
            WI(16) =  .6003390803336D-1 
            WI(17) =  .5585282133268D-1 
            WI(18) =  .5083313531552D-1 
            WI(19) =  .4503586702456D-1 
            WI(20) =  .3740291895436D-1 
            WI(21) =  .2466440325160D-1 
            WI(22) =  .9744678312737D-2 
            WI(23) =  .2070623134384D-2 
            WI(24) =  .2174351376926D-3 
            GO TO 10 
          END IF 
      END IF 
      IF(NCASE.EQ.0) THEN 
          IF(N.EQ.1) THEN 
            XI(1) =  .1021589505488D0 
            WI(1) = -.1162208914753D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) =  .8407503043439D-1 
            XI(2) =  .1800671378671D0 
            WI(1) = -.9432607769038D-3 
            WI(2) = -.2189481378490D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) =  .7589114747476D-1 
            XI(2) =  .1410236389321D0 
            XI(3) =  .2736505742436D0 
            WI(1) = -.7285010897828D-3 
            WI(2) = -.4165152110159D-3 
            WI(3) = -.1719261395408D-4 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) =  .7111046129975D-1 
            XI(2) =  .1214130193630D0 
            XI(3) =  .2141024748807D0 
            XI(4) =  .3759288207714D0 
            WI(1) = -.5705471849468D-3 
            WI(2) = -.5250812216357D-3 
            WI(3) = -.6564480498667D-4 
            WI(4) = -.9357031835974D-6 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) =  .6794441088072D-1 
            XI(2) =  .1092621250134D0 
            XI(3) =  .1819604384409D0 
            XI(4) =  .2980101350362D0 
            XI(5) =  .4833773591336D0 
            WI(1) = -.4566574209351D-3 
            WI(2) = -.5699626929461D-3 
            WI(3) = -.1291446258307D-3 
            WI(4) = -.6402385627304D-5 
            WI(5) = -.4178941349624D-7 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) =  .6568137871527D-1 
            XI(2) =  .1008924701570D0 
            XI(3) =  .1611934042275D0 
            XI(4) =  .2536382183260D0 
            XI(5) =  .3890976421620D0 
            XI(6) =  .5942916230654D0 
            WI(1) = -.3731204569292D-3 
            WI(2) = -.5772886267199D-3 
            WI(3) = -.1926185534613D-3 
            WI(4) = -.1870424242957D-4 
            WI(5) = -.4753922941905D-6 
            WI(6) = -.1642918535157D-8 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) =  .6397785846566D-1 
            XI(2) =  .9473642752008D-1 
            XI(3) =  .1464851806735D0 
            XI(4) =  .2239487507919D0 
            XI(5) =  .3332972942952D0 
            XI(6) =  .4852442338751D0 
            XI(7) =  .7077237644642D0 
            WI(1) = -.3104089140185D-3 
            WI(2) = -.5640269911863D-3 
            WI(3) = -.2483008577601D-3 
            WI(4) = -.3745831296862D-4 
            WI(5) = -.1984323571825D-5 
            WI(6) = -.2945625559097D-7 
            WI(7) = -.5899188111549D-10 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) =  .6264639326881D-1 
            XI(2) =  .8999949609733D-1 
            XI(3) =  .1354477044954D0 
            XI(4) =  .2023955398790D0 
            XI(5) =  .2949298474489D0 
            XI(6) =  .4188147460788D0 
            XI(7) =  .5851843646069D0 
            XI(8) =  .8230829800884D0 
            WI(1) = -.2622588392148D-3 
            WI(2) = -.5401955916053D-3 
            WI(3) = -.2935158049191D-3 
            WI(4) = -.6084522180839D-4 
            WI(5) = -.5221203069393D-5 
            WI(6) = -.1706526567708D-6 
            WI(7) = -.1599500843881D-8 
            WI(8) = -.1978255320642D-11 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) =  .6157552120515D-1 
            XI(2) =  .8623174319813D-1 
            XI(3) =  .1268236487425D0 
            XI(4) =  .1859231185148D0 
            XI(5) =  .2664957952010D0 
            XI(6) =  .3721681088661D0 
            XI(7) =  .5088333481385D0 
            XI(8) =  .6880965434964D0 
            XI(9) =  .9399724845798D0 
            WI(1) = -.2245311506833D-3 
            WI(2) = -.5115351999107D-3 
            WI(3) = -.3282311058291D-3 
            WI(4) = -.8676605881308D-4 
            WI(5) = -.1055279410174D-4 
            WI(6) = -.5799648044579D-6 
            WI(7) = -.1256208909247D-7 
            WI(8) = -.7845838042681D-10 
            WI(9) = 0.D0 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) =  .6069460471460D-1 
            XI(2) =  .8315764753572D-1 
            XI(3) =  .1198802070163D0 
            XI(4) =  .1728707637813D0 
            XI(5) =  .2444098755278D0 
            XI(6) =  .3370229788789D0 
            XI(7) =  .4543153625599D0 
            XI(8) =  .6024428910596D0 
            XI(9) =  .7934113895484D0 
            XI(10) =  .1058111388730D1 
            WI(1) = -.1944360312097D-3 
            WI(2) = -.4812850685489D-3 
            WI(3) = -.3535558726335D-3 
            WI(4) = -.1134108921338D-3 
            WI(5) = -.1802255770270D-4 
            WI(6) = -.1443312793091D-5 
            WI(7) = -.5435739142044D-7 
            WI(8) = -.8187893807816D-9 
            WI(9) = -.3548382944316D-11 
            WI(10) = 0.D0 
            GO TO 10 
          END IF 
      END IF 
               
      WRITE (6,FMT=*) 'CASE N=',N,' IS NOT PROVIDED' 
      STOP 'GAUFDT2' 
               
 10   CONTINUE 
      RETURN   
               
      END 
c     ************************************************ 
      SUBROUTINE GAUFDT3(XI,WI,XSCALE,N,NCASE) 
c     ************************************************ 
C     .. Scalar Arguments .. 
      INTEGER N,NCASE 
      DOUBLE PRECISION XSCALE 
C     .. 
C     .. Array Arguments .. 
      DOUBLE PRECISION WI(*),XI(*) 
C     .. 
C     XZERO = LOG(0.5D0*(1.D0+SQRT(33.D0))) 
      XSCALE = 13.D0 * LOG(10.D0) 
      IF(NCASE.LT.0) RETURN 
      IF(NCASE.EQ.1) THEN 
          IF(N.EQ.1) THEN 
            XI(1) = -.4995533681254D0 
            WI(1) =  .1000761987062D1 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) = -.7883224548166D0 
            XI(2) = -.2102819290488D0 
            WI(1) =  .5008158550582D0 
            WI(2) =  .4999461320039D0 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) = -.8869539433051D0 
            XI(2) = -.4985246780636D0 
            XI(3) = -.1108245947517D0 
            WI(1) =  .2786217843736D0 
            WI(2) =  .4456497254143D0 
            WI(3) =  .2764904772742D0 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) = -.9302178825562D0 
            XI(2) = -.6683413826957D0 
            XI(3) = -.3267748861729D0 
            XI(4) = -.6632872582360D-1 
            WI(1) =  .1748032651591D0 
            WI(2) =  .3276748673990D0 
            WI(3) =  .3274369294384D0 
            WI(4) =  .1708469250656D0 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) = -.9527407996346D0 
            XI(2) = -.7675230060642D0 
            XI(3) = -.4963258887656D0 
            XI(4) = -.2252747907540D0 
            XI(5) = -.4216664075665D-1 
            WI(1) =  .1193444792156D0 
            WI(2) =  .2410800778804D0 
            WI(3) =  .2864792667402D0 
            WI(4) =  .2407706809735D0 
            WI(5) =  .1130874822523D0 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) = -.9658969580830D0 
            XI(2) = -.8289125044364D0 
            XI(3) = -.6155201820064D0 
            XI(4) = -.3745783773295D0 
            XI(5) = -.1613998250627D0 
            XI(6) = -.2700792447703D-1 
            WI(1) =  .8651896935133D-1 
            WI(2) =  .1821788582659D0 
            WI(3) =  .2362658406337D0 
            WI(4) =  .2361847747100D0 
            WI(5) =  .1819123947716D0 
            WI(6) =  .7770114932957D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) = -.9742352867528D0 
            XI(2) = -.8691483838048D0 
            XI(3) = -.6992113619995D0 
            XI(4) = -.4937759777189D0 
            XI(5) = -.2883977118401D0 
            XI(6) = -.1187513190379D0 
            XI(7) = -.1637168296200D-1 
            WI(1) =  .6555313912636D-1 
            WI(2) =  .1416009165601D0 
            WI(3) =  .1932911105270D0 
            WI(4) =  .2115492091186D0 
            WI(5) =  .1931705548580D0 
            WI(6) =  .1413970157967D0 
            WI(7) =  .5420004107533D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) = -.9798498206001D0 
            XI(2) = -.8968228486795D0 
            XI(3) = -.7592452215445D0 
            XI(4) = -.5856688968206D0 
            XI(5) = -.3995489056231D0 
            XI(6) = -.2260656108595D0 
            XI(7) = -.8892349762228D-1 
            XI(8) = -.8246922084314D-2 
            WI(1) =  .5136647571727D-1 
            WI(2) =  .1128413697565D0 
            WI(3) =  .1591766497299D0 
            WI(4) =  .1840130556861D0 
            WI(5) =  .1839733136119D0 
            WI(6) =  .1590368253873D0 
            WI(7) =  .1125479507908D0 
            WI(8) =  .3780634638229D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) = -.9838101424051D0 
            XI(2) = -.9166256465659D0 
            XI(3) = -.8034104659177D0 
            XI(4) = -.6564087364531D0 
            XI(5) = -.4915555505365D0 
            XI(6) = -.3267346293442D0 
            XI(7) = -.1798601708258D0 
            XI(8) = -.6732338663147D-1 
            XI(9) = -.1659902506918D-2 
            WI(1) =  .4132629265948D-1 
            WI(2) =  .9185474698102D-1 
            WI(3) =  .1325103135755D0 
            WI(4) =  .1588080108884D0 
            WI(5) =  .1678851961947D0 
            WI(6) =  .1587404871326D0 
            WI(7) =  .1323707165312D0 
            WI(8) =  .9112302999266D-1 
            WI(9) =  .2614319310647D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) = -.9867080097201D0 
            XI(2) = -.9312636851631D0 
            XI(3) = -.8366935400727D0 
            XI(4) = -.7113797737451D0 
            XI(5) = -.5664592035960D0 
            XI(6) = -.4148184465354D0 
            XI(7) = -.2699548299190D0 
            XI(8) = -.1448105316849D0 
            XI(9) = -.5124638050307D-1 
            XI(10) =  .3876985754771D-2 
            WI(1) =  .3396228971432D-1 
            WI(2) =  .7612965821898D-1 
            WI(3) =  .1115992600415D0 
            WI(4) =  .1371553532831D0 
            WI(5) =  .1505187831840D0 
            WI(6) =  .1504933086840D0 
            WI(7) =  .1370657174448D0 
            WI(8) =  .1114791562589D0 
            WI(9) =  .7447769699637D-1 
            WI(10) =  .1788076323624D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.11) THEN 
            XI(1) = -.9888922456674D0 
            XI(2) = -.9423794513789D0 
            XI(3) = -.8623247380363D0 
            XI(4) = -.7546476815945D0 
            XI(5) = -.6273349354398D0 
            XI(6) = -.4898341620377D0 
            XI(7) = -.3523552657755D0 
            XI(8) = -.2251217530272D0 
            XI(9) = -.1176716550598D0 
            XI(10) = -.3898097853369D-1 
            XI(11) =  .8611970343671D-2 
            WI(1) =  .2840212129474D-1 
            WI(2) =  .6407069796335D-1 
            WI(3) =  .9504326288205D-1 
            WI(4) =  .1189696838132D0 
            WI(5) =  .1340694298872D0 
            WI(6) =  .1392178412673D0 
            WI(7) =  .1340240894737D0 
            WI(8) =  .1188652290630D0 
            WI(9) =  .9494995206872D-1 
            WI(10) =  .6101392221102D-1 
            WI(11) =  .1213575713787D-1 
            GO TO 10 
          END IF 
          IF(N.EQ.12) THEN 
            XI(1) = -.9905793303603D0 
            XI(2) = -.9510136568551D0 
            XI(3) = -.8824442827205D0 
            XI(4) = -.7891644373516D0 
            XI(5) = -.6770352779590D0 
            XI(6) = -.5531055430578D0 
            XI(7) = -.4251690449875D0 
            XI(8) = -.3012785234489D0 
            XI(9) = -.1892505824663D0 
            XI(10) = -.9627688919320D-1 
            XI(11) = -.2938497440044D-1 
            XI(12) =  .1267997469169D-1 
            WI(1) =  .2410185782712D-1 
            WI(2) =  .5463493466581D-1 
            WI(3) =  .8178248984596D-1 
            WI(4) =  .1037940744194D0 
            WI(5) =  .1192820178175D0 
            WI(6) =  .1272703026722D0 
            WI(7) =  .1272522460132D0 
            WI(8) =  .1192202358657D0 
            WI(9) =  .1036860181614D0 
            WI(10) =  .8169312537476D-1 
            WI(11) =  .4981621437285D-1 
            WI(12) =  .8228470026286D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.13) THEN 
            XI(1) = -.9919094337205D0 
            XI(2) = -.9578506464976D0 
            XI(3) = -.8985054087579D0 
            XI(4) = -.8170598176291D0 
            XI(5) = -.7179043086802D0 
            XI(6) = -.6063864062590D0 
            XI(7) = -.4885222243421D0 
            XI(8) = -.3706737602115D0 
            XI(9) = -.2592100658733D0 
            XI(10) = -.1601793312505D0 
            XI(11) = -.7914539056046D-1 
            XI(12) = -.2167628537590D-1 
            XI(13) =  .1616921840130D-1 
            WI(1) =  .2070806582161D-1 
            WI(2) =  .4712104777684D-1 
            WI(3) =  .7103445551404D-1 
            WI(4) =  .9112094610770D-1 
            WI(5) =  .1062939776652D0 
            WI(6) =  .1157337194988D0 
            WI(7) =  .1189286197331D0 
            WI(8) =  .1157014147159D0 
            WI(9) =  .1062194587894D0 
            WI(10) =  .9102316746100D-1 
            WI(11) =  .7087861832629D-1 
            WI(12) =  .4038294761725D-1 
            WI(13) =  .5615548034994D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.14) THEN 
            XI(1) = -.9929765853682D0 
            XI(2) = -.9633549540853D0 
            XI(3) = -.9115184675497D0 
            XI(4) = -.8398793762977D0 
            XI(5) = -.7517869444306D0 
            XI(6) = -.6513615772170D0 
            XI(7) = -.5433017097774D0 
            XI(8) = -.4326648794127D0 
            XI(9) = -.3246331980937D0 
            XI(10) = -.2242757804908D0 
            XI(11) = -.1363355258018D0 
            XI(12) = -.6524071432834D-1 
            XI(13) = -.1531780250814D-1 
            XI(14) =  .1915227138769D-1 
            WI(1) =  .1798301075845D-1 
            WI(2) =  .4104498532019D-1 
            WI(3) =  .6222312482936D-1 
            WI(4) =  .8049416058229D-1 
            WI(5) =  .9500073408697D-1 
            WI(6) =  .1050632539944D0 
            WI(7) =  .1102096214263D0 
            WI(8) =  .1101963893298D0 
            WI(9) =  .1050191702400D0 
            WI(10) =  .9491855259004D-1 
            WI(11) =  .8041957796351D-1 
            WI(12) =  .6187348087075D-1 
            WI(13) =  .3243461739307D-1 
            WI(14) =  .3881307676994D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.15) THEN 
            XI(1) = -.9938457782642D0 
            XI(2) = -.9678506470704D0 
            XI(3) = -.9222013725832D0 
            XI(4) = -.8587566453353D0 
            XI(5) = -.7801133128415D0 
            XI(6) = -.6894918955640D0 
            XI(7) = -.5906042148820D0 
            XI(8) = -.4875017941585D0 
            XI(9) = -.3844109811466D0 
            XI(10) = -.2855620724009D0 
            XI(11) = -.1950219230146D0 
            XI(12) = -.1165645182362D0 
            XI(13) = -.5381890743374D-1 
            XI(14) = -.9943821567788D-2 
            XI(15) =  .2169651898369D-1 
            WI(1) =  .1576202189773D-1 
            WI(2) =  .3606472521727D-1 
            WI(3) =  .5492195351381D-1 
            WI(4) =  .7153285846324D-1 
            WI(5) =  .8521483879293D-1 
            WI(6) =  .9540681741943D-1 
            WI(7) =  .1016905250154D0 
            WI(8) =  .1038070636272D0 
            WI(9) =  .1016668163536D0 
            WI(10) =  .9535314775269D-1 
            WI(11) =  .8513176947055D-1 
            WI(12) =  .7148848608289D-1 
            WI(13) =  .5419847864068D-1 
            WI(14) =  .2579449296820D-1 
            WI(15) =  .2727991846477D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.16) THEN 
            XI(1) = -.9945631015081D0 
            XI(2) = -.9715692358410D0 
            XI(3) = -.9310744550952D0 
            XI(4) = -.8745332952801D0 
            XI(5) = -.8039886809703D0 
            XI(6) = -.7219908444891D0 
            XI(7) = -.6315046196188D0 
            XI(8) = -.5358024494972D0 
            XI(9) = -.4383466406198D0 
            XI(10) = -.3426653298852D0 
            XI(11) = -.2522273179161D0 
            XI(12) = -.1703245046996D0 
            XI(13) = -.1000082404954D0 
            XI(14) = -.4433176669004D-1 
            XI(15) = -.5308069629211D-2 
            XI(16) =  .2386524136408D-1 
            WI(1) =  .1392812610340D-1 
            WI(2) =  .3193347056698D-1 
            WI(3) =  .4881207021139D-1 
            WI(4) =  .6392847838595D-1 
            WI(5) =  .7673416277277D-1 
            WI(6) =  .8676557905619D-1 
            WI(7) =  .9365949043373D-1 
            WI(8) =  .9716572089012D-1 
            WI(9) =  .9715579849125D-1 
            WI(10) =  .9362711002680D-1 
            WI(11) =  .8670487158445D-1 
            WI(12) =  .7665833027354D-1 
            WI(13) =  .6390995206287D-1 
            WI(14) =  .4750132810316D-1 
            WI(15) =  .2032375390711D-1 
            WI(16) =  .1953744192376D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.17) THEN 
            XI(1) = -.9951619737956D0 
            XI(2) = -.9746795926219D0 
            XI(3) = -.9385216819701D0 
            XI(4) = -.8878423020107D0 
            XI(5) = -.8242697065723D0 
            XI(6) = -.7498475481206D0 
            XI(7) = -.6669686770646D0 
            XI(8) = -.5782982960423D0 
            XI(9) = -.4866885911817D0 
            XI(10) = -.3950876608305D0 
            XI(11) = -.3064459461404D0 
            XI(12) = -.2236240689427D0 
            XI(13) = -.1493120816663D0 
            XI(14) = -.8602068876563D-1 
            XI(15) = -.3636537468612D-1 
            XI(16) = -.1245851608908D-2 
            XI(17) =  .2571607745805D-1 
            WI(1) =  .1239639031248D-1 
            WI(2) =  .2846974964904D-1 
            WI(3) =  .4365246949293D-1 
            WI(4) =  .5743399681537D-1 
            WI(5) =  .6936948082606D-1 
            WI(6) =  .7907474136863D-1 
            WI(7) =  .8623737073738D-1 
            WI(8) =  .9062650185079D-1 
            WI(9) =  .9209999586224D-1 
            WI(10) =  .9060866246350D-1 
            WI(11) =  .8619777151910D-1 
            WI(12) =  .7901028566166D-1 
            WI(13) =  .6930937127817D-1 
            WI(14) =  .5742026466147D-1 
            WI(15) =  .4153812140107D-1 
            WI(16) =  .1589032542299D-1 
            WI(17) =  .1426487739228D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.18) THEN 
            XI(1) = -.9956671041283D0 
            XI(2) = -.9773071928330D0 
            XI(3) = -.9448310503637D0 
            XI(4) = -.8991656930087D0 
            XI(5) = -.8416242208629D0 
            XI(6) = -.7738622466544D0 
            XI(7) = -.6978297806686D0 
            XI(8) = -.6157151409295D0 
            XI(9) = -.5298821665387D0 
            XI(10) = -.4428025561182D0 
            XI(11) = -.3569854098982D0 
            XI(12) = -.2749063129255D0 
            XI(13) = -.1989393697985D0 
            XI(14) = -.1313051209642D0 
            XI(15) = -.7410895180024D-1 
            XI(16) = -.2960133588207D-1 
            XI(17) =  .2353336183123D-2 
            XI(18) =  .2729974999930D-1 
            WI(1) =  .1110394557917D-1 
            WI(2) =  .2553784297504D-1 
            WI(3) =  .3925892222823D-1 
            WI(4) =  .5185227275979D-1 
            WI(5) =  .6295393440167D-1 
            WI(6) =  .7224407399041D-1 
            WI(7) =  .7945508748142D-1 
            WI(8) =  .8437908627695D-1 
            WI(9) =  .8687373447062D-1 
            WI(10) =  .8686613488126D-1 
            WI(11) =  .8435466991809D-1 
            WI(12) =  .7940974907962D-1 
            WI(13) =  .7218005073605D-1 
            WI(14) =  .6291689551036D-1 
            WI(15) =  .5180128889647D-1 
            WI(16) =  .3615534093992D-1 
            WI(17) =  .1235766504196D-1 
            WI(18) =  .1061291895073D-2 
            GO TO 10 
          END IF 
          IF(N.EQ.19) THEN 
            XI(1) = -.9960970800047D0 
            XI(2) = -.9795468131263D0 
            XI(3) = -.9502217968080D0 
            XI(4) = -.9088749745069D0 
            XI(5) = -.8565766212974D0 
            XI(6) = -.7946813987734D0 
            XI(7) = -.7247928175086D0 
            XI(8) = -.6487216865496D0 
            XI(9) = -.5684393114784D0 
            XI(10) = -.4860266410463D0 
            XI(11) = -.4036207440845D0 
            XI(12) = -.3233601473604D0 
            XI(13) = -.2473308301527D0 
            XI(14) = -.1775164609448D0 
            XI(15) = -.1157699190471D0 
            XI(16) = -.6389171686342D-1 
            XI(17) = -.2379218917788D-1 
            XI(18) =  .5565172017514D-2 
            XI(19) =  .2865970486611D-1 
            WI(1) =  .1000345351336D-1 
            WI(2) =  .2303468274990D-1 
            WI(3) =  .3548901936582D-1 
            WI(4) =  .4702574846822D-1 
            WI(5) =  .5734458973567D-1 
            WI(6) =  .6617790883695D-1 
            WI(7) =  .7329666925006D-1 
            WI(8) =  .7851619092186D-1 
            WI(9) =  .8170083919595D-1 
            WI(10) =  .8276741384112D-1 
            WI(11) =  .8168711566715D-1 
            WI(12) =  .7848623261868D-1 
            WI(13) =  .7324735844250D-1 
            WI(14) =  .6611928660392D-1 
            WI(15) =  .5733498817501D-1 
            WI(16) =  .4687219857820D-1 
            WI(17) =  .3126896437873D-1 
            WI(18) =  .9585574522969D-2 
            WI(19) =  .8037521960421D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.20) THEN 
            XI(1) = -.9964661008992D0 
            XI(2) = -.9814711055748D0 
            XI(3) = -.9548630774035D0 
            XI(4) = -.9172597346551D0 
            XI(5) = -.8695419096347D0 
            XI(6) = -.8128281683983D0 
            XI(7) = -.7484481735637D0 
            XI(8) = -.6779114825157D0 
            XI(9) = -.6028722210703D0 
            XI(10) = -.5250904360996D0 
            XI(11) = -.4463910641311D0 
            XI(12) = -.3686215524308D0 
            XI(13) = -.2936092806376D0 
            XI(14) = -.2231203110272D0 
            XI(15) = -.1588238195721D0 
            XI(16) = -.1022840052412D0 
            XI(17) = -.5506978975286D-1 
            XI(18) = -.1874520625923D-1 
            XI(19) =  .8443709661071D-2 
            XI(20) =  .2983247601395D-1 
            WI(1) =  .9058722199130D-2 
            WI(2) =  .2088090502780D-1 
            WI(3) =  .3223149127518D-1 
            WI(4) =  .4282801285269D-1 
            WI(5) =  .5242076697774D-1 
            WI(6) =  .6078458362845D-1 
            WI(7) =  .6772322241627D-1 
            WI(8) =  .7307382569644D-1 
            WI(9) =  .7671066987793D-1 
            WI(10) =  .7854803656808D-1 
            WI(11) =  .7854210360156D-1 
            WI(12) =  .7669183255542D-1 
            WI(13) =  .7303928737690D-1 
            WI(14) =  .6767215852967D-1 
            WI(15) =  .6073677174762D-1 
            WI(16) =  .5243796625556D-1 
            WI(17) =  .4248408089805D-1 
            WI(18) =  .2684164155094D-1 
            WI(19) =  .7437073003140D-2 
            WI(20) =  .6188350235119D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.21) THEN 
            XI(1) = -.9967851649852D0 
            XI(2) = -.9831365075493D0 
            XI(3) = -.9588870206807D0 
            XI(4) = -.9245481166140D0 
            XI(5) = -.8808511408287D0 
            XI(6) = -.8287274762025D0 
            XI(7) = -.7692883051236D0 
            XI(8) = -.7038008833352D0 
            XI(9) = -.6336615627123D0 
            XI(10) = -.5603661084185D0 
            XI(11) = -.4854779571098D0 
            XI(12) = -.4105951329087D0 
            XI(13) = -.3373166042294D0 
            XI(14) = -.2672089778006D0 
            XI(15) = -.2017750372394D0 
            XI(16) = -.1424297283126D0 
            XI(17) = -.9051022670652D-1 
            XI(18) = -.4740502491774D-1 
            XI(19) = -.1431089564281D-1 
            XI(20) =  .1102992355745D-1 
            XI(21) =  .3084843863494D-1 
            WI(1) =  .8241700220721D-2 
            WI(2) =  .1901462386284D-1 
            WI(3) =  .2939850960882D-1 
            WI(4) =  .3915711904265D-1 
            WI(5) =  .4808126636962D-1 
            WI(6) =  .5598046949792D-1 
            WI(7) =  .6268621241258D-1 
            WI(8) =  .6805540959665D-1 
            WI(9) =  .7197340677520D-1 
            WI(10) =  .7435637631207D-1 
            WI(11) =  .7515303125471D-1 
            WI(12) =  .7434561180886D-1 
            WI(13) =  .7195022472656D-1 
            WI(14) =  .6801734584682D-1 
            WI(15) =  .6263613386203D-1 
            WI(16) =  .5594876372368D-1 
            WI(17) =  .4811758308393D-1 
            WI(18) =  .3851692142134D-1 
            WI(19) =  .2286126320249D-1 
            WI(20) =  .5786280317806D-2 
            WI(21) =  .4837341147850D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.22) THEN 
            XI(1) = -.9970628976474D0 
            XI(2) = -.9845873915028D0 
            XI(3) = -.9623979982516D0 
            XI(4) = -.9309216470102D0 
            XI(5) = -.8907704902381D0 
            XI(6) = -.8427260517410D0 
            XI(7) = -.7877236576255D0 
            XI(8) = -.7268341847192D0 
            XI(9) = -.6612432344523D0 
            XI(10) = -.5922281065856D0 
            XI(11) = -.5211330260545D0 
            XI(12) = -.4493431271160D0 
            XI(13) = -.3782577438868D0 
            XI(14) = -.3092636089128D0 
            XI(15) = -.2437087146701D0 
            XI(16) = -.1828785511301D0 
            XI(17) = -.1279819125816D0 
            XI(18) = -.8017727786813D-1 
            XI(19) = -.4070526355618D-1 
            XI(20) = -.1037412938785D-1 
            XI(21) =  .1335689708536D-1 
            XI(22) =  .3173269299216D-1 
            WI(1) =  .7530375265603D-2 
            WI(2) =  .1738702512348D-1 
            WI(3) =  .2692008166844D-1 
            WI(4) =  .3593037068353D-1 
            WI(5) =  .4424143110682D-1 
            WI(6) =  .5169125346169D-1 
            WI(7) =  .5813470727063D-1 
            WI(8) =  .6344625344919D-1 
            WI(9) =  .6752234972692D-1 
            WI(10) =  .7028343381036D-1 
            WI(11) =  .7167542664921D-1 
            WI(12) =  .7167071496519D-1 
            WI(13) =  .7026861000894D-1 
            WI(14) =  .6749549807665D-1 
            WI(15) =  .6340593211832D-1 
            WI(16) =  .5808883320740D-1 
            WI(17) =  .5168012421561D-1 
            WI(18) =  .4428025502278D-1 
            WI(19) =  .3487795825750D-1 
            WI(20) =  .1932393041695D-1 
            WI(21) =  .4524029480596D-2 
            WI(22) =  .3833930763003D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.23) THEN 
            XI(1) = -.9973061402275D0 
            XI(2) = -.9858590283297D0 
            XI(3) = -.9654793090894D0 
            XI(4) = -.9365261022364D0 
            XI(5) = -.8995156085099D0 
            XI(6) = -.8551082865294D0 
            XI(7) = -.8040967362825D0 
            XI(8) = -.7473915078203D0 
            XI(9) = -.6860048604970D0 
            XI(10) = -.6210327318830D0 
            XI(11) = -.5536352384748D0 
            XI(12) = -.4850160682911D0 
            XI(13) = -.4164011573913D0 
            XI(14) = -.3490170745811D0 
            XI(15) = -.2840695915616D0 
            XI(16) = -.2227231518557D0 
            XI(17) = -.1660833561664D0 
            XI(18) = -.1151914327516D0 
            XI(19) = -.7106498760489D-1 
            XI(20) = -.3481361009995D-1 
            XI(21) = -.6846873679018D-2 
            XI(22) =  .1545275772982D-1 
            XI(23) =  .3250592934764D-1 
            WI(1) =  .6907279176068D-2 
            WI(2) =  .1595919907190D-1 
            WI(3) =  .2473992517590D-1 
            WI(4) =  .3308027767654D-1 
            WI(5) =  .4083043310284D-1 
            WI(6) =  .4785187401184D-1 
            WI(7) =  .5401919681278D-1 
            WI(8) =  .5922224892950D-1 
            WI(9) =  .6336806361893D-1 
            WI(10) =  .6638249654894D-1 
            WI(11) =  .6821151981100D-1 
            WI(12) =  .6882214132519D-1 
            WI(13) =  .6820293093990D-1 
            WI(14) =  .6636420075220D-1 
            WI(15) =  .6333823043415D-1 
            WI(16) =  .5918122929146D-1 
            WI(17) =  .5398108578956D-1 
            WI(18) =  .4786414616860D-1 
            WI(19) =  .4084474900605D-1 
            WI(20) =  .3150023082023D-1 
            WI(21) =  .1622267422010D-1 
            WI(22) =  .3560138783270D-2 
            WI(23) =  .3077155951781D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.24) THEN 
            XI(1) = -.9975203751421D0 
            XI(2) = -.9869797361798D0 
            XI(3) = -.9681980664831D0 
            XI(4) = -.9414795617082D0 
            XI(5) = -.9072624967742D0 
            XI(6) = -.8661087049282D0 
            XI(7) = -.8186940463209D0 
            XI(8) = -.7657972636522D0 
            XI(9) = -.7082871974312D0 
            XI(10) = -.6471085411825D0 
            XI(11) = -.5832663683531D0 
            XI(12) = -.5178096916087D0 
            XI(13) = -.4518143386191D0 
            XI(14) = -.3863654505847D0 
            XI(15) = -.3225399362506D0 
            XI(16) = -.2613892825368D0 
            XI(17) = -.2039234869192D0 
            XI(18) = -.1510987922779D0 
            XI(19) = -.1038200030779D0 
            XI(20) = -.6299312959427D-1 
            XI(21) = -.2960082342791D-1 
            XI(22) = -.3662099454240D-2 
            XI(23) =  .1734222645710D-1 
            XI(24) =  .3318520360174D-1 
            WI(1) =  .6358407435823D-2 
            WI(2) =  .1469983053915D-1 
            WI(3) =  .2281240053445D-1 
            WI(4) =  .3055142135617D-1 
            WI(5) =  .3778891245077D-1 
            WI(6) =  .4440584196419D-1 
            WI(7) =  .5029347126377D-1 
            WI(8) =  .5535504935231D-1 
            WI(9) =  .5950737566703D-1 
            WI(10) =  .6268214993547D-1 
            WI(11) =  .6482707438937D-1 
            WI(12) =  .6590668375261D-1 
            WI(13) =  .6590288458704D-1 
            WI(14) =  .6481520669811D-1 
            WI(15) =  .6266088091640D-1 
            WI(16) =  .5947534380340D-1 
            WI(17) =  .5531521451131D-1 
            WI(18) =  .5026678792674D-1 
            WI(19) =  .4444180941173D-1 
            WI(20) =  .3774050527736D-1 
            WI(21) =  .2834034016257D-1 
            WI(22) =  .1354170374828D-1 
            WI(23) =  .2822870053565D-2 
            WI(24) =  .2498213245034D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.25) THEN 
            XI(1) = -.9977100373536D0 
            XI(2) = -.9879724605948D0 
            XI(3) = -.9706088171686D0 
            XI(4) = -.9458784410673D0 
            XI(5) = -.9141558036027D0 
            XI(6) = -.8759217665954D0 
            XI(7) = -.8317560090128D0 
            XI(8) = -.7823281942520D0 
            XI(9) = -.7283878172238D0 
            XI(10) = -.6707528569988D0 
            XI(11) = -.6102974035856D0 
            XI(12) = -.5479384498963D0 
            XI(13) = -.4846220574983D0 
            XI(14) = -.4213091207369D0 
            XI(15) = -.3589609705027D0 
            XI(16) = -.2985250855155D0 
            XI(17) = -.2409212799096D0 
            XI(18) = -.1870292694828D0 
            XI(19) = -.1376809693468D0 
            XI(20) = -.9367001253362D-1 
            XI(21) = -.5581286614835D-1 
            XI(22) = -.2495989225599D-1 
            XI(23) = -.7687174405023D-3 
            XI(24) =  .1904739329867D-1 
            XI(25) =  .3378460022423D-1 
            WI(1) =  .5872426843688D-2 
            WI(2) =  .1358349209854D-1 
            WI(3) =  .2110020616374D-1 
            WI(4) =  .2829801937269D-1 
            WI(5) =  .3506698731052D-1 
            WI(6) =  .4130432624987D-1 
            WI(7) =  .4691540880859D-1 
            WI(8) =  .5181511536616D-1 
            WI(9) =  .5592910217125D-1 
            WI(10) =  .5919491599600D-1 
            WI(11) =  .6156292767667D-1 
            WI(12) =  .6299706561933D-1 
            WI(13) =  .6347533412837D-1 
            WI(14) =  .6299010889188D-1 
            WI(15) =  .6154823872376D-1 
            WI(16) =  .5917115467285D-1 
            WI(17) =  .5589581043642D-1 
            WI(18) =  .5177865604179D-1 
            WI(19) =  .4690361518390D-1 
            WI(20) =  .4136084345832D-1 
            WI(21) =  .3490661005830D-1 
            WI(22) =  .2537492111058D-1 
            WI(23) =  .1125505052707D-1 
            WI(24) =  .2256704413381D-2 
            WI(25) =  .2049457384408D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.26) THEN 
            XI(1) = -.9978787456321D0 
            XI(2) = -.9888559514846D0 
            XI(3) = -.9727562536373D0 
            XI(4) = -.9498020585951D0 
            XI(5) = -.9203152015375D0 
            XI(6) = -.8847095845132D0 
            XI(7) = -.8434851022034D0 
            XI(8) = -.7972205809148D0 
            XI(9) = -.7465656492906D0 
            XI(10) = -.6922316282846D0 
            XI(11) = -.6349815641969D0 
            XI(12) = -.5756195463635D0 
            XI(13) = -.5149794644814D0 
            XI(14) = -.4539133724176D0 
            XI(15) = -.3932796370067D0 
            XI(16) = -.3339310644402D0 
            XI(17) = -.2767032299031D0 
            XI(18) = -.2224033881725D0 
            XI(19) = -.1718010780566D0 
            XI(20) = -.1256245028078D0 
            XI(21) = -.8457665361068D-1 
            XI(22) = -.4940018604596D-1 
            XI(23) = -.2080208514240D-1 
            XI(24) =  .1872533101082D-2 
            XI(25) =  .2058810506331D-1 
            XI(26) =  .3431578327418D-1 
            WI(1) =  .5440086790047D-2 
            WI(2) =  .1258936734749D-1 
            WI(3) =  .1957263009058D-1 
            WI(4) =  .2628202685296D-1 
            WI(5) =  .3262260389452D-1 
            WI(6) =  .3850519667546D-1 
            WI(7) =  .4384716117258D-1 
            WI(8) =  .4857345783151D-1 
            WI(9) =  .5261768547154D-1 
            WI(10) =  .5592300375739D-1 
            WI(11) =  .5844292180110D-1 
            WI(12) =  .6014193817620D-1 
            WI(13) =  .6099602042241D-1 
            WI(14) =  .6099291528211D-1 
            WI(15) =  .6013229350125D-1 
            WI(16) =  .5842579326759D-1 
            WI(17) =  .5589726848049D-1 
            WI(18) =  .5258426859446D-1 
            WI(19) =  .4854280785737D-1 
            WI(20) =  .4385315019153D-1 
            WI(21) =  .3857479795187D-1 
            WI(22) =  .3229131740839D-1 
            WI(23) =  .2259592141781D-1 
            WI(24) =  .9328143821698D-2 
            WI(25) =  .1819471247603D-2 
            WI(26) =  .1697377561173D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.27) THEN 
            XI(1) = -.9980294763731D0 
            XI(2) = -.9896456502301D0 
            XI(3) = -.9746772694043D0 
            XI(4) = -.9533161251235D0 
            XI(5) = -.9258403241032D0 
            XI(6) = -.8926080212559D0 
            XI(7) = -.8540525047533D0 
            XI(8) = -.8106765058582D0 
            XI(9) = -.7630456397652D0 
            XI(10) = -.7117810379579D0 
            XI(11) = -.6575512637797D0 
            XI(12) = -.6010636172202D0 
            XI(13) = -.5430549453049D0 
            XI(14) = -.4842820834762D0 
            XI(15) = -.4255120618787D0 
            XI(16) = -.3675122195605D0 
            XI(17) = -.3110403837243D0 
            XI(18) = -.2568353180326D0 
            XI(19) = -.2056078652997D0 
            XI(20) = -.1580341647736D0 
            XI(21) = -.1147557663797D0 
            XI(22) = -.7640168939581D-1 
            XI(23) = -.4365087112060D-1 
            XI(24) = -.1705395284911D-1 
            XI(25) =  .4292311518709D-2 
            XI(26) =  .2198218280711D-1 
            XI(27) =  .3478844699902D-1 
            WI(1) =  .5053776482567D-2 
            WI(2) =  .1170028507766D-1 
            WI(3) =  .1820421110362D-1 
            WI(4) =  .2447164861933D-1 
            WI(5) =  .3042018049054D-1 
            WI(6) =  .3597211570191D-1 
            WI(7) =  .4105501927395D-1 
            WI(8) =  .4560258630100D-1 
            WI(9) =  .4955548897793D-1 
            WI(10) =  .5286214221941D-1 
            WI(11) =  .5547936936981D-1 
            WI(12) =  .5737295641234D-1 
            WI(13) =  .5851808535741D-1 
            WI(14) =  .5889963917698D-1 
            WI(15) =  .5851237528707D-1 
            WI(16) =  .5736098596292D-1 
            WI(17) =  .5546015503866D-1 
            WI(18) =  .5283503267129D-1 
            WI(19) =  .4952329446361D-1 
            WI(20) =  .4558031366082D-1 
            WI(21) =  .4108070639850D-1 
            WI(22) =  .3604222903069D-1 
            WI(23) =  .2985193298902D-1 
            WI(24) =  .2000525889454D-1 
            WI(25) =  .7720885190616D-2 
            WI(26) =  .1479507802426D-2 
            WI(27) =  .1418051074982D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.28) THEN 
            XI(1) = -.9981646957167D0 
            XI(2) = -.9903543656416D0 
            XI(3) = -.9764025324682D0 
            XI(4) = -.9564754357674D0 
            XI(5) = -.9308146163646D0 
            XI(6) = -.8997315240376D0 
            XI(7) = -.8636035075075D0 
            XI(8) = -.8228691943894D0 
            XI(9) = -.7780231606877D0 
            XI(10) = -.7296099310444D0 
            XI(11) = -.6782173781688D0 
            XI(12) = -.6244696015413D0 
            XI(13) = -.5690193736608D0 
            XI(14) = -.5125402490416D0 
            XI(15) = -.4557184375832D0 
            XI(16) = -.3992445503809D0 
            XI(17) = -.3438053338965D0 
            XI(18) = -.2900755254492D0 
            XI(19) = -.2387100318122D0 
            XI(20) = -.1903369386122D0 
            XI(21) = -.1455530359225D0 
            XI(22) = -.1049273946759D0 
            XI(23) = -.6902849536803D-1 
            XI(24) = -.3847664441132D-1 
            XI(25) = -.1365492454756D-1 
            XI(26) =  .6515115199890D-2 
            XI(27) =  .2324557805003D-1 
            XI(28) =  .3521068094907D-1 
            WI(1) =  .4707188199662D-2 
            WI(2) =  .1090198094342D-1 
            WI(3) =  .1697370439949D-1 
            WI(4) =  .2284016724569D-1 
            WI(5) =  .2842949782545D-1 
            WI(6) =  .3367371843268D-1 
            WI(7) =  .3850912106474D-1 
            WI(8) =  .4287697567502D-1 
            WI(9) =  .4672422742078D-1 
            WI(10) =  .5000413418462D-1 
            WI(11) =  .5267682887581D-1 
            WI(12) =  .5470979726607D-1 
            WI(13) =  .5607826404133D-1 
            WI(14) =  .5676548085924D-1 
            WI(15) =  .5676291212803D-1 
            WI(16) =  .5607032205661D-1 
            WI(17) =  .5469580057415D-1 
            WI(18) =  .5265589253944D-1 
            WI(19) =  .4997636459951D-1 
            WI(20) =  .4669480789283D-1 
            WI(21) =  .4286563529586D-1 
            WI(22) =  .3855500784858D-1 
            WI(23) =  .3372605859845D-1 
            WI(24) =  .2755481943422D-1 
            WI(25) =  .1760962695684D-1 
            WI(26) =  .6391052566346D-2 
            WI(27) =  .1213184375178D-2 
            WI(28) =  .1194157620547D-3 
            GO TO 10 
          END IF 
      END IF 
      IF(NCASE.EQ.0) THEN 
          IF(N.EQ.1) THEN 
            XI(1) =  .8658765285328D-1 
            WI(1) = -.7619870621007D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.2) THEN 
            XI(1) =  .7018635586069D-1 
            XI(2) =  .1648420028631D0 
            WI(1) = -.6299550438768D-3 
            WI(2) = -.1320320182239D-3 
            GO TO 10 
          END IF 
          IF(N.EQ.3) THEN 
            XI(1) =  .6287547487759D-1 
            XI(2) =  .1263743732614D0 
            XI(3) =  .2593759997905D0 
            WI(1) = -.4982067606548D-3 
            WI(2) = -.2538660296513D-3 
            WI(3) = -.9914271794602D-5 
            GO TO 10 
          END IF 
          IF(N.EQ.4) THEN 
            XI(1) =  .5862019673633D-1 
            XI(2) =  .1072723348447D0 
            XI(3) =  .2000815208147D0 
            XI(4) =  .3621773739436D0 
            WI(1) = -.3988321542456D-3 
            WI(2) = -.3248325101772D-3 
            WI(3) = -.3779242448766D-4 
            WI(4) = -.5299731901904D-6 
            GO TO 10 
          END IF 
          IF(N.EQ.5) THEN 
            XI(1) =  .5580149442192D-1 
            XI(2) =  .9554929637825D-1 
            XI(3) =  .1681377964080D0 
            XI(4) =  .2844692409692D0 
            XI(5) =  .4699534159045D0 
            WI(1) = -.3253981806946D-3 
            WI(2) = -.3583602728554D-3 
            WI(3) = -.7459792098562D-4 
            WI(4) = -.3607259481678D-5 
            WI(5) = -.2342808337208D-7 
            GO TO 10 
          END IF 
          IF(N.EQ.6) THEN 
            XI(1) =  .5378327422155D-1 
            XI(2) =  .8753419090927D-1 
            XI(3) =  .1475486013330D0 
            XI(4) =  .2402487243581D0 
            XI(5) =  .3758503900280D0 
            XI(6) =  .5811026397676D0 
            WI(1) = -.2703300036972D-3 
            WI(2) = -.3688980274340D-3 
            WI(3) = -.1119817564992D-3 
            WI(4) = -.1051118301031D-4 
            WI(5) = -.2651769841178D-6 
            WI(6) = -.9144758837307D-9 
            GO TO 10 
          END IF 
          IF(N.EQ.7) THEN 
            XI(1) =  .5226050178454D-1 
            XI(2) =  .8167201793277D-1 
            XI(3) =  .1330077430515D0 
            XI(4) =  .2106754928942D0 
            XI(5) =  .3201864352788D0 
            XI(6) =  .4722018206809D0 
            XI(7) =  .6947158658171D0 
            WI(1) = -.2281728860198D-3 
            WI(2) = -.3660876714974D-3 
            WI(3) = -.1455661669608D-3 
            WI(4) = -.2104116674274D-4 
            WI(5) = -.1102812564599D-5 
            WI(6) = -.1632565898516D-7 
            WI(7) = -.3265635912326D-10 
            GO TO 10 
          END IF 
          IF(N.EQ.8) THEN 
            XI(1) =  .5106728052492D-1 
            XI(2) =  .7718021960064D-1 
            XI(3) =  .1221289767460D0 
            XI(4) =  .1892182750840D0 
            XI(5) =  .2819262314361D0 
            XI(6) =  .4058914349362D0 
            XI(7) =  .5722993373830D0 
            XI(8) =  .8102207662384D0 
            WI(1) = -.1952400591912D-3 
            WI(2) = -.3558179868265D-3 
            WI(3) = -.1737167800734D-3 
            WI(4) = -.3422270858288D-4 
            WI(5) = -.2894387465101D-5 
            WI(6) = -.9425663627888D-7 
            WI(7) = -.8822349721881D-9 
            WI(8) = -.1090309811291D-11 
            GO TO 10 
          END IF 
          IF(N.EQ.9) THEN 
            XI(1) =  .5010512229952D-1 
            XI(2) =  .7361870488332D-1 
            XI(3) =  .1136546981758D0 
            XI(4) =  .1728303171172D0 
            XI(5) =  .2535784114466D0 
            XI(6) =  .3593428369478D0 
            XI(7) =  .4960513470331D0 
            XI(8) =  .6753388967268D0 
            XI(9) =  .9272309938135D0 
            WI(1) = -.1690371806203D-3 
            WI(2) = -.3416084306686D-3 
            WI(3) = -.1962399978254D-3 
            WI(4) = -.4893420518833D-4 
            WI(5) = -.5840847627894D-5 
            WI(6) = -.3194491217583D-6 
            WI(7) = -.6907905644352D-8 
            WI(8) = -.4310823258919D-10 
            WI(9) = 0.D0 
            GO TO 10 
          END IF 
          IF(N.EQ.10) THEN 
            XI(1) =  .4931163905232D-1 
            XI(2) =  .7071970245523D-1 
            XI(3) =  .1068516437886D0 
            XI(4) =  .1598557585797D0 
            XI(5) =  .2315636522361D0 
            XI(6) =  .3242793858035D0 
            XI(7) =  .4416207272297D0 
            XI(8) =  .5897748390558D0 
            XI(9) =  .7807601609886D0 
            XI(10) =  .1045472162940D1 
            WI(1) = -.1478477544962D-3 
            WI(2) = -.3255507058100D-3 
            WI(3) = -.2135905657586D-3 
            WI(4) = -.6420612938099D-4 
            WI(5) = -.9968393584411D-5 
            WI(6) = -.7932473564127D-6 
            WI(7) = -.2981508455441D-7 
            WI(8) = -.4486851140015D-9 
            WI(9) = -.1943368334675D-11 
            WI(10) = 0.D0 
            GO TO 10 
          END IF 
      END IF 
               
      WRITE (6,FMT=*) 'CASE N=',N,' IS NOT PROVIDED' 
               
 10   CONTINUE 
      RETURN   
               
      END 
