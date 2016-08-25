c 20.09.95 ***************************************************************
      SUBROUTINE EMESHT(EZ,DF,NPNT,EBOT,EMU,TK,NPOL,NPNT1,NPNT2,NPNT3)
c ************************************************************************
C This subroutine provides the energy mesh in array EZ and the
C appropriate integration weights in array DF.
C
C Poles of the Fermi function C (Matsubara frequencies) and
C a contour in the complex energy are used as described in (????).
C
C The contour consists of three straight lines with
C NPNT1, NPNT2, and NPNT3 integration points and is determined by
C the input arguments: EBOT, EMU, TK, and NPOL.
C
C            TK   = temperature in K
C            EMU  = chemical potential in Ry
C            EBOT = bottom of contour in Ry
C            NPOL = number of Matsubara frequencies
C
C The three lines are defined by:
C
C  1. the line from EBOT to EBOT+2*NPOL*pi*i*k*TK
C              with NPNT1 integration points (Gauss-Legendre rule)
C
C  2. the line from EBOT+2*NPOL*pi*i*k*TK to EMU+(2*NPOL*pi*i-30)*k*TK
C              with NPNT2 integration points (Gauss-Legendre rule)
C
C  3. the line from EMU+(2*NPOL*pi*i-30)*k*TK to infinity
C              with NPNT3 integration points (Gauss-Fermi-Dirac rule)
C
C  The total number of integration points is given by:
C              NPNT=NPNT1+NPNT2+NPNT3+NPOL
C
C  The integration points and weights on three lines are chosen
C  according to Gauss integration rules. Only in third interval
C  the Fermi function matters since exp(x) < 10**(-10) for x < -25.
C
C  There are two special cases determined by NPOL = 0 and NPOL < 0.
C
C  a) NPOL = 0 leads to density-of-states calculations
C  with constant integration weights and equally distributed points
C  between EBOT - pi*i*k*TK and EMU - pi*i*k*TK.
C
C  The total number of integration points is given by:
C              NPNT=NPNT2
C
C
C  b) NPOL < 0 is meant for calculations where the Fermi-Dirac function
C  is replaced by a step function with step at EMU. When this option is
C  used no poles of the Fermi-Dirac function are used and the contour
C  consists of the three straight lines:
C
C  1. the line from EBOT to EBOT-2*NPOL*pi*i*k*TK
C              with NPNT1 integration points (Gauss-Legendre rule)
C
C  2. the line from EBOT-2*NPOL*pi*i*k*TK to EMU-2*NPOL*pi*i*k*TK
C              with NPNT2 integration points (Gauss-Legendre rule)
C
C  3. the line from EMU-2*NPOL*pi*i*k*TK to EMU
C              with NPNT3 integration points (Gauss-Legendre rule)
C
C
C  The total number of integration points is given by:
C              NPNT=NPNT1+NPNT2+NPNT3
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      DOUBLE PRECISION EBOT,EMU,TK
      INTEGER NPNT,NPNT1,NPNT2,NPNT3,NPOL
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
      DATA PI /3.14159265358979312d0/
      DATA KB /0.6333659D-5/
      DATA RYD/13.6058D0/
C     ..
      write(6,*) '>>> EMESHT: generates a complex E contour'
      ETK = PI*KB*TK
      IF (NPOL.EQ.0) THEN
        DE = (EMU-EBOT)
        IF (NPNT2.GT.1) DE = DE/(NPNT2-1)
        NPNT = 0
        DO 10 I = 1,NPNT2
          NPNT = NPNT + 1
          ER = EBOT + (I-1)*DE
          EZ(NPNT) = DCMPLX(ER,ETK)
          DF(NPNT) = DE
   10   CONTINUE
        WRITE (6,FMT=9000) NPNT,ETK,ETK*RYD
c 9000   format('density-of-states calculation',/,
c     +       'for',I4,' energy points with broadening',E12.4,'Ry'
 9000   format(' density-of-states calculation',/,
     +         ' for',I4,' energy points with broadening',
     +           3p,f10.3,' mRy = ',f10.3,' meV')

      ELSE IF (NPOL.GT.0) THEN
        CALL GAULEG(XI,WI,NPNT1)
        DE = NPOL*DCMPLX(0.0D0,ETK)
        NPNT = 0
        DO 20 I = 1,NPNT1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT
          DF(NPNT) = WI(I)*DE
   20   CONTINUE
        CALL GAULEG(XI,WI,NPNT2)
        DE = (EMU-30*KB*TK-EBOT)*0.5D0
        DO 30 I = 1,NPNT2
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT + 2*NPOL*DCMPLX(0.0D0,ETK)
          DF(NPNT) = WI(I)*DE
   30   CONTINUE
        CALL GAUFD(XI,WI,NPNT3)
        DE = 30*KB*TK
        DO 40 I = 1,NPNT3
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + EMU + 2*NPOL*DCMPLX(0.0D0,ETK)
          DF(NPNT) = WI(I)*DE
   40   CONTINUE
        DO 50 I = NPOL,1,-1
          NPNT = NPNT + 1
          EZ(NPNT) = EMU + (2*I-1)*DCMPLX(0.0D0,ETK)
          DF(NPNT) = -2*DCMPLX(0.0D0,ETK)
   50   CONTINUE

      ELSE
        IF (NPNT1.GT.0) CALL GAULEG(XI,WI,NPNT1)
        DE = -NPOL*DCMPLX(0.0D0,ETK)
        NPNT = 0
        DO 60 I = 1,NPNT1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT
          DF(NPNT) = WI(I)*DE
   60   CONTINUE
        CALL GAULEG(XI,WI,NPNT2)
        DE = (EMU-EBOT)*0.5D0
        DO 70 I = 1,NPNT2
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EBOT - 2*NPOL*DCMPLX(0.0D0,ETK)
          IF (OPT('GF-EF   ')) EZ(NPNT) = EMU + NPOL*DCMPLX(0.0D0,ETK)
         DF(NPNT) = WI(I)*DE
   70   CONTINUE
        IF (NPNT3.GT.0) CALL GAULEG(XI,WI,NPNT3)
        DE = -NPOL*DCMPLX(0.0D0,ETK)
        DO 80 I = NPNT3,1,-1
          NPNT = NPNT + 1
          EZ(NPNT) = XI(I)*DE + DE + EMU
          DF(NPNT) = -WI(I)*DE
   80   CONTINUE
      END IF

      RETURN

      END
c***********************************************************************
      SUBROUTINE EMESHT_SEMICORE(EZ,DF,IEMXD,EFERMI,NPNT,IESEMICORE,
     &   EBOTSEMI,EMUSEMI,TKSEMI,NPOLSEMI,NPNT1SEMI,NPNT2SEMI,NPNT3SEMI,
     &   EBOTVAL, EMUVAL, TKVAL, NPOLVAL, NPNT1VAL, NPNT2VAL, NPNT3VAL)
c Calls the routine EMESHT once for the semicore contour and once for
c the valence contour. In the semicore contour, -NPOLSEMI is used
c to create rectangular contour.
      implicit none
      INTEGER IEMXD
      COMPLEX*16 EZ(*),DF(*),EZSEMI(IEMXD),DFSEMI(IEMXD)
      COMPLEX*16 EZVAL(IEMXD),DFVAL(IEMXD)
      REAL*8 EBOTSEMI,EMUSEMI,TKSEMI,EBOTVAL,EMUVAL,TKVAL,EFERMI
      INTEGER NPOLSEMI,NPNT1SEMI,NPNT2SEMI,NPNT3SEMI
      INTEGER NPOLVAL, NPNT1VAL, NPNT2VAL, NPNT3VAL
      INTEGER IESEMICORE,NPNT,NPNTSEMI,NPNTVAL
      INTEGER IE,JE

      IESEMICORE = 0 
      IF (NPNT1SEMI+NPNT2SEMI+NPNT3SEMI.GT.0) THEN
         write(*,*) 'entering emesht for semicore' ! test fivos
      CALL EMESHT(EZSEMI,DFSEMI,NPNTSEMI,EBOTSEMI,EMUSEMI,TKSEMI
     &            ,-NPOLSEMI,NPNT1SEMI,NPNT2SEMI,NPNT3SEMI)
      IESEMICORE = NPNTSEMI
      ENDIF

         write(*,*) 'entering emesht for valence' ! test fivos
      CALL EMESHT(EZVAL,DFVAL,NPNTVAL,EBOTVAL,EMUVAL,TKVAL,
     &            NPOLVAL,NPNT1VAL,NPNT2VAL,NPNT3VAL)

      NPNT = IESEMICORE + NPNTVAL

      DO IE = 1,IESEMICORE
         EZ(IE) = EZSEMI(IE)
         DF(IE) = DFSEMI(IE)
c        write(*,*) IE, EZ(IE)! test fivos
      ENDDO

      DO IE = IESEMICORE+1,NPNT
         JE = IE - IESEMICORE
         EZ(IE) = EZVAL(JE)
         DF(IE) = DFVAL(JE)
c        write(*,*) IE,JE,EZ(IE)! test fivos
      ENDDO

      RETURN
      END

