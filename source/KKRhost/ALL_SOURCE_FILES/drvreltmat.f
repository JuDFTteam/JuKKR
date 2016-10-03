C*==drvreltmat.f    processed by SPAG 6.05Rc at 11:48 on 10 May 2004
      SUBROUTINE DRVRELTMAT(ERYD,TMATLL,VT,BT,R,DRDI,R2DRDI,ZAT,JWS,
     &                      SOLVER,SOCTL,CTL,LMMAXD,LMAXD,IRMD)
C   ********************************************************************
C   *                                                                  *
C   * driving routine to call relativistic < SSITE > routine           *
C   * only to calculate the single-site t matrix                       *
C   * v.popescu, munich, may 2004                                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C PARAMETER definitions
C
      INTEGER NRMAX
      PARAMETER ( NRMAX=900 )
      INTEGER NLAMAX,NQMAX,NTMAX,NMMAX
      PARAMETER (NLAMAX=1,NQMAX=1,NTMAX=1,NMMAX=1)
      INTEGER NLMAX,NKMMAX,NMUEMAX,NKMPMAX,NKMAX,LINMAX
      PARAMETER ( NLMAX = 5 ) ! this should be >= LMAXD + 1
      PARAMETER ( NKMMAX = 2*NLMAX**2, NKMAX = 2*NLMAX-1 )
      PARAMETER ( NKMPMAX = NKMMAX+2*NLMAX, NMUEMAX = 2*NLMAX)
      PARAMETER ( LINMAX = 2*NLMAX*(2*NLMAX-1) )
C
C Dummy arguments
C
      INTEGER LMAXD,LMMAXD,IRMD
      INTEGER ZAT(NTMAX),JWS(NMMAX)
      DOUBLE COMPLEX TMATLL(LMMAXD,LMMAXD)
!       DOUBLE PRECISION SOCTL(NTMAX,NLMAX)
!       DOUBLE PRECISION CTL(NTMAX,NLMAX)
      DOUBLE PRECISION SOCTL(NLMAX)
      DOUBLE PRECISION CTL(NLMAX)
!       DOUBLE PRECISION VT(NRMAX,NTMAX),BT(NRMAX,NTMAX)
      DOUBLE PRECISION VT(NRMAX),BT(NRMAX)
!             DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL))
!       DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL))
!       DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL)),
!      &                 R2DRDIREL(IRMD*KREL+(1-KREL)),
!      &                 RMREL(IRMD*KREL+(1-KREL))
!       
      DOUBLE PRECISION R(NRMAX,NMMAX),R2DRDI(NRMAX,NMMAX)
      DOUBLE PRECISION DRDI(NRMAX,NMMAX)
C
C Local variables
C
      REAL*8 AMEOPO(NKMMAX,NKMMAX,NLAMAX,3),AT(NRMAX,NLAMAX,3,NTMAX)
      COMPLEX*16 BZJ(LINMAX,NTMAX),BZZ(LINMAX,NTMAX),
     &           DZJ(LINMAX,NTMAX),
     &           DZZ(LINMAX,NTMAX),ERYD,
     &           MSST(NKMMAX,NKMMAX,NTMAX),
     &           OZJ(LINMAX,NTMAX),
     &           OZZ(LINMAX,NTMAX),P,
     &           QZJ(LINMAX,NTMAX),QZZ(LINMAX,NTMAX),
     &           SZJ(LINMAX,NTMAX),SZZ(LINMAX,NTMAX)
      COMPLEX*16 TSST(NKMMAX,NKMMAX,NTMAX),TSSTLIN(LINMAX,NTMAX),
     &           TZJ(LINMAX,NTMAX),TZZ(LINMAX,NTMAX)
      COMPLEX*16 OZZS(LINMAX,NTMAX,2),OZJS(LINMAX,NTMAX,2)
      LOGICAL CALCINT,GETIRRSOL
      REAL*8 CGC(NKMPMAX,2)
      INTEGER I,IHYPER,IKM1LIN(LINMAX),IKM2LIN(LINMAX),IL,
     &        IMT(NTMAX),IMUE,IPRINT,IQ,IQAT(NQMAX,NTMAX),
     &        IT,IWRIRRWF,IWRREGWF,J,LOPT(NTMAX),
     &        MMAX,NKM,NKMQ(NQMAX),NL,
     &        NLINQ(NQMAX),NLQ(NQMAX),NT,NUCLEUS
      INTEGER NFILCBWF
      INTEGER NSOLLM(NLMAX,NMUEMAX),LTAB(NMUEMAX),
     &        KAPTAB(NMUEMAX),NMUETAB(NMUEMAX)
      CHARACTER*10 SOLVER
      INTEGER ICALL
C
      DATA ICALL / 0 /
C
      SAVE ICALL,IKM1LIN,IKM2LIN,LOPT,NLQ,NKMQ,
     &     IQAT,IMT,
     &     NKM,IHYPER,IPRINT,IT,NT,NUCLEUS,
     &     IWRREGWF,IWRIRRWF,CALCINT,GETIRRSOL,NFILCBWF
C
      ICALL = ICALL + 1
C
C=======================================================================
C       initialise relativistic and dummy variables and SAVE them
C=======================================================================
      IF ( ICALL.EQ.1 ) THEN
C
         IF ( LMAXD.GT.NLMAX-1) THEN
            WRITE(6,*) ' LMAXD = ',LMAXD, ' > NLMAX-1 = ',NLMAX - 1
            STOP  ' Increase NLMAX in < DRVRELTMAT > '
         END IF
C         
         IF ( IRMD.GT.NRMAX ) THEN
            WRITE (6,*) ' IRMD = ',IRMD,' > NRMAX = ',NRMAX
            WRITE (6,*) ' Increase NRMAX in < sprkkr_rmesh.dim > '
            STOP ' and in < DRVRELTMAT > '
         END IF
C
         NL = LMAXD+1           ! no need to save, used only here once
         IPRINT = 0
C     
         DO I = 1,NMUEMAX
            LTAB(I) = I/2
            IF( 2*LTAB(I).EQ.I ) THEN
               KAPTAB(I) =  LTAB(I)
            ELSE
               KAPTAB(I) = -LTAB(I) - 1
            END IF
            NMUETAB(I) = 2*ABS(KAPTAB(I))
         END DO
C          
         DO IL = 1,NLMAX
            MMAX = 2*IL
            DO IMUE = 1,MMAX
               IF ( (IMUE.EQ.1) .OR. (IMUE.EQ.MMAX) ) THEN
                  NSOLLM(IL,IMUE) = 1
               ELSE
                  NSOLLM(IL,IMUE) = 2
               END IF
            END DO
         END DO
C     
         CALL IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,
     &               LINMAX,NLMAX)
C     
         CALL CALCCGC(LTAB,KAPTAB,NMUETAB,CGC,NKMAX,NMUEMAX,NKMPMAX)
C     
         DO IT = 1,NTMAX
            IMT(IT) = 1
            LOPT(IT) = -1       ! this should change for Brooks' OP
         END DO
C          
         DO IQ = 1,NQMAX
            NLQ(IQ) = NL
            NKMQ(IQ) = LMMAXD
            NLINQ(IQ) = 2*NLQ(IQ)*(2*NLQ(IQ)-1)
            IQAT(IQ,1) = 1
         END DO
C     
         NKM = LMMAXD
         IHYPER = 0
         NT = 1
         IT = 1
         NUCLEUS = 0
C     
         IWRREGWF = 0
         IWRIRRWF = 0
         CALCINT = .FALSE.
         GETIRRSOL = .FALSE.
         NFILCBWF = 0
      END IF                    ! ICALL.EQ.1
C=======================================================================
C     
      CALL SSITE(IWRREGWF,IWRIRRWF,NFILCBWF,CALCINT,GETIRRSOL,SOCTL,CTL,
     &     ERYD,P,IHYPER,IPRINT,IKM1LIN,IKM2LIN,NLQ,NKMQ,NLINQ,NT,
     &     NKM,IQAT,TSST,MSST,TSSTLIN,DZZ,DZJ,SZZ,SZJ,OZZ,OZJ,BZZ,
     &     BZJ,QZZ,QZJ,TZZ,TZJ,VT,BT,AT,ZAT,NUCLEUS,R,DRDI,R2DRDI,
     &     JWS,IMT,AMEOPO,LOPT,SOLVER,CGC,OZZS,OZJS,NLMAX,NQMAX,
     &     LINMAX,NRMAX,NMMAX,NTMAX,NKMMAX,NKMPMAX,NLAMAX)
C     
      DO J = 1,NKM
          CALL ZCOPY(NKM,TSST(1,J,IT),1,TMATLL(1,J),1)
      END DO
C     
      RETURN
      END
