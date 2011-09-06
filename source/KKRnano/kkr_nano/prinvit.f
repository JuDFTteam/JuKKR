      SUBROUTINE PRINVIT(
     >                   INFO,IAT,
     >                   NUMN0,INDN0,
     >                   TMATLL,GLLH1,X0,PRSC,SPRS,
     <                   GLLKE1)
C
      IMPLICIT NONE
C ------------------------------------------------------------------------
C Preconditioning of iterative inversion performed in invit
C Options:
C a) sc prec
C    store solution obtained at the last self-consistency iteration
C                                                         A. Thiess Nov'09
C ------------------------------------------------------------------------
C     .. parameters ..
      include 'inc.p'
      include 'inc.cls'
C
C
      INTEGER          LMMAXD
      PARAMETER        (LMMAXD= (LMAXD+1)**2)
      INTEGER          ALM
      PARAMETER        (ALM   = NAEZD*LMMAXD)
      INTEGER          NGTBD
      PARAMETER        (NGTBD = NACLSD*LMMAXD)
      REAL             CUT
      PARAMETER        (CUT   = 1.0D-5)
      DOUBLE COMPLEX   CZERO
      PARAMETER        (CZERO = ( 0.0D0,0.0D0 ))
      DOUBLE COMPLEX   CONE
      PARAMETER        (CONE  = ( 1.0D0,0.0D0 ))
c     ..
c     .. GLOBAL SCALER ARGUMENTS ..
      CHARACTER        INFO
      INTEGER          IAT
c     ..
c     .. GLOBAL ARRAY ARGUMENTS ..
      DOUBLE COMPLEX   GLLH1(LMMAXD,NGTBD,NAEZD),
     +                 TMATLL(LMMAXD,LMMAXD,NAEZD),
     +                 TMATP(ALM,LMMAXD),
     +                 GLLKE1(ALM,LMMAXD)
      COMPLEX          PRSC(NGUESSD*LMMAXD)
      INTEGER          NUMN0(NAEZD),
     +                 INDN0(NAEZD,NACLSD),
     +                 SPRS(NGUESSD*LMMAXD+1)
C     ..
C     .. LOCAL SCALARS ..
      INTEGER          I1,LM1,LM2,JLM,IL1,JSP,ISP
      LOGICAL          LSAME
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE COMPLEX   X0(ALM,LMMAXD)
      LOGICAL          DONE(LMMAXD)
C     ..
C     .. EXTERNAL FUNCTIONS ..
      EXTERNAL         LSAME
C     ..
c     ..
c-----------------------------------------------------------------------
c
C ================================================================
      IF (LSAME(INFO,'I')) THEN
C ================================================================
C Ia)  translate sparse format of PRSC back to X0 ..
C
        CALL CINIT(ALM*LMMAXD,X0)
C
        DO JSP = 1, NGUESSD*LMMAXD
C
          IF (SPRS(JSP).EQ.(NAEZD*LMMAXD*LMMAXD+9999))
     +      GOTO 99
C
          JLM = INT((SPRS(JSP)-1)/LMMAXD) + 1
          LM2 = MOD((SPRS(JSP)-1),LMMAXD) + 1
C
          X0(JLM,LM2) =
     +    DCMPLX(REAL(PRSC(JSP)),AIMAG(PRSC(JSP)))
C
        ENDDO
C
   99   CONTINUE
C ..
C ================================================================

C ================================================================
C Ib)  obtain b' = b - Ax
C                        0            ..
C
        DO LM1=1,LMMAXD
          DONE(LM1) = .FALSE.
        ENDDO
C
C---------  initialize TMATLL for the next IE by GLLKE0  ----------

        CALL SPRSZMM(                               
     >               IAT,GLLH1,NUMN0,INDN0,X0,DONE,
     >               -CONE,CZERO,                    
     <               TMATP)

C----------   \Delta t' = \Delta t - X0 + GLLHX0  -------------------
C
        DO I1=1,NAEZD
          DO LM1=1,LMMAXD
            IL1=LMMAXD*(I1-1)+LM1
            DO LM2=1,LMMAXD

              IF (I1.EQ.IAT) THEN
                TMATLL(LM1,LM2,I1) = TMATP(IL1,LM2) + TMATLL(LM1,LM2,I1)
              ELSE
                TMATLL(LM1,LM2,I1) = TMATP(IL1,LM2)
              ENDIF

            ENDDO
          ENDDO
        ENDDO

      ENDIF
C ..
C ================================================================



      IF (LSAME(INFO,'F')) THEN
C ================================================================
C Fa) determine true solution by adding initial guess ..
C
        DO I1=1,NAEZD
          DO LM1=1,LMMAXD
            IL1=LMMAXD*(I1-1)+LM1
            DO LM2=1,LMMAXD
              GLLKE1(IL1,LM2) = GLLKE1(IL1,LM2) + X0(IL1,LM2)
            ENDDO
          ENDDO
        ENDDO
C ..
C ================================================================

C ================================================================
C Fb) store new result as initial guess for the next self-consistency
C     iteration in sparse format ..
C
        ISP = 1
C
        DO I1=1,NAEZD
          DO LM1=1,LMMAXD
            JLM=LMMAXD*(I1-1)+LM1
            DO LM2=1,LMMAXD

C sparse >>
            IF ((ABS(DREAL(GLLKE1(JLM,LM2))).GT.CUT).OR.
     +          (ABS(DIMAG(GLLKE1(JLM,LM2))).GT.CUT)) THEN
C
              SPRS(ISP)=  LMMAXD*(JLM-1) + LM2
C
              PRSC(ISP) =
     +        CMPLX(DREAL(GLLKE1(JLM,LM2)),DIMAG(GLLKE1(JLM,LM2)))
C
              ISP  =  ISP + 1
C
            ENDIF
C sparse <<

            ENDDO
          ENDDO
        ENDDO
C
        SPRS(ISP) = NAEZD*LMMAXD*LMMAXD + 9999 !STOP signature
C ..
C ================================================================
C
      ENDIF
C
      RETURN
C
      END
