C*==relpotcvt.f    processed by SPAG 6.05Rc at 11:31 on 10 May 2004
      SUBROUTINE RELPOTCVT(ICALL,VM2Z,ZIN,RIN,DRDIIN,IRCUT,VTREL,BTREL,
     &                     ZREL,RMREL,JWSREL,DRDIREL,R2DRDIREL,IRSHIFT,
     &                     IPAND,IRMD,NPOTD,NATYPD)
C   ********************************************************************
C   *                                                                  *
C   * driving routine to convert the TB-KKR potential from the non-    *
C   * relativistic representation VM2Z(IRMD,NPOTD), with IPOTD the     *
C   * combined index for ATOM and SPIN to the relativistic one         *
C   *      VTREL =  (VUP+VDN)/2.0D0                                    *
C   *      BTREL =  (VUP-VDN)/2.0D0                                    *
C   *                                                                  *
C   * IMPORTANT 1:  because this routine is called only IF KREL.EQ.0,  *
C   *               the number of spins in VM2Z is always 2!           *
C   * IMPORTANT 2:  so far, only SPHERICAL part implemented            *
C   *                                                                  *
C   * Additionally, for compatibility with the relativistic routines   *
C   * included in the package, VTREL includes the Coulomb term, and    *
C   * the auxiliary arrays                                             *
C   *      ZREL, RMREL, JWSREL, DRDI, R2DRDI and IRSHIFT               *
C   * are created. IRSHIFT(NATYPD) accounts for the shift in the       *
C   * radial mesh, since the first point (sometimes first two points)  *
C   * of VM2Z ( = 0D0 ) are skipped. The relativistic routines require *
C   * an ODD number of radial points (Simpson integration routine)     *
C   *                                                                  *
C   * v.popescu, munich, may 2004                                      *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C PARAMETER definitions
C
      INTEGER NSPINPOT
      PARAMETER (NSPINPOT=2)
C ..
C .. Scalar arguments
      INTEGER ICALL,IPAND,IRMD,NPOTD,NATYPD
C ..
C .. Array arguments
      DOUBLE PRECISION VM2Z(IRMD,NPOTD)
      DOUBLE PRECISION ZIN(NATYPD),RIN(IRMD,NATYPD)
      DOUBLE PRECISION DRDIIN(IRMD,NATYPD)
      INTEGER IRCUT(0:IPAND,NATYPD)
C
      DOUBLE PRECISION VTREL(IRMD,NATYPD),BTREL(IRMD,NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD,NATYPD),R2DRDIREL(IRMD,NATYPD)
      DOUBLE PRECISION RMREL(IRMD,NATYPD)
      INTEGER IRSHIFT(NATYPD),JWSREL(NATYPD),ZREL(NATYPD)
C ..
C .. Local scalars
      DOUBLE PRECISION VDN, VUP
      INTEGER IT,IR,IP,IPOT,ISHIFT,JR
C ..
C .. Intrinsic Functions
      INTRINSIC NINT
C ..
C .. External Subroutines
      EXTERNAL RINIT
C
C ------------------------------------------------------- INITIALISATION
      IF ( ICALL.EQ.1 ) THEN
         CALL RINIT(IRMD*NATYPD,RMREL)
         CALL RINIT(IRMD*NATYPD,DRDIREL)
         CALL RINIT(IRMD*NATYPD,R2DRDIREL)
         DO IT = 1,NATYPD
            JWSREL(IT) = 0
            IRSHIFT(IT) = 0
            ZREL(IT) = 0
         END DO
      END IF
      CALL RINIT(IRMD*NATYPD,VTREL)
      CALL RINIT(IRMD*NATYPD,BTREL)
C ------------------------------------------------------- INITIALISATION
C
C *************************************************************** NATYPD
      DO IT = 1,NATYPD
C ================================================================ ICALL
C                                       variables require init only once
         IF ( ICALL.EQ.1 ) THEN
C
C skip first mesh point and also the second if IRCUT(1,IT) = WS-rad odd,
C since JWSREL(IT) must be odd
C
            ISHIFT = 1
            IF ( MOD(IRCUT(1,IT),2).EQ.1 ) ISHIFT = 2
            IR = 0
C ----------------------------------------------------------------------
            DO JR = 1 + ISHIFT,IRCUT(1,IT)
               IR = IR + 1
               RMREL(IR,IT) = RIN(JR,IT)
               DRDIREL(IR,IT) = DRDIIN(JR,IT)
               R2DRDIREL(IR,IT) = RMREL(IR,IT)*RMREL(IR,IT)
     &                            *DRDIREL(IR,IT)
            END DO
C ----------------------------------------------------------------------
            JWSREL(IT) = IR
            IRSHIFT(IT) = ISHIFT
            ZREL(IT) = NINT(ZIN(IT))
         END IF
C ================================================================ ICALL
         IPOT = (IT-1)*NSPINPOT + 1
         ISHIFT = IRSHIFT(IT)
C ----------------------------------------------------------------------
         DO IR = 1,JWSREL(IT)
            IP = IR + ISHIFT
            VDN = -2D0*ZIN(IT)/RIN(IP,IT) + VM2Z(IP,IPOT)
            VUP = -2D0*ZIN(IT)/RIN(IP,IT) + VM2Z(IP,IPOT+1)
            VTREL(IR,IT) = (VUP+VDN)/2.0D0
            BTREL(IR,IT) = (VUP-VDN)/2.0D0
         END DO
C ----------------------------------------------------------------------
      END DO
C *************************************************************** NATYPD
      END
