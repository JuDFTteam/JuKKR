      SUBROUTINE CALCROTMAT( NK,IREL, ALFDEG,BETDEG,GAMDEG,
     &                   ROT, FACT, NKMMAX )
C   ********************************************************************
C   *                                                                  *
C   *   SETS UP THE ROTATION-MATRICES FOR THE EULER ANGLES             *
C   *           ( ALFDEG, BETDEG, GAMDEG )                             *
C   *                                                                  *
C   *   SEE:     E.M. ROSE  ELEMENTARY THEORY OF ANGULAR MOMENTUM      *
C   *            EQS. (4.8), (4.12) AND (4.13)                         *
C   *                                                                  *
C   *   for IREL=0,1   NK == NL           non-relativistic (l,m_l)     *
C   *       IREL=3     NK == odd          relativistic (kappa,mue)     *
C   *                                                                  *
C   *   12/11/96  HE  deal with beta = 0                               *
C   ********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMPLEX*16 CI, C0
      PARAMETER ( CI = (0.0D0,1.0D0), C0 = (0.0D0,0.0D0) )
      REAL*8 PI
      PARAMETER ( PI = 3.141592653589793238462643D0 )   
C
      REAL*8     NUM, MSB05, MSB05SQ, MSB05PW, J,M1,M2, RFAC,X
      REAL*8     FACT(0:100)

      INTEGER    S, SLOW, SHIGH, OFF
      COMPLEX*16 EMIM2A, EMIM1G, ROT(NKMMAX,NKMMAX) 
C                       
C INLINE FUNCTION    FACTORIAL FOR REAL ARGUMENT
      RFAC(X) = FACT( NINT(X) )
C
      IF( IREL .EQ. 2 ) CALL ERRORTRAP('CALCROTMAT',12,1) 
      IF( IREL .EQ. 3 .AND. MOD(NK,2).EQ.0) 
     &            CALL ERRORTRAP('CALCROTMAT',13,1) 

      DO 20 I2=1,NKMMAX 
      DO 20 I1=1,NKMMAX 
20    ROT(I1,I2) = C0
C
       CB05   =   DCOS( BETDEG*0.5D0*PI/180.0D0 )
       CB05SQ =   CB05 *  CB05
      MSB05   = - DSIN( BETDEG*0.5D0*PI/180.0D0 )
      MSB05SQ =  MSB05 * MSB05
C     
      OFF = 0
      DO 100 K=1,NK 
      IF( IREL .LT. 2 ) THEN
         L = K - 1
         J = L
      ELSE
         L = K/2
         IF( L*2 .EQ. K ) THEN
            J = L - 0.5D0
         ELSE 
            J = L + 0.5D0
         END IF         
      END IF

      NMUE = NINT( 2*J + 1 )
C
         DO 90 IM2 = 1, NMUE
         M2 = - J + (IM2-1.0D0)
         EMIM2A = CDEXP( -CI*M2*ALFDEG*PI/180.0D0 )
C
            DO 80 IM1 = 1, NMUE
            M1 = - J + (IM1-1.0D0)
            EMIM1G = CDEXP( -CI*M1*GAMDEG*PI/180.0D0 )
C
            IF( DABS(BETDEG) .LT. 1D-8 ) THEN
               IF( IM1 .EQ. IM2 ) THEN 
                  SUM = 1.0D0
               ELSE
                  SUM = 0.0D0
               END IF
            ELSE
               SLOW   = MAX(          0, NINT(M1-M2) )
               SHIGH  = MIN( NINT(J-M2), NINT( J+M1) )
                CB05PW =  CB05**NINT(2*J+M1-M2-2*SLOW    +2)
               MSB05PW = MSB05**NINT(    M2-M1+2*SLOW    -2)
               DOM = (-1.0D0)**(SLOW-1) *
     &          DSQRT( RFAC(J+M1)*RFAC(J-M1)*RFAC(J+M2)*RFAC(J-M2) )
               SUM = 0.0D0
C
               DO S=SLOW,SHIGH
                  DOM = -DOM
                  NUM =    FACT(S) * RFAC(J-M2-S)
     &                             * RFAC(J+M1-S) * RFAC(M2-M1+S)
                   CB05PW =  CB05PW /  CB05SQ
                  MSB05PW = MSB05PW * MSB05SQ
                  SUM = SUM + (DOM/NUM) * CB05PW * MSB05PW
               END DO
            END IF
C
80          ROT(OFF+IM2,OFF+IM1) = EMIM1G * SUM * EMIM2A
C
90       CONTINUE

      OFF = OFF + NMUE
100   CONTINUE
C
      RETURN
      END
