      SUBROUTINE IKMLIN(IPRINT,NSOLLM,IKM1LIN,IKM2LIN,NLMAX,NMUEMAX,
     &                  LINMAX,NL)
C
C   ********************************************************************
C   *                                                                  *
C   * SETUP TABLE OF INDICES    IKM(INT)                               *
C   *                                                                  *
C   *  IKM IS STANDARD INDEX IN  (KAPPA,MUE)-REPRESENTATION            *
C   *  IKM = 2*L*(J+1/2) + J + MUE + 1                                 *
C   *                                                                  *
C   *  INT NUMBERS LINEARLY ONLY NON-VANISHING ELEMENTS OF M-SS        *
C   *  USED TO CALCULATE DOS ...                                       *
C   *                                                                  *
C   ********************************************************************
C
      use mod_types, only: t_inc
      IMPLICIT NONE
C
C
C Dummy arguments
C
      INTEGER IPRINT,LINMAX,NL,NLMAX,NMUEMAX
      INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),NSOLLM(NLMAX,NMUEMAX)
C
C Local variables
C
      INTEGER I,IL,IMUE,K1,K2,KAP(2),L,LIN,MUEM05,NSOL
      INTEGER IKAPMUE
C
      LIN = 0
C
      DO IL = 1,NL
         L = IL - 1
         MUEM05 = -IL - 1
         KAP(1) = -L - 1
         KAP(2) = +L
C
         DO IMUE = 1,2*IL
            MUEM05 = MUEM05 + 1
            NSOL = NSOLLM(IL,IMUE)
C
            DO K2 = 1,NSOL
               DO K1 = 1,NSOL
                  LIN = LIN + 1
                  IKM1LIN(LIN) = IKAPMUE(KAP(K1),MUEM05)
                  IKM2LIN(LIN) = IKAPMUE(KAP(K2),MUEM05)
               END DO
            END DO
C
         END DO
      END DO
C
      IF ( IPRINT.LT.2 ) RETURN
      if(t_inc%i_write>0) then
      WRITE (1337,FMT='('' INT='',I3,''  IKM=('',I3,'','',I3,'')'')')
     &       (I,IKM1LIN(I),IKM2LIN(I),I=1,LIN)
      endif
      END
