C*==projtau.f    processed by SPAG 6.05Rc at 18:34 on 11 Nov 2000
      SUBROUTINE PROJTAU(ICPAFLAG,CPACHNG,KMROT,WRTAU,WRTAUMQ,IFILTAU,
     &                   ERYD,NT,NQ,NKMQ,MSST,MSSQ,NLINQ,IQAT,CONC,TAUQ,
     &                   TAUT,TAUTLIN,IKM1LIN,IKM2LIN,DROTQ,NTMAX,NQMAX,
     &                   NKMMAX,LINMAX)
C   ********************************************************************
C   *                                                                  *
C   *   calculate the component projected TAU - matrices               *
C   *                                                                  *
C   *      TAU(IT) =  TAU(IQ) * ( 1 + (m(t)-m(c))*TAU(IQ) )**(-1)      *
C   *                                                                  *
C   *   NOTE: it is assumed that all equivalent sites  IQ  have the    *
C   *   same TAU-matrix  TAUQ(IQ). To get  TAU(IT)  the first site IQ  *
C   *   occupied by type IT is taken to be representative for          *
C   *   all other (NAT(IT)-1) sites occupied by IT                     *
C   *                                                                  *
C   *   allows an atom type IT to have different orientation of        *
C   *   its moment on different but equivalent sites  IQ               *
C   *                                                                  *
C   * 01/11/00                                                         *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C PARAMETER definitions
C
      DOUBLE PRECISION TOL
      PARAMETER (TOL=1.0D-6)
      DOUBLE COMPLEX C0,C1
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      DOUBLE PRECISION CPACHNG
      DOUBLE COMPLEX ERYD
      INTEGER ICPAFLAG,IFILTAU,KMROT,LINMAX,NKMMAX,NQ,NQMAX,NT,NTMAX
      LOGICAL WRTAU,WRTAUMQ
      DOUBLE PRECISION CONC(NTMAX)
      DOUBLE COMPLEX DROTQ(NKMMAX,NKMMAX,NQMAX),
     &           MSSQ(NKMMAX,NKMMAX,NQMAX),
     &           MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX),
     &           TAUT(NKMMAX,NKMMAX,NTMAX),TAUTLIN(LINMAX,NTMAX)
      INTEGER IKM1LIN(LINMAX),IKM2LIN(LINMAX),IQAT(NTMAX),
     &        NKMQ(NQMAX),NLINQ(NQMAX)
C
C Local variables
C
      DOUBLE PRECISION CPAC
      DOUBLE COMPLEX DMAMC(NKMMAX,NKMMAX),DMATTG(NKMMAX,NKMMAX),
     &           DTILTG(NKMMAX,NKMMAX),RMSS,RTAU,W1(NKMMAX,NKMMAX)
      INTEGER I,ICPAF,IQ,IT,J,LIN,M,N
C
C*** End of declarations rewritten by SPAG
C
C
      DO IT = 1,NT
C
C ---------- pick first site IQ occupied by type IT to be representative
C ----------- all other (NAT(IT)-1) occupied sites have to be equivalent
C
         IQ = IQAT(IT)
         M = NKMMAX
         N = NKMQ(IQ)
C
         IF ( CONC(IT).LT.0.995 ) THEN
C
C ------------------------- rotate the single site m-matrix if necessary
            IF ( KMROT.NE.0 ) THEN
C
               CALL ROTATE(MSST(1,1,IT),'L->G',W1,N,DROTQ(1,1,IQ),M)
C
               CALL GETDMAT(TAUQ(1,1,IQ),DMATTG,DTILTG,DMAMC,N,
     &                      MSSQ(1,1,IQ),W1,M)
C
            ELSE
C
               CALL GETDMAT(TAUQ(1,1,IQ),DMATTG,DTILTG,DMAMC,N,
     &                      MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
            END IF
C
C     -------------------------------------------
C              TAU(t) = TAU * D~(t)
C     ----------------------------------------
            CALL ZGEMM('N','N',N,N,N,C1,TAUQ(1,1,IQ),M,DTILTG,M,C0,
     &                 TAUT(1,1,IT),M)
C
            ICPAF = ICPAFLAG
            CPAC = CPACHNG
C
         ELSE
C
C     CONC > 0.995:  COPY TAU TO TAUTLIN
C
            DO J = 1,N
               CALL ZCOPY(N,TAUQ(1,J,IQ),1,TAUT(1,J,IT),1)
            END DO
C
            ICPAF = 0
            CPAC = 0D0
C
         END IF
C
C     -------------------------------------------
C            rotate  TAU(t)  if required
C     -------------------------------------------
C
         IF ( KMROT.NE.0 ) THEN
C
            DO J = 1,N
               CALL ZCOPY(N,TAUT(1,J,IT),1,W1(1,J),1)
            END DO
C
            CALL ROTATE(W1,'G->L',TAUT(1,1,IT),N,DROTQ(1,1,IQ),M)
C
         END IF
C
C     -------------------------------------------
C        STORE TAU(t) IN LINEAR ARRAY TAUTLIN
C     -------------------------------------------
C
         DO LIN = 1,NLINQ(IQ)
            TAUTLIN(LIN,IT) = TAUT(IKM1LIN(LIN),IKM2LIN(LIN),IT)
         END DO
C
         IF ( WRTAU ) THEN
            WRITE (IFILTAU,99001) ERYD,IT,IQ,ICPAF,CPAC
            DO I = 1,N
               DO J = 1,N
                  IF ( I.EQ.J ) THEN
                     WRITE (IFILTAU,99003) I,J,TAUT(I,J,IT)
                  ELSE
                     IF ( CDABS(TAUT(I,J,IT)/TAUT(I,I,IT)).GT.TOL )
     &                    WRITE (IFILTAU,99003) I,J,TAUT(I,J,IT)
                  END IF
               END DO
            END DO
         END IF
C
      END DO
C================================================================= IT ==
C
      IF ( WRTAUMQ ) THEN
         DO IQ = 1,NQ
            WRITE (IFILTAU,99002) ERYD,IQ,ICPAFLAG,CPACHNG
            DO I = 1,N
               DO J = 1,N
                  IF ( I.EQ.J ) THEN
                     WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),MSSQ(I,J,IQ)
                  ELSE
                     RTAU = TAUQ(I,J,IQ)/TAUQ(I,I,IQ)
                     RMSS = MSSQ(I,J,IQ)/MSSQ(I,I,IQ)
                     IF ( (CDABS(RTAU).GT.TOL) .OR. (CDABS(RMSS).GT.TOL)
     &                    ) WRITE (IFILTAU,99003) I,J,TAUQ(I,J,IQ),
     &                             MSSQ(I,J,IQ)
                  END IF
C
               END DO
            END DO
C
         END DO
      END IF
C--------------------------------------------------------- FORMAT IFMT=2
99001 FORMAT (/,80('*')/,2F21.15,' RYD   TAU FOR IT=',I2,'  IQ=',I2,:,
     &        '  CPA:',I2,F15.6)
99002 FORMAT (/,80('*')/,2F21.15,' RYD   TAU-C M-C  FOR IQ=',I2,:,
     &        '  CPA:',I2,F15.6)
99003 FORMAT (2I5,1P,4E22.14)
C
      END
