C*==cpamillsx.f    processed by SPAG 6.05Rc at 10:24 on 15 Dec 2003
      SUBROUTINE CPAMILLSX(ITCPA,CPAERR,CPACORR,CPACHNG,IPRINT,ICPA,NQ,
     &                     NKMQ,NOQ,ITOQ,CONC,MSSQ,MSST,TAUQ,DMSSQ,
     &                     KMROT,DROTQ,NTMAX,NQMAX,NKMMAX)
C   ********************************************************************
C   *                                                                  *
C   *   perform  CPA-iteration according    MILLS's  algorithm         *
C   *                                                                  *
C   *   the CPA - iteration step for site IQ is omitted if             *
C   *   ICPA(IQ) = 0 ( set in <INITALL> )                              *
C   *                                                                  *
C   *   only the projection matrix  DTILT(G) is needed                 *
C   *   this is set up with respect to the global frame                *
C   *   for this reason MSST has to be rotated prior calling <GETDMAT> *
C   *                                                                  *
C   *   allows an atom type IT to have different orientation of        *
C   *   its moment on different but equivalent sites  IQ               *
C   *                                                                  *
C   * 15/12/03                                                         *
C   ********************************************************************
      use mod_types, only: t_inc
      IMPLICIT COMPLEX*16(A-H,O-Z)
C
C PARAMETER definitions
C
      REAL*8 TOL,SCLSTD
      PARAMETER (TOL=10D0,SCLSTD=1D0)
      COMPLEX*16 C0,C1
      PARAMETER (C0=(0.0D0,0.0D0),C1=(1.0D0,0.0D0))
C
C Dummy arguments
C
      REAL*8 CPACHNG,CPACORR,CPAERR
      INTEGER IPRINT,ITCPA,KMROT,NKMMAX,NQ,NQMAX,NTMAX
      REAL*8 CONC(NTMAX)
      COMPLEX*16 DROTQ(NKMMAX,NKMMAX,NQMAX),MSSQ(NKMMAX,NKMMAX,NQMAX),
     &     DMSSQ(NKMMAX,NKMMAX,NQMAX),
     &     MSST(NKMMAX,NKMMAX,NTMAX),TAUQ(NKMMAX,NKMMAX,NQMAX)
      INTEGER ICPA(NQMAX),ITOQ(NTMAX,NQMAX),NKMQ(NQMAX),NOQ(NQMAX)
C
C Local variables
C
      LOGICAL CHECK
      REAL*8 CPACHNGL,CPACORRL,SCL
      COMPLEX*16 CSUM
      COMPLEX*16 DMAMC(NKMMAX,NKMMAX),DMATTG(NKMMAX,NKMMAX),
     &     DQ(NKMMAX,NQMAX),
     &     DTILTG(NKMMAX,NKMMAX),ERR(NKMMAX,NKMMAX),
     &     W1(NKMMAX,NKMMAX),W2(NKMMAX,NKMMAX)
      INTEGER I,ICPARUN,INFO,IO
      INTEGER IPIV(NKMMAX),IQ,IT,IW,IW0,J,M,N
      DOUBLE PRECISION P1,P2
      SAVE CPACHNGL,CPACORRL,SCL
C
      DATA ICPARUN/0/
C
      CPAERR = 0.0D0
      CPACORR = 0.0D0
      CPACHNG = 0.0D0
      CHECK = .TRUE.
      CHECK = .FALSE.
C
      IF ( ITCPA.EQ.1 ) THEN
         SCL = SCLSTD
         CPACHNGL = 1D+20
         CPACORRL = 1D+20
         ICPARUN = ICPARUN + 1
      END IF
C
      DO IQ = 1,NQ
         IF ( ICPA(IQ).NE.0 ) THEN
C
            M = NKMMAX
            N = NKMQ(IQ)
C
            DO J = 1,N
               DQ(J,IQ) = -C1
               CALL CINIT(N,ERR(1,J))
            END DO
C
C================================================================= IT ==
            DO IO = 1,NOQ(IQ)
               IT = ITOQ(IO,IQ)
C
C ------------------------- rotate the single site m-matrix if necessary
               IF ( KMROT.NE.0 ) THEN
C
                  CALL ROTATE(MSST(1,1,IT),'L->G',W1,N,DROTQ(1,1,IQ),M)
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMATTG,DTILTG,DMAMC,N,
     &                         MSSQ(1,1,IQ),W1,M)
C
               ELSE
C
                  CALL GETDMAT(TAUQ(1,1,IQ),DMATTG,DTILTG,DMAMC,N,
     &                         MSSQ(1,1,IQ),MSST(1,1,IT),M)
C
               END IF
C
               DO I = 1,N
                  DQ(I,IQ) = DQ(I,IQ) + CONC(IT)*DTILTG(I,I)
               END DO
C
C     -------------------------------------------
C            - E[a] = D~[a] * ( m[a] - m[c] )
C     -------------------------------------------
               CALL ZGEMM('N','N',N,N,N,C1,DTILTG(1,1),M,DMAMC,M,C0,W1,
     &                    M)
C
C     -------------------------------------------
C            E = SUM[a]  c[a] *  E[a]
C     -------------------------------------------
               DO J = 1,N
                  DO I = 1,N
                     ERR(I,J) = ERR(I,J) - CONC(IT)*W1(I,J)
                  END DO
               END DO
C
            END DO
C================================================================= IT ==
C
C     -------------------------------------------
C                   E * TAU
C     -------------------------------------------
C
            CALL ZGEMM('N','N',N,N,N,C1,ERR,M,TAUQ(1,1,IQ),M,C0,W2,M)
C
C     -------------------------------------------
C                1 + E * TAU
C     -------------------------------------------
            DO I = 1,N
               W2(I,I) = C1 + W2(I,I)
            END DO
C
C     -------------------------------------------
C               ( 1 + E * TAU )**(-1)
C     -------------------------------------------
C
            CALL ZGETRF(N,N,W2,M,IPIV,INFO)
            CALL ZGETRI(N,W2,M,IPIV,W1,M*M,INFO)
C
C     -------------------------------------------
C           ( 1 + E * TAU )**(-1) * E
C     -------------------------------------------
C
            CALL ZGEMM('N','N',N,N,N,C1,W2,M,ERR,M,C0,W1,M)
C
C     -------------------------------------------
C           APPLY CORRECTION  TO  MEFF
C     m{n+1} = m{n} -  ( 1 + E * TAU )**(-1) * E
C     -------------------------------------------
            DO J = 1,N
               CPAERR = CPAERR + ABS(DREAL(DQ(J,IQ)))
     &                  + ABS(DIMAG(DQ(J,IQ)))
               CPACORR = CPACORR + ABS(DREAL(W1(J,J)))
     &                   + ABS(DIMAG(W1(J,J)))
               CPACHNG = MAX(CPACHNG,ABS(W1(J,J)/MSSQ(J,J,IQ)))
            END DO
            CPACHNG = SCL*CPACHNG
            CPACORR = SCL*CPACORR
C
            IF ( CPACHNG.GT.TOL*CPACHNGL .OR. CPACORR.GT.CPACORRL ) THEN
               WRITE (*,*) '############### CPA step back'
C
               P1 = 0.5D0   ! P1 = 0.05D0
               P2 = 1D0 - P1
               CPACHNG = P1*CPACHNGL
               CPACORR = P1*CPACORRL
               CPACHNG = CPACHNGL
               CPACORR = CPACORRL
               SCL = P1
               DO J = 1,N
                  DO I = 1,N
                     MSSQ(I,J,IQ) = MSSQ(I,J,IQ) + P2*DMSSQ(I,J,IQ)
                     DMSSQ(I,J,IQ) = DMSSQ(I,J,IQ)*P1
                  END DO
               END DO
            ELSE
               DO J = 1,N
                  DO I = 1,N
                     W1(I,J) = SCL*W1(I,J)
                     MSSQ(I,J,IQ) = MSSQ(I,J,IQ) - W1(I,J)
                     DMSSQ(I,J,IQ) = W1(I,J)
                  END DO
               END DO
            END IF
C
            CPAERR = CPACHNG
C
            IF ( IPRINT.GE.2 .OR. CHECK ) THEN
               CSUM = C0
               DO I = 1,N
                  CSUM = CSUM + MSSQ(I,I,IQ)
               END DO
               if(t_inc%i_write>0) then
                  WRITE (1337,99001) IQ,CPAERR,CPACORR,CSUM
               endif
            END IF
C-----------------------------------------------------------------------
            IF ( CHECK ) THEN
               IF ( ICPARUN.EQ.2 ) STOP 'CPA-iter written to for...'
               IW0 = 100*IQ
               DO I = 1,N
                  IW = IW0 + I
                  WRITE (IW,'(I4,4E14.4)') ITCPA,MSSQ(I,I,IQ),W1(I,I)
               END DO
            END IF
C-----------------------------------------------------------------------
         END IF
C
      END DO
C================================================================= IQ ==
C
      CPACHNGL = CPACHNG
      CPACORRL = CPACORR
C
99001 FORMAT (' CPA:  IQ',I3,'  ERR',F12.5,'  CORR',F13.5,'  M',
     &        18(1X,2(1PE14.6)))
      END
