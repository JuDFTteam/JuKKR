      MODULE MOD_RHOOUT
      CONTAINS
!-------------------------------------------------------------------------
!> Summary: Valence charge density outside MT-sphere (non-spherical part)
!> Category: physical-observables, KKRimp
!>     
!>     calculates the charge density from r(irmin) to r(irc)
!>      in case of a non spherical input potential .
!>
!>     fills the array cden for the complex density of states
!>
!>     attention : the gaunt coeffients are stored in index array
!>                   (see subroutine gaunt)
!>
!>     the structured part of the greens-function (gmat) is symmetric in
!>       its lm-indices , therefore only one half of the matrix is
!>       calculated in the subroutine for the back-symmetrisation .
!>       the gaunt coeffients are symmetric too (since the are calculated
!>       using the real spherical harmonics) . that is why the lm2- and
!>       the lm02- loops are only only going up to lm1 or lm01 and the
!>       summands are multiplied by a factor of 2 in the case of lm1 .ne.
!>       lm2 or lm01 .ne. lm02 .
!>
!>             (see notes by b.drittler)
!>
!>                        b.drittler   aug. 1988
!-------------------------------------------------------------------------
      SUBROUTINE RHOOUT(CDEN,DF,GMAT,EK,PNS,QNS,RHO2NS,THETAS,IFUNM,
     +                    IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND
     +                   ,CDENLM,  ! lm-dos
     +                   NCLEB,LMAXD,LMMAXD,LMPOTD,IRMD,IRMIND,IRID,
     +                   NFUND)
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
C     .. Parameters ..
!       INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      IMPLICIT NONE
      INTEGER LMAXD
      INTEGER LMMAXD
!       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NCLEB
      INTEGER IRMD
      INTEGER IRMIND
      INTEGER IRID
      INTEGER NFUND
!       PARAMETER (IRMIND=IRMD-IRNSD)
C     ..
C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IMT1,IPAN1,NSRA
C     ..
C     .. Array Arguments ..
      DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),
     +               GMAT(LMMAXD,LMMAXD),
     +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
     +              QNSI(LMMAXD,LMMAXD),
     +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
     +              ,CDENLM(IRMD,LMMAXD) ! lm-dos
      DOUBLE PRECISION CLEB(NCLEB),RHO2NS(IRMD,LMPOTD),
     +                 THETAS(IRID,NFUND)
      INTEGER ICLEB(NCLEB,4),IFUNM(*),LMSP(*)
C     ..
C     .. Local Scalars ..
      DOUBLE COMPLEX CLTDF,CONE,CZERO
      DOUBLE PRECISION C0LL
      INTEGER I,IFUN,IR,J,L1,LM1,LM2,LM3,M1
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX WR(LMMAXD,LMMAXD,IRMIND:IRMD)
C     ..
C     .. External Subroutines ..
      EXTERNAL ZGEMM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
C     ..
C     .. Save statement ..
      SAVE
C     ..
C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
C     ..
c
C     C0LL = 1/sqrt(4*pi)
      C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
c
c
c---> initialize array for complex charge density
c
      DO L1 = 0,LMAXD
        DO I = 1,IRMD
          CDEN(I,L1) = CZERO
      ENDDO    
      ENDDO    
C------------------------------------------------------------------
c
c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
c                                      summed over lm3
c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
c                                               summed over lm3
      DO 50 IR = IRMIND + 1,IRMD
      DO LM1=1,LMMAXD
      DO LM2=1,LMMAXD
      QNSI(LM1,LM2)=QNS(LM1,LM2,IR,1)
      ENDDO
      ENDDO
        CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1),
     +             LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
        CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1),
     +             LMMAXD,QNSI,LMMAXD,CZERO,WR(1,1,IR),LMMAXD)
        IF (NSRA.EQ.2) THEN
      DO LM1=1,LMMAXD
      DO LM2=1,LMMAXD
      QNSI(LM1,LM2)=QNS(LM1,LM2,IR,2)
      ENDDO
      ENDDO
          CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2),
     +               LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
          CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2),
     +               LMMAXD,QNSI,LMMAXD,CONE,WR(1,1,IR),LMMAXD)
        END IF

        DO LM1 = 1,LMMAXD
          DO LM2 = 1,LM1 - 1
            WR(LM1,LM2,IR) = WR(LM1,LM2,IR) + WR(LM2,LM1,IR)
      ENDDO
      ENDDO
 50   CONTINUe
c
c---> first calculate only the spherically symmetric contribution
c
      DO 100 L1 = 0,LMAXD
        DO 70 M1 = -L1,L1
          LM1 = L1* (L1+1) + M1 + 1
          DO 60 IR = IRMIND + 1,IRMD
c
c---> fill array for complex density of states
c
            CDEN(IR,L1) = CDEN(IR,L1) + WR(LM1,LM1,IR)
            CDENLM(IR,LM1) = WR(LM1,LM1,IR) ! lm-dos
   60     CONTINUE
   70   CONTINUE
c
c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
c
        DO 80 IR = IRMIND + 1,IRMD
          RHO2NS(IR,1) = RHO2NS(IR,1) + C0LL*DIMAG(CDEN(IR,L1)*DF)
   80   CONTINUE
c
c
c
        IF (IPAN1.GT.1) THEN
          DO 90 I = IMT1 + 1,IRMD
            CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1)*C0LL

            DO M1 = -L1,L1                                         ! lm-dos
               LM1 = L1* (L1+1) + M1 + 1                           ! lm-dos
               CDENLM(I,LM1) = CDENLM(I,LM1)*THETAS(I-IMT1,1)*C0LL ! lm-dos
            ENDDO                                                  ! lm-dos

   90     CONTINUE
        END IF

  100 CONTINUE
c
      IF (IPAN1.GT.1) THEN
        DO 110 I = 1,IRMD
          CDENNS(I) = 0.0D0
  110   CONTINUE
      END IF

      DO 140 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        CLTDF = DF*CLEB(J)
c
c---> calculate the non spherically symmetric contribution
c
        DO 120 IR = IRMIND + 1,IRMD
          RHO2NS(IR,LM3) = RHO2NS(IR,LM3) + DIMAG(CLTDF*WR(LM1,LM2,IR))
  120   CONTINUE
c
        IF (IPAN1.GT.1 .AND. LMSP(LM3).GT.0) THEN
c       IF (IPAN1.GT.1) THEN
          IFUN = IFUNM(LM3)
          DO 130 I = IMT1 + 1,IRMD
            CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
     +                  THETAS(I-IMT1,IFUN)
  130     CONTINUE

        END IF

  140 CONTINUE


      END SUBROUTINE
      END MODULE MOD_RHOOUT 
