      SUBROUTINE RHOOUT(CDEN,DF,GMAT,EK,PNS,QNS,RHO2NS,THETAS,IFUNM, &
                 IPAN1,IMT1,IRMIN,IRMAX,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND &      ! Added IRMIN,IRMAX 1.7.2014
                         ,CDENLM,CWR)   ! lm-dos
!-----------------------------------------------------------------------
!
!     calculates the charge density from r(irmin) to r(irc)
!      in case of a non spherical input potential .
!
!     fills the array cden for the complex density of states
!
!     attention : the gaunt coeffients are stored in index array
!                   (see subroutine gaunt)
!
!     the structured part of the greens-function (gmat) is symmetric in
!       its lm-indices , therefore only one half of the matrix is
!       calculated in the subroutine for the back-symmetrisation .
!       the gaunt coeffients are symmetric too (since the are calculated
!       using the real spherical harmonics) . that is why the lm2- and
!       the lm02- loops are only only going up to lm1 or lm01 and the
!       summands are multiplied by a factor of 2 in the case of lm1 .ne.
!       lm2 or lm01 .ne. lm02 .
!
!             (see notes by b.drittler)
!
!                               b.drittler   aug. 1988
!-----------------------------------------------------------------------
!     .. Parameters ..
      IMPLICIT NONE
      INCLUDE 'inc.p'
!
! *********************************************************************
! * For KREL = 1 (relativistic mode)                                  *
! *                                                                   *
! *  NPOTD = 2 * NATYPD                                               *
! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! *  NSPIND = 1                                                       *
! *                                                                   *
! *********************************************************************
!
      INTEGER LMMAXD
      parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
!     ..
!     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IMT1,IPAN1,NSRA,IRMIN,IRMAX
!     ..
!     .. Array Arguments ..
DOUBLE COMPLEX CDEN(IRMD,0:*),CDENNS(*),GMAT(LMMAXD,LMMAXD), &
               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2), &
              QNSI(LMMAXD,LMMAXD), &
               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2) &
              ,CDENLM(IRMD,*),CWR(IRMD,LMMAXD,LMMAXD) ! lm-dos
DOUBLE PRECISION CLEB(*),RHO2NS(IRMD,LMPOTD),THETAS(IRID,NFUND)
INTEGER ICLEB(NCLEB,4),IFUNM(*),LMSP(*)
!     ..
!     .. Local Scalars ..
      DOUBLE COMPLEX CLTDF,CONE,CZERO
      DOUBLE PRECISION C0LL
      INTEGER I,IFUN,IR,J,L1,LM1,LM2,LM3,M1
!     ..
!     .. Local Arrays ..
DOUBLE COMPLEX WR(LMMAXD,LMMAXD,IRMIND:IRMD), &
               WR2(LMMAXD,LMMAXD,IRMIND:IRMD)
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
!     ..
!     .. Save statement ..
      SAVE
!     ..
!     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
      LOGICAL OPT
!     ..
!
!     C0LL = 1/sqrt(4*pi)
      C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
!
!
!---> initialize array for complex charge density
!
      CDEN(1:IRMD,0:LMAXD) = CZERO
      CWR(:,:,:) = CZERO
!------------------------------------------------------------------
!
!---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
!                                      summed over lm3
!---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
!                                               summed over lm3
DO 50 IR = IRMIN + 1,IRMAX
DO LM1=1,LMMAXD
DO LM2=1,LMMAXD
QNSI(LM1,LM2)=QNS(LM1,LM2,IR,1)
ENDDO
ENDDO
  CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1), &
             LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
  CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1), &
             LMMAXD,QNSI,LMMAXD,CZERO,WR(1,1,IR),LMMAXD)
  IF (NSRA.EQ.2) THEN
DO LM1=1,LMMAXD
DO LM2=1,LMMAXD
QNSI(LM1,LM2)=QNS(LM1,LM2,IR,2)
ENDDO
ENDDO
    CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2), &
               LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
    CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2), &
               LMMAXD,QNSI,LMMAXD,CONE,WR(1,1,IR),LMMAXD)
  END IF

  DO LM1 = 1,LMMAXD
    DO LM2 = 1,LM1 - 1
      WR(LM1,LM2,IR) = WR(LM1,LM2,IR) + WR(LM2,LM1,IR)
    ENDDO ! LM2
    DO LM2 = 1,LMMAXD
       WR2(LM1,LM2,IR) = WR(LM1,LM2,IR)
    ENDDO ! LM2
  ENDDO ! LM1
 50   CONTINUE  ! IR
!
!---> first calculate only the spherically symmetric contribution
!
DO 100 L1 = 0,LMAXD
  DO 70 M1 = -L1,L1
    LM1 = L1* (L1+1) + M1 + 1
    DO 60 IR = IRMIN + 1,IRMAX
!
!---> fill array for complex density of states
!
      CDEN(IR,L1) = CDEN(IR,L1) + WR(LM1,LM1,IR)
      CDENLM(IR,LM1) = WR(LM1,LM1,IR) ! lm-dos
      DO 61 LM2 = 1,LMMAXD                    ! lmlm-dos
         CWR(IR,LM1,LM2) = WR2(LM1,LM2,IR)    ! lmlm-dos
   61       CONTINUE                                ! lmlm-dos
   60     CONTINUE ! IR
   70   CONTINUE ! M1
!
!---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
!
        DO 80 IR = IRMIN + 1,IRMAX
          RHO2NS(IR,1) = RHO2NS(IR,1) + C0LL*DIMAG(CDEN(IR,L1)*DF)   ! Implicit integration over energies
   80   CONTINUE
!
!
!
        IF (IPAN1.GT.1) THEN
          DO 90 I = IMT1 + 1,IRMAX
            CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1)*C0LL

            DO M1 = -L1,L1                                                   ! lm-dos
               LM1 = L1* (L1+1) + M1 + 1                                     ! lm-dos
               CDENLM(I,LM1) = CDENLM(I,LM1)*THETAS(I-IMT1,1)*C0LL           ! lm-dos
               DO LM2 = 1,LMMAXD                                             ! lmlm-dos
                  CWR(I,LM1,LM2) = CWR(I,LM1,LM2)*THETAS(I-IMT1,1)*C0LL      ! lmlm-dos
! if LDAU, integrate up to MT
                IF (OPT('LDA+U   ')) THEN         ! LDAU
                 CWR(I,LM1,LM2) = CZERO           ! LDAU
                ENDIF                             ! LDAU
               ENDDO                                                         ! lmlm-dos
            ENDDO                                                            ! lm-dos

   90     CONTINUE
        END IF

  100 CONTINUE ! L1
!
      IF (IPAN1.GT.1) THEN
         CDENNS(1:IRMD) = 0.0D0
      END IF

      DO 140 J = 1,IEND
        LM1 = ICLEB(J,1)
        LM2 = ICLEB(J,2)
        LM3 = ICLEB(J,3)
        CLTDF = DF*CLEB(J)
!
!---> calculate the non spherically symmetric contribution
!
        DO 120 IR = IRMIN + 1,IRMAX
          RHO2NS(IR,LM3) = RHO2NS(IR,LM3) + DIMAG(CLTDF*WR(LM1,LM2,IR))
  120   CONTINUE
!
        IF (IPAN1.GT.1 .AND. LMSP(LM3).GT.0) THEN
!       IF (IPAN1.GT.1) THEN
          IFUN = IFUNM(LM3)
          DO 130 I = IMT1 + 1,IRMAX
            CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)* &
                  THETAS(I-IMT1,IFUN)
  130     CONTINUE

        END IF

  140 CONTINUE


      END
