      MODULE MOD_RHOOUTNEW
       CONTAINS




      SUBROUTINE RHOOUTNEW(gauntcoeff,DF,GMAT,EK,cellnew,wavefunction,RHO2NSC, &
                          NSRA, &
                      LMAXD,LMMAXD,LMSIZE,LMSIZE2,LMPOTD,IRMD,irmind,&
                         ISPIN,NSPINDEN)
use type_gauntcoeff
use type_cellnew
use type_wavefunction

!       CALL RHOOUTNEW(density%den(:,ispin,ie),DF,GMATll,EK,RLL(:,:,:,1),SLL(:,:,:,1),density%RHO2NS(:,:,ispin),shapefun%thetas,shapefun%lmused,cell%npan, &
!                     1,shapefun%lm2index,CDENNS,config%NSRA,gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND &
!                    ,CDENLM, & ! lm-dos
!                     gauntcoeff%NCLEB,LMAXD,LMMAXATOM,(2*LMAXD+1)**2,cellnew%nrmaxnew,1,IRID, &
!                     shapefun%nlmshaped)
! c-----------------------------------------------------------------------
! c
! c     calculates the charge density from r(irmin) to r(irc)
! c      in case of a non spherical input potential .
! c
! c     fills the array cden for the complex density of states
! c
! c     attention : the gaunt coeffients are stored in index array
! c                   (see subroutine gaunt)
! c
! c     the structured part of the greens-function (gmat) is symmetric in
! c       its lm-indices , therefore only one half of the matrix is
! c       calculated in the subroutine for the back-symmetrisation .
! c       the gaunt coeffients are symmetric too (since the are calculated
! c       using the real spherical harmonics) . that is why the lm2- and
! c       the lm02- loops are only only going up to lm1 or lm01 and the
! c       summands are multiplied by a factor of 2 in the case of lm1 .ne.
! c       lm2 or lm01 .ne. lm02 .
! c
! c             (see notes by b.drittler)
! c
! c                               b.drittler   aug. 1988
! c-----------------------------------------------------------------------
! C     .. Parameters ..
! !       INCLUDE 'inc.p'
! C
! C *********************************************************************
! C * For KREL = 1 (relativistic mode)                                  *
! C *                                                                   *
! C *  NPOTD = 2 * NATYPD                                               *
! C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! C *  NSPIND = 1                                                       *
! C *                                                                   *
! C *********************************************************************
! C
      IMPLICIT NONE
      TYPE(GAUNTCOEFF_TYPE)  :: GAUNTCOEFF
      TYPE(CELL_TYPENEW)  :: CELLNEW
      TYPE(WAVEFUNCTION_TYPE)  :: WAVEFUNCTION
      INTEGER LMAXD
      INTEGER LMMAXD
      INTEGER LMSIZE,LMSIZE2
!       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER NCLEB
      INTEGER IRMD
      INTEGER IRMIND
      INTEGER IRID
      INTEGER NFUND
      INTEGER NSPIN
!       PARAMETER (IRMIND=IRMD-IRNSD)
! C     ..
! C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IEND,IMT1,IPAN1,NSRA
! C     ..
! C     .. Array Arguments ..
      DOUBLE COMPLEX CDEN(NSPINDEN,IRMD,0:LMAXD),CDENNS(IRMD), &
                     GMAT(LMSIZE,LMSIZE), &
!                      PNS(LMSIZE,LMSIZE,IRMD,2), &
                    QNSI(LMSIZE,LMSIZE), &
!                      QNS(LMSIZE,LMSIZE,IRMD,2) &
                    CDENLM(NSPINDEN,IRMD,LMMAXD), & ! lm-dos
                    RHO2NSC(IRMD,LMPOTD,NSPINDEN)
!       DOUBLE PRECISION  &
!                        THETAS(IRID,NFUND)
!       INTEGER IFUNM(*),LMSP(*)
! C     ..
! C     .. Local Scalars ..
      DOUBLE COMPLEX CLTDF,CONE,CZERO
      DOUBLE PRECISION C0LL
      INTEGER I,IFUN,IR,J,L1,LM1,LM2,LM3,M1
! C     ..
! C     .. Local Arrays ..
      DOUBLE COMPLEX WR(LMSIZE,LMSIZE,IRMD)
      INTEGER          :: ISPIN,JSPIN

      INTEGER     :: SPININDEX1(4) !=(/1,2,1,2 /)
      INTEGER     :: SPININDEX2(4) !=(/1,2,2,1 /)
      INTEGER               :: LMSHIFT1(4)
      INTEGER               :: LMSHIFT2(4)
      INTEGER               :: NSPINSTART,NSPINSTOP,NSPINDEN
! C     ..
! C     .. External Subroutines ..
      EXTERNAL ZGEMM
! C     ..
! C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
! C     ..
! C     .. Save statement ..
      SAVE
! C     ..
! C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
! C     ..
! c


      IF (NSPINDEN==4) THEN
        NSPINSTART=1
        NSPINSTOP=NSPINDEN
        SPININDEX1   =(/1,2,1,2 /)
        SPININDEX2   =(/1,2,2,1 /)
        LMSHIFT1=(LMAXD+1)**2*(SPININDEX1-1)
        LMSHIFT2=(LMAXD+1)**2*(SPININDEX2-1)
      ELSE
        NSPINSTART=ISPIN
        NSPINSTOP=ISPIN
        SPININDEX1   =(/1,1,0,0 /)
        SPININDEX2   =(/1,1,0,0 /)
        LMSHIFT1=(LMAXD+1)**2*(SPININDEX1-1)
        LMSHIFT2=(LMAXD+1)**2*(SPININDEX2-1)
      END IF
! print *,LMSHIFT1
! print *,LMSHIFT2

!    write(*,*) 'lmsize',lmsize,nspinstart,nspinstop


! C     C0LL = 1/sqrt(4*pi)
      C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
! c
! c
! c---> initialize array for complex charge density
! c

      DO L1 = 0,LMAXD
        DO I = 1,IRMD
          CDEN(:,I,L1) = CZERO
      ENDDO    
      ENDDO    
! C------------------------------------------------------------------
! c
! c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
! c                                      summed over lm3
! c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
! c                                               summed over lm3

!          GMAT=(1.0D0,1.0D0)
!          CELLNEW%SLL(:,:,:,1)=(1.0D0,1.0D0)
!          CELLNEW%RLL(:,:,:,1)=(1.0D0,1.0D0)

      DO    IR = 1,IRMD
        DO LM1=1,LMSIZE
          DO LM2=1,LMSIZE
            QNSI(LM1,LM2)=WAVEFUNCTION%SLL(LM1,LM2,IR,1)
          ENDDO
        ENDDO
        CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,WAVEFUNCTION%RLL(1,1,IR,1), &
                  LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)
        CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,WAVEFUNCTION%RLL(1,1,IR,1), &
                    LMSIZE,QNSI,LMSIZE,CZERO,WR(1,1,IR),LMSIZE)
        IF (NSRA.EQ.2) THEN
          DO LM1=1,LMSIZE
            DO LM2=1,LMSIZE
            QNSI(LM1,LM2)=WAVEFUNCTION%SLL(LM1,LM2,IR,2)
            ENDDO
          ENDDO
          CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,WAVEFUNCTION%RLL(1,1,IR,2), &
                    LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)
          CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,WAVEFUNCTION%RLL(1,1,IR,2), &
                    LMSIZE,QNSI,LMSIZE,CONE,WR(1,1,IR),LMSIZE)
        END IF

!           if (ir==IRMD) then
!           DO LM1=1,LMSIZE
!            write(5432,'(5000E)') WR(lm1,:,IR)
!           end do
!           stop
!           end if
        DO JSPIN=NSPINSTART,NSPINSTOP
          DO LM1 = 1,LMMAXD
            DO LM2 = 1,LM1 - 1
              WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR) = &
               WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR) + WR(LM2+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
            ENDDO
          ENDDO
        END DO !JSPIN



      END DO !IR


! c
! c---> first calculate only the spherically symmetric contribution
! c

      DO L1 = 0,LMAXD
        DO M1 = -L1,L1
          LM1 = L1* (L1+1) + M1 + 1
          DO IR = 1,IRMD
! c
! c---> fill array for complex density of states

! c
            DO JSPIN=NSPINSTART,NSPINSTOP
              CDEN(JSPIN,IR,L1) = CDEN(JSPIN,IR,L1)  &
                     + WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
!               CDENLM(JSPIN,IR,LM1) = &
!                        WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR) ! lm-dos
            END DO

          END DO !IR
        END DO !M1

! c
! c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
! c


        DO JSPIN=NSPINSTART,NSPINSTOP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! c
! c---> check this routine!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111! c





!             write(1337,*) '[WARNING TO BE CHECKED]'



            DO IR = 1,IRMD
!             DO IR = irmind,IRMD


!               RHO2NSC(JSPIN,IR,1) = RHO2NSC(JSPIN,IR,1) 
!      +                          + C0LL*DIMAG(CDEN(JSPIN,IR,L1)*DF)
              RHO2NSC(IR,1,JSPIN) = RHO2NSC(IR,1,JSPIN)  &
                                + C0LL*(CDEN(JSPIN,IR,L1)*DF)
            END DO
        END DO




!         IF (IPAN1.GT.1) THEN
!           DO 90 I = IMT1 + 1,IRMD
!             CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1)*C0LL
! 
!             DO M1 = -L1,L1                                         ! lm-dos
!                LM1 = L1* (L1+1) + M1 + 1                           ! lm-dos
!                CDENLM(I,LM1) = CDENLM(I,LM1)*THETAS(I-IMT1,1)*C0LL ! lm-dos
!             ENDDO                                                  ! lm-dos
! 
!    90     CONTINUE
!         END IF

      END DO !L1

!       IF (IPAN1.GT.1) THEN
!         DO 110 I = 1,IRMD
!           CDENNS(I) = 0.0D0
!   110   CONTINUE
!       END IF

      DO J = 1,gauntcoeff%IEND
        LM1 = gauntcoeff%ICLEB(J,1)
        LM2 = gauntcoeff%ICLEB(J,2)
        LM3 = gauntcoeff%ICLEB(J,3)
        CLTDF = DF*gauntcoeff%CLEB(J,1)
!         print *, lm3
! c
! c---> calculate the non spherically symmetric contribution
! c
!         write(1337,*) '[WARNING TO BE CHECKED]'
        DO JSPIN=NSPINSTART,NSPINSTOP
            DO IR = 1,IRMD
!             DO IR = IRMIND,IRMD

!               RHO2NSC(JSPIN,IR,LM3) = RHO2NSC(JSPIN,IR,LM3) 
!      +                  + DIMAG(CLTDF*WR(LM1+LMSHIFT1(JSPIN),
!      +                                   LM2+LMSHIFT2(JSPIN),IR))
              RHO2NSC(IR,LM3,JSPIN) = RHO2NSC(IR,LM3,JSPIN)  &
                        + (CLTDF*WR(LM1+LMSHIFT1(JSPIN), &
                                         LM2+LMSHIFT2(JSPIN),IR)) !&
!                         + (CLTDF*WR(LM2+LMSHIFT2(JSPIN), &
!                                          LM1+LMSHIFT1(JSPIN),IR)) 
            END DO
        END DO
! c
!         IF (IPAN1.GT.1 .AND. LMSP(LM3).GT.0) THEN
! c       IF (IPAN1.GT.1) THEN
!           IFUN = IFUNM(LM3)
!           DO 130 I = IMT1 + 1,IRMD
!             CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
!      +                  THETAS(I-IMT1,IFUN)
!   130     CONTINUE

!         END IF

      END DO !J
!          do lm1=1,(2*LMAXD+1)**2
!            write(9013,'(5000F)') rho2nsc(:,lm1,nspinstart)
!          end do 
! stop
      END SUBROUTINE

! 
!       SUBROUTINE RHOOUTNEW_7Feb2011(DF,GMAT,EK,PNS,QNS,RHO2NS,THETAS,
!      +                   IFUNM,
!      +                    IPAN1,IMT1,LMSP,CDENNS,NSRA,CLEB,ICLEB,IEND
!      +                   ,CDENLM,  ! lm-dos
!      +                   NCLEB,LMAXD,LMMAXD,LMSIZE,LMPOTD,IRMD,IRMIND,
!      +                   IRID, NFUND)
! 
! !       CALL RHOOUTNEW(density%den(:,ispin,ie),DF,GMATll,EK,RLL(:,:,:,1),SLL(:,:,:,1),density%RHO2NS(:,:,ispin),shapefun%thetas,shapefun%lmused,cell%npan, &
! !                     1,shapefun%lm2index,CDENNS,config%NSRA,gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND &
! !                    ,CDENLM, & ! lm-dos
! !                     gauntcoeff%NCLEB,LMAXD,LMMAXATOM,(2*LMAXD+1)**2,cellnew%nrmaxnew,1,IRID, &
! !                     shapefun%nlmshaped)
! c-----------------------------------------------------------------------
! c
! c     calculates the charge density from r(irmin) to r(irc)
! c      in case of a non spherical input potential .
! c
! c     fills the array cden for the complex density of states
! c
! c     attention : the gaunt coeffients are stored in index array
! c                   (see subroutine gaunt)
! c
! c     the structured part of the greens-function (gmat) is symmetric in
! c       its lm-indices , therefore only one half of the matrix is
! c       calculated in the subroutine for the back-symmetrisation .
! c       the gaunt coeffients are symmetric too (since the are calculated
! c       using the real spherical harmonics) . that is why the lm2- and
! c       the lm02- loops are only only going up to lm1 or lm01 and the
! c       summands are multiplied by a factor of 2 in the case of lm1 .ne.
! c       lm2 or lm01 .ne. lm02 .
! c
! c             (see notes by b.drittler)
! c
! c                               b.drittler   aug. 1988
! c-----------------------------------------------------------------------
! C     .. Parameters ..
! !       INCLUDE 'inc.p'
! C
! C *********************************************************************
! C * For KREL = 1 (relativistic mode)                                  *
! C *                                                                   *
! C *  NPOTD = 2 * NATYPD                                               *
! C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! C *  NSPIND = 1                                                       *
! C *                                                                   *
! C *********************************************************************
! C
!       IMPLICIT NONE
!       INTEGER LMAXD
!       INTEGER LMMAXD
!       INTEGER LMSIZE
! !       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
!       INTEGER LMPOTD
! !       PARAMETER (LMPOTD= (LPOTD+1)**2)
!       INTEGER NCLEB
!       INTEGER IRMD
!       INTEGER IRMIND
!       INTEGER IRID
!       INTEGER NFUND
! !       PARAMETER (IRMIND=IRMD-IRNSD)
! C     ..
! C     .. Scalar Arguments ..
!       DOUBLE COMPLEX DF,EK
!       INTEGER IEND,IMT1,IPAN1,NSRA
! C     ..
! C     .. Array Arguments ..
!       DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),
!      +               GMAT(LMSIZE,LMSIZE),
!      +               PNS(LMSIZE,LMSIZE,IRMIND:IRMD,2),
!      +              QNSI(LMSIZE,LMSIZE),
!      +               QNS(LMSIZE,LMSIZE,IRMIND:IRMD,2)
!      +              ,CDENLM(IRMD,LMMAXD) ! lm-dos
!       DOUBLE PRECISION CLEB(NCLEB),RHO2NS(IRMD,LMPOTD),
!      +                 THETAS(IRID,NFUND)
!       INTEGER ICLEB(NCLEB,4),IFUNM(*),LMSP(*)
! C     ..
! C     .. Local Scalars ..
!       DOUBLE COMPLEX CLTDF,CONE,CZERO
!       DOUBLE PRECISION C0LL
!       INTEGER I,IFUN,IR,J,L1,LM1,LM2,LM3,M1
! C     ..
! C     .. Local Arrays ..
!       DOUBLE COMPLEX WR(LMSIZE,LMSIZE,IRMIND:IRMD)
! C     ..
! C     .. External Subroutines ..
!       EXTERNAL ZGEMM
! C     ..
! C     .. Intrinsic Functions ..
!       INTRINSIC ATAN,DIMAG,SQRT
! C     ..
! C     .. Save statement ..
!       SAVE
! C     ..
! C     .. Data statements ..
!       DATA CZERO/ (0.0D0,0.0D0)/
!       DATA CONE/ (1.0D0,0.0D0)/
! C     ..
! c
! C     C0LL = 1/sqrt(4*pi)
!       C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
! c
! c
! c---> initialize array for complex charge density
! c
! 
!       DO L1 = 0,LMAXD
!         DO I = 1,IRMD
!           CDEN(I,L1) = CZERO
!       ENDDO    
!       ENDDO    
! C------------------------------------------------------------------
! c
! c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
! c                                      summed over lm3
! c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
! c                                               summed over lm3
!       DO 50 IR = IRMIND + 1,IRMD
!       DO LM1=1,LMSIZE
!       DO LM2=1,LMSIZE
!       QNSI(LM1,LM2)=QNS(LM1,LM2,IR,1)
!       ENDDO
!       ENDDO
!         CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,PNS(1,1,IR,1),
!      +             LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)
!         CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,PNS(1,1,IR,1),
!      +             LMSIZE,QNSI,LMSIZE,CZERO,WR(1,1,IR),LMSIZE)
!         IF (NSRA.EQ.2) THEN
!       DO LM1=1,LMSIZE
!       DO LM2=1,LMSIZE
!       QNSI(LM1,LM2)=QNS(LM1,LM2,IR,2)
!       ENDDO
!       ENDDO
!           CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,PNS(1,1,IR,2),
!      +               LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)
!           CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,PNS(1,1,IR,2),
!      +               LMSIZE,QNSI,LMSIZE,CONE,WR(1,1,IR),LMSIZE)
!         END IF
! 
!         DO LM1 = 1,LMSIZE
!           DO LM2 = 1,LM1 - 1
!             WR(LM1,LM2,IR) = WR(LM1,LM2,IR) + WR(LM2,LM1,IR)
!       ENDDO
!       ENDDO
!  50   CONTINUe
! c
! c---> first calculate only the spherically symmetric contribution
! c
!       DO 100 L1 = 0,LMAXD
!         DO 70 M1 = -L1,L1
!           LM1 = L1* (L1+1) + M1 + 1
!           DO 60 IR = IRMIND + 1,IRMD
! c
! c---> fill array for complex density of states
! c
!             CDEN(IR,L1) = CDEN(IR,L1) + WR(LM1,LM1,IR)
!             CDENLM(IR,LM1) = WR(LM1,LM1,IR) ! lm-dos
!    60     CONTINUE
!    70   CONTINUE
! c
! c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
! c
!         DO 80 IR = IRMIND + 1,IRMD
!           RHO2NS(IR,1) = RHO2NS(IR,1) + C0LL*DIMAG(CDEN(IR,L1)*DF)
!    80   CONTINUE
! c
! c
! c
! !         IF (IPAN1.GT.1) THEN
! !           DO 90 I = IMT1 + 1,IRMD
! !             CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1)*C0LL
! ! 
! !             DO M1 = -L1,L1                                         ! lm-dos
! !                LM1 = L1* (L1+1) + M1 + 1                           ! lm-dos
! !                CDENLM(I,LM1) = CDENLM(I,LM1)*THETAS(I-IMT1,1)*C0LL ! lm-dos
! !             ENDDO                                                  ! lm-dos
! ! 
! !    90     CONTINUE
! !         END IF
! 
!   100 CONTINUE
! c
! !       IF (IPAN1.GT.1) THEN
! !         DO 110 I = 1,IRMD
! !           CDENNS(I) = 0.0D0
! !   110   CONTINUE
! !       END IF
! 
!       DO 140 J = 1,IEND
!         LM1 = ICLEB(J,1)
!         LM2 = ICLEB(J,2)
!         LM3 = ICLEB(J,3)
!         CLTDF = DF*CLEB(J)
! c
! c---> calculate the non spherically symmetric contribution
! c
!         DO 120 IR = IRMIND + 1,IRMD
!           RHO2NS(IR,LM3) = RHO2NS(IR,LM3) + DIMAG(CLTDF*WR(LM1,LM2,IR))
!   120   CONTINUE
! c
! !         IF (IPAN1.GT.1 .AND. LMSP(LM3).GT.0) THEN
! ! c       IF (IPAN1.GT.1) THEN
! !           IFUN = IFUNM(LM3)
! !           DO 130 I = IMT1 + 1,IRMD
! !             CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
! !      +                  THETAS(I-IMT1,IFUN)
! !   130     CONTINUE
! 
! !         END IF
! 
!   140 CONTINUE
! !          do lm1=1,(2*LMAXD+1)**2
! !            write(9013,'(5000F)') rho2ns(:,lm1)
! !          end do 
! 
!       END SUBROUTINE
! 
! 
!       SUBROUTINE RHOOUTNEW_notUSED(DF,GMAT,EK,PNS,QNS,
!      +                    NSRA,
!      +                   LMAXD,LMMAXD,LMPOTD,IRMD,IRMIND,GAUNTCOEFF)
! 
! !       CALL RHOOUTNEW(density%den(:,ispin,ie),DF,GMATll,EK,RLL(:,:,:,1),SLL(:,:,:,1),density%RHO2NS(:,:,ispin),shapefun%thetas,shapefun%lmused,cell%npan, &
! !                     1,shapefun%lm2index,CDENNS,config%NSRA,gauntcoeff%CLEB,gauntcoeff%ICLEB,gauntcoeff%IEND &
! !                    ,CDENLM, & ! lm-dos
! !                     gauntcoeff%NCLEB,LMAXD,LMMAXATOM,(2*LMAXD+1)**2,cellnew%nrmaxnew,1,IRID, &
! !                     shapefun%nlmshaped)
! c-----------------------------------------------------------------------
! c
! c     calculates the charge density from r(irmin) to r(irc)
! c      in case of a non spherical input potential .
! c
! c     fills the array cden for the complex density of states
! c
! c     attention : the gaunt coeffients are stored in index array
! c                   (see subroutine gaunt)
! c
! c     the structured part of the greens-function (gmat) is symmetric in
! c       its lm-indices , therefore only one half of the matrix is
! c       calculated in the subroutine for the back-symmetrisation .
! c       the gaunt coeffients are symmetric too (since the are calculated
! c       using the real spherical harmonics) . that is why the lm2- and
! c       the lm02- loops are only only going up to lm1 or lm01 and the
! c       summands are multiplied by a factor of 2 in the case of lm1 .ne.
! c       lm2 or lm01 .ne. lm02 .
! c
! c             (see notes by b.drittler)
! c
! c                               b.drittler   aug. 1988
! c-----------------------------------------------------------------------
! C     .. Parameters ..
! !       INCLUDE 'inc.p'
! C
! C *********************************************************************
! C * For KREL = 1 (relativistic mode)                                  *
! C *                                                                   *
! C *  NPOTD = 2 * NATYPD                                               *
! C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
! C *  NSPIND = 1                                                       *
! C *                                                                   *
! C *********************************************************************
! C
!       use type_gauntcoeff
!       IMPLICIT NONE
!       TYPE(GAUNTCOEFF_TYPE) :: GAUNTCOEFF
!       INTEGER LMAXD
!       INTEGER LMMAXD
! !       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
!       INTEGER LMPOTD
! !       PARAMETER (LMPOTD= (LPOTD+1)**2)
! !       INTEGER NCLEB
!       INTEGER IRMD
!       INTEGER IRMIND
! !       INTEGER IRID
!       INTEGER NFUND
! !       PARAMETER (IRMIND=IRMD-IRNSD)
! C     ..
! C     .. Scalar Arguments ..
!       DOUBLE COMPLEX DF,EK
!       INTEGER IEND,IMT1,IPAN1,NSRA
! C     ..
! C     .. Array Arguments ..
!       DOUBLE COMPLEX CDEN(IRMD,0:LMAXD),CDENNS(IRMD),
!      +               GMAT(LMMAXD,LMMAXD),
!      +               PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
!      +              QNSI(LMMAXD,LMMAXD),
!      +               QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
!      +              ,CDENLM(IRMD,LMMAXD) ! lm-dos
!       DOUBLE PRECISION RHO2NS(IRMD,LMPOTD)
! !      +                 THETAS(IRID,NFUND)
! !       INTEGER ICLEB(NCLEB,4),IFUNM(*),LMSP(*)
! C     ..
! C     .. Local Scalars ..
!       DOUBLE COMPLEX CLTDF,CONE,CZERO
!       DOUBLE PRECISION C0LL
!       INTEGER I,IFUN,IR,J,L1,LM1,LM2,LM3,M1
! C     ..
! C     .. Local Arrays ..
!       DOUBLE COMPLEX WR(LMMAXD,LMMAXD,IRMIND:IRMD)
! C     ..
! C     .. External Subroutines ..
!       EXTERNAL ZGEMM
! C     ..
! C     .. Intrinsic Functions ..
!       INTRINSIC ATAN,DIMAG,SQRT
! C     ..
! C     .. Save statement ..
!       SAVE
! C     ..
! C     .. Data statements ..
!       DATA CZERO/ (0.0D0,0.0D0)/
!       DATA CONE/ (1.0D0,0.0D0)/
! C     ..
! c
! C     C0LL = 1/sqrt(4*pi)
!       C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
!       RHO2NS=0.0
! c
! c
! c---> initialize array for complex charge density
! c
! 
!       DO L1 = 0,LMAXD
!         DO I = 1,IRMD
!           CDEN(I,L1) = CZERO
!       ENDDO    
!       ENDDO    
! C------------------------------------------------------------------
! c
! c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
! c                                      summed over lm3
! c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
! c                                               summed over lm3
!       DO 50 IR = IRMIND + 1,IRMD
!       DO LM1=1,LMMAXD
!       DO LM2=1,LMMAXD
!       QNSI(LM1,LM2)=QNS(LM1,LM2,IR,1)
!       ENDDO
!       ENDDO
!         CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1),
!      +             LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
!         CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,1),
!      +             LMMAXD,QNSI,LMMAXD,CZERO,WR(1,1,IR),LMMAXD)
!         IF (NSRA.EQ.2) THEN
!       DO LM1=1,LMMAXD
!       DO LM2=1,LMMAXD
!       QNSI(LM1,LM2)=QNS(LM1,LM2,IR,2)
!       ENDDO
!       ENDDO
!           CALL ZGEMM('N','N',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2),
!      +               LMMAXD,GMAT,LMMAXD,EK,QNSI,LMMAXD)
!           CALL ZGEMM('N','T',LMMAXD,LMMAXD,LMMAXD,CONE,PNS(1,1,IR,2),
!      +               LMMAXD,QNSI,LMMAXD,CONE,WR(1,1,IR),LMMAXD)
!         END IF
! 
!         DO LM1 = 1,LMMAXD
!           DO LM2 = 1,LM1 - 1
!             WR(LM1,LM2,IR) = WR(LM1,LM2,IR) + WR(LM2,LM1,IR)
!       ENDDO
!       ENDDO
!  50   CONTINUe
! c
! c---> first calculate only the spherically symmetric contribution
! c
!       DO 100 L1 = 0,LMAXD
!         DO 70 M1 = -L1,L1
!           LM1 = L1* (L1+1) + M1 + 1
!           DO 60 IR = IRMIND + 1,IRMD
! c
! c---> fill array for complex density of states
! c
!             CDEN(IR,L1) = CDEN(IR,L1) + WR(LM1,LM1,IR)
!             CDENLM(IR,LM1) = WR(LM1,LM1,IR) ! lm-dos
!    60     CONTINUE
!    70   CONTINUE
! c
! c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
! c
!         DO 80 IR = IRMIND + 1,IRMD
! !           write(*,*) CDEN(IR,L1)
!           RHO2NS(IR,1) = RHO2NS(IR,1) + C0LL*DIMAG(CDEN(IR,L1)*DF)
!    80   CONTINUE
! !           stop
! c
! c
! c
! !         IF (IPAN1.GT.1) THEN
! !           DO 90 I = IMT1 + 1,IRMD
! !             CDEN(I,L1) = CDEN(I,L1)*THETAS(I-IMT1,1)*C0LL
! ! 
! !             DO M1 = -L1,L1                                         ! lm-dos
! !                LM1 = L1* (L1+1) + M1 + 1                           ! lm-dos
! !                CDENLM(I,LM1) = CDENLM(I,LM1)*THETAS(I-IMT1,1)*C0LL ! lm-dos
! !             ENDDO                                                  ! lm-dos
! ! 
! !    90     CONTINUE
! !         END IF
! 
!   100 CONTINUE
! c
! !       IF (IPAN1.GT.1) THEN
! !         DO 110 I = 1,IRMD
! !           CDENNS(I) = 0.0D0
! !   110   CONTINUE
! !       END IF
! 
!       DO 140 J = 1,GAUNTCOEFF%IEND
!         LM1 = GAUNTCOEFF%ICLEB(J,1)
!         LM2 = GAUNTCOEFF%ICLEB(J,2)
!         LM3 = GAUNTCOEFF%ICLEB(J,3)
!         CLTDF = DF*GAUNTCOEFF%CLEB(J,1)
! c
! c---> calculate the non spherically symmetric contribution
! c
!         DO 120 IR = IRMIND + 1,IRMD
!           RHO2NS(IR,LM3) = RHO2NS(IR,LM3) + DIMAG(CLTDF*WR(LM1,LM2,IR))
!   120   CONTINUE
! c
! !         IF (IPAN1.GT.1 .AND. LMSP(LM3).GT.0) THEN
! ! c       IF (IPAN1.GT.1) THEN
! !           IFUN = IFUNM(LM3)
! !           DO 130 I = IMT1 + 1,IRMD
! !             CDENNS(I) = CDENNS(I) + CLEB(J)*WR(LM1,LM2,I)*
! !      +                  THETAS(I-IMT1,IFUN)
! !   130     CONTINUE
! 
! !         END IF
! 
!   140 CONTINUE
! 
! !          do lm1=1,(2*LMAXD+1)**2
! !            write(9013,'(5000F)') rho2ns(:,lm1)
! !          end do 
! 
!       END SUBROUTINE





      END MODULE MOD_RHOOUTNEW
