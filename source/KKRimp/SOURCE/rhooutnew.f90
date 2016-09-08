      MODULE MOD_RHOOUTNEW
       CONTAINS




      SUBROUTINE RHOOUTNEW(gauntcoeff,DF,GMATIN,EK,cellnew,wavefunction,RHO2NSC, &
                          NSRA, &
                      LMAXD,LMAXATOM,LMMAXATOM,LMSIZE,LMSIZE2,LMPOTD,IRMD,&
                         ISPIN,NSPINDEN,imt1,cden,cdenlm,cdenns,shapefun,corbital, &
                        gflle_part )                                              ! lda+u

use type_gauntcoeff
use type_cellnew
use type_wavefunction
USE mod_physic_params,only: cvlight
use mod_basistransform
use mod_mathtools, only: transpose
use mod_config, only: config_testflag
use type_shapefun
use mod_orbitalmoment, only: calc_orbitalmoment
use mod_chebyshev   ! lda+u

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
      INTEGER LMAXD,LMAXATOM
      INTEGER LMMAXATOM
      INTEGER LMSIZE,LMSIZE2
!       parameter (lmmaxd= (krel+1) * (lmaxd+1)**2)
      INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER IRMD
! C     ..
! C     .. Scalar Arguments ..
      DOUBLE COMPLEX DF,EK
      INTEGER IMT1,NSRA
      type(shapefun_type),intent(in)            :: shapefun
      integer        :: corbital

! C     ..
! C     .. Array Arguments ..
      DOUBLE COMPLEX :: CDEN(IRMD,0:LMAXD,NSPINDEN)
      DOUBLE COMPLEX CDENNS(IRMD,NSPINDEN), &
                     GMAT(LMSIZE,LMSIZE), &
                     GMATIN(LMSIZE,LMSIZE), &
!                      PNS(LMSIZE,LMSIZE,IRMD,2), &
                    QNSI(LMSIZE,LMSIZE),RLLTEMP(LMSIZE,LMSIZE), &
!                      QNS(LMSIZE,LMSIZE,IRMD,2) &
                    CDENLM(IRMD,LMMAXATOM,NSPINDEN), & ! lm-dos
                    RHO2NSC(IRMD,LMPOTD,NSPINDEN), &
                    gflle_part(lmsize,lmsize) ! lda+u

!       DOUBLE PRECISION  &
!                        THETAS(IRID,NFUND)
!       INTEGER IFUNM(*),LMSP(*)
! C     ..
! C     .. Local Scalars ..
      DOUBLE COMPLEX CLTDF,CONE,CZERO
      DOUBLE PRECISION C0LL
      INTEGER IFUN,IR,J,L1,LM1,LM2,LM3,M1
! C     ..
! C     .. Local Arrays ..
!       DOUBLE COMPLEX WR(LMSIZE,LMSIZE,IRMD)
      DOUBLE COMPLEX,allocatable ::  WR(:,:,:), &
                                     cwr(:) ! lda+u

      INTEGER          :: ISPIN,JSPIN

      INTEGER     :: SPININDEX1(4) !=(/1,2,1,2 /)
      INTEGER     :: SPININDEX2(4) !=(/1,2,2,1 /)
      INTEGER     :: LMSHIFT1(4)
      INTEGER     :: LMSHIFT2(4)
      INTEGER     :: NSPINSTART,NSPINSTOP,NSPINDEN

      DOUBLE COMPLEX :: Loperator(lmsize,lmsize,3)

! C     ..
! C     .. External Subroutines ..
      EXTERNAL ZGEMM
! C     ..
! C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DIMAG,SQRT
! C     ..
! C     .. Save statement ..
!       SAVE
! C     ..
! C     .. Data statements ..
      DATA CZERO/ (0.0D0,0.0D0)/
      DATA CONE/ (1.0D0,0.0D0)/
! C     ..
! c

      ALLOCATE ( WR(LMSIZE,LMSIZE,IRMD), &
                 cwr(irmd) )                    ! lda+u

      GMAT=GMATIN

      IF (NSPINDEN==4) THEN
        NSPINSTART=1
        NSPINSTOP=NSPINDEN
        SPININDEX1   =(/1,2,1,2 /)
        SPININDEX2   =(/1,2,2,1 /)
        LMSHIFT1=LMMAXATOM*(SPININDEX1-1) 
        LMSHIFT2=LMMAXATOM*(SPININDEX2-1) 
      ELSE
        NSPINSTART=ISPIN
        NSPINSTOP=ISPIN
        SPININDEX1   =(/1,1,0,0 /)
        SPININDEX2   =(/1,1,0,0 /)
        LMSHIFT1=LMMAXATOM*(SPININDEX1-1)
        LMSHIFT2=LMMAXATOM*(SPININDEX2-1) 
     END IF
     
     if (corbital/=0) then
        CALL calc_orbitalmoment(lmaxatom,Loperator)
     end if

!    write(*,*) 'lmsize',lmsize,nspinstart,nspinstop


! C     C0LL = 1/sqrt(4*pi)
     C0LL = 1.0d0/SQRT(16.0D0*ATAN(1.0D0))
! c
! c
! c---> initialize array for complex charge density
! c


     DO JSPIN=NSPINSTART,NSPINSTOP
        CDEN(:,:,JSPIN) = CZERO
        CDENLM(:,:,JSPIN) = CZERO
     END DO !JSPIN

! C------------------------------------------------------------------
! c
! c---> set up array ek*qns(lm1,lm2) + { gmat(lm3,lm2)*pns(lm1,lm3) }
! c                                      summed over lm3
! c---> set up of wr(lm1,lm2) = { pns(lm1,lm3)*qns(lm2,lm3) }
! c                                               summed over lm3


     DO    IR = 1,IRMD


! ########################################################3
!
! WATCH OUT CHECK IF A FACTOR OF M_0 needs to be added into the Greensfunction
!
! #########################################################3
        IF (allocated(WAVEFUNCTION%SLLleft)) then
           QNSI(:,:)=WAVEFUNCTION%SLLleft(1:LMSIZE,1:LMSIZE,IR,1)
        else
           QNSI(:,:)=WAVEFUNCTION%SLL(1:LMSIZE,1:LMSIZE,IR,1)
        end if

        IF (allocated(WAVEFUNCTION%SLLleft)) then
           RLLTEMP=WAVEFUNCTION%RLLleft(1:LMSIZE,1:LMSIZE,IR,1)
        else
           RLLTEMP=WAVEFUNCTION%RLL(1:LMSIZE,1:LMSIZE,IR,1)
        end if

        ! changed the second mode to transpose - bauer
        CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,RLLTEMP, &
                  LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)

        RLLTEMP(:,:)=WAVEFUNCTION%RLL(1:LMSIZE,1:LMSIZE,IR,1)

        CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,RLLTEMP, &
                    LMSIZE,QNSI,LMSIZE,CZERO,WR(1,1,IR),LMSIZE)



        IF (NSRA.EQ.2 .and. (.not. config_testflag('nosmallcomp')) ) THEN

           if (allocated(WAVEFUNCTION%SLLleft)) then
              QNSI(:,:)= - WAVEFUNCTION%SLLleft(LMSIZE+1:2*LMSIZE,1:LMSIZE,IR,1) ! Attention to the
                                                                                 ! additional minus sign
              ! ##########################################################################################
              ! Drittler assumes that for the left solution, is given by the right solution with an
              ! additional minus sign. This minus sign is contained inside the equations to calculate
              ! the electronic density. While calculating the left solution, the minus sign is already 
              ! included in the left solution. To make calculations consistant a factor of -1 is included
              ! which cancels out by the routines of Drittler
              ! ##########################################################################################
           else
              QNSI(:,:)=WAVEFUNCTION%SLL(LMSIZE+1:2*LMSIZE,1:LMSIZE,IR,1)
           end if
           if (allocated(WAVEFUNCTION%RLLleft)) then
              RLLTEMP= - WAVEFUNCTION%RLLleft(LMSIZE+1:2*LMSIZE,1:LMSIZE,IR,1) ! Attention to the
                                                                               ! additional minus sign
              ! ##########################################################################################
              ! Drittler assumes that for the left solution, is given by the right solution with an
              ! additional minus sign. This minus sign is contained inside the equations to calculate
              ! the electronic density. While calculating the left solution, the minus sign is already 
              ! included in the left solution. To make calculations consistant a factor of -1 is included
              ! which cancels out by the routines of Drittler
              ! ##########################################################################################
           else
              RLLTEMP=WAVEFUNCTION%RLL(LMSIZE+1:2*LMSIZE,1:LMSIZE,IR,1)
           end if





!          changed the second mode to transpose - bauer
           CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,RLLTEMP, &
                      LMSIZE,GMAT,LMSIZE,EK,QNSI,LMSIZE)

           RLLTEMP=WAVEFUNCTION%RLL(LMSIZE+1:2*LMSIZE,1:LMSIZE,IR,1)!/CVLIGHT

           CALL ZGEMM('N','T',LMSIZE,LMSIZE,LMSIZE,CONE,RLLTEMP, &
                       LMSIZE,QNSI,LMSIZE,CONE,WR(1,1,IR),LMSIZE)
        END IF


        if (corbital/=0) then

           CALL ZGEMM('N','N',LMSIZE,LMSIZE,LMSIZE,CONE,Loperator(:,:,corbital), &
                      LMSIZE,WR(:,:,IR),LMSIZE,CZERO,RLLTEMP,LMSIZE)
           WR(:,:,IR)=RLLTEMP
        end if


! Phivos lda+u: Place here r-integration of wr(lms1,lms2,ir) to obtain cnll(lms1,lms2).    ! lda+u
! Integrate only up to muffin-tin radius.                                                  ! lda+u
        gflle_part(:,:) = czero                                                            ! lda+u
        do lm2 = 1,lmsize                                                                  ! lda+u
           do lm1 = 1,lmsize                                                               ! lda+u
              cwr(1:imt1) = wr(lm1,lm2,1:imt1)                                             ! lda+u
              cwr(imt1+1:irmd) = czero                                                     ! lda+u
              call intcheb_cell(cwr,cellnew,gflle_part(lm1,lm2))                           ! lda+u
           enddo                                                                           ! lda+u
        enddo                                                                              ! lda+u


! Change by Phivos 12.6.2012: this part was within previous IR-loop, now moved here        ! lda+u
! in order to calculate the complex density just before this averaging.                    ! lda+u


        DO JSPIN=NSPINSTART,NSPINSTOP
           DO LM1 = 1,LMMAXATOM
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

     DO L1 = 0,LMAXATOM

        DO M1 = -L1,L1
           LM1 = L1* (L1+1) + M1 + 1
           DO IR = 1,IRMD
              ! c
              ! c---> fill array for complex density of states
              ! c
              DO JSPIN=NSPINSTART,NSPINSTOP
                 CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)  &
                      + WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR)
                 CDENLM(IR,LM1,JSPIN) = &
                      WR(LM1+LMSHIFT1(JSPIN),LM1+LMSHIFT2(JSPIN),IR) ! lm-dos
              END DO

           END DO !IR
        END DO !M1

        ! c
        ! c---> remember that the gaunt coeffients for that case are 1/sqrt(4 pi)
        ! c

        DO JSPIN=NSPINSTART,NSPINSTOP

           DO IR = 1,IRMD
              RHO2NSC(IR,1,JSPIN) = RHO2NSC(IR,1,JSPIN)  &
                                + C0LL*(CDEN(IR,L1,JSPIN)*DF)
           END DO

           DO IR = IMT1 + 1,IRMD
              CDEN(IR,L1,JSPIN) = CDEN(IR,L1,JSPIN)*cellnew%shapefun(IR,1)*C0LL

              DO M1 = -L1,L1                                         ! lm-dos
                 LM1 = L1* (L1+1) + M1 + 1                           ! lm-dos
                 CDENLM(IR,LM1,JSPIN) = CDENLM(IR,LM1,JSPIN)*cellnew%shapefun(IR,1)*C0LL ! lm-dos
              ENDDO                                                  ! lm-dos
           END DO

        END DO


     END DO ! L1 = 0,LMAXATOM



     DO JSPIN=NSPINSTART,NSPINSTOP
        CDENNS(:,JSPIN) = 0.0D0
     END DO




     DO J = 1,gauntcoeff%IEND
        LM1 = gauntcoeff%ICLEB(J,1)
        LM2 = gauntcoeff%ICLEB(J,2)
        LM3 = gauntcoeff%ICLEB(J,3)
        CLTDF = DF*gauntcoeff%CLEB(J,1)
! c
! c---> calculate the non spherically symmetric contribution
! c
        DO JSPIN=NSPINSTART,NSPINSTOP
           DO IR = 1,IRMD
              RHO2NSC(IR,LM3,JSPIN) = RHO2NSC(IR,LM3,JSPIN)  &
                        + (CLTDF*WR(LM1+LMSHIFT1(JSPIN), &
                                         LM2+LMSHIFT2(JSPIN),IR)) !&
           END DO
! c
           IF (shapefun%LMused(lm3)==1) THEN
              IFUN = shapefun%lm2index(lm3) !IFUNM(LM3)
              DO IR = IMT1 + 1,IRMD
                 CDENNS(IR,JSPIN) = CDENNS(IR,JSPIN) &
                      + gauntcoeff%CLEB(J,1)*WR(LM1+LMSHIFT1(JSPIN),LM2+LMSHIFT2(JSPIN),IR)* &
                      cellnew%shapefun(IR,IFUN)
              END DO
           END IF
        END DO


      END DO !J

      deallocate(wr,cwr)


      END SUBROUTINE

      END MODULE MOD_RHOOUTNEW


