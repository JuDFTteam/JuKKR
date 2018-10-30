      MODULE mod_vllmat
      CONTAINS
      SUBROUTINE VLLMAT(VNSPLL,VINS,LMAXATOM,LMMAXATOM,LMPOTATOM,
     +                 IRMIND,IRMD,GAUNTCOEFF,ZATOM,RMESH,
     +                 VLLLMMAX,use_fullgmat,NSPIN,JSPIN,CMODE )
      USE TYPE_GAUNTCOEFF
      IMPLICIT NONE
c-----------------------------------------------------------------------
c     to determine the non - spherical wavefunctions the potential
c         has to be lm1 and lm2 dependent . the potential is stored
c         only as lm dependent , therefore a transformation in the
c         following way has to be done :
c
c        vnsll(r,lm1,lm2)   =   {  c(lm1,lm2,lm3) *vins(r,lm3)  }
c                                  (summed over lm3 at the right site )
c        where c(lm1,lm2,lm3) are the gaunt coeffients .
c
c             (see notes by b.drittler)
c
c     attention : the gaunt coeffients are stored in an index array
c                  only for lm1.gt.lm2
c                 (see subroutine gaunt)
c
c                               b.drittler   july 1988
c-----------------------------------------------------------------------
c                          modified by R. Zeller Sep. 2000
c-----------------------------------------------------------------------
C     .. Parameters ..
!       INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXATOM = 2 * (LMAXATOM+1)^2                                         *
C *  NSPIND = 1                                                       *
C *                                                                   *
C *********************************************************************
C
      INTEGER LMAXATOM,LMMAXATOM,LMPOTATOM,IRMIND,IRMD 
      TYPE(GAUNTCOEFF_TYPE) GAUNTCOEFF
      INTEGER VLLLMMAX
      CHARACTER(LEN=*) CMODE
      DOUBLE PRECISION VINS(IRMD,LMPOTATOM,NSPIN),
     +                 RMESH(IRMD),ZATOM
      DOUBLE COMPLEX   VNSPLL(VLLLMMAX,VLLLMMAX,IRMD)
C     .. Local Scalars ..
      INTEGER IR,J,LM1,LM2,LM3,NSPIN
      INTEGER use_fullgmat,LMMAXSHIFT,NSPINORBIT,ISPIN,ISPIN2,jspin

      IF     (CMODE=='NS') THEN
      ELSEIF (CMODE=='SPH') THEN
      ELSEIF (CMODE=='NOSPH') THEN
      ELSE
        STOP '[VLLMAT_NEW] error in cmode'
      END IF


!************************************************************************************
! integrity check
!************************************************************************************
      IF (use_fullgmat==1 .and. 
     +    VLLLMMAX/=2*LMMAXATOM) stop 'Error in vllmat:VLLLMMAX'

!************************************************************************************
! In case of a SPINORBIT calculation use a twice as big VLL matrix:
!                 ( Vdndn   Vdnup )
!      VLL =      (               )
!                 ( Vupdn   Vupup )
!************************************************************************************
      VNSPLL(:,:,:) = (0.0D0,0.0D0)
      IF (use_fullgmat==1) THEN
        NSPINORBIT=2
      ELSE
        NSPINORBIT=1
      END IF

      DO ISPIN=1,NSPINORBIT
        IF (ISPIN==1) THEN 
          LMMAXSHIFT=0
        ELSE
          LMMAXSHIFT=LMMAXATOM
        END IF
!************************************************************************************
! If Spinorbit is used for NSPIN=1 calculation use the same potential for
! both spin blocks
!************************************************************************************
        ISPIN2=ISPIN
        IF (NSPINORBIT==1 .and. JSPIN==2) ISPIN2=2
        IF (NSPINORBIT==2 .and. JSPIN==2) stop 'error in vllmatnew'

!         write(*,*) '***************',ispin,nspin,'***********'

        IF (.not. CMODE=='SPH') THEN
          DO J = 1,GAUNTCOEFF%IEND
            LM1 = GAUNTCOEFF%ICLEB(J,1)
            LM2 = GAUNTCOEFF%ICLEB(J,2)
            LM3 = GAUNTCOEFF%ICLEB(J,3)
            IF (LM1 <= LMMAXATOM .AND. LM2 <= LMMAXATOM ) THEN
              DO IR = IRMIND,IRMD
                VNSPLL(LM1+LMMAXSHIFT,LM2+LMMAXSHIFT,IR) = 
     +                     VNSPLL(LM1+LMMAXSHIFT,LM2+LMMAXSHIFT,IR)
     +                   + GAUNTCOEFF%CLEB(J,1)*VINS(IR,LM3,ISPIN2)
              END DO
            END IF
          END DO
  !
  !---> use symmetry of the gaunt coef.
  !
          DO LM1 = 1,LMMAXATOM
            DO LM2 = 1,LM1 - 1
                DO IR = IRMIND,IRMD
                  VNSPLL(LMMAXSHIFT+LM2,LMMAXSHIFT+LM1,IR) =
     +            VNSPLL(LMMAXSHIFT+LM1,LMMAXSHIFT+LM2,IR)
                END DO
            END DO
          END DO
        END IF


        IF (.not. CMODE=='NOSPH') THEN
          DO LM1 = 1,LMMAXATOM
            DO IR = 1,IRMD
              VNSPLL(LMMAXSHIFT+LM1,LMMAXSHIFT+LM1,IR) = 
     +         VNSPLL(LMMAXSHIFT+LM1,LMMAXSHIFT+LM1,IR) 
     +         + (VINS(IR,1,ISPIN2)-2.0D0*ZATOM/rmesh(ir))
            END DO
          END DO
        END IF !(CMODE=='NS')
      END DO !NSPINORBIT

      END SUBROUTINE
      END MODULE mod_vllmat
