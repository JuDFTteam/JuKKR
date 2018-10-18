!-------------------------------------------------------------------------------
!> Summary: Convolutes potentials with shape functions
!> Author:  
!> 
!> Calculate convolution of potential with shapefunction.
!-------------------------------------------------------------------------------
      MODULE MOD_CONVOL
      CONTAINS
!-------------------------------------------------------------------------------
!> Summary: Driver routine for the convolution module
!> Author:  
!> Date: 
!> Category: KKRimp, shape-functions, potential 
!> Deprecated: False 
!>
!> Driver routine for the convolution of potential with shapefunction.
!-------------------------------------------------------------------------------
      SUBROUTINE CONVOLDRV(nrmaxd,LMAXD,NSPIN,NATOM,LMAXATOM,cell,
     +                     gauntshape,shapefun,ZAT, vons)
      use type_cell
      use type_gauntshape
      use type_shapefun
      IMPLICIT NONE
      INTEGER :: NRMAXD
      INTEGER :: LMAXD
      INTEGER :: NSPIN
      INTEGER :: NATOM
      INTEGER :: LMAXATOM(NATOM)
      TYPE(CELL_TYPE)      :: cell(NATOM)
      TYPE(gauntshape_type)      :: gauntshape(lmaxd)
      TYPE(shapefun_type)      :: shapefun(NATOM)
      REAL(KIND=8)            :: ZAT(NATOM)
      real*8                ::  vons(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),

      REAL(KIND=8),PARAMETER  :: RFPI=3.5449077018110318
!       real(kind=dp),allocatable        ::  THESME(:,:)              ! shape function
!local
      INTEGER              :: ISPIN,IATOM,LMPOT,LMAX
!       INTEGER              :: NTCELL(NATOM)
!       allocate(ntcell(natom))
!        ntcell=1
!           ntcell=1
    
          DO ISPIN = 1,NSPIN
              DO IATOM = 1,NATOM
                LMAX  = LMAXATOM(IATOM)
                LMPOT = (2*LMAXATOM(IATOM)+1)**2
!                 write(*,*) ispin,iatom,natom
!                 write(*,*) '*******************************'
!                 write(*,*) CELL(IATOM)%NRCUT(1)
!                 write(*,*) CELL(IATOM)%NRCUT(CELL(IATOM)%NPAN)
! !                 write(*,*) NTCELL(IATOM)
!                 write(*,*) gauntshape%IMAXSH(LMPOT)
!                 write(*,*) lbound(gauntshape%ILM),ubound(gauntshape%ILM)
!                 write(*,*) ''
!                 write(*,*) '*******************************'
! 
! 
!                 write(101010+iatom,* ) CELL(IATOM)%NRCUT(1),
!      &                   CELL(IATOM)%NRCUT(CELL(IATOM)%NPAN),
!      &                   gauntshape%IMAXSH(LMPOT),
!      &                   gauntshape%ILM,
!      &                   shapefun(iatom)%lm2index,
!      &                   (2*lmaxd+1)**2,
! !      &                   'GSH',
!      &                   gauntshape%GSH,
! !      &                   'THETAS',
!      &                   shapefun(iatom)%THETAS,
!      &                   ZAT(IATOM),RFPI,
!      &                   CELL(IATOM)%RMESH(:),
!      &                   VONS(:,:,ISPIN,IATOM),
!      &                   shapefun(iatom)%lmused,
!      &                   shapefun(iatom)%NRSHAPE,
!      &                   shapefun(iatom)%NLMSHAPE,
!      &                   NRMAXD,
!      &                   GAUNTSHAPE%NGSHD
!                 write(*,*) '*******************************'


            CALL CONVOL( CELL(IATOM)%NRCUT(1),
     &                   CELL(IATOM)%NRCUT(CELL(IATOM)%NPAN),
     &                   gauntshape(LMAX)%IMAXSH(LMPOT),
     &                   gauntshape(LMAX)%ILM,
     &                   shapefun(iatom)%lm2index,
     &                   (2*lmaxd+1)**2,
     &                   gauntshape(LMAX)%GSH,
     &                   shapefun(iatom)%THETAS,
     &                   ZAT(IATOM),RFPI,
     &                   CELL(IATOM)%RMESH(:),
     &                   VONS(:,:,ISPIN,IATOM),
     &                   shapefun(iatom)%lmused,
     &                   shapefun(iatom)%NRSHAPED,
     &                   shapefun(iatom)%NLMSHAPE,
     &                   NRMAXD,
     &                   GAUNTSHAPE(LMAX)%NGSHD)
    !              CALL CONVOL(IRCUT(1,IATOM),IRC(IATOM),NTCELL(IATOM),
    !      &                   IMAXSH(LMPOT),ILM,IFUNM,LMPOT,GSH,
    !      &                   THETAS,THESME,ZAT(IATOM),RFPI,
    !      &                   R(1,IATOM),VONS(1,1,IPOT),VSPSMDUM(1,1),LMSP)
              END DO
          END DO

!         stop

      END SUBROUTINE CONVOLDRV



!-------------------------------------------------------------------------------
!> Summary: Convolutes potentials with shape functions 
!> Author:  
!> Date: 
!> Category: KKRimp, shape-fucntions, potential
!> Deprecated: False
!> 
!> Calculate convolution of potential with shapefunction.
!-------------------------------------------------------------------------------
      SUBROUTINE CONVOL(IMT1,IRC1,IMAXSH,ILM,IFUNM,LMPOT,GSH,
     +                  THETAS,Z,RFPI,R,VONS,LMSP,IRID,NFUND,IRMD,NGSHD)
c ************************************************************************
C     .. Parameters ..
!       include 'inc.p'
!       INTEGER LMPOTD
!       PARAMETER (LMPOTD= (LPOTD+1)**2)
C     ..
C     .. Scalar Arguments ..
      INTEGER IRID,NFUND,IRMD,NGSHD
      DOUBLE PRECISION RFPI,Z
      INTEGER IMAXSH,IMT1,IRC1,LMPOT
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION GSH(*),R(IRMD),THETAS(IRID,NFUND),
     +                 VONS(:,:)

      INTEGER IFUNM(:),ILM(NGSHD,3),LMSP(:)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ZZOR
      INTEGER I,IFUN,IR,IRH,LM,LM1,LM2,LM3
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION VSTORE(IRID,LMPOT) !,VSTSME(IRID,LMPOTD)
C     ..
      DO 20 LM = 1,LMPOT
        DO 10 IR = 1,IRC1 - IMT1
          VSTORE(IR,LM) = 0.0D0
!           VSTSME(IR,LM) = 0.0D0
   10   CONTINUE
   20 CONTINUE

      DO 30 IR = IMT1 + 1,IRC1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VONS(IR,1) - ZZOR
   30 CONTINUE

      DO 50 I = 1,IMAXSH
        LM1 = ILM(I,1)
        LM2 = ILM(I,2)
        LM3 = ILM(I,3)
        IF (LMSP(LM3).GT.0) THEN    !LMSP(ICELL,LM3).GT.0
          IFUN = IFUNM(LM3) !IFUNM(ICELL,LM3)
          DO 40 IR = IMT1 + 1,IRC1
            IRH = IR - IMT1
            VSTORE(IRH,LM1) = VSTORE(IRH,LM1) +
     +                        GSH(I)*VONS(IR,LM2)*THETAS(IRH,IFUN) !THETAS(IRH,IFUN,ICELL)
!             VSTSME(IRH,LM1) = VSTSME(IRH,LM1) +
!      +           GSH(I)*VONS(IR,LM2)*THESME(IRH,IFUN,ICELL)
   40     CONTINUE
        END IF
   50 CONTINUE

      DO 60 IR = IMT1 + 1,IRC1
        IRH = IR - IMT1
        ZZOR = 2.0D0*Z/R(IR)*RFPI
        VONS(IR,1) = VSTORE(IRH,1) + ZZOR
!         VSPSMO(IR) = (VSTSME(IRH,1) + ZZOR) /RFPI
   60 CONTINUE

C     COPY THE PART INSIDE THE MT SPHERE
!       DO IR = 1,IMT1
!         VSPSMO(IR) = VONS(IR,1)/RFPI
!       END DO

      DO 80 LM = 2,LMPOT
        DO 70 IR = IMT1 + 1,IRC1
          IRH = IR - IMT1
          VONS(IR,LM) = VSTORE(IRH,LM)
   70   CONTINUE
   80 CONTINUE

      RETURN

      END SUBROUTINE
      END MODULE MOD_CONVOL
