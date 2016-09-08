MODULE MOD_VINTRAS
  CONTAINS

SUBROUTINE VINTRAS(NATOM, NSPIN, NRMAXD,LMAXD, LMAXATOM,CELL, VPOT_OUT, SHAPEFUN, GAUNTSHAPE, DENSITY,CMOM,CMOM_INTERST,ins)
 Use mod_sinwk
 Use mod_soutk
 use type_cell, only: CELL_TYPE
 use type_shapefun, only: shapefun_TYPE
 use type_gauntshape, only: gauntshape_TYPE
 use type_density, only: density_type
 Implicit None
!-----------------------------------------------------------------------
!     calculate the electron-intracell-potentials and the charge-
!     moments of given charge densities . ( for each spin-direc-
!     tion the potential is the same in the polarized case . )
!     initialize the potential v with the electron-intracell-potentials
!    the intracell-potential is expanded into spherical harmonics .
!     the lm-term of the intracell-potential of the representive atom i
!     is given by
!                    8pi        r      r'** l
!      v(r,lm,i) =  ----- *  (  s dr' --------   rho2ns(r',lm,i,1)
!                   2*l+1       0     r **(l+1)
!
!                                 rcut    r ** l
!                               +  s dr' ---------   rho2ns(r',lm,i,1) )
!                                  r     r' **(l+1)
!
!     the lm contribution of the charge moment of the representive
!     atom i is given by
!
!                             rcut
!              cmom(lm,i) =    s dr' r'** l rho2ns(r',lm,i,1)
!                              0
!
!             (see notes by b.drittler and u.klemradt)
!
!              rcut is muffin tin or wigner seitz sphere radius,
!              depending on kshape turned on or off
!
!     attention : rho2ns(...,1) is the real charge density times r**2
!                 developed into spherical harmonics . (see deck rholm)
!
!                               b.drittler   may 1987
!-----------------------------------------------------------------------
integer :: natom
integer :: nspin
integer :: nrmaxd
integer :: lmaxd
integer :: lmaxatom(natom)
real*8,allocatable  :: cmom(:,:)
real*8,allocatable  :: CMOM_INTERST(:,:)
real*8                ::  vpot_out(nrmaxd,(2*lmaxd+1)**2,nspin,natom) !thetas(iri,nfund,*),

type(cell_type)                     :: cell(natom)
type(shapefun_type)                 :: shapefun(natom)
type(gauntshape_type)               :: gauntshape(lmaxd)
type(density_type)                  :: density(natom)

REAL*8 FAC,PI,RL


INTEGER I,IATOM,IEND,IFUN,IRC1,IRS1,ISTART,J,LVAL,LM,LM2,LM3,MVAL,LMAX
REAL*8 V1(NRMaxD),V2(NRMaxD),VINT1(NRMaxD),VINT2(NRMaxD)
INTEGER ins
integer ipand
INTEGER,allocatable :: IRCUTM(:)
!       INTEGER IFUNM(NATYPD,*)

!       INTEGER ins,LMAX,NEND,NSPIN,NSTART, irm, iri
!C     ..
!C     .. Array Arguments ..
!       REAL*8 CMINST(LMPOTD,*),CMOM(LMPOTD,*)
! ,DRDI(IRM,*),
!      +                 GSH(*),R(IRM,*),RHO2NS(IRM,LMPOTD,NATYPD,*),
!     REAL*8                  V(IRM,LMPOTD,NSPIN) !THETAS(IRI,NFUND,*),
!       INTEGER IFUNM(NATYPD,*)
! ,ILM(NGSHD,3),IMAXSH(0:LMPOTD),IPAN(*),
!      +        IRCUT(0:IPAND,*),IRWS(*),LMSP(NATYPD,*),NTCELL(*)
!C     ..
!C     .. Local Scalars ..
!       REAL*8 FAC,PI,RL

!C     ..
!C     .. Local Arrays ..
!       REAL*8 V1(CELL(1)%NRMaxD),V2(CELL(1)%NRMaxD),VINT1(CELL(1)%NRMaxD),VINT2(CELL(1)%NRMaxD)
!C     ..
!C     .. External Subroutines ..
!       EXTERNAL SINWK,SOUTK
!C     ..
!C     .. Intrinsic Functions ..
!       INTRINSIC DATAN,REAL
!       SAVE
!C     ..

! ######################################################
! calculate the maximum number of panels
! ######################################################
ipand=0
do iatom=1,natom
  if (CELL(IATOM)%NPAN> ipand) ipand=CELL(IATOM)%NPAN
end do !natom

! write(*,*) 'ipand',ipand

! allocate( cmom((2*lmaxd+1)**2,natom)  )



allocate(IRCUTM(0:IPAND))

PI = 4.D0*DATAN(1.D0)

DO IATOM = 1,NATOM
!   write(*,*) 'ins',ins
  IF (ins.NE.0) THEN
    IRS1 = CELL(IATOM)%NRCUT(      1           ) !IRCUT(1,IATOM)
!     write(*,*) 'irs1',irs1
    IRC1 = CELL(IATOM)%NRCUT( CELL(IATOM)%NPAN ) !IRCUT(IPAN(IATOM),IATOM)
!     write(*,*) 'irc1',irc1
!     ICELL = IATOM
!cccc ERROR IN HERE
!cccc          DO 10 I = 0,IPAN(ICELL)
!cccc Next line is correct one (T.Korhonen, Nov 1997)
    IRCUTM=0
    DO I = 0, CELL(IATOM)%NPAN  !IPAN(IATOM)
      IRCUTM(I) = CELL(IATOM)%NRCUT(I) !IRCUT(I,IATOM)
    END DO
!     write(*,*) 'ircutm',ircutm
  ELSE
      IRS1 = CELL(IATOM)%NRMAX !IRWS(IATOM) ????
      IRC1 = IRS1
      IRCUTM(0) = CELL(IATOM)%NRCUT(0) !IRCUT(0,IATOM)
      IRCUTM(1) = IRC1
  END IF
!c---> determine the right potential numbers
!         IPOT = NSPIN*IATOM

  LMAX=LMAXATOM(IATOM)
  DO LVAL = 0,2*LMAXATOM(IATOM)
    V1=0.0D0
    V2=0.0D0
    VINT1=0.0D0
    VINT2=0.0D0

    FAC = 8.0D0*PI/REAL(2*LVAL+1)
    DO MVAL = -LVAL,LVAL
      LM = LVAL*LVAL + LVAL + MVAL + 1
!c
!c---> set up of the integrands v1 and v2
!c
      V1(1) = 0.0D0
      V2(1) = 0.0D0
      DO I = 2,IRS1
        RL = CELL(IATOM)%RMESH(I)**LVAL
        V1(I) = DENSITY(IATOM)%RHO2NS(I,LM,1)*RL*CELL(IATOM)%DRMESHDI(I)
        V2(I) = DENSITY(IATOM)%RHO2NS(I,LM,1)/CELL(IATOM)%RMESH(I)/RL*CELL(IATOM)%DRMESHDI(I)
      END DO ! I
!c
!c---> convolute charge density of interstial with shape function
!c        if ins.gt.0
!c
      IF (ins.NE.0) THEN
        DO I = IRS1 + 1,IRC1
          V1(I) = 0.0D0
        END DO !I
        ISTART = GAUNTSHAPE(LMAX)%IMAXSH(LM-1) + 1
        IEND   = GAUNTSHAPE(LMAX)%IMAXSH(LM)
        DO J = ISTART,IEND
          LM2 = GAUNTSHAPE(LMAX)%ILM(J,2)
          LM3 = GAUNTSHAPE(LMAX)%ILM(J,3)
          IF (SHAPEFUN(IATOM)%lmused(lm3).GT.0) THEN !LMSP(ICELL,LM3).GT.0
            IFUN = SHAPEFUN(IATOM)%LM2INDEX(LM3) !IFUNM(ICELL,LM3)
            DO I = IRS1 + 1,IRC1
              V1(I) = V1(I) + GAUNTSHAPE(LMAX)%GSH(J)*DENSITY(IATOM)%RHO2NS(I,LM2,1)* &
                              SHAPEFUN(IATOM)%THETAS(I-IRS1,IFUN)
            END DO !I
          END IF
        END DO !J

        DO I = IRS1 + 1,IRC1
          RL = CELL(IATOM)%RMESH(I)**LVAL
          V2(I) = V1(I)/CELL(IATOM)%RMESH(I)/RL*CELL(IATOM)%DRMESHDI(I)
          V1(I) = V1(I)*RL*CELL(IATOM)%DRMESHDI(I)
        END DO !I
      END IF
!c
!c---> now integrate v1 and v2
!c
!             write(*,*) 'call soutk'
!             write(*,*) ircutm
!             write(*,*) CELL(IATOM)%NRCUT
      CALL SOUTK(V1,VINT1,CELL(IATOM)%NPAN,IRCUTM)
!             write(*,*) 'call sinwk'
      CALL SINWK(V2,VINT2,CELL(IATOM)%NPAN,IRCUTM)
!c
!c---> gather all parts
!c
      IF (LM.EQ.1) THEN
        VPOT_OUT(1,LM,1,IATOM) = FAC*VINT2(1)
      ELSE
        VPOT_OUT(1,LM,1,IATOM) = 0.0D0
      END IF

      DO I = 2,IRC1
        RL = CELL(IATOM)%RMESH(I)**LVAL
        VPOT_OUT(I,LM,1,IATOM) = FAC* (VINT1(I)/CELL(IATOM)%RMESH(I)/RL+VINT2(I)*RL)
      END DO
!c
!c---> store charge moment - in case of kshape.gt.0 this is the moment
!c      of the charge in the muffin tin sphere
!c            CMOM(LM,IATYP) = VINT1(IRS1)

!c
!c---> store charge moment of interstial in case of kshape.gt.0
!c
!             IF (KSHAPE.NE.0) CMINST(LM,IATOM) = VINT1(IRC1) - VINT1(IRS1)
      IF (ins.NE.0) THEN
         CMOM_INTERST(LM,IATOM) = VINT1(IRC1) - VINT1(IRS1)
      END IF
!       ELSE
        CMOM(LM,IATOM) = VINT1(IRC1)
!       END IF
!c
      IF (NSPIN.EQ.2) THEN
        DO I = 1,IRC1
          VPOT_OUT(I,LM,2,IATOM) = VPOT_OUT(I,LM,1,IATOM)
        END DO  
      END IF

    END DO !M

  END DO !L

END DO !IATOM

END SUBROUTINE
END MODULE MOD_VINTRAS
