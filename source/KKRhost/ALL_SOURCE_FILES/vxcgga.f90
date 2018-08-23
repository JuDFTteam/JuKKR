module mod_vxcgga

contains

SUBROUTINE VXCGGA(EXC,KTE,KXC,LMAX,NSPIN,IATYP,RHO2NS,V,R,DRDI,A, &
                  IRWS,IRCUT,IPAN,KSHAPE,GSH,ILM,IMAXSH, &
                  IFUNM,THETAS,WTYR,IJEND,LMSP,THET,YLM,DYLMT1, &
                  DYLMT2,DYLMF1,DYLMF2,DYLMTF)
!-----------------------------------------------------------------------
! add the exchange-correlation-potential to the given potential
! and if total energies should be calculated (kte=1) the exchange-
! correlation-energies are calculated .
! use as input the charge density times r**2 (rho2ns(...,1)) and
! in the spin-polarized case (nspin=2) the spin density times r**2
! (rho2ns(...,2)) .
! the density times 4 pi is generated at an angular mesh .
! the exchange-correlation potential and the exchange-correlation
! energy are calculated at those mesh points with a subroutine .
! in the non-spin-polarized case the "spin-density" is
! set equal zero .
! after that the exchange-correlation potential and in the case of
! total energies (kte=1) the exchange-correlation energy are
! expanded into spherical harmonics .
! the ex.-cor. potential is added to the given potential .
! the expansion into spherical harmonics uses the orthogonality
! of these harmonics . - therefore a gauss-legendre integration
! for "theta" and a gauss-tschebyscheff integration for "phi"
! is used .
! all needed values for the angular mesh and angular integration
! are generate in the subroutine sphere .
!
! the ex.-cor. potential is extrapolated to the origin only
! for the lm=1 value .
!
!                           b.drittler   june 1987
!
! modified for shape functions
!                                   b. drittler oct. 1989
! simplified and modified for Paragon X/PS
!                                   R. Zeller Nov. 1993
!-----------------------------------------------------------------------
!INCLUDE 'inc.p'
! Parameters ..
!INTEGER LMPOTD
!PARAMETER (LMPOTD= (LPOTD+1)**2)
!INTEGER LMXSPD
!PARAMETER (LMXSPD= (2*LPOTD+1)**2)
use global_variables
   use mod_mkxcpe2
   use mod_mkxcpe
   use mod_gradrl
   use mod_simpk
implicit none

! Scalar Arguments ..
DOUBLE PRECISION A
INTEGER IATYP,IJEND,IPAN,IRWS,KSHAPE,KTE,KXC,LMAX,NSPIN

! Array Arguments ..
DOUBLE PRECISION DRDI(IRMD),R(IRMD), &
                 DYLMF1(IJEND,LMPOTD),DYLMF2(IJEND,LMPOTD), &
                 DYLMT1(IJEND,LMPOTD),DYLMT2(IJEND,LMPOTD), &
                 DYLMTF(IJEND,LMPOTD),EXC(0:LPOTD,*),GSH(*), &
                 RHO2NS(IRMD,LMPOTD,2),THET(IJEND),WTYR(IJEND,*), &
                 THETAS(IRID,NFUND),V(IRMD,LMPOTD,2), &
                 YLM(IJEND,LMPOTD)
INTEGER IFUNM(LMXSPD)
INTEGER ILM(NGSHD,3),IMAXSH(0:LMPOTD),IRCUT(0:IPAND), &
        LMSP(LMXSPD)

! Local Scalars ..
DOUBLE PRECISION CHGDEN,DX,ELMXC,FPI,R1,R2,RPOINT,SPIDEN,VLMXC, &
                 VXC1,VXC2,VXC3,ZERO,ZERO1
INTEGER IFUN,IPAN1,IPOT,IR,IRC0,IRC1,IRH,IRS1,ISPIN,J,L,L1MAX,LM, &
        LM2,LMMAX,M,MESH,NSPIN2

! Local Arrays ..
DOUBLE PRECISION DDRRL(IRMD,LMPOTD),DDRRUL(IRMD,LMPOTD), &
                 DRRL(IRMD,LMPOTD),DRRUL(IRMD,LMPOTD), &
                 ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJEND), &
                 RHOL(IRMD,2,LMPOTD),RHOLM(LMPOTD,2),VXC(IJEND,2), &
                 VXCR(2:3,2)

! External Functions ..
DOUBLE PRECISION DDOT
EXTERNAL DDOT

! External Subroutines ..
EXTERNAL GRADRL,MKXCPE,SIMP3,SIMPK,MKXCPE2

! Intrinsic Functions ..
INTRINSIC ABS,ATAN,MOD

! Data statements ..
DATA ZERO,ZERO1/0.d0,1.d-12/

WRITE (1337,FMT=*) ' GGA CALCULATION '
FPI = 16.0D0*ATAN(1.0D0)
LMMAX = (LMAX+1)* (LMAX+1)

! loop over given representive atoms

IF (KSHAPE.NE.0) THEN
  IPAN1 = IPAN
  IRC1 = IRCUT(IPAN)
  IRS1 = IRCUT(1)
  IRC0 = 2
  IF (KREL.EQ.1) STOP ' REL + FULL POTENTIAL N/A '
ELSE

  IRC1 = IRWS
  IRS1 = IRC1
  IPAN1 = 1
  IRC0 = 2
  IF (KREL.EQ.1) IRC0 = 2 + MOD(IRCUT(1),2)
END IF

DO ISPIN = 1,NSPIN
  VXCR(2,ISPIN) = 0.0D0
  VXCR(3,ISPIN) = 0.0D0
end do

! initialize for ex.-cor. energy

IF (KTE.EQ.1) THEN
  DO L = 0,LMAX
    EXC(L,IATYP) = 0.0D0
    DO IR = 1,IRC1
      ER(IR,L) = 0.0D0
    END DO
  END DO

  DO LM = 1,LMMAX
    DO IR = 1,IRC1
      ESTOR(IR,LM) = 0.0D0
    END DO
  END DO
END IF

L1MAX = LMAX + 1
MESH = IRWS
DX = A

IF (NSPIN.EQ.2) THEN
  DO LM = 1,LMMAX
    DO IR = 2,MESH
      R1 = R(IR)
      R2 = R1*R1
      CHGDEN = RHO2NS(IR,LM,1)/R2
      SPIDEN = RHO2NS(IR,LM,2)/R2
      IF (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
      IF (ABS(SPIDEN).LE.ZERO1) SPIDEN = ZERO
      RHOL(IR,2,LM) = (CHGDEN+SPIDEN)/2.d0
      RHOL(IR,1,LM) = (CHGDEN-SPIDEN)/2.d0
    END DO

!   extrapolate

    RHOL(1,1,LM) = RHOL(2,1,LM)
    RHOL(1,2,LM) = RHOL(2,2,LM)
  END DO

ELSE

  DO LM = 1,LMMAX
    DO IR = 2,MESH
      R1 = R(IR)
      R2 = R1*R1

      CHGDEN = RHO2NS(IR,LM,1)/R2
      IF (ABS(CHGDEN).LE.ZERO1) CHGDEN = ZERO
      RHOL(IR,1,LM) = CHGDEN/2.d0
      RHOL(IR,2,LM) = CHGDEN/2.d0
    END DO

!   extrapolate
    RHOL(1,1,LM) = RHOL(2,1,LM)
    RHOL(1,2,LM) = RHOL(2,2,LM)
  END DO
END IF


CALL GRADRL(NSPIN,MESH,L1MAX,DX,RHOL,R,DRDI,IPAN1,IPAND,IRCUT, &
            DRRL,DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)


! loop over radial mesh


DO IR = IRC0,IRC1
  RPOINT = R(IR)

! calculate the ex.-cor. potential

  NSPIN2 = 2

  DO ISPIN = 1,NSPIN2
    DO LM = 1,LMMAX
      RHOLM(LM,ISPIN) = RHOL(IR,ISPIN,LM)
    END DO
  END DO

! only for spin-polarized

  ! PW91 functional
  IF(KXC.EQ.3)THEN
     CALL MKXCPE(NSPIN2,IR,IJEND,L1MAX,RPOINT,RHOLM,VXC,EXCIJ, &
                 THET,YLM,DYLMT1,DYLMT2,DYLMF1,DYLMF2,DYLMTF,DRRL, &
                 DDRRL,DRRUL,DDRRUL,IRMD,LMPOTD)
  ! PBE functional
  ELSEIF(KXC.EQ.4)THEN
     CALL MKXCPE2(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1, &
                  DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL, &
                  IRMD,LMPOTD,LMMAX,.false.)
  ! PBEsol functional
  ELSEIF(KXC.EQ.5)THEN
     CALL MKXCPE2(IR,IJEND,RPOINT,RHOLM,VXC,EXCIJ,YLM,DYLMT1, &
                  DYLMF1,DYLMF2,DYLMTF,DRRL,DDRRL,DRRUL,DDRRUL, &
                  IRMD,LMPOTD,LMMAX,.true.)
  ELSE
     WRITE(1337,*) ' KXC ???'
     STOP
  ENDIF




! expand the ex.-cor. potential into spherical harmonics ,
!   using the orthogonality

  DO ISPIN = 1,NSPIN

! determine the corresponding potential number

    IPOT = ISPIN
    DO LM = 1,LMMAX
      VLMXC = DDOT(IJEND,VXC(1,ISPIN),1,WTYR(1,LM),1)
      V(IR,LM,IPOT) = V(IR,LM,IPOT) + VLMXC

! store the ex.-c. potential of ir=2 and =3 for the extrapolation

      IF (LM.EQ.1 .AND. (IR.EQ.2.OR.IR.EQ.3)) VXCR(IR, &
          ISPIN) = VLMXC
    END DO
  END DO

! file er in case of total energies

  IF (KTE.EQ.1) THEN

! expand ex.-cor. energy into spherical harmonics
!   using the orthogonality

    DO L = 0,LMAX
      DO M = -L,L
        LM = L*L + L + M + 1
        ELMXC = DDOT(IJEND,EXCIJ,1,WTYR(1,LM),1)

! multiply the lm-component of the ex.-cor. energy with the same
! lm-component of the charge density times r**2 and sum over lm
! this corresponds to a integration over the angular .

        IF ((KSHAPE.NE.0) .AND. (IR.GT.IRS1)) THEN
          ESTOR(IR,LM) = ELMXC

        ELSE

          ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*ELMXC
        END IF

      END DO

    END DO

  END IF

END DO

! integrate er in case of total energies to get exc

IF (KTE.EQ.1) THEN
  IF (KSHAPE.EQ.0) THEN
    DO L = 0,LMAX
      CALL SIMP3(ER(1,L),EXC(L,IATYP),1,IRS1,DRDI)
    END DO

  ELSE

    DO L = 0,LMAX
      DO M = -L,L
        LM = L*L + L + M + 1

! convolute with shape function

        DO J = IMAXSH(LM-1) + 1,IMAXSH(LM)
          LM2 = ILM(J,2)
          IF (LMSP(ILM(J,3)).GT.0) THEN
            IFUN = IFUNM(ILM(J,3))
            DO IR = IRS1 + 1,IRC1
              IRH = IR - IRS1
              ER(IR,L) = ER(IR,L) + RHO2NS(IR,LM,1)*GSH(J)* &
                         THETAS(IRH,IFUN)*ESTOR(IR,LM2)
            END DO
          END IF
        END DO
      END DO
      CALL SIMPK(ER(1,L),EXC(L,IATYP),IPAN1,IRCUT,DRDI)
    END DO
  END IF

END IF

! extrapolate ex.-cor potential to the origin only for lm=1

DO ISPIN = 1,NSPIN
  IPOT = ISPIN

  VXC2 = VXCR(2,ISPIN)
  VXC3 = VXCR(3,ISPIN)
  VXC1 = VXC2 - R(2)* (VXC3-VXC2)/ (R(3)-R(2))

  V(1,1,IPOT) = V(1,1,IPOT) + VXC1
END DO

END

end module mod_vxcgga
