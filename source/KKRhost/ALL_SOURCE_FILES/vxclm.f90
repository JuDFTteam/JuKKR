module mod_vxclm

contains

SUBROUTINE VXCLM(EXC,KTE,KXC,LMAX,NSPIN,IATYP,RHO2NS,V,R,DRDI, &
                 IRWS,IRCUT,IPAN,KSHAPE,GSH,ILM,IMAXSH, &
                 IFUNM,THETAS,YR,WTYR,IJEND,LMSP)
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
! in the paramagnetic case the "spin-density" is set equal zero .
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
!                        cor error 23/6/1996
!-----------------------------------------------------------------------
!INCLUDE 'inc.p'
!.. Parameters ..
!integer LMPOTD
!parameter (LMPOTD= (LPOTD+1)**2)
!integer LMXSPD
!parameter (LMXSPD= (2*LPOTD+1)**2)
use mod_DataTypes, only: dp
use global_variables
   use mod_vosko
   use mod_vxcspo
  use mod_simpk
  use mod_simp3
implicit none
!..
!.. Scalar Arguments ..
integer IATYP,IJEND,IPAN,IRWS,KSHAPE,KTE,KXC,LMAX,NSPIN
!..
!.. Array Arguments ..
real (kind=dp) DRDI(IRMD),EXC(0:LPOTD,*),GSH(*),R(IRMD), &
                 RHO2NS(IRMD,LMPOTD,2),THETAS(IRID,NFUND), &
                 V(IRMD,LMPOTD,2),WTYR(IJEND,*),YR(IJEND,*)
integer IFUNM(LMXSPD)
integer ILM(NGSHD,3),IMAXSH(0:LMPOTD),IRCUT(0:IPAND), &
        LMSP(LMXSPD)
!..
!.. Local Scalars ..
real (kind=dp) ELMXC,FPI,FPIPR2,VLMXC,VXC1,VXC2,VXC3,factor
integer IFUN,IJ,IPOT,IR,IRC1,IRH,IRS1,IS,ISPIN,J,L,LM,LM2,LMMAX,M
!..
!.. Local Arrays ..
real (kind=dp) ER(IRMD,0:LPOTD),ESTOR(IRMD,LMPOTD),EXCIJ(IJEND), &
                 FPRHO(IJEND,2),VXC(IJEND,2),VXCR(2:3,2)
!..
!.. External Functions ..
real (kind=dp) DDOT
external DDOT

WRITE(1337,*) 'Including cutoff of vxc for small density'
FPI = 16.0D0*ATAN(1.0D0)
LMMAX = (LMAX+1)* (LMAX+1)


! loop over given representive atoms

IF (KSHAPE.NE.0) THEN
  IRC1 = IRCUT(IPAN)
  IRS1 = IRCUT(1)

ELSE

  IRC1 = IRWS
  IRS1 = IRC1
END IF

DO ISPIN = 1,NSPIN
  VXCR(2,ISPIN) = 0.0D0
  VXCR(3,ISPIN) = 0.0D0
END DO

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

! loop over radial mesh


DO IR = 2,IRC1


! generate the densities on an angular mesh

  DO IS = 1,2
    DO IJ = 1,IJEND
      FPRHO(IJ,IS) = 0.D0
    END DO
  END DO

  FPIPR2 = FPI/R(IR)**2
  DO ISPIN = 1,NSPIN
    DO LM = 1,LMMAX
      CALL DAXPY(IJEND,RHO2NS(IR,LM,ISPIN)*FPIPR2,YR(1,LM),1, &
                 FPRHO(1,ISPIN),1)
    END DO
  END DO

! calculate the ex.-cor. potential

  IF (KXC.LE.1) THEN
    CALL VXCSPO(EXCIJ,FPRHO,VXC,KXC,IJEND,IJEND)
  ELSE
    CALL VOSKO(EXCIJ,FPRHO,VXC,IJEND,IJEND)
  END IF

    do ij=1,ijend
    factor = (1.d0-exp(-abs(fprho(ij,1))*1000.d0))
    do ispin=1,nspin
    vxc(ij,ispin) = &
      vxc(ij,ispin) * factor  !cutoff
    enddo
    enddo


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
      CALL SIMPK(ER(1,L),EXC(L,IATYP),IPAN,IRCUT,DRDI)
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

END SUBROUTINE VXCLM

end module mod_vxclm
