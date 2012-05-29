subroutine CALCDTMAT &
(LDAU,NLDAU,ICST, &
NSRA,EZ,DZ, &
DRDI,R,VINS, &
VISP,ZAT,IPAN, &
IRCUT,CLEB,LOFLM,ICLEB,IEND, &
DTDE,TR_ALPH,LMAX, &
LLDAU,WMLDAU_ISPIN, &
!                        new input parameters after inc.p replace
ncleb, ipand, irmd, irnsd)

  implicit none

  integer ncleb
  integer ipand
  integer irmd
  integer irnsd

  !      PARAMETER          (LMMAXD= (LMAXD+1)**2)
  !      PARAMETER          (LMAXD1 = LMAXD + 1)
  !      PARAMETER          (MMAXD=2*LMAXD+1)
  !      PARAMETER          (LM2D= (2*LMAXD+1)**2)
  !      PARAMETER          (LMPOTD= (LPOTD+1)**2)
  !                                = (2*LMAX+1)**2)
  !      PARAMETER          (IRMIND=IRMD-IRNSD)

  integer            LMAX,LM1,LM2,IEND
  !     ..

  double complex     DZ
  double complex     EZ,EZ1,EZ2
  double precision   CLEB(NCLEB,2)

  double precision   VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2), &
  VISP(IRMD), &
  WMLDAU_ISPIN(2*LMAX+1, 2*LMAX+1, LMAX + 1)
  !     ..
  !     DOUBLE COMPLEX     DTDE(LMMAXD,LMMAXD)
  double complex     DTDE((LMAX+1)**2,(LMAX+1)**2)
  double complex     TR_ALPH,TR_ALPH1,TR_ALPH2

  !     DOUBLE COMPLEX     TMATN1(LMMAXD,LMMAXD),TMATN2(LMMAXD,LMMAXD)
  double complex     TMATN1((LMAX+1)**2,(LMAX+1)**2)
  double complex     TMATN2((LMAX+1)**2,(LMAX+1)**2)
  double precision   DRDI(IRMD)
  double precision   R(IRMD)
  double precision   ZAT
  ! ----------------------------------------------------------------------
  integer            IPAN,NLDAU

  !     INTEGER            IRCUT(0:IPAND),LLDAU(LMAXD1)
  integer            IRCUT(0:IPAND),LLDAU(LMAX + 1)
  !     INTEGER            ICLEB(NCLEB,3),LOFLM(LM2D)
  integer            ICLEB(NCLEB,3),LOFLM((2*LMAX+1)**2)

  integer            ICST,NSRA
  logical            LDAU
  !     ..

  integer             LMMAXD

  LMMAXD= (LMAX+1)**2

  ! Interpolation points for difference quotient
  !.. perpendicular to the contour
  EZ1 = EZ + DZ
  !.. parallel to the contour
  EZ2 = EZ - DZ

  call CALCTMAT(LDAU,NLDAU,ICST, &
  NSRA,EZ1, &
  DRDI,R,VINS, &
  VISP,ZAT,IPAN, &
  IRCUT,CLEB,LOFLM,ICLEB,IEND, &
  TMATN1,TR_ALPH1,LMAX, &
  LLDAU,WMLDAU_ISPIN, &
  ncleb, ipand, irmd, irnsd)

  call CALCTMAT(LDAU,NLDAU,ICST, &
  NSRA,EZ2, &
  DRDI,R,VINS, &
  VISP,ZAT,IPAN, &
  IRCUT,CLEB,LOFLM,ICLEB,IEND, &
  TMATN2,TR_ALPH2,LMAX, &
  LLDAU,WMLDAU_ISPIN, &
  ncleb, ipand, irmd, irnsd)
  !====================================================================

  ! dT(E)
  ! -----
  !  dE
  do LM2 = 1,LMMAXD
    do LM1 = 1,LMMAXD
      DTDE(LM1,LM2) = (TMATN1(LM1,LM2)-TMATN2(LM1,LM2))*0.5D0/DZ
    enddo
  enddo

  TR_ALPH = -(TR_ALPH1-TR_ALPH2)*0.5D0/DZ
  return
end

!------------------------------------------------------------------------------
!> Get \Delta E(z) for difference quotient to calculate the derivative of
!> \Delta T_ref.
!> @param[out] DZ Delta E(z)
subroutine calcdtmat_DeltaEz(DZ, IE, NPNT1, NPNT2, NPNT3, TK)
  implicit none
  double complex :: DZ
  integer :: IE
  integer :: NPNT1
  integer :: NPNT2
  integer :: NPNT3
  double precision :: TK

  !------------------------------------
  double precision :: KB
  double precision :: PI

  PI = 4.0D0*ATAN(1.0D0)
  KB = 0.6333659D-5

  DZ = DCMPLX(0.0D0, 0.0D0)

  if (IE.le.NPNT1 .or. IE.gt.(NPNT1+NPNT2+NPNT3)) then
    DZ = DCMPLX(0.01D0*PI*KB*TK,0.0D0)
  else
    DZ = DCMPLX(0.0D0,0.01D0*PI*KB*TK)
  end if
end subroutine
