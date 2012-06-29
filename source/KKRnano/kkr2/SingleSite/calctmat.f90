subroutine CALCTMAT(LDAU,NLDAU,ICST, &
NSRA,EZ, &
DRDI,R,VINS,VISP,ZAT,IPAN, &
IRCUT,CLEB,LOFLM,ICLEB,IEND, &
TMATN,TR_ALPH,LMAX, &
LLDAU,WMLDAU_ISPIN, &
!                        new input parameters after inc.p removal
ncleb, ipand, irmd, irnsd)
  implicit none

  !     .. Parameters ..

  integer ncleb
  integer ipand
  integer irmd
  integer irnsd

  !      PARAMETER          (LMMAXD= (LMAXD+1)**2)
  !      PARAMETER          (LMAXD1 = LMAXD + 1)
  !      PARAMETER          (MMAXD=2*LMAXD+1)
  !      PARAMETER          (LMPOTD= (LPOTD+1)**2)
  !                                = (2*LMAX+1)**2
  !      PARAMETER          (IRMIND=IRMD-IRNSD)
  !      PARAMETER          (LM2D= (2*LMAXD+1)**2)

  double precision    CVLIGHT
  double complex, parameter :: CZERO = (0.0d0, 0.0d0)

  !     ..
  !     .. Scalar Arguments ..
  double precision   ZAT
  integer            ICST,IEND,IPAN,NSRA,LMAX,NLDAU
  double complex     TR_ALPH,EZ
  logical            LDAU
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX     TMATN(LMMAXD,LMMAXD)
  double complex     TMATN((LMAX+1)**2,(LMAX+1)**2)

  double precision   CLEB(NCLEB,2),DRDI(IRMD),R(IRMD)

  !     DOUBLE PRECISION   VINS(IRMIND:IRMD,LMPOTD)
  double precision   VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2)

  double precision   VISP(IRMD)

  double precision   WMLDAU_ISPIN(2*LMAX+1, 2*LMAX+1, LMAX + 1)

  integer            ICLEB(NCLEB,3),IRCUT(0:IPAND)

  !     INTEGER            LOFLM(LM2D)
  integer            LOFLM((2*LMAX+1)**2)

  !     INTEGER            LLDAU(LMAXD1)
  integer            LLDAU(LMAX + 1)

  !     ..
  !     .. Local Scalars ..
  double complex     ERYD,EK,DET
  double precision   PI
  integer            LM1,LM2,L, &
  LMLO,LMHI,MMAX,IM,ILDAU
  !     ..
  !     .. Local Arrays ..
  !      DOUBLE COMPLEX     ALPHA(0:LMAXD),
  !     +                   FZ(IRMD,0:LMAXD),
  !     +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !     +                   PZ(IRMD,0:LMAXD),
  !     +                   QZ(IRMD,0:LMAXD),
  !     +                   SZ(IRMD,0:LMAXD),TMAT(0:LMAXD)
  !      DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
  !     +                   LDAUCUT(IRMD),
  !     +                   WMLDAUAV(LMAXD1)

  double complex    ALPHA(0:LMAX), &
  FZ(IRMD,0:LMAX)

  !     DOUBLE COMPLEX    PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2)
  double complex    PNS((LMAX+1)**2,(LMAX+1)**2,(IRMD-IRNSD):IRMD,2)
  double complex    PZ(IRMD,0:LMAX), &
  QZ(IRMD,0:LMAX), &
  SZ(IRMD,0:LMAX),TMAT(0:LMAX)
  double precision  RS(IRMD,0:LMAX),S(0:LMAX), &
  LDAUCUT(IRMD), &
  WMLDAUAV(LMAX + 1)

  !     ..
  !     .. External Subroutines ..
  external CRADWF,PNSTMAT,WFMESH

  integer             LMMAXD

  LMMAXD= (LMAX + 1)**2

  PI = 4.D0*ATAN(1.D0)
  CVLIGHT=274.0720442D0  ! must not be parameter?

!------------------------------------------------------------------------------
! Initialisation of local variables to be on the safe side
  EK = CZERO
  DET = CZERO
  ALPHA = CZERO
  FZ = CZERO
  PNS = CZERO
  PZ = CZERO
  QZ = CZERO
  SZ = CZERO
  TMAT = CZERO
  RS = 0.0d0
  S = 0.0d0
  LDAUCUT = 0.0d0
  WMLDAUAV = 0.0d0
! -----------------------------------------------------------------------------

  !LDAU

  if (LDAU) then

    do ILDAU=1,NLDAU

      WMLDAUAV(ILDAU) = 0.0D0
      LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
      LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
      MMAX = LMHI - LMLO + 1
      do IM = 1,MMAX
        WMLDAUAV(ILDAU)=WMLDAUAV(ILDAU)+WMLDAU_ISPIN(IM,IM,ILDAU)
      enddo
      WMLDAUAV(ILDAU) = WMLDAUAV(ILDAU)/DBLE(MMAX)

    enddo

  ! -> Note: Application if WLDAU makes the potential discontinuous.
  !    A cutoff can be used if necessary to make the potential continuous
  !    for example (array bounds should be adjusted):

  !       IF(TEST('CUTOFF  ')) THEN
  !         DO IR = 1,IRMD
  !           LDAUCUT(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) *
  !    &                    ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
  !           LDAUCUT(IR) = 1D0/LDAUCUT(IR)
  !         ENDDO
  !       ELSE
  !         DO IR = 1,IRMD
  !           LDAUCUT(IR) = 1.D0
  !         ENDDO
  !       ENDIF

  endif

  !LDAU


  do LM2 = 1,LMMAXD
    do LM1 = 1,LMMAXD
      TMATN(LM1,LM2) = (0.0D0,0.0D0)
    end do
  end do

  ERYD = EZ

  call WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN), &
  IRMD,LMAX)

  call CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S, &
  PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT, &
  LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT, &
  lmax, irmd, ipand)

  !-----------------------------------------------------------------------
  call PNSTMAT(DRDI,EK,ICST,PZ,QZ,FZ,SZ,PNS, &
  TMATN, &
  VINS,IPAN,IRCUT,NSRA,CLEB,ICLEB,IEND,LOFLM, &
  TMAT,DET,LMAX, &
  LDAU,NLDAU,LLDAU, &
  WMLDAU_ISPIN,WMLDAUAV,LDAUCUT, &
  lmax, irmd, irnsd, ipand, ncleb)

  TR_ALPH = LOG(DET)
  do L=0,LMAX
    TR_ALPH = TR_ALPH + (L+L+1)*LOG(ALPHA(L))
  enddo


  return
end
