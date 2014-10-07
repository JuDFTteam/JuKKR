subroutine RHOVAL(LDORHOEF,ICST,IELAST,NSRA, &
                  ISPIN,NSPIN, &
                  EZ,WEZ,DRDI,R,IRMIN, &
                  VINS,VISP,ZAT,IPAN,IRCUT, &
                  THETAS,IFUNM,LMSP,RHO2NS,R2NEF,DEN, &
                  ESPV,CLEB,LOFLM,ICLEB,IEND,JEND, &
                  GMATN, &                                 ! input
                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &       ! input
                  DMATLDAU, &                              ! output
                  ! new parameters after inc.p removal
                  iemxd, &
                  lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

  implicit none

  integer :: iemxd
  integer :: lmaxd
  integer :: irmd
  integer :: ncleb
  integer :: irnsd
  integer :: irid
  integer :: ipand
  integer :: nfund

  double precision ::    CVLIGHT
  parameter          (CVLIGHT=274.0720442D0)
  double complex      CONE
  parameter          ( CONE=(1D0,0D0) )
  double complex, parameter :: CZERO = (0.0d0, 0.0d0)
  !     ..
  !     .. Scalar Arguments ..
  double precision ::   ZAT
  integer ::            ICST,IELAST,IEND,IPAN,ISPIN,NSPIN,NSRA, &
  IRMIN,NLDAU
  logical ::            LDORHOEF,LDAU
  !     ..
  !     .. Array Arguments ..
  !     DOUBLE COMPLEX     DEN(0:LMAXD1,IEMXD),EZ(IEMXD),
  !    +                   WEZ(IEMXD),
  !    +                   PHILDAU(IRMD,LMAXD1),
  !    +                   DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  !     DOUBLE PRECISION   CLEB(NCLEB,2),DRDI(IRMD),
  !    +                   ESPV(0:LMAXD1,1),
  !    +                   R(IRMD),RHO2NS(IRMD,LMPOTD,2),
  !    +                   R2NEF(IRMD,LMPOTD,2),   ! at fermi energy
  !    +                   THETAS(IRID,NFUND),VINS(IRMIND:IRMD,LMPOTD),
  !    +                   VISP(IRMD),
  !    +                   WMLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  !     INTEGER            ICLEB(NCLEB,3),IFUNM(LMXSPD),IRCUT(0:IPAND),
  !    +                   JEND(LMPOTD,0:LMAXD,0:LMAXD),
  !    +                   LMSP(LMXSPD),LOFLM(LM2D),
  !    +                   LLDAU(LMAXD1)

  double complex     DEN(0:LMAXD+1,IEMXD)
  double complex     EZ(IEMXD)
  double complex     WEZ(IEMXD)
  double complex     PHILDAU(IRMD,LMAXD+1)
  double complex     DMATLDAU(2*LMAXD+1,2*LMAXD+1,NSPIN,LMAXD+1)

  double complex    GMATN((LMAXD+1)**2,(LMAXD+1)**2,IEMXD,NSPIN)

  double precision ::   CLEB(NCLEB,2)
  double precision ::   DRDI(IRMD)
  double precision ::   ESPV(0:LMAXD+1,1)
  double precision ::   R(IRMD)
  double precision ::   RHO2NS(IRMD,(2*LMAXD+1)**2,2)
  double precision ::   R2NEF(IRMD,(2*LMAXD+1)**2,2)
  double precision ::   THETAS(IRID,NFUND)
  double precision ::   VINS(IRMD-IRNSD:IRMD,(2*LMAXD+1)**2)
  double precision ::   VISP(IRMD)
  double precision ::   WMLDAU(2*LMAXD+1,2*LMAXD+1,LMAXD+1,NSPIN)

  integer ::            ICLEB(NCLEB,3)
  integer ::            IFUNM((4*LMAXD+1)**2)
  integer ::            IRCUT(0:IPAND)
  integer ::            JEND((2*LMAXD+1)**2,0:LMAXD,0:LMAXD)
  integer ::            LMSP((4*LMAXD+1)**2)
  integer ::            LOFLM((2*LMAXD+1)**2)
  integer ::            LLDAU(LMAXD+1)

  !     .. Local Scalars ..
  double complex     DF,ERYD,EK
  integer ::            IDIM,IE,IR,L,LM1,LM2, &
  LMLO,LMHI,MMAX,IM,ILDAU
  !     ..
  !     .. Local Arrays ..
  !     DOUBLE COMPLEX     ALPHA(0:LMAXD),AR(LMMAXD,LMMAXD),
  !    +                   DR(LMMAXD,LMMAXD),
  !    +                   CR(LMMAXD,LMMAXD),
  !    +                   EKL(0:LMAXD),FZ(IRMD,0:LMAXD),
  !    +                   GMATLL(LMMAXD,LMMAXD),
  !    +                   GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND),
  !    +                   PNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !    +                   PZ(IRMD,0:LMAXD),
  !    +                   QNS(LMMAXD,LMMAXD,IRMIND:IRMD,2),
  !    +                   QZ(IRMD,0:LMAXD),SZ(IRMD,0:LMAXD),
  !    +                   TMAT(0:LMAXD)
  !     DOUBLE PRECISION   RS(IRMD,0:LMAXD),S(0:LMAXD),
  !    +                   LDAUCUT(IRMD),
  !    +                   WMLDAUAV(LMAXD1)
  !     DOUBLE COMPLEX     DENDUM(0:LMAXD1)

  ! The following arrays are local
  double complex    ALPHA(0:LMAXD)
  double complex    AR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    DR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    CR((LMAXD+1)**2,(LMAXD+1)**2)
  double complex    EKL(0:LMAXD)
  double complex    FZ(IRMD,0:LMAXD)
  double complex    GMATLL((LMAXD+1)**2,(LMAXD+1)**2)

  double complex    QZ(IRMD,0:LMAXD)
  double complex    SZ(IRMD,0:LMAXD)
  double complex    TMAT(0:LMAXD)
  double complex    PZ(IRMD,0:LMAXD)

  double precision ::   RS(IRMD,0:LMAXD)
  double precision ::   S(0:LMAXD)
  double precision ::   LDAUCUT(IRMD)
  double precision ::   WMLDAUAV(LMAXD+1)

  double complex     DENDUM(0:LMAXD+1)

  ! dynamically allocate large arrays
  ! DOUBLE COMPLEX    PNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
  ! DOUBLE COMPLEX    QNS((LMAXD+1)**2,(LMAXD+1)**2,IRMD-IRNSD:IRMD,2)
  double complex, dimension(:,:,:,:), allocatable :: PNS
  double complex, dimension(:,:,:,:), allocatable :: QNS

  !     ..
  !     .. External Subroutines ..
  external DAXPY,DSCAL,CRADWF,PNSQNS,RHONS,WFMESH
  !     ..
  !     .. Intrinsic Functions ..
  intrinsic ATAN,DBLE,DIMAG,SQRT
  !     ..

  integer :: memory_stat
  logical :: memory_fail

  integer ::             LMMAXD
  integer ::             LMAXD1
  integer ::             NSPIND
  integer ::             LMPOTD

  LMPOTD = (2*LMAXD+1)**2
  LMMAXD= (LMAXD+1)**2
  LMAXD1= LMAXD+1
  NSPIND = NSPIN

  !------------------- Array allocations ---------------------------------
  memory_stat = 0
  memory_fail = .false.

  allocate(PNS(LMMAXD,LMMAXD,IRMD-IRNSD:IRMD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(QNS(LMMAXD,LMMAXD,IRMD-IRNSD:IRMD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "RHOVAL: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  end if

  !-----------------------------------------------------------------------
  ! Initialise local variables to be on the safe side
  EK = CZERO
  PNS = CZERO
  QNS = CZERO
  ALPHA = CZERO
  AR = CZERO
  DR = CZERO
  CR = CZERO
  EKL = CZERO
  FZ = CZERO
  !GMATLL = CZERO initialised further down
  QZ = CZERO
  SZ = CZERO
  TMAT = CZERO
  PZ = CZERO
  RS = 0.0d0
  S = 0.0d0
  LDAUCUT = 0.0d0
  WMLDAUAV = 0.0d0
  DENDUM = CZERO

  !-----------------------------------------------------------------------
  ! LDAU

  if (LDAU) then
    do ILDAU=1,NLDAU
      WMLDAUAV(ILDAU) = 0.D0
      LMLO = LLDAU(ILDAU)*LLDAU(ILDAU) + 1
      LMHI = (LLDAU(ILDAU)+1)*(LLDAU(ILDAU)+1)
      MMAX = LMHI - LMLO + 1
      do IM = 1,MMAX
        WMLDAUAV(ILDAU)=WMLDAUAV(ILDAU)+WMLDAU(IM,IM,ILDAU,ISPIN)
      enddo
      WMLDAUAV(ILDAU) = WMLDAUAV(ILDAU)/DBLE(MMAX)
    enddo
    
    ! -> Note: Application if WLDAU makes the potential discontinuous.
    !    A cutoff can be used if necessary to make the potential continuous
    !    for example (array bounds should be adjusted):
    
!    if(TEST('CUTOFF  ')) then
!      do IR = 1,IRMD
!        LDAUCUT(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) * &
!        ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
!        LDAUCUT(IR) = 1D0/LDAUCUT(IR)
!      enddo
!    else
!      do IR = 1,IRMD
!        LDAUCUT(IR) = 1.D0
!      enddo
!    endif
  endif

  ! LDAU
  !-----------------------------------------------------------------------



  do LM1 = 1,LMPOTD
    do IR = 1,IRMD
      RHO2NS(IR,LM1,ISPIN) = 0.0D0
      R2NEF(IR,LM1,ISPIN) = 0.0D0
    end do
  end do

  ESPV = 0.0D0

  do IE = 1,IELAST

    do LM2 = 1,LMMAXD
      do LM1 = 1,LMMAXD
        GMATLL(LM1,LM2) = GMATN(LM1,LM2,IE,ISPIN)
      end do
    end do

    ERYD = EZ(IE)
    DF = WEZ(IE)/DBLE(NSPIN)
    
    !=======================================================================
    call WFMESH(ERYD,EK,CVLIGHT,NSRA,ZAT,R,S,RS,IRCUT(IPAN),IRMD,LMAXD)

    call CRADWF(ERYD,EK,NSRA,ALPHA,IPAN,IRCUT,CVLIGHT,RS,S, &
                PZ,FZ,QZ,SZ,TMAT,VISP,DRDI,R,ZAT, &
                LDAU,NLDAU,LLDAU,WMLDAUAV,LDAUCUT, &
                lmaxd, irmd, ipand)
    !-----------------------------------------------------------------------
    ! non-spherical
    
    call PNSQNS(AR,CR,DR,DRDI,EK,ICST,PZ,QZ,FZ,SZ, &
                PNS,QNS,NSRA,VINS,IPAN,IRCUT, &
                CLEB,ICLEB,IEND,LOFLM,LMAXD,ISPIN, &
                LDAU,NLDAU,LLDAU, &
                WMLDAU,WMLDAUAV,LDAUCUT, &
                lmaxd, nspind, irmd, irnsd, ipand, ncleb)


    do L = 0,LMAXD
      EKL(L) = EK*DBLE(2*L+1)
    end do
    !-----------------------------------------------------------------------
    call RHONS(DEN(0,IE),DF,DRDI,GMATLL,EK, &
               RHO2NS(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP, &
               NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB, &
               JEND,IEND,EKL, &
               lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

    !-----------------------------------------------------------------------


    do L = 0,LMAXD1
      ESPV(L,1) = ESPV(L,1) + DIMAG(ERYD*DEN(L,IE)*DF)
    end do
    
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !     get the charge at the Fermi energy (IELAST)
    !     call with the energy weight CONE --> not overwrite DF
    !          with the dummy DENDUM       --> not overwrite DEN
    
    if ( (IE == IELAST) .and. (LDORHOEF) ) then
      call RHONS(DENDUM,CONE,DRDI,GMATLL,EK, &
                 R2NEF(1,1,ISPIN),IPAN,IRCUT,THETAS,IFUNM,LMSP, &
                 NSRA,QNS,PNS,AR,CR,PZ,FZ,QZ,SZ,CLEB(1,1),ICLEB, &
                 JEND,IEND,EKL, &
                 lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)
    end if

    
    if (LDAU .and. NLDAU >= 1) then
      call LDAUDMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL, &
                    IPAN,IRCUT,DRDI,EK, &
                    IRMIN,LLDAU,PHILDAU,NLDAU, &
                    DMATLDAU,ISPIN, &
                    lmaxd, nspind, irmd, irnsd, ipand)
        
    endif


  end do

  ! this should really be separated into another routine
  if (ISPIN == 2) then
    IDIM = IRMD*LMPOTD
    call DSCAL(IDIM,2.D0,RHO2NS(1,1,1),1)
    call DAXPY(IDIM,-0.5D0,RHO2NS(1,1,1),1,RHO2NS(1,1,2),1)
    call DAXPY(IDIM,1.0D0,RHO2NS(1,1,2),1,RHO2NS(1,1,1),1)
    
    ! --> do the same at the Fermi energy
    
    call DSCAL(IDIM,2.D0,R2NEF(1,1,1),1)
    call DAXPY(IDIM,-0.5D0,R2NEF(1,1,1),1,R2NEF(1,1,2),1)
    call DAXPY(IDIM,1.0D0,R2NEF(1,1,2),1,R2NEF(1,1,1),1)
  end if

  !------------------- Array deallocations ---------------------------------
  deallocate(PNS)
  deallocate(QNS)
!-----------------------------------------------------------------------

end subroutine RHOVAL
