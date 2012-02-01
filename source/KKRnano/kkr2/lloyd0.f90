!> This routine uses the previously calculated Lloyd's formula terms
!> to calculate the correction to the DOS and renormalized weights for
!> the energy integration.
!> the routine uses communication
!> TODO: split, simplify, extract communication
subroutine LLOYD0(EZ,WEZ,CLEB,DRDI,R,IRMIN, &
                  VINS,VISP,THETAS,ZAT,ICLEB, &
                  IFUNM1,IPAN,IRCUT,LMSP1,JEND,LOFLM,NTCELL,ICST, &
                  IELAST,IEND,NAEZ,NSPIN,NSRA, &
                  WEZRN,RNORM, &                                     ! <
                  GMATN, &                                           ! >
                  LLY_GRDT, &                                        ! >
                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &                 ! >
                  DMATLDAU, &                                        ! <
                  LMPIC,MYLRANK, &                                   ! >
                  LCOMM,LSIZE, &                                     ! >
                  !   new input parameters after inc.p removal
                  prod_lmpid_smpid_empid, lmax, irmd, irnsd, iemxd, &
                  irid, nfund, ncelld, ipand, ncleb)

  !use mpi
  implicit none
  include 'mpif.h'

  integer :: lmax
  integer :: irmd
  integer :: irnsd
  integer :: iemxd
  integer :: prod_lmpid_smpid_empid
  integer :: irid
  integer :: nfund
  integer :: ncelld
  integer :: ipand
  integer :: ncleb


!  !     .. Parameters ..
!  INTEGER             LMMAXD,LMPOTD
!  PARAMETER          (LMPOTD= (LPOTD+1)**2)
!  PARAMETER          (LMMAXD= (LMAXD+1)**2)
!  INTEGER             LMAXD1
!  PARAMETER          (LMAXD1=LMAXD+1)
!  INTEGER             MMAXD
!  PARAMETER          (MMAXD=2*LMAXD+1)
!  INTEGER             LM2D
!  PARAMETER          (LM2D= (2*LMAXD+1)**2)
!  INTEGER             LMXSPD
!  PARAMETER          (LMXSPD= (2*LPOTD+1)**2)
!  INTEGER             IRMIND
!  PARAMETER          (IRMIND=IRMD-IRNSD)
  DOUBLE COMPLEX      CZERO
  PARAMETER          (CZERO=(0.0D0,0.0D0))
  !     ..
  !     .. Scalars ..
  integer::I1
  integer::ICELL
  integer::ICST
  integer::IE
  integer::IELAST
  integer::IEND
  integer::ISPIN
  integer::L
  integer::NAEZ
  integer::NSPIN
  integer::NSRA

  double complex :: WEZRN(IEMXD,2)
  double complex :: EZ(IEMXD)
  double complex :: WEZ(IEMXD)
  !double complex :: GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND)
  double complex :: GMATN((LMAX+1)**2, (LMAX+1)**2, IEMXD, NSPIN)
  double complex :: LLY_GRDT(IEMXD,NSPIN)


  double precision::DLOYD
  double precision::DLOYDINT
  double precision::DLOC
  double precision::D0LOC
  double precision::D0LOCINT
  double precision::D1LOC
  double precision::D1LOCINT
  double precision::RNORM(IEMXD,2)
  double precision::CLEB(NCLEB,2)
  double precision::DRDI(IRMD,NAEZ)
  double precision::R(IRMD,NAEZ)

  double precision::VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2, 2)
  double precision::VISP(IRMD,2)
  double precision::THETAS(IRID,NFUND,NCELLD)
  double precision::ZAT(NAEZ)

  integer::ICLEB(NCLEB,3)
  !integer::IFUNM1(LMXSPD,NAEZ) (2*LPOTD+1)**2
  integer::IFUNM1((4*LMAX+1)**2,NAEZ)
  integer::IPAN(NAEZ)
  integer::IRCUT(0:IPAND,NAEZ)
  integer::IRMIN(NAEZ)
  !integer::LMSP1(LMXSPD,NAEZD)
  integer::LMSP1((4*LMAX+1)**2, NAEZ)
  integer::JEND((LMAX+1)**2, 0:LMAX, 0:LMAX)
  integer::LOFLM((2*LMAX+1)**2)
  integer::NTCELL(NAEZ)


  ! -------------- LDA+U --------------------------------------------
  !double complex :: DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  double complex :: DMATLDAU(2*LMAX+1,2*LMAX+1,NSPIN,LMAX+1)
  double complex :: PHILDAU(IRMD,LMAX+1)
  double precision::WMLDAU(2*LMAX+1,2*LMAX+1,NSPIN,LMAX+1)
  integer::NLDAU
  logical::LDAU
  integer::LLDAU(LMAX+1)
  ! -----------------------------------------------------------------

  !     ..
  !     .. MPI ..
  !     .. N-MPI
  integer:: IERR
  integer::MAPBLOCK

  integer::MYLRANK(prod_lmpid_smpid_empid)
  integer::LCOMM(prod_lmpid_smpid_empid)
  integer::LSIZE(prod_lmpid_smpid_empid)
  integer::LMPIC

  !     Dynamically allocated arrays
  !double complex :: WORK1(4*IEMXD)
  !double complex :: WORK2(4*IEMXD)
  !double complex :: DOS(IEMXD,2) ! local
  !double complex :: DOS0(IEMXD)  ! loc
  !double complex :: DOS1(IEMXD)  ! loc
  !double complex :: DEN0(0:LMAX+1,IEMXD,NSPIN) ! loc
  !double precision::RHO2N1(IRMD,LMPOTD,2)
  !double precision::RHO2N1(IRMD,(2*LMAX+1)**2, 2)   ! loc
  !double precision::RHO2N2(IRMD,LMPOTD,2)
  !double precision::RHO2N2(IRMD,(2*LMAX+1)**2, 2)   ! loc
  !double precision::ESPV(0:LMAX+1,NSPIN) ! loc

  double complex, dimension(:),       allocatable :: WORK1
  double complex, dimension(:),       allocatable :: WORK2
  double complex, dimension(:,:),     allocatable :: DOS
  double complex, dimension(:),       allocatable :: DOS0
  double complex, dimension(:),       allocatable :: DOS1
  double complex, dimension(:,:,:),   allocatable :: DEN0
  double precision, dimension(:,:,:), allocatable :: RHO2N1
  double precision, dimension(:,:,:), allocatable :: RHO2N2
  double precision, dimension(:,:),   allocatable :: ESPV

  !     ..
  !     .. External Functions ..
  logical:: TEST
  external TEST
  !     ..

  integer :: lmaxd
  integer :: nspind
  integer :: irmind
  integer :: lmaxd1
  integer :: lmpotd

  integer :: memory_stat
  logical :: memory_fail

  lmaxd = lmax
  nspind = nspin
  irmind = irmd - irnsd
  lmaxd1 = lmax + 1
  lmpotd = (2*LMAX+1)**2

  ! ----------- allocate work arrays ----------------------------------
  memory_fail = .false.

  allocate(WORK1(4*IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(WORK2(4*IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DOS(IEMXD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DOS0(IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DOS1(IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(DEN0(0:LMAX+1,IEMXD,NSPIN), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(RHO2N1(IRMD,LMPOTD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(RHO2N2(IRMD,LMPOTD,2), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(ESPV(0:LMAX+1,NSPIN), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "LLOYD0: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  end if
  ! -------------------------------------------------------------------

  do IE=1,IELAST
    RNORM(IE,1)=1.D0
    RNORM(IE,2)=1.D0
    DOS0(IE)       =CZERO
    DOS1(IE)       =CZERO
    do ISPIN = 1, NSPIN
      DOS(IE,ISPIN)=CZERO
    enddo
  enddo

  !=================================================================
  !==  calculate DOS  ==============================================
  !=================================================================

  call CINIT(IEMXD*(LMAXD+2)*NSPIND,DEN0)

  do I1 = 1,NAEZ
    if (MYLRANK(LMPIC) == MAPBLOCK(I1,1,NAEZ,1,0,LSIZE(LMPIC)-1)) then
      do ISPIN = 1,NSPIN
        ICELL = NTCELL(I1)

        call RHOVAL(.false.,ICST,IELAST,NSRA, &
                    ISPIN,NSPIN, &
                    EZ,WEZ,DRDI(1,I1),R(1,I1),IRMIN(I1), &
                    VINS(IRMIND,1,ISPIN),VISP(1,ISPIN), &
                    ZAT(I1),IPAN(I1),IRCUT(0,I1), &
                    THETAS(1,1,ICELL),IFUNM1(1,ICELL), &
                    LMSP1(1,ICELL), &
                    RHO2N1,RHO2N2,DEN0(0,1,ISPIN),ESPV(0,ISPIN), &
                    CLEB,LOFLM,ICLEB,IEND,JEND, &
                    GMATN, &
                    LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                    DMATLDAU, &
                    iemxd, &
                    lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb)

        do IE = 1,IELAST
          do L = 0,LMAXD1
            DOS(IE,ISPIN) = DOS(IE,ISPIN) + WEZ(IE)*DEN0(L,IE,ISPIN)
          end do
        end do

      end do

      do IE = 1,IELAST
        call RHOVAL0(EZ(IE),WEZ(IE),DRDI(1,I1),R(1,I1), &
                     IPAN(I1),IRCUT(0,I1), &
                     THETAS(1,1,ICELL), &
                     DOS0(IE),DOS1(IE), &
                     lmaxd, irmd, irid, ipand, nfund)
      end do

    endif
  end do


  !==  allreduce DOS  ==============================================

  call ZCOPY(IEMXD,DOS0,1,WORK1,1)
  call ZCOPY(IEMXD,DOS1,1,WORK1(IEMXD+1),1)
  call ZCOPY(2*IEMXD,DOS,1,WORK1(2*IEMXD+1),1)

  call MPI_ALLREDUCE(WORK1,WORK2,4*IEMXD, &
  MPI_DOUBLE_COMPLEX,MPI_SUM,LCOMM(LMPIC), &
  IERR)

  call ZCOPY(IEMXD,WORK2,1,DOS0,1)
  call ZCOPY(IEMXD,WORK2(IEMXD+1),1,DOS1,1)
  call ZCOPY(2*IEMXD,WORK2(2*IEMXD+1),1,DOS,1)

  ! ================================================================ NAEZ



  ! ======================================================================

  do ISPIN=1,NSPIN

    D0LOCINT=0.0D0
    D1LOCINT=0.0D0
    DLOYDINT=0.0D0
    do IE=1,IELAST

      DLOYD=DIMAG(WEZ(IE)*LLY_GRDT(IE,ISPIN))

      DLOC=DIMAG(DOS(IE,ISPIN))
      D0LOC=DIMAG(DOS0(IE))
      D1LOC=DIMAG(DOS1(IE))

      RNORM(IE,ISPIN) = (DLOYD+D0LOC)/DLOC

      D0LOCINT=D0LOCINT+D0LOC
      D1LOCINT=D1LOCINT+D1LOC
      DLOYDINT=DLOYDINT+DLOYD

      WEZRN(IE,ISPIN) = WEZ(IE)*RNORM(IE,ISPIN)

    enddo

    do IE=1,IELAST
      WEZRN(IE,ISPIN) = WEZRN(IE,ISPIN)* &
      (DLOYDINT+D0LOCINT-D1LOCINT)/(DLOYDINT+D0LOCINT)
      RNORM(IE,ISPIN) = RNORM(IE,ISPIN)* &
      (DLOYDINT+D0LOCINT-D1LOCINT)/(DLOYDINT+D0LOCINT)
    enddo

  enddo

  ! ----------- deallocate work arrays ----------------------------------
  deallocate(WORK1)
  deallocate(WORK2)
  deallocate(DOS)
  deallocate(DOS0)
  deallocate(DOS1)
  deallocate(DEN0)
  deallocate(RHO2N1)
  deallocate(RHO2N2)
  deallocate(ESPV)
  ! -------------------------------------------------------------------

end subroutine LLOYD0
