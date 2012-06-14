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
                  MYLRANK, &                                         ! >
                  communicator, comm_size, &                         ! >
                  !   new input parameters after inc.p removal
                  lmax, irmd, irnsd, iemxd, &
                  irid, nfund, ncelld, ipand, ncleb)
  implicit none

  integer :: lmax
  integer :: irmd
  integer :: irnsd
  integer :: iemxd
  integer :: irid
  integer :: nfund
  integer :: ncelld
  integer :: ipand
  integer :: ncleb

  double complex :: EZ(IEMXD)  ! in
  double complex :: WEZ(IEMXD) ! in
  double precision::CLEB(NCLEB,2)   !in
  double precision::DRDI(IRMD,NAEZ) !in
  double precision::R(IRMD,NAEZ)    !in
  integer::IRMIN(NAEZ) !in
  double precision::VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2, 2)  !in?
  double precision::VISP(IRMD,2) ! in?
  double precision::THETAS(IRID,NFUND,NCELLD) !in
  double precision::ZAT(NAEZ) !in
  integer::ICLEB(NCLEB,3)
  !integer::IFUNM1(LMXSPD,NAEZ) (2*LPOTD+1)**2
  integer::IFUNM1((4*LMAX+1)**2,NAEZ) !in?
  integer::IPAN(NAEZ) !in
  integer::IRCUT(0:IPAND,NAEZ) !in
  !integer::LMSP1(LMXSPD,NAEZD)
  integer::LMSP1((4*LMAX+1)**2, NAEZ) !in
  integer::JEND((LMAX+1)**2, 0:LMAX, 0:LMAX) !in
  integer::LOFLM((2*LMAX+1)**2) !in
  integer::NTCELL(NAEZ) !in
  integer::ICST
  integer::IELAST
  integer::IEND
  integer::NAEZ  !in
  integer::NSPIN !in
  integer::NSRA  !in
  double complex :: WEZRN(IEMXD,2)  ! out
  double precision::RNORM(IEMXD,2)  !out
  double complex :: GMATN((LMAX+1)**2, (LMAX+1)**2, IEMXD, NSPIN) ! inout?
  double complex :: LLY_GRDT(IEMXD,NSPIN) ! in

  ! -------------- LDA+U --------------------------------------------
  !double complex :: DMATLDAU(MMAXD,MMAXD,NSPIND,LMAXD1)
  double complex :: DMATLDAU(2*LMAX+1,2*LMAX+1,NSPIN,LMAX+1)  ! out
  double complex :: PHILDAU(IRMD,LMAX+1) !in?
  double precision::WMLDAU(2*LMAX+1,2*LMAX+1,LMAX+1,NSPIN) !in?
  integer::NLDAU !in
  logical::LDAU  !in
  integer::LLDAU(LMAX+1) !in?
  !------------------------------------------------------------------

  integer, intent(in) :: MYLRANK
  integer, intent(in) :: communicator
  integer, intent(in) :: comm_size

!---------------- Local variables -----------------------------------

!  !     .. Parameters ..

  DOUBLE COMPLEX      CZERO
  PARAMETER          (CZERO=(0.0D0,0.0D0))
  !     ..
  !     .. Scalars ..
  integer::I1
  integer::ICELL
  integer::IE

  integer::ISPIN
  integer::L
  !     ..
  !     .. MPI ..
  !     .. N-MPI
  integer::MAPBLOCK

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

  double complex, dimension(:,:),     allocatable :: DOS
  double complex, dimension(:),       allocatable :: DOS0
  double complex, dimension(:),       allocatable :: DOS1
  double complex, dimension(:,:,:),   allocatable :: DEN0
  double precision, dimension(:,:,:), allocatable :: RHO2N1
  double precision, dimension(:,:,:), allocatable :: RHO2N2
  double precision, dimension(:,:),   allocatable :: ESPV

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

  ! TODO: get rid of this "loop" somehow - we are already atom-parallel
  do I1 = 1,NAEZ
    if (MYLRANK == MAPBLOCK(I1,1,NAEZ,1,0,comm_size-1)) then
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

        ! result: DEN0

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

  ! communicate the DOS results
  call lloyd_communicate(DOS, DOS0, DOS1, iemxd, communicator)
  call lloyd_calcRenormalisation(DOS, DOS0, DOS1, LLY_GRDT, RNORM, WEZ, WEZRN, NSPIN, IELAST, iemxd)

  ! ----------- deallocate work arrays ----------------------------------
  deallocate(DOS)
  deallocate(DOS0)
  deallocate(DOS1)
  deallocate(DEN0)
  deallocate(RHO2N1)
  deallocate(RHO2N2)
  deallocate(ESPV)
  ! -------------------------------------------------------------------

end subroutine LLOYD0

!------------------------------------------------------------------------------
subroutine lloyd_calcRenormalisation(DOS, DOS0, DOS1, LLY_GRDT, RNORM, WEZ, WEZRN, NSPIN, IELAST, iemxd)
  implicit none
  integer :: iemxd
  integer :: NSPIN

  double complex :: DOS(IEMXD,2)
  double complex :: DOS0(IEMXD)
  double complex :: DOS1(IEMXD)

  doublecomplex :: LLY_GRDT(IEMXD,NSPIN)
  double precision :: RNORM(IEMXD,2)
  double complex :: WEZ(IEMXD)
  double complex :: WEZRN(IEMXD,2)

  integer :: IE
  integer :: IELAST
  integer :: ISPIN

  double precision :: D0LOC
  double precision :: D0LOCINT
  double precision :: D1LOC
  double precision :: D1LOCINT
  double precision :: DLOC
  double precision :: DLOYD
  double precision :: DLOYDINT
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
end subroutine


!------------------------------------------------------------------------------
subroutine lloyd_communicate(DOS, DOS0, DOS1, iemxd, communicator)
  implicit none
  include 'mpif.h'

  double complex, intent(inout) :: DOS(IEMXD,2)
  double complex, intent(inout) :: DOS0(IEMXD)
  double complex, intent(inout) :: DOS1(IEMXD)
  integer, intent(in) :: iemxd
  integer, intent(in) :: communicator

  integer :: IERR
  double complex, dimension(:), allocatable :: WORK1
  double complex, dimension(:), allocatable :: WORK2

  integer :: memory_stat
  logical :: memory_fail

  memory_fail = .false.

  allocate(WORK1(4*IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.
  allocate(WORK2(4*IEMXD), stat = memory_stat)
  if (memory_stat /= 0) memory_fail = .true.

  if (memory_fail .eqv. .true.) then
    write(*,*) "LLOYD0: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  end if

  !==  allreduce DOS  ==============================================

  call ZCOPY(IEMXD,DOS0,1,WORK1,1)
  call ZCOPY(IEMXD,DOS1,1,WORK1(IEMXD+1),1)
  call ZCOPY(2*IEMXD,DOS,1,WORK1(2*IEMXD+1),1)

  call MPI_ALLREDUCE(WORK1,WORK2,4*IEMXD, &
  MPI_DOUBLE_COMPLEX,MPI_SUM,communicator, &
  IERR)

  call ZCOPY(IEMXD,WORK2,1,DOS0,1)
  call ZCOPY(IEMXD,WORK2(IEMXD+1),1,DOS1,1)
  call ZCOPY(2*IEMXD,WORK2(2*IEMXD+1),1,DOS,1)

  deallocate(WORK1)
  deallocate(WORK2)
end subroutine
