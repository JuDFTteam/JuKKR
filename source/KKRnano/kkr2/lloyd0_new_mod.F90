#include "DebugHelpers/test_macros.h"

module lloyd0_new_mod
implicit none
  private
  public :: lloyd0_wrapper_com

  contains

!----------------------------------------------------------------------------
!> Lloyd's formula.
!> @param[in,out] emesh Energy mesh. Renormalized energy weights are updated
subroutine lloyd0_wrapper_com(atomdata, mpi_comm, LLY_GRDT, emesh, RNORM, LLY, ICST, NSRA, GMATN, gaunts, ldau_data, method)
  use BasisAtom_mod, only: BasisAtom
  use GauntCoefficients_mod, only: GauntCoefficients
  use EnergyMesh_mod, only: EnergyMesh
  use LDAUData_mod, only: LDAUData

  integer, intent(in) :: lly  !< use lloyd 0/1
  integer, intent(in) :: icst !< num. born iterations
  integer, intent(in) :: nsra !< flag scalar relativistic
  double complex, intent(in) :: lly_grdt(:,:) ! in

  double complex, intent(in) :: gmatn(:,:,:,:) !in
  type(basisatom), intent(inout) :: atomdata !in or inout?
  type(gauntcoefficients), intent(in) :: gaunts
  type(energymesh), intent(inout) :: emesh !inout or in?
  type(ldaudata), intent(inout) :: ldau_data ! inout?
  double precision, intent(inout) :: rnorm(:,:)  !out
  integer, intent(in) :: mpi_comm
  integer, intent(in) :: method !< single site method

  integer :: nspind, irmind, irnsd, lmaxd, ielast, ie

  nspind = atomdata%nspin

#define  mesh atomdata%mesh_ptr
#define  cell atomdata%cell_ptr
  
  irmind = atomdata%potential%irmind
  irnsd = atomdata%potential%irmd - atomdata%potential%irmind

  lmaxd = atomdata%potential%lpot / 2

  CHECKASSERT (lmaxd*2 == atomdata%potential%lpot)

  ielast = emesh%ielast

  if (lly==1) then
    ! get wezrn and rnorm, the important input from previous
    ! calculations is lly_grdt_all

    call lloyd0_new(emesh%ez,emesh%wez,gaunts%cleb,mesh%drdi,mesh%r,mesh%irmin, &
                    atomdata%potential%vins,atomdata%potential%visp,cell%shdata%theta,atomdata%z_nuclear,gaunts%icleb, &
                    cell%shdata%ifunm,mesh%ipan,mesh%ircut,cell%shdata%lmsp, &
                    gaunts%jend,gaunts%loflm,icst,ielast,gaunts%iend,nspind,nsra, &
                    emesh%wezrn,rnorm, &
                    gmatn, &
                    lly_grdt, &
                    ldau_data%ldau,ldau_data%nldau,ldau_data%lldau,ldau_data%phildau,ldau_data%wmldau,ldau_data%dmatldau, &
                    mpi_comm, &
                    lmaxd, mesh%irmd, irnsd, ielast, &
                    cell%shdata%irid, cell%shdata%nfund, mesh%ipand, gaunts%ncleb, method)

  else ! no lloyd

    do ie=1,ielast
      emesh%wezrn(ie,1) = emesh%wez(ie)
      emesh%wezrn(ie,2) = emesh%wez(ie)
    enddo
  endif

#undef  mesh
#undef  cell

endsubroutine

!> This routine uses the previously calculated Lloyd's formula terms
!> to calculate the correction to the DOS and renormalized weights for
!> the energy integration.
!>
!> the routine uses communication
!> TODO: split, simplify, extract communication
subroutine LLOYD0_NEW(EZ,WEZ,CLEB,DRDI,R,IRMIN, &
                  VINS,VISP,THETAS,ZAT,ICLEB, &
                  IFUNM1,IPAN,IRCUT,LMSP1,JEND,LOFLM,ICST, &
                  IELAST,IEND,NSPIN,NSRA, &
                  WEZRN,RNORM, &                                     ! <
                  GMATN, &                                           ! >
                  LLY_GRDT, &                                        ! >
                  LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &                 ! >
                  DMATLDAU, &                                        ! <
                  communicator, &                         ! >
                  lmax, irmd, irnsd, iemxd, &
                  irid, nfund, ipand, ncleb, method)
  use ValenceDensity_mod, only: RHOVAL, RHOVAL0

  !USE_ARRAYTEST_MOD

  integer :: lmax
  integer :: irmd
  integer :: irnsd
  integer :: iemxd
  integer :: irid
  integer :: nfund
  integer :: ipand
  integer :: ncleb

  double complex :: EZ(IEMXD)  ! in
  double complex :: WEZ(IEMXD) ! in
  double precision::CLEB(NCLEB,2)   !in
  double precision::DRDI(IRMD) !in
  double precision::R(IRMD)    !in
  integer::IRMIN !in
  double precision::VINS((IRMD-IRNSD):IRMD,(2*LMAX+1)**2, 2)  !in?
  double precision::VISP(IRMD,2) ! in?
  double precision::THETAS(IRID,NFUND) !in
  double precision::ZAT !in
  integer::ICLEB(NCLEB,3)
  !integer::IFUNM1(LMXSPD,NAEZ) (2*LPOTD+1)**2
  integer::IFUNM1((4*LMAX+1)**2) !in?
  integer::IPAN !in
  integer::IRCUT(0:IPAND) !in
  !integer::LMSP1(LMXSPD,NAEZD)
  integer::LMSP1((4*LMAX+1)**2) !in
  integer::JEND((LMAX+1)**2, 0:LMAX, 0:LMAX) !in
  integer::LOFLM((2*LMAX+1)**2) !in
  integer::ICST
  integer::IELAST
  integer::IEND
  integer::NSPIN !in
  integer::NSRA  !in
!  integer::FRED  !in
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

  integer, intent(in) :: communicator
  integer, intent(in) :: method

!---------------- Local variables -----------------------------------

!  !     .. Parameters ..

  DOUBLE COMPLEX      CZERO
  PARAMETER          (CZERO=(0.0D0,0.0D0))
  !     ..
  !     .. Scalars ..
  integer::IE

  integer::ISPIN
  integer::L
  !     ..

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

  double complex, allocatable :: DOS(:,:)
  double complex, allocatable :: DOS0(:)
  double complex, allocatable :: DOS1(:)
  double complex, allocatable :: DEN0(:,:,:)
  double precision, allocatable :: RHO2N1(:,:,:)
  double precision, allocatable :: RHO2N2(:,:,:)
  double precision, allocatable :: ESPV(:,:)

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

  if (ielast /= iemxd) then !assertion
    write(*,*) "lloyd0_new: ielast /= iemxd"
    stop
  endif

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
  endif
  ! -------------------------------------------------------------------

  ! initialise
  RNORM = 1.D0
  DOS0 = CZERO
  DOS1 = CZERO
  DOS  = CZERO
  DEN0 = CZERO

  !=================================================================
  !==  calculate DOS  ==============================================
  !=================================================================

      do ISPIN = 1,NSPIN

        call RHOVAL(.false.,ICST,IELAST,NSRA, &
                    ISPIN,NSPIN, &
                    EZ,WEZ,DRDI,R,IRMIN, &
                    VINS(IRMIND,1,ISPIN),VISP(1,ISPIN), &
                    ZAT,IPAN,IRCUT, &
                    THETAS,IFUNM1, &
                    LMSP1, &
                    RHO2N1,RHO2N2,DEN0(0,1,ISPIN),ESPV(0,ISPIN), &
                    CLEB,LOFLM,ICLEB,IEND,JEND, &
                    GMATN, &
                    LDAU,NLDAU,LLDAU,PHILDAU,WMLDAU, &
                    DMATLDAU, &
                    iemxd, &
                    lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb, method)

        ! result: DEN0

        do IE = 1,IELAST
          do L = 0,LMAXD1
            DOS(IE,ISPIN) = DOS(IE,ISPIN) + WEZ(IE)*DEN0(L,IE,ISPIN)
          enddo
        enddo

      enddo

      do IE = 1,IELAST
        call RHOVAL0(EZ(IE),WEZ(IE),DRDI,R, &
                     IPAN,IRCUT, &
                     THETAS, &
                     DOS0(IE),DOS1(IE), &
                     lmaxd, irmd, irid, ipand, nfund)
      enddo

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

endsubroutine LLOYD0_NEW

!------------------------------------------------------------------------------
subroutine lloyd_calcRenormalisation(DOS, DOS0, DOS1, LLY_GRDT, RNORM, WEZ, WEZRN, NSPIN, IELAST, iemxd)
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
endsubroutine


!------------------------------------------------------------------------------
subroutine lloyd_communicate(DOS, DOS0, DOS1, iemxd, communicator)
  include 'mpif.h'

  double complex, intent(inout) :: DOS(IEMXD,2)
  double complex, intent(inout) :: DOS0(IEMXD)
  double complex, intent(inout) :: DOS1(IEMXD)
  integer, intent(in) :: iemxd
  integer, intent(in) :: communicator

  integer :: IERR
  double complex, allocatable :: WORK1(:), WORK2(:)
  double complex, parameter :: CZERO = (0.0d0, 0.0d0)

  integer :: memory_stat

  allocate(WORK1(4*IEMXD), stat = memory_stat)
  
  if (memory_stat /= 0) then
    write(*,*) "LLOYD0: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  endif

  allocate(WORK2(4*IEMXD), stat = memory_stat)

  if (memory_stat /= 0) then
    write(*,*) "LLOYD0: FATAL Error, failure to allocate memory."
    write(*,*) "        Probably out of memory."
    stop
  endif

  WORK1 = CZERO
  WORK2 = CZERO

  !==  allreduce DOS  ==============================================

  !call ZCOPY(IEMXD,DOS0,1,WORK1,1)
  !call ZCOPY(IEMXD,DOS1,1,WORK1(IEMXD+1),1)
  !call ZCOPY(2*IEMXD,DOS,1,WORK1(2*IEMXD+1),1)

  WORK1(1:IEMXD) = DOS0
  WORK1(IEMXD+1:2*IEMXD) = DOS1
  WORK1(2*IEMXD+1:3*IEMXD) = DOS(:,1)
  WORK1(3*IEMXD+1:4*IEMXD) = DOS(:,2)

  call MPI_ALLREDUCE(WORK1,WORK2,4*IEMXD, &
  MPI_DOUBLE_COMPLEX,MPI_SUM,communicator, &
  IERR)

  DOS0 = WORK2(1:IEMXD)
  DOS1 = WORK2(IEMXD+1:2*IEMXD)
  DOS(:,1) = WORK2(2*IEMXD+1:3*IEMXD)
  DOS(:,2) = WORK2(3*IEMXD+1:4*IEMXD)

  !call ZCOPY(IEMXD,WORK2,1,DOS0,1)
  !call ZCOPY(IEMXD,WORK2(IEMXD+1),1,DOS1,1)
  !call ZCOPY(2*IEMXD,WORK2(2*IEMXD+1),1,DOS,1)

  deallocate(WORK1)
  deallocate(WORK2)
endsubroutine

endmodule
