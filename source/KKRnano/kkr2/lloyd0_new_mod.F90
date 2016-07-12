#include "DebugHelpers/test_macros.h"

module lloyd0_new_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
implicit none
  private
  public :: lloyd0_wrapper_com

  contains

!----------------------------------------------------------------------------
!> Lloyd's formula.
!> @param[in,out] emesh Energy mesh. Renormalized energy weights are updated
subroutine lloyd0_wrapper_com(atomdata, communicator, lly_grdt, emesh, rnorm, lly, icst, nsra, gmatn, gaunts, ldau_data, method)
  use BasisAtom_mod, only: BasisAtom
  use GauntCoefficients_mod, only: GauntCoefficients
  use EnergyMesh_mod, only: EnergyMesh
  use LDAUData_mod, only: LDAUData

  integer, intent(in) :: lly  !< use lloyd 0/1
  integer, intent(in) :: icst !< num. born iterations
  integer, intent(in) :: nsra !< flag scalar relativistic
  double complex, intent(in) :: lly_grdt(:,:) ! in

  double complex, intent(in) :: gmatn(:,:,:,:) !in
  type(BasisAtom), intent(inout) :: atomdata !in or inout?
  type(GauntCoefficients), intent(in) :: gaunts
  type(EnergyMesh), intent(inout) :: emesh !inout or in?
  type(LDAUData), intent(inout) :: ldau_data ! inout?
  double precision, intent(inout) :: rnorm(:,:)  !out
  integer, intent(in) :: communicator
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

  if (lly == 1) then
    ! get wezrn and rnorm, the important input from previous
    ! calculations is lly_grdt_all

    call lloyd0_new(emesh%ez,emesh%wez,gaunts%cleb,mesh%drdi,mesh%r,mesh%irmin, &
                    atomdata%potential%vins,atomdata%potential%visp,cell%theta,atomdata%z_nuclear,gaunts%icleb, &
                    cell%ifunm,mesh%ipan,mesh%ircut,cell%lmsp, &
                    gaunts%jend,gaunts%loflm,icst,ielast,gaunts%iend,nspind,nsra, &
                    emesh%wezrn,rnorm, &
                    gmatn, &
                    lly_grdt, &
                    ldau_data%ldau,ldau_data%nldau,ldau_data%lldau,ldau_data%phildau,ldau_data%wmldau,ldau_data%dmatldau, &
                    communicator, &
                    lmaxd, mesh%irmd, irnsd, ielast, &
                    cell%irid, cell%nfund, mesh%ipand, gaunts%ncleb, method)

  else ! no lloyd

    do ie = 1, ielast
      emesh%wezrn(ie,1) = emesh%wez(ie)
      emesh%wezrn(ie,2) = emesh%wez(ie)
    enddo ! ie
    
  endif

#undef  mesh
#undef  cell
endsubroutine

!> this routine uses the previously calculated lloyd's formula terms
!> to calculate the correction to the dos and renormalized weights for
!> the energy integration.
!>
!> the routine uses communication
!> todo: split, simplify, extract communication
subroutine lloyd0_new(ez,wez,cleb,drdi,r,irmin, &
                  vins,visp,thetas,zat,icleb, &
                  ifunm1,ipan,ircut,lmsp1,jend,loflm,icst, &
                  ielast,iend,nspin,nsra, &
                  wezrn,rnorm, &                                     ! <
                  gmatn, &                                           ! >
                  lly_grdt, &                                        ! >
                  ldau,nldau,lldau,phildau,wmldau, &                 ! >
                  dmatldau, &                                        ! <
                  communicator, &                         ! >
                  lmax, irmd, irnsd, iemxd, &
                  irid, nfund, ipand, ncleb, method)
  use ValenceDensity_mod, only: rhoval, rhoval0

  !use_arraytest_mod

  integer :: lmax
  integer :: irmd
  integer :: irnsd
  integer :: iemxd
  integer :: irid
  integer :: nfund
  integer :: ipand
  integer :: ncleb

  double complex :: ez(iemxd)  ! in
  double complex :: wez(iemxd) ! in
  double precision::cleb(ncleb,2)   !in
  double precision::drdi(irmd) !in
  double precision::r(irmd)    !in
  integer::irmin !in
  double precision::vins((irmd-irnsd):irmd,(2*lmax+1)**2, 2)  !in?
  double precision::visp(irmd,2) ! in?
  double precision::thetas(irid,nfund) !in
  double precision::zat !in
  integer::icleb(ncleb,3)
  !integer::ifunm1(lmxspd,naez) (2*lpotd+1)**2
  integer::ifunm1((4*lmax+1)**2) !in?
  integer::ipan !in
  integer::ircut(0:ipand) !in
  !integer::lmsp1(lmxspd,naezd)
  integer::lmsp1((4*lmax+1)**2) !in
  integer::jend((lmax+1)**2, 0:lmax, 0:lmax) !in
  integer::loflm((2*lmax+1)**2) !in
  integer::icst
  integer::ielast
  integer::iend
  integer::nspin !in
  integer::nsra  !in
!  integer::fred  !in
  double complex :: wezrn(iemxd,2)  ! out
  double precision::rnorm(iemxd,2)  !out
  double complex :: gmatn((lmax+1)**2, (lmax+1)**2, iemxd, nspin) ! inout?
  double complex :: lly_grdt(iemxd,nspin) ! in

  ! -------------- lda+u --------------------------------------------
  !double complex :: dmatldau(mmaxd,mmaxd,nspind,lmaxd1)
  double complex :: dmatldau(2*lmax+1,2*lmax+1,nspin,lmax+1)  ! out
  double complex :: phildau(irmd,lmax+1) !in?
  double precision::wmldau(2*lmax+1,2*lmax+1,lmax+1,nspin) !in?
  integer::nldau !in
  logical::ldau  !in
  integer::lldau(lmax+1) !in?
  !------------------------------------------------------------------

  integer, intent(in) :: communicator
  integer, intent(in) :: method

!---------------- local variables -----------------------------------

  double complex, parameter :: czero=(0.d0, 0.d0)
  integer :: ie, ispin, l

  !     dynamically allocated arrays
  !double complex :: work1(4*iemxd)
  !double complex :: work2(4*iemxd)
  !double complex :: dos(iemxd,2) ! local
  !double complex :: dos0(iemxd)  ! loc
  !double complex :: dos1(iemxd)  ! loc
  !double complex :: den0(0:lmax+1,iemxd,nspin) ! loc
  !double precision::rho2n1(irmd,lmpotd,2)
  !double precision::rho2n1(irmd,(2*lmax+1)**2, 2)   ! loc
  !double precision::rho2n2(irmd,lmpotd,2)
  !double precision::rho2n2(irmd,(2*lmax+1)**2, 2)   ! loc
  !double precision::espv(0:lmax+1,nspin) ! loc

  double complex, allocatable :: dos(:,:)
  double complex, allocatable :: dos0(:)
  double complex, allocatable :: dos1(:)
  double complex, allocatable :: den0(:,:,:)
  double precision, allocatable :: rho2n1(:,:,:)
  double precision, allocatable :: rho2n2(:,:,:)
  double precision, allocatable :: espv(:,:)

  integer :: lmaxd
  integer :: nspind
  integer :: irmind
  integer :: lmaxd1
  integer :: lmpotd

  integer :: ist

  lmaxd = lmax
  nspind = nspin
  irmind = irmd - irnsd
  lmaxd1 = lmax + 1
  lmpotd = (2*lmax+1)**2

  if (ielast /= iemxd) die_here("ielast /= iemxd")

  allocate(dos(iemxd,2), dos0(iemxd), dos1(iemxd), den0(0:lmax+1,iemxd,nspin), &
           rho2n1(irmd,lmpotd,2), rho2n2(irmd,lmpotd,2), espv(0:lmax+1,nspin), stat=ist)
  if (ist /= 0) die_here("fatal error, failure to allocate memory. Probably out of memory.")

  ! initialise
  rnorm = 1.d0
  dos0 = czero
  dos1 = czero
  dos  = czero
  den0 = czero

  !=================================================================
  !==  calculate dos  ==============================================
  !=================================================================

  do ispin = 1, nspin

    call rhoval(.false.,icst,ielast,nsra, &
                ispin,nspin, &
                ez,wez,drdi,r,irmin, &
                vins(irmind,1,ispin),visp(1,ispin), &
                zat,ipan,ircut, &
                thetas,ifunm1, &
                lmsp1, &
                rho2n1,rho2n2,den0(0,1,ispin),espv(0,ispin), &
                cleb,loflm,icleb,iend,jend, &
                gmatn, &
                ldau,nldau,lldau,phildau,wmldau, &
                dmatldau, &
                iemxd, &
                lmaxd, irmd, irnsd, irid, ipand, nfund, ncleb, method)

    ! result: den0

    do ie = 1, ielast
      do l = 0, lmaxd1
        dos(ie,ispin) = dos(ie,ispin) + wez(ie)*den0(l,ie,ispin)
      enddo ! l
    enddo ! ie

  enddo ! ispin

  do ie = 1, ielast
    call rhoval0(ez(ie),wez(ie),drdi,r, ipan,ircut, thetas, dos0(ie),dos1(ie), lmaxd, irmd, irid, ipand, nfund)
  enddo ! ie

  ! communicate the dos results
  call lloyd_communicate(dos, dos0, dos1, iemxd, communicator)
  call lloyd_calcrenormalisation(dos, dos0, dos1, lly_grdt, rnorm, wez, wezrn, nspin, ielast, iemxd)

  deallocate(dos, dos0, dos1, den0, rho2n1, rho2n2, espv, stat=ist) ! deallocate work arrays 

endsubroutine ! lloyd0_new

!------------------------------------------------------------------------------
subroutine lloyd_calcrenormalisation(dos, dos0, dos1, lly_grdt, rnorm, wez, wezrn, nspin, ielast, iemxd)

  double complex, intent(in) :: dos(iemxd,2)
  double complex, intent(in) :: dos0(iemxd)
  double complex, intent(in) :: dos1(iemxd)
  double complex, intent(in) :: lly_grdt(iemxd,nspin)
  double precision, intent(out) :: rnorm(iemxd,2)
  double complex, intent(in) :: wez(iemxd)
  double complex, intent(out) :: wezrn(iemxd,2)
  integer, intent(in) :: nspin, ielast, iemxd

  integer :: ie, ispin

  double precision :: d0loc, d0locint
  double precision :: d1loc, d1locint
  double precision :: dloyd, dloydint
  double precision :: dloc
  ! ======================================================================

  do ispin = 1, nspin

    d0locint = 0.d0
    d1locint = 0.d0
    dloydint = 0.d0
    do ie = 1, ielast

      dloyd = dimag(wez(ie)*lly_grdt(ie,ispin))

      dloc  = dimag(dos(ie,ispin))
      d0loc = dimag(dos0(ie))
      d1loc = dimag(dos1(ie))

      rnorm(ie,ispin) = (dloyd + d0loc)/dloc

      d0locint = d0locint + d0loc
      d1locint = d1locint + d1loc
      dloydint = dloydint + dloyd

      wezrn(ie,ispin) = wez(ie)*rnorm(ie,ispin)

    enddo ! ie

    do ie = 1, ielast
      wezrn(ie,ispin) = wezrn(ie,ispin)*(dloydint + d0locint - d1locint)/(dloydint + d0locint)
      rnorm(ie,ispin) = rnorm(ie,ispin)*(dloydint + d0locint - d1locint)/(dloydint + d0locint)
    enddo ! ie

  enddo
endsubroutine


!------------------------------------------------------------------------------
subroutine lloyd_communicate(dos, dos0, dos1, iemxd, communicator)
  include 'mpif.h'

  double complex, intent(inout) :: dos(iemxd,2)
  double complex, intent(inout) :: dos0(iemxd)
  double complex, intent(inout) :: dos1(iemxd)
  integer, intent(in) :: iemxd
  integer, intent(in) :: communicator

  integer :: ierr
  double complex, allocatable :: work1(:), work2(:)
  double complex, parameter :: czero = (0.d0, 0.d0)

  integer :: ist

  allocate(work1(4*iemxd), work2(4*iemxd), stat=ist)
  
  if (ist /= 0) stop "lloyd0: fatal error, failure to allocate memory. probably out of memory."

  work1 = czero
  work2 = czero

  !==  allreduce dos  ==============================================

  ! call zcopy(iemxd,dos0,1,work1,1)
  ! call zcopy(iemxd,dos1,1,work1(iemxd+1),1)
  ! call zcopy(2*iemxd,dos,1,work1(2*iemxd+1),1)

  ! pack
  work1(1:iemxd) = dos0
  work1(iemxd+1:2*iemxd) = dos1
  work1(2*iemxd+1:3*iemxd) = dos(:,1)
  work1(3*iemxd+1:4*iemxd) = dos(:,2)

  call mpi_allreduce(work1, work2, 4*iemxd, mpi_double_complex, mpi_sum, communicator, ierr)

  ! unpack
  dos0 = work2(1:iemxd)
  dos1 = work2(iemxd+1:2*iemxd)
  dos(:,1) = work2(2*iemxd+1:3*iemxd)
  dos(:,2) = work2(3*iemxd+1:4*iemxd)

  ! call zcopy(iemxd,work2,1,dos0,1)
  ! call zcopy(iemxd,work2(iemxd+1),1,dos1,1)
  ! call zcopy(2*iemxd,work2(2*iemxd+1),1,dos,1)

  deallocate(work1, work2)
endsubroutine

endmodule ! lloyd
