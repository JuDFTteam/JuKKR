#include "../DebugHelpers/test_macros.h"

!> Define the block-circulant preconditioning matrix

! SEVERE restriction of preconditioner code:
! The number of non-zero blocks must be the same in each row
! this is not the case for all crystal structures (e.g. perovskite)
! or amorphous structures

module BCPOperator_mod
  use ClusterInfo_mod, only: ClusterInfo
  implicit none
  private
  public :: BCPOperator, destroy, create, calc, multiply
  
  !> Represents the Block-Circulant preconditioning matrix
  type :: BCPOperator
    ! private
    double complex, allocatable :: GLLHBLCK(:,:)
    logical :: active = .false.
    integer :: xyzdim(3) = 1
    integer :: natbld = 1
    type(ClusterInfo), pointer :: cluster_info
    integer :: lmmaxd
  endtype

  interface create
    module procedure create_BCPOperator
  endinterface

  interface calc
    module procedure calc_BCPOperator
  endinterface

  interface multiply
    module procedure multiply_BCPOperator
  endinterface
  
  interface destroy
    module procedure destroy_BCPOperator
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  !> Setup of BCP preconditioner
  subroutine create_BCPOperator(self, natbld, xyzdim, cluster_info, lmmaxd)
    type(BCPOperator) :: self
    integer, intent(in) :: natbld, xyzdim(3)
    type(ClusterInfo), target  :: cluster_info
    integer, intent(in) :: lmmaxd

    integer :: naezd, blocks_per_row

    naezd = cluster_info%naez_trc

    CHECKASSERT(naezd == natbld*product(xyzdim(1:3)))

    allocate(self%gllhblck(natbld*lmmaxd,naezd*lmmaxd))

    blocks_per_row = cluster_info%numn0_trc(1)

    ! SEVERE restriction of preconditioner code:
    ! The number of non-zero blocks must be the same in each row
    ! this is not the case for all crystal structures (e.g. perovskite)
    CHECKASSERT(all(cluster_info%numn0_trc == blocks_per_row ))

    self%natbld = natbld
    self%xyzdim = xyzdim
    self%cluster_info => cluster_info
    self%lmmaxd = lmmaxd
    self%active = .true.

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Calculation of BCP preconditioner.
  !>
  !> uses bcpwupper
  subroutine calc_BCPOperator(self, GLLH)
    type(BCPOperator) :: self
    double complex, intent(in) :: GLLH(:)

    integer :: naezd, blocks_per_row
    
    if (.not. self%active) return
    
    naezd = self%cluster_info%naez_trc
    blocks_per_row = self%cluster_info%numn0_trc(1)

    call bcpwupper(GLLH, self%GLLHBLCK, naezd, self%cluster_info%numn0_trc, self%cluster_info%indn0_trc, &
                   self%lmmaxd, self%natbld, &
                   self%xyzdim(1), self%xyzdim(2), self%xyzdim(3), &
                   blocks_per_row)

  endsubroutine ! calc

  !----------------------------------------------------------------------------
  !> Applies Preconditioner/Operator on mat_X and returns result in mat_AX.
  subroutine multiply_BCPOperator(self, mat_X, mat_AX)
    type(BCPOperator) :: self
    double complex, intent(in)  :: mat_X(:,:)
    double complex, intent(out) :: mat_AX(:,:)

    integer :: num_columns, naez, natbld, xyzdim(3)

    natbld = self%natbld
    xyzdim(:) = self%xyzdim(1:3)

    num_columns = size(mat_X, 2)
    naez = self%cluster_info%naez_trc

    ASSERT(naez == natbld*product(xyzdim(1:3)))

    call appblckcirc(mat_X, mat_AX, self%GLLHBLCK, naez, self%lmmaxd, natbld, xyzdim(1), xyzdim(2), xyzdim(3), num_columns)

  endsubroutine ! apply

  !----------------------------------------------------------------------------
  elemental subroutine destroy_BCPOperator(self)
    type(BCPOperator), intent(inout) :: self

    if (allocated(self%GLLHBLCK)) deallocate(self%GLLHBLCK)
    nullify(self%cluster_info)
  endsubroutine ! destroy

  
  subroutine bcpwupper(gllh, gllhblck, naez, numn0, indn0, lmmaxd, natbld, xdim, ydim, zdim, blocks_per_row)
    integer, intent(in) :: naez !< total number of atoms
    integer, intent(in) :: xdim, ydim, zdim !< number of preconditioning blocks in each direction
    integer, intent(in) :: lmmaxd !< block dimension
    integer, intent(in) :: blocks_per_row !< number of non-zero blocks per row
    integer, intent(in) :: natbld !< number of atoms per preconditioning block

    integer, parameter :: nsymaxd = 48
    integer, parameter :: nintactd = 19 ! without corners of the interactions-cube, would be 27 if corners were included 

    double complex :: gllh(lmmaxd,blocks_per_row*lmmaxd,naez)
    double complex :: gllhblck(natbld*lmmaxd,natbld*xdim*ydim*zdim*lmmaxd)

    integer(kind=2), intent(in) :: indn0(blocks_per_row,naez)
    integer, intent(in) :: numn0(naez)

    
    double complex, parameter :: cone=(1.d0,0.d0), imone=(0.d0,1.d0), czero=(0.d0,0.d0)
    double complex :: blav(natbld*lmmaxd,natbld*lmmaxd,nintactd)
    double complex :: tmpblck(natbld*lmmaxd,natbld*lmmaxd)
    integer :: nblckd !< nblckd = xdim*ydim*zdim
    integer :: ipiv(natbld*lmmaxd)
    double precision :: pi, invnblckd
    integer :: i1,ix,iy,iz,ndim,ix2,iy2,iz2,cen,lef,rig,dow,upp,fow,bac,ledo,leup,leba,lefo,rido,riup,riba,rifo,doba,dofo,upba,upfo, info

    nblckd = xdim*ydim*zdim
    pi = 4.d0*atan(1.d0)
    invnblckd = 1.d0/dble(nblckd)
    
    gllhblck(:,:) = czero
    blav(:,:,:) = czero
    
    do ix = 1, xdim     ; ix2 = modulo(ix - 2, xdim)
      do iy = 1, ydim   ; iy2 = modulo(iy - 2, ydim)
        do iz = 1, zdim ; iz2 = modulo(iz - 2, zdim)

          cen = mod(ix-1,xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          lef = mod(ix2, xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          rig = mod(ix,  xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          dow = mod(ix-1,xdim) + 1 + mod(iy2, ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          upp = mod(ix-1,xdim) + 1 + mod(iy,  ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          bac = mod(ix-1,xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz2, zdim)*xdim*ydim
          fow = mod(ix-1,xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz,  zdim)*xdim*ydim

          call genblav(gllh, blav(1,1, 1), cen, cen, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 2), cen, lef, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 3), cen, rig, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 4), cen, dow, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 5), cen, upp, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 6), cen, bac, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 7), cen, fow, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)

          ledo = mod(ix2, xdim) + 1 + mod(iy2, ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          leup = mod(ix2, xdim) + 1 + mod(iy,  ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          leba = mod(ix2, xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz2, zdim)*xdim*ydim
          lefo = mod(ix2, xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz,  zdim)*xdim*ydim
          rido = mod(ix,  xdim) + 1 + mod(iy2, ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          riup = mod(ix,  xdim) + 1 + mod(iy,  ydim)*xdim + mod(iz-1,zdim)*xdim*ydim
          riba = mod(ix,  xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz2, zdim)*xdim*ydim
          rifo = mod(ix,  xdim) + 1 + mod(iy-1,ydim)*xdim + mod(iz,  zdim)*xdim*ydim
          doba = mod(ix-1,xdim) + 1 + mod(iy2, ydim)*xdim + mod(iz2, zdim)*xdim*ydim
          dofo = mod(ix-1,xdim) + 1 + mod(iy2, ydim)*xdim + mod(iz,  zdim)*xdim*ydim
          upba = mod(ix-1,xdim) + 1 + mod(iy,  ydim)*xdim + mod(iz2, zdim)*xdim*ydim
          upfo = mod(ix-1,xdim) + 1 + mod(iy,  ydim)*xdim + mod(iz,  zdim)*xdim*ydim
  
          call genblav(gllh, blav(1,1, 8), cen, ledo, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1, 9), cen, leup, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,10), cen, leba, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,11), cen, lefo, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,12), cen, rido, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,13), cen, riup, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,14), cen, riba, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,15), cen, rifo, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,16), cen, doba, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,17), cen, dofo, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,18), cen, upba, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
          call genblav(gllh, blav(1,1,19), cen, upfo, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)

        enddo ! iz
      enddo ! iy
    enddo ! ix

! normalize ...
    blav(:,:,:) = blav(:,:,:)*invnblckd

    ndim = (natbld*lmmaxd)*(natbld*lmmaxd)
!     simplified openmp code: e. rabel

!$omp parallel private (ix,iy,iz,tmpblck,info,ipiv)
!$omp do collapse(3)
    do ix = 1, xdim
      do iy = 1, ydim
        do iz = 1, zdim

          tmpblck(:,:) = czero

! .. constant times a vector plus a vector >> tmpblck = tmpblck + fac*blav ..
          call zaxpy(ndim, cone, blav(1,1, 1),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(ix-1)/xdim), blav(1,1, 2),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(ix-1)/xdim), blav(1,1, 3),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(iy-1)/ydim), blav(1,1, 4),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(iy-1)/ydim), blav(1,1, 5),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(iz-1)/zdim), blav(1,1, 6),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(iz-1)/zdim), blav(1,1, 7),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(ix-1)/xdim+2*pi*imone*(iy-1)/ydim), blav(1,1, 8),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(ix-1)/xdim-2*pi*imone*(iy-1)/ydim), blav(1,1, 9),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(ix-1)/xdim+2*pi*imone*(iz-1)/zdim), blav(1,1,10),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(ix-1)/xdim-2*pi*imone*(iz-1)/zdim), blav(1,1,11),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(ix-1)/xdim+2*pi*imone*(iy-1)/ydim), blav(1,1,12),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(ix-1)/xdim-2*pi*imone*(iy-1)/ydim), blav(1,1,13),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(ix-1)/xdim+2*pi*imone*(iz-1)/zdim), blav(1,1,14),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(ix-1)/xdim-2*pi*imone*(iz-1)/zdim), blav(1,1,15),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(iy-1)/ydim+2*pi*imone*(iz-1)/zdim), blav(1,1,16),1,tmpblck,1)
          call zaxpy(ndim, exp( 2*pi*imone*(iy-1)/ydim-2*pi*imone*(iz-1)/zdim), blav(1,1,17),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(iy-1)/ydim+2*pi*imone*(iz-1)/zdim), blav(1,1,18),1,tmpblck,1)
          call zaxpy(ndim, exp(-2*pi*imone*(iy-1)/ydim-2*pi*imone*(iz-1)/zdim), blav(1,1,19),1,tmpblck,1)
          
! .. initialize gllhblck as identity-matrix
          do i1 = 1, natbld*lmmaxd
            gllhblck(i1,natbld*lmmaxd*((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)+i1) = cone
          enddo ! i1

! ..        solve a*x = 1
! ..        gllhblck = tmpblck**-1
!           gllhblck contains inverse blocks
!           of first column of block-circulant
!           preconditioning matrix  (commented: e.r.)
          call zgesv(natbld*lmmaxd,natbld*lmmaxd,tmpblck,natbld*lmmaxd,ipiv,gllhblck(1,natbld*lmmaxd*((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)+1),natbld*lmmaxd,info)

        enddo ! iz
      enddo ! iy
    enddo ! ix
!$omp end do
!$omp end parallel

  endsubroutine ! bcpwupper

  subroutine genblav(gllh, blav, i, j, indn0, numn0, naez, lmmaxd, natbld, blocks_per_row)
    integer, intent(in) :: naez
    integer, intent(in) :: lmmaxd
    integer, intent(in) :: natbld
    integer, intent(in) :: blocks_per_row
    integer(kind=2), intent(in) :: indn0(blocks_per_row,naez)
    integer, intent(in) :: numn0(naez)
    double complex, intent(in) :: gllh(lmmaxd,lmmaxd*blocks_per_row,naez)
    double complex, intent(inout) :: blav(lmmaxd*natbld,lmmaxd*natbld)

    integer :: i, j, i1, i2, lm1, lm2, il1, il2, il2b, isrh

    do i1 = 1, natbld
      do i2 = 1, natbld

        il2b = 0

        do isrh = 1, numn0((i-1)*natbld+i1)
          if (indn0(isrh,(i-1)*natbld+i1) == (j-1)*natbld+i2) il2b = isrh
        enddo ! isrh
!
! also a pointer array pointto(naezd,naezd) could come up for this job
!          il2b = pointto((j-1)*natbld+i2,(i-1)*natbld+i1)
!
        if (il2b /= 0) then
          do lm1 = 1, lmmaxd
            il1 = lmmaxd*(i1-1)+lm1
            do lm2 = 1, lmmaxd
              il2 = lmmaxd*(i2-1)+lm2
              blav(il1,il2) = blav(il1,il2) + gllh(lm1,lmmaxd*(il2b-1)+lm2,(i-1)*natbld+i1)
            enddo ! lm2
          enddo ! lm1
        endif
        
      enddo ! i2
    enddo ! i1

  endsubroutine genblav
  
  !> apply block circulant preconditioner on vecs_in and put result into vecs_out
  subroutine appblckcirc(vecs_in, vecs_out, gllhblck, naez,lmmaxd, natbld, xdim, ydim, zdim, num_columns)
    integer, intent(in) :: naez !< number of all atoms
    integer, intent(in) :: lmmaxd !< number of atoms per preconditioning block
    integer, intent(in) :: natbld !< number of preconditioning blocks in each direction
    integer, intent(in) :: xdim, ydim, zdim
    integer, intent(in) :: num_columns
    double complex, intent(in)  :: vecs_in(naez*lmmaxd,num_columns)
    double complex, intent(out) :: vecs_out(naez*lmmaxd,num_columns)
    double complex, intent(in)  :: gllhblck(natbld*lmmaxd,xdim*ydim*zdim*natbld*lmmaxd)

    include 'fftw3.f'
    
    double complex, parameter :: cone  = (1.d0, 0.d0)
    double complex, parameter :: czero = (0.d0, 0.d0)

    !     local arrays - large, stack based!!!
    double complex :: tblck(natbld*lmmaxd,natbld*lmmaxd)
    double complex :: txk(natbld*lmmaxd)
    double complex :: tyk(natbld*lmmaxd)

    !     local arrays - large, stack based!!!
    double complex :: x(xdim,ydim,zdim)
    double complex :: xk(natbld*lmmaxd,xdim,ydim,zdim)
    double complex :: y(xdim,ydim,zdim)
    double complex :: yk(natbld*lmmaxd,xdim,ydim,zdim)

  !ibm* align(32, xk, yk)

    ! ..
    ! local scalars ..
    double precision :: fac
    integer :: lm1, lmatbl, ix, iy, iz, j, num
    integer(kind=8) :: FFTwplan_fwd, FFTwplan_bwd

    num = lmmaxd*natbld

    fac = 1.d0/dble(xdim*ydim*zdim)

    ! all threads use the same FFTw3-plans
    ! - possible according to FFTw3 documentation
    call dFFTw_plan_dft_3d(FFTwplan_bwd, xdim, ydim, zdim, x, x, FFTw_backward, FFTw_estimate)
    call dFFTw_plan_dft_3d(FFTwplan_fwd, xdim, ydim, zdim, y, y, FFTw_forward,  FFTw_estimate)

    !$omp parallel private(x,y,tblck,txk,tyk,xk,yk)
    !$omp do
    do lm1 = 1, num_columns

      !=======================================================================

      !-----------------------------------------------------------------------
      ! perform fourier-transform backward              begin
      !-----------------------------------------------------------------------

      do lmatbl = 1, num

        do iz = 1, zdim
          do iy = 1, ydim
            do ix = 1, xdim

              x(ix,iy,iz) = vecs_in(((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)*num+lmatbl,lm1)

            enddo ! ix
          enddo ! iy
        enddo ! iz

        call dFFTw_execute_dft(FFTwplan_bwd, x, x)

        do iz = 1, zdim
          do iy = 1, ydim
            do ix = 1, xdim
              xk(lmatbl,ix,iy,iz) = x(ix,iy,iz)*fac
            enddo ! ix
          enddo ! iy
        enddo ! iz

      enddo ! lmatbl

      !-----------------------------------------------------------------------
      ! perform fourier-transform backward              end
      !-----------------------------------------------------------------------

      do iz = 1, zdim
        do iy = 1, ydim
          do ix = 1, xdim

            txk(1:num) = xk(1:num,ix,iy,iz)

            do j = 1, num
              tblck(1:num,j) = gllhblck(1:num,num*((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)+j)
            enddo ! j

            call zgemv('n',num, num, cone, tblck, num, txk, 1, czero, tyk, 1)

            yk(1:num,ix,iy,iz) = tyk(1:num)

          enddo ! ix
        enddo ! iy
      enddo ! iz

      !-----------------------------------------------------------------------
      ! perform fourier-transform forward               begin
      !-----------------------------------------------------------------------

      do lmatbl = 1, num

        do iz = 1, zdim
          do iy = 1, ydim
            do ix = 1, xdim
              y(ix,iy,iz) = yk(lmatbl,ix,iy,iz)
            enddo ! ix
          enddo ! iy
        enddo ! iz

        call dFFTw_execute_dft(FFTwplan_fwd, y, y)

        do iz = 1, zdim
          do iy = 1, ydim
            do ix = 1, xdim

              vecs_out(((ix-1)+(iy-1)*xdim+(iz-1)*xdim*ydim)*num+lmatbl,lm1) = y(ix,iy,iz)

            enddo ! ix
          enddo ! iy
        enddo ! iz

      enddo ! lmatbl

    !-----------------------------------------------------------------------
    ! perform fourier-transform forward               end
    !-----------------------------------------------------------------------

    enddo ! lm1
    !$omp enddo
    !$omp endparallel

    call dFFTw_destroy_plan(FFTwplan_fwd)
    call dFFTw_destroy_plan(FFTwplan_bwd)

  endsubroutine ! appblckcirc
  
endmodule ! BCPOperator_mod
