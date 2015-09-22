!------------------------------------------------------------------------------
!> Module to save solutions obtained in self consistency step for use as
!> initial guesses in next self consistency step.
!> @author: Elias Rabel
module InitialGuess_mod
  implicit none
  private
  public :: InitialGuess, create, destroy
  public :: update_ekm, iguess_set_energy_ind, iguess_set_spin_ind, iguess_set_k_ind, iguess_load, iguess_save
  public :: iguess_init, iguess_destroy ! deprecated

  type InitialGuess
    !PRIVATE
    complex, allocatable :: PRSC(:,:,:)
    integer, allocatable :: ek_indices(:)
    integer :: ekm
    integer :: e_ind
    integer :: k_ind
    integer :: spin_ind
    integer :: iguess
  endtype

  interface create
    module procedure iguess_init
  endinterface
  
  interface destroy
    module procedure iguess_destroy
  endinterface

  contains

!------------------------------------------------------------------------------
subroutine iguess_init(self, nofks, num_spin, blocksize, iguess)
  type(InitialGuess) :: self
  integer, intent(in) :: nofks(:)
  integer, intent(in) :: num_spin
  integer, intent(in) :: blocksize
  integer, intent(in) :: iguess
  integer :: ii
  integer :: ekmd
  double complex, parameter :: czero=(0.d0, 0.d0)

  allocate(self%ek_indices(size(nofks) + 1))

  self%ek_indices(1) = 1
  do ii = 1, size(nofks)
    self%ek_indices(ii+1) = nofks(ii) + self%ek_indices(ii)
  enddo

  ekmd = self%ek_indices(size(self%ek_indices))

  ! only allocate storage if iguess flag is set
  if (iguess == 1) then
    allocate(self%prsc(blocksize, ekmd, num_spin))
    self%prsc = CZERO
  endif

  self%iguess = iguess
  self%ekm = 1
  self%k_ind = 1
  self%e_ind = 1
  self%spin_ind = 1

endsubroutine ! init

!------------------------------------------------------------------------------
subroutine iguess_set_energy_ind(self, e_ind)
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: e_ind
  self%e_ind = e_ind
  call update_ekm(self)
endsubroutine ! set

!------------------------------------------------------------------------------
subroutine iguess_set_spin_ind(self, spin_ind)
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: spin_ind
  self%spin_ind = spin_ind
endsubroutine ! set

!------------------------------------------------------------------------------
subroutine iguess_set_k_ind(self, k_ind)
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: k_ind
  self%k_ind = k_ind
  call update_ekm(self)
endsubroutine ! set

!------------------------------------------------------------------------------
subroutine iguess_load(self, x0)
  type(InitialGuess), intent(inout) :: self
  double complex, intent(out) :: x0(:,:)

  double complex, parameter :: czero=(0.d0,0.d0)
  if (self%iguess == 1) then
    call initialGuess_load_impl(x0, self%prsc(:,self%ekm,self%spin_ind))
  else
    x0 = czero
  endif
endsubroutine ! load

!------------------------------------------------------------------------------
subroutine iguess_save(self, solution)
  type(InitialGuess), intent(inout) :: self
  double complex, intent(in) :: solution(:,:)

  if (self%iguess == 1) call initialGuess_save_impl(self%prsc(:,self%ekm,self%spin_ind), solution)
endsubroutine ! save

!------------------------------------------------------------------------------
subroutine update_ekm(self)
  type(InitialGuess), intent(inout) :: self
  self%ekm = self%ek_indices(self%e_ind) + self%k_ind - 1
endsubroutine ! update

!------------------------------------------------------------------------------
subroutine iguess_destroy(self)
  type(InitialGuess), intent(inout) :: self
  if (allocated(self%prsc)) then
    deallocate(self%prsc)
  endif
endsubroutine ! destroy

!=================== private helper routines ===============================

subroutine initialGuess_load_impl(x0, prsc)
  ! ------------------------------------------------------------------------
  ! Initial guess optimization: Using the result of previous
  ! self-consistency step as starting vector for iterative solver
  ! Options:
  ! a) sc prec
  !    store solution obtained at the last self-consistency iteration
  !                                                         A. Thiess Nov'09
  ! simplified: Elias Rabel, 2012
  ! ------------------------------------------------------------------------
  complex, intent(in) :: prsc(:)
  double complex, intent(inout) :: x0(:,:)

  double complex, parameter :: czero=(0.d0,0.d0)
  integer :: lm2, site_lm_index, ind, lmmaxd

  !-----------------------------------------------------------------------

  ! translate sparse format of prsc back to x0 ..

  lmmaxd = size(x0, 2)
  x0 = czero

  do ind = 1, size(prsc)

    site_lm_index = int((ind-1)/lmmaxd) + 1
    lm2 = mod((ind-1),lmmaxd) + 1

    x0(site_lm_index,lm2) = dcmplx(real(prsc(ind)),aimag(prsc(ind)))

  enddo

endsubroutine ! load

!------------------------------------------------------------------------------
subroutine initialGuess_save_impl(prsc, gllke1)
  ! ------------------------------------------------------------------------
  ! Initial guess optimization: Using the result of previous
  ! self-consistency step as starting vector for iterative solver
  ! Options:
  ! a) sc prec
  !    store solution obtained at the last self-consistency iteration
  !                                                         A. Thiess Nov'09
  ! simplified: Elias Rabel, 2012
  ! ------------------------------------------------------------------------
  double complex, intent(in) :: gllke1(:,:)
  complex, intent(inout) :: prsc(:)

  integer :: ind, lm2, site_lm_index

  ! ================================================================
  ! fb) store new result as initial guess for the next self-consistency
  !     iteration in

  ind = 1

    !do site_index=1,naez
    !do lm1=1,lmmaxd
    do site_lm_index = 1, size(gllke1, 1)
      ! site_lm_index=lmmaxd*(site_index-1)+lm1  ! use a combined site and lm-index
      do lm2=1, size(gllke1, 2)

        ! convert to single precision
        prsc(ind) = cmplx(dreal(gllke1(site_lm_index,lm2)),dimag(gllke1(site_lm_index,lm2)))

        ind =  ind + 1

      enddo ! lm2
    enddo ! site_lm_index

! ================================================================

endsubroutine ! save


endmodule InitialGuess_mod

#ifdef TEST_INITIALGUESS_MOD
program test
  use InitialGuess_mod
  implicit none
  integer, dimension(6) :: nofks = (/1, 2, 3, 4, 5, 6 /)
  double complex, dimension(2,2) :: block
  double complex, dimension(2,2) :: block_read
  type(InitialGuess) :: ig
  integer :: e, spin, k
  integer :: ekm

  call iguess_init(ig, nofks, 2, 4, 1)
 
  ekm = 1
  do e = 1, 6
    do k = 1, nofks(e)
      do spin = 1, 2
        block = dcmplx(ekm, spin)
        call iguess_set_energy_ind(ig, e)
        call iguess_set_k_ind(ig, k)
        call iguess_set_spin_ind(ig, spin) 
        call iguess_save(ig, block)
      enddo
      ekm = ekm + 1
    enddo
  enddo

  do e = 1, 6
    do k = 1, nofks(e)
      do spin = 1, 2
        block_read = -3424
        call iguess_set_energy_ind(ig, e)
        call iguess_set_k_ind(ig, k)
        call iguess_set_spin_ind(ig, spin)
        call iguess_load(ig, block_read)
        write(*,*) block_read
        write(*,*) "-----------------------------"
      enddo
    enddo
  enddo
 
  
  call iguess_destroy(ig)

endprogram
#endif
