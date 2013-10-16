!------------------------------------------------------------------------------
!> Module to save solutions obtained in self consistency step for use as
!> initial guesses in next self consistency step.
!> @author: Elias Rabel
module InitialGuess_mod

type InitialGuess
  !PRIVATE
  complex, allocatable, dimension(:,:,:) :: PRSC
  integer, allocatable :: ek_indices(:)
  integer :: ekm
  integer :: e_ind
  integer :: k_ind
  integer :: spin_ind
  integer :: iguess
end type

private :: update_ekm
private :: initialGuess_load_impl
private :: initialGuess_save_impl

CONTAINS

!------------------------------------------------------------------------------
subroutine iguess_init(self, nofks, num_spin, blocksize, iguess)
  implicit none
  type(InitialGuess) :: self
  integer, intent(in) :: nofks(:)
  integer, intent(in) :: num_spin
  integer, intent(in) :: blocksize
  integer, intent(in) :: iguess
  integer :: ii
  integer :: ekmd
  double complex, parameter :: CZERO = (0.0d0, 0.0d0)

  allocate(self%ek_indices(size(nofks) + 1))

  self%ek_indices(1) = 1
  do ii = 1, size(nofks)
    self%ek_indices(ii+1) = nofks(ii) + self%ek_indices(ii)
  end do

  ekmd = self%ek_indices(size(self%ek_indices))

  ! only allocate storage if iguess flag is set
  if (iguess == 1) then
    allocate(self%prsc(blocksize, ekmd, num_spin))
    self%prsc = CZERO
  end if

  self%iguess = iguess
  self%ekm = 1
  self%k_ind = 1
  self%e_ind = 1
  self%spin_ind = 1

end subroutine

!------------------------------------------------------------------------------
subroutine iguess_set_energy_ind(self, e_ind)
  implicit none
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: e_ind
  self%e_ind = e_ind
  call update_ekm(self)
end subroutine

!------------------------------------------------------------------------------
subroutine iguess_set_spin_ind(self, spin_ind)
  implicit none
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: spin_ind
  self%spin_ind = spin_ind
end subroutine

!------------------------------------------------------------------------------
subroutine iguess_set_k_ind(self, k_ind)
  implicit none
  type(InitialGuess), intent(inout) :: self
  integer, intent(in) :: k_ind
  self%k_ind = k_ind
  call update_ekm(self)
end subroutine

!------------------------------------------------------------------------------
subroutine iguess_load(self, x0)
  implicit none
  type(InitialGuess), intent(inout) :: self
  double complex, intent(out) :: x0(:,:)

  double complex, parameter :: CZERO = (0.0d0, 0.0d0)
  if (self%iguess == 1) then
    call initialGuess_load_impl(x0, self%prsc(:, self%ekm, self%spin_ind))
  else
    x0 = CZERO
  end if
end subroutine

!------------------------------------------------------------------------------
subroutine iguess_save(self, solution)
  implicit none
  type(InitialGuess), intent(inout) :: self
  double complex, intent(in) :: solution(:,:)

  if (self%iguess == 1) then
    call initialGuess_save_impl( &
         self%prsc(:, self%ekm, self%spin_ind), solution)
  end if
end subroutine

!------------------------------------------------------------------------------
subroutine update_ekm(self)
  implicit none
  type(InitialGuess), intent(inout) :: self
  self%ekm = self%ek_indices(self%e_ind) + self%k_ind - 1
end subroutine

!------------------------------------------------------------------------------
subroutine iguess_destroy(self)
  implicit none
  type(InitialGuess), intent(inout) :: self
  if (allocated(self%prsc)) then
    deallocate(self%prsc)
  end if
end subroutine

!=================== private helper routines ===============================

subroutine initialGuess_load_impl(X0,PRSC)

  implicit none
  ! ------------------------------------------------------------------------
  ! Initial guess optimization: Using the result of previous
  ! self-consistency step as starting vector for iterative solver
  ! Options:
  ! a) sc prec
  !    store solution obtained at the last self-consistency iteration
  !                                                         A. Thiess Nov'09
  ! simplified: Elias Rabel, 2012
  ! ------------------------------------------------------------------------


  double complex :: CZERO
  parameter        (CZERO = ( 0.0D0,0.0D0 ))

  !     .. ARRAY ARGUMENTS ..
  complex, intent(in)         :: PRSC(:)

  double complex, intent(inout)     :: X0(:,:)

  !     .. LOCAL SCALARS ..
  integer::LM2
  integer::site_lm_index

  integer::ind
  integer :: lmmaxd

  !-----------------------------------------------------------------------

  ! translate sparse format of PRSC back to X0 ..

  lmmaxd = size(X0, 2)
  X0 = CZERO

  do ind = 1, size(PRSC)

    site_lm_index = INT((ind-1)/LMMAXD) + 1
    LM2 = MOD((ind-1),LMMAXD) + 1

    X0(site_lm_index,LM2) = &
    DCMPLX(REAL(PRSC(ind)),AIMAG(PRSC(ind)))

  enddo

end subroutine

!------------------------------------------------------------------------------
subroutine initialGuess_save_impl(PRSC, GLLKE1)

  implicit none
  ! ------------------------------------------------------------------------
  ! Initial guess optimization: Using the result of previous
  ! self-consistency step as starting vector for iterative solver
  ! Options:
  ! a) sc prec
  !    store solution obtained at the last self-consistency iteration
  !                                                         A. Thiess Nov'09
  ! simplified: Elias Rabel, 2012
  ! ------------------------------------------------------------------------

  !     .. ARRAY ARGUMENTS ..

  double complex, intent(in) :: GLLKE1(:,:)
  complex, intent(inout)     :: PRSC(:)

  !     .. LOCAL SCALARS ..
  integer::ind
  integer::LM2
  integer::site_lm_index

  ! ================================================================
  ! Fb) store new result as initial guess for the next self-consistency
  !     iteration in

  ind = 1

    !do site_index=1,NAEZ
    !do LM1=1,LMMAXD
    do site_lm_index = 1, size(GLLKE1, 1)
      ! site_lm_index=LMMAXD*(site_index-1)+LM1  ! use a combined site and lm-index
      do LM2=1, size(GLLKE1, 2)

        ! Convert to single precision
        PRSC(ind) = &
        CMPLX(DREAL(GLLKE1(site_lm_index,LM2)),DIMAG(GLLKE1(site_lm_index,LM2)))

        ind =  ind + 1

      enddo
    enddo

! ================================================================

end subroutine


end module

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
      end do
      ekm = ekm + 1
    end do
  end do

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
      end do
    end do
  end do
 
  
  call iguess_destroy(ig)

end program
#endif
