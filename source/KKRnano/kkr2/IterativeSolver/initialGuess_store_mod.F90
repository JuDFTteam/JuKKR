module initialGuess_store_mod

type InitialGuess
  complex, pointer, dimension(:,:,:) :: PRSC => null()
  integer :: ekmd
  integer :: nspin
  integer :: ekm
  integer :: ispin
end type

CONTAINS

!subroutine createInitialGuess(self, ekmd, nspin)
!  implicit none
!  type (InitialGuess), intent(inout) :: self
!  self%ekmd = ekmd
!  self%nspin = nspin
!  self%ekm = 1
!  self%ispin = 1
!end subroutine
!
!subroutine storeInitialGuess(self)
!  implicit none
!  type (InitialGuess), intent(inout) :: self
!
!end subroutine
!
!subroutine loadInitialGuess(self)
!  implicit none
!  type (InitialGuess), intent(in) :: self
!
!end subroutine
!
!subroutine destroyInitialGuess(self)
!  implicit none
!  type (InitialGuess), intent(inout) :: self
!
!end subroutine
!
!subroutine setEkmIndex(self, ekm)
!  implicit none
!  type (InitialGuess), intent(inout) :: self
!
!  self%ekm = ekm
!
!end subroutine
!
!subroutine setIspinIndex(self, ispin)
!  implicit none
!  type (InitialGuess), intent(inout) :: self
!
!  self%ispin = ispin
!
!end subroutine


!------------------------------------------------------------------------------
! For compatibility with old solver:
!------------------------------------------------------------------------------
subroutine initialGuess_load(X0,PRSC)

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
  !complex::         PRSC(NGUESSD*LMMAXD)
  complex, intent(in)         :: PRSC(:)

  !double complex :: X0(NAEZ*LMMAXD,LMMAXD)
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
subroutine initialGuess_save(PRSC, GLLKE1)

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

  !double complex, intent(in) :: GLLKE1(NAEZ*LMMAXD,LMMAXD)
  !complex, intent(inout)     :: PRSC(NGUESSD*LMMAXD)
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
