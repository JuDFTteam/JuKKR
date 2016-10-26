!------------------------------------------------------------------------------
!> Module to save solutions obtained in self consistency step for use as
!> initial guesses in next self consistency step.
!> @author: Elias Rabel, Paul F. Baumeister
module InitialGuess_mod
  implicit none
  private
  public :: InitialGuess, create, destroy, load, store

  type InitialGuess
    integer :: prec ! precision: 0:no saving, 1:single precision, 2:double complex
    integer, allocatable :: ek_offsets(:)
           complex, allocatable :: prsc(:,:,:) ! if prec == 1
    double complex, allocatable :: prsz(:,:,:) ! if prec == 2
  endtype ! InitialGuess

  interface create
    module procedure createInitialGuess
  endinterface
  
  interface destroy
    module procedure destroyInitialGuess
  endinterface

  interface load
    module procedure iguess_load
  endinterface

  interface store
    module procedure iguess_save
  endinterface

  
  contains

  subroutine createInitialGuess(self, nofks, nspin, datasize, prec)
    type(InitialGuess) :: self
    integer, intent(in) :: nofks(:)
    integer, intent(in) :: nspin ! number of collinear spins, must bin in [1,2]
    integer, intent(in) :: datasize
    integer, intent(in) :: prec ! must be in [0,1,2]
    
    integer :: ik, nk, ekmd, ist

    nk = size(nofks)
    allocate(self%ek_offsets(nk+1))

    self%ek_offsets(1) = 0
    do ik = 1, nk
      self%ek_offsets(ik+1) = self%ek_offsets(ik) + nofks(ik) 
    enddo ! ik
    ekmd = self%ek_offsets(nk+1) ! == sum of all

    ! only allocate storage if prec flag is set
    self%prec = 0
    if (prec == 1) then
      allocate(self%prsc(datasize,ekmd,nspin), stat=ist)
      self%prsc = 0.
      self%prec = 1
    else if(prec == 2) then
      allocate(self%prsz(datasize,ekmd,nspin), stat=ist)
      self%prsz = 0.d0
      self%prec = 2
    endif

  endsubroutine ! create

#define EKM  (self%ek_offsets(ie) + ik)

  subroutine iguess_load(self, startval, ik, is, ie)
    type(InitialGuess), intent(inout) :: self
    double complex, intent(out) :: startval(:,:)
    integer, intent(in) :: ik, is, ie
    
    if (self%prec == 1) then
      startval = reshape(dcmplx(self%prsc(:,EKM,is)), shape(startval))
    else if (self%prec == 2) then
      startval = reshape(       self%prsz(:,EKM,is) , shape(startval))
    else
      startval = 0.d0
    endif
  endsubroutine ! load
  
  subroutine iguess_save(self, solution, ik, is, ie)
    type(InitialGuess), intent(inout) :: self
    double complex, intent(in) :: solution(:,:)
    integer, intent(in) :: ik, is, ie

    if (self%prec == 1) then
      self%prsc(:,EKM,is) = reshape(solution, [size(solution)])
    else if (self%prec == 2) then
      self%prsz(:,EKM,is) = reshape(solution, [size(solution)])
    endif
  endsubroutine ! save

  subroutine destroyInitialGuess(self)
    type(InitialGuess), intent(inout) :: self
    
    integer :: ist
    deallocate(self%prsc, self%prsz, self%ek_offsets, stat=ist) ! ignore status
  endsubroutine ! destroy

endmodule ! InitialGuess_mod

#ifdef TEST_INITIALGUESS_MOD
program test
  use InitialGuess_mod
  implicit none
  double complex :: block(2,2,1), block_read(2,2,1)
  type(InitialGuess) :: ig
  integer :: ie, is, ik, nofks(6) = [1, 2, 3, 4, 5, 6]!, ekm

  call create(ig, nofks, 2, 4, prec=2) ! igues=1 means storage in single precision complex, 2:double complex

  do ie = 1, 6
    do ik = 1, nofks(ie)
      do is = 1, 2
        block = dcmplx(ik+10000.*ie, is)
        call store(ig, block, ik, is, ie)
      enddo ! is
    enddo ! ik
  enddo ! ie

  do ie = 1, 6
    do ik = 1, nofks(ie)
      do is = 1, 2
        block_read = -3424
        call load(ig, block_read, ik, is, ie)
        write(*,*) block_read
        write(*,*) "-----------------------------"
      enddo ! is
    enddo ! ik
  enddo ! ie

  call destroy(ig)

endprogram
#endif
