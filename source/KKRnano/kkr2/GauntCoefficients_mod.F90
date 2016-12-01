!------------------------------------------------------------------------------
!> Module to calculate the Gaunt coefficients.

!> Declares a data structure GauntCoefficients that contains
!> Gaunt related data needed for a KKR calculation
!> It wraps the previously used routines gaunt and gaunt2
!> @author Elias Rabel and authors of gaunt and gaunt2

! Some macros for checked allocation/deallocation
! they need an integer variable named ist declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=ist); CHECKALLOC(ist)

module GauntCoefficients_mod
  implicit none
  private
  public :: GauntCoefficients, create, destroy

  type GauntCoefficients
    !> contains the gaunt coefficients
    !> cleb(:,1) and cleb(:,2) contain gaunts but with different prefactors
    double precision, allocatable :: cleb(:,:)
    integer, allocatable :: icleb(:,:)
    integer, allocatable :: jend(:,:,:)
    integer, allocatable :: loflm(:) ! gives l from lm index
    integer :: iend
    integer :: ncleb
    integer :: lmax
  endtype

  interface create
    module procedure createGauntCoefficients
  endinterface
  
  interface destroy
    module procedure destroyGauntCoefficients
  endinterface
  
  ! ncleb = (2*lmaxd+1)**2 * (lmaxd+1)**2

  contains

  !----------------------------------------------------------------------------
  subroutine createGauntCoefficients(self, lmax)
    use Harmonics_mod, only: Gaunt, Gaunt2 ! initialization of wg and yrg
    type(GauntCoefficients), intent(inout) :: self
    integer, intent(in) :: lmax
    
    integer :: ist
    integer :: lassld, lpot, lm2d, ncleb, lmpotd
    double precision, allocatable :: wg(:), yrg(:,:,:)

    ncleb = (lmax*2+1)**2 * (lmax+1)**2
    lpot = 2*lmax
    lassld = 4*lmax
    lm2d = (2*lmax+1)**2
    lmpotd = (lpot+1) ** 2

    self%lmax = lmax
    self%ncleb = ncleb

    ALLOCATECHECK(self%cleb(ncleb,2))
    ALLOCATECHECK(self%icleb(ncleb,3))
    ALLOCATECHECK(self%jend(lmpotd,0:lmax,0:lmax))
    ALLOCATECHECK(self%loflm(lm2d))
    
    self%cleb = 0.d0
    self%icleb = -1
    self%jend = -1
    self%loflm = -1

    ALLOCATECHECK(wg(lassld))
    ALLOCATECHECK(yrg(lassld,0:lassld,0:lassld))

    call Gaunt2(wg, yrg, lmax)
    call Gaunt(lmax, lpot, wg, yrg, self%cleb, self%loflm, self%icleb, self%iend, self%jend, self%ncleb)

    deallocate(wg, yrg, stat=ist)
  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyGauntCoefficients(self)
    type(GauntCoefficients), intent(inout) :: self
    integer :: ist ! ignore status
    deallocate(self%cleb, self%icleb, self%jend, self%loflm, stat=ist)
  endsubroutine ! destroy

endmodule ! GauntCoefficients_mod
