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
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=ist); CHECKALLOC(ist)
#define DEALLOCATECHECK(X) deallocate(X, stat=ist); CHECKDEALLOC(ist)

module GauntCoefficients_mod
  implicit none
  private
  public :: GauntCoefficients, create, destroy
  public :: createGauntCoefficients, destroyGauntCoefficients ! deprecated

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
  
  !ncleb = (2*lmaxd+1)**2 * (lmaxd+1)**2

  contains

  !----------------------------------------------------------------------------
  subroutine createGauntCoefficients(coeff, lmax)
    use Harmonics_mod, only: Gaunt, Gaunt2 ! initialization of wg and yrg
    type(GauntCoefficients), intent(inout) :: coeff
    integer, intent(in) :: lmax
    
    integer :: ist
    integer :: lassld, lpot, lm2d, ncleb, lmpotd
    double precision, allocatable :: wg(:), yrg(:,:,:)

    ncleb = (lmax*2+1)**2 * (lmax+1)**2
    lpot = 2*lmax
    lassld = 4*lmax
    lm2d = (2*lmax+1)**2
    lmpotd = (lpot+1) ** 2

    coeff%lmax = lmax
    coeff%ncleb = ncleb

    ALLOCATECHECK(coeff%cleb(ncleb,2))
    ALLOCATECHECK(coeff%icleb(ncleb,3))
    ALLOCATECHECK(coeff%jend(lmpotd,0:lmax,0:lmax))
    ALLOCATECHECK(coeff%loflm(lm2d))
    
    coeff%cleb = 0.d0
    coeff%icleb = -1
    coeff%jend = -1
    coeff%loflm = -1

    ALLOCATECHECK(wg(lassld))
    ALLOCATECHECK(yrg(lassld,0:lassld,0:lassld))

    call Gaunt2(wg, yrg, lmax)
    call Gaunt(lmax, lpot, wg, yrg, coeff%cleb, coeff%loflm, coeff%icleb, coeff%iend, coeff%jend, coeff%ncleb)

    DEALLOCATECHECK(wg)
    DEALLOCATECHECK(yrg)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyGauntCoefficients(coeff)
    type(GauntCoefficients), intent(inout) :: coeff
    integer :: ist
    
    DEALLOCATECHECK(coeff%CLEB)
    DEALLOCATECHECK(coeff%ICLEB)
    DEALLOCATECHECK(coeff%JEND)
    DEALLOCATECHECK(coeff%LOFLM)
  endsubroutine ! destroy

endmodule GauntCoefficients_mod
