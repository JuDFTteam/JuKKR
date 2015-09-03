!------------------------------------------------------------------------------
!> Module to calculate the Gaunt coefficients.

!> Declares a data structure GauntCoefficients that contains
!> Gaunt related data needed for a KKR calculation
!> It wraps the previously used routines gaunt and gaunt2
!> @author Elias Rabel and authors of gaunt and gaunt2

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module GauntCoefficients_mod
  implicit none
  private
  public :: GauntCoefficients, create, destroy
  public :: createGauntCoefficients, destroyGauntCoefficients ! deprecated

  type GauntCoefficients
    !> Contains the Gaunt coefficients
    !> CLEB(:,1) and CLEB(:,2) contain Gaunts but with different prefactors
    double precision, allocatable :: CLEB(:,:)
    integer, allocatable :: ICLEB(:,:)
    integer, allocatable :: JEND(:,:,:)
    integer, allocatable :: LOFLM(:) ! gives l from LM index
    integer :: IEND
    integer :: ncleb
    integer :: lmax
  endtype

  interface create
    module procedure createGauntCoefficients
  endinterface
  
  interface destroy
    module procedure destroyGauntCoefficients
  endinterface
  
  !NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2

  contains

  !----------------------------------------------------------------------------
  subroutine createGauntCoefficients(coeff, lmax)
    type(GauntCoefficients), intent(inout) :: coeff
    integer, intent(in) :: lmax
    !---------------------------
    integer :: memory_stat
    integer :: LASSLD
    integer :: LPOT
    integer :: LM2D
    integer :: NCLEB
    integer :: LMPOTD
    double precision, allocatable :: WG(:)
    double precision, allocatable :: YRG(:,:,:)

    NCLEB = (lmax*2+1)**2 * (lmax+1)**2
    LPOT = 2*lmax
    LASSLD = 4*lmax
    LM2D = (2*lmax+1)**2
    LMPOTD = (LPOT+1) ** 2

    coeff%lmax = lmax
    coeff%NCLEB = NCLEB

    ALLOCATECHECK(coeff%CLEB(NCLEB,2))
    ALLOCATECHECK(coeff%ICLEB(NCLEB,3))
    ALLOCATECHECK(coeff%JEND(LMPOTD,0:LMAX,0:LMAX))
    ALLOCATECHECK(coeff%LOFLM(LM2D))
    coeff%CLEB = 0.d0
    coeff%ICLEB = -1
    coeff%Jend = -1
    coeff%LOFLM = -1

    ALLOCATECHECK(WG(LASSLD))
    ALLOCATECHECK(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG, YRG, lmax)
    call GAUNT(lmax, LPOT, WG, YRG, coeff%CLEB, coeff%LOFLM, coeff%ICLEB, coeff%IEND, coeff%JEND, coeff%NCLEB)

    DEALLOCATECHECK(WG)
    DEALLOCATECHECK(YRG)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyGauntCoefficients(coeff)
    type(GauntCoefficients), intent(inout) :: coeff
    integer :: memory_stat
    DEALLOCATECHECK(coeff%CLEB)
    DEALLOCATECHECK(coeff%ICLEB)
    DEALLOCATECHECK(coeff%JEND)
    DEALLOCATECHECK(coeff%LOFLM)
  endsubroutine ! destroy

endmodule GauntCoefficients_mod
