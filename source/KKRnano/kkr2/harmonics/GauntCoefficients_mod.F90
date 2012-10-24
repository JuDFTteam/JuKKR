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

  type GauntCoefficients
    !> Contains the Gaunt coefficients
    !> CLEB(:,1) and CLEB(:,2) contain Gaunts but with different prefactors
    double precision, dimension(:,:), allocatable :: CLEB
    integer, dimension(:,:), allocatable :: ICLEB
    integer, dimension(:,:,:), allocatable :: JEND
    integer, dimension(:), allocatable :: LOFLM ! gives l from LM index
    integer :: IEND
    integer :: ncleb
    integer :: lmax
  end type

  !NCLEB = (LMAXD*2+1)**2 * (LMAXD+1)**2

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createGauntCoefficients(coeff, lmax)
    implicit none
    type (GauntCoefficients), intent(inout) :: coeff
    integer, intent(in) :: lmax
    !---------------------------
    integer :: memory_stat
    integer :: LASSLD
    integer :: LPOT
    integer :: LM2D
    integer :: NCLEB
    integer :: LMPOTD
    double precision, dimension(:),     allocatable :: WG
    double precision, dimension(:,:,:), allocatable :: YRG

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
    coeff%CLEB = 0.0d0
    coeff%ICLEB = -1
    coeff%JEND = -1
    coeff%LOFLM = -1

    ALLOCATECHECK(WG(LASSLD))
    ALLOCATECHECK(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG,YRG,lmax)
    call GAUNT(lmax,LPOT,WG,YRG,coeff%CLEB,coeff%LOFLM, &
               coeff%ICLEB,coeff%IEND,coeff%JEND,coeff%NCLEB)

    DEALLOCATECHECK(WG)
    DEALLOCATECHECK(YRG)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyGauntCoefficients(coeff)
    implicit none
    type (GauntCoefficients), intent(inout) :: coeff
    integer :: memory_stat
    DEALLOCATECHECK(coeff%CLEB)
    DEALLOCATECHECK(coeff%ICLEB)
    DEALLOCATECHECK(coeff%JEND)
    DEALLOCATECHECK(coeff%LOFLM)
  end subroutine

end module GauntCoefficients_mod
