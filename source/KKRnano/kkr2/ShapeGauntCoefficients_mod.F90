!------------------------------------------------------------------------------
!> Module to calculate the Gaunt coefficients for Shape functions.

!> This separated datastructure was necessary since it seems that a different
!> convention for the Gaunt coefficients was used in shape function related
!> code.
!> It wraps the previously used routines gaunt and shape
!> @author Elias Rabel

! Some macros for checked allocation/deallocation
! they need an integer variable named memory_stat declared in each routine
! they are used.

#define CHECKALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Allocation error. ", __FILE__, __LINE__; STOP; endif;
#define CHECKDEALLOC(STAT) if( (STAT) /= 0) then; write(*,*) "Deallocation error. ", __FILE__, __LINE__; STOP; endif;
#define ALLOCATECHECK(X) allocate(X, stat=memory_stat); CHECKALLOC(memory_stat)
#define DEALLOCATECHECK(X) deallocate(X, stat=memory_stat); CHECKDEALLOC(memory_stat)

module ShapeGauntCoefficients_mod
  implicit none
  private
  public :: ShapeGauntCoefficients, create, destroy
  public :: createShapeGauntCoefficients, destroyShapeGauntCoefficients ! deprecated

  type ShapeGauntCoefficients
    double precision, allocatable :: GSH(:)
    integer, allocatable :: ILM(:,:)
    integer, allocatable :: IMAXSH(:)
    integer :: NGSHD
    integer :: lmax
    integer :: lmpotd
  endtype
  
  interface create
    module procedure createShapeGauntCoefficients
  endinterface
  
  interface destroy
    module procedure destroyShapeGauntCoefficients
  endinterface

  contains

  !----------------------------------------------------------------------------
  subroutine createShapeGauntCoefficients(coeff, lmax)
    type(ShapeGauntCoefficients), intent(inout) :: coeff
    integer, intent(in) :: lmax
    !---------------------------

    integer :: memory_stat
    integer :: LASSLD
    integer :: LPOT
    integer :: LMPOTD
    integer :: NGSHD
    double precision, allocatable :: WG(:) 
    double precision, allocatable :: YRG(:,:,:)

    LPOT = 2*lmax
    LASSLD = 4*lmax

    LMPOTD = (LPOT+1) ** 2

    coeff%lmax = lmax
    coeff%lmpotd = LMPOTD

    ALLOCATECHECK(WG(LASSLD))
    ALLOCATECHECK(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG,YRG,lmax)
    call SHAPEG_count(LPOT,WG,YRG,LMAX, NGSHD) ! determine number of coefficients

    ALLOCATECHECK(coeff%GSH(NGSHD))
    ALLOCATECHECK(coeff%ILM(NGSHD,3))
    ALLOCATECHECK(coeff%IMAXSH(0:LMPOTD))

    coeff%GSH = 0.0d0
    coeff%ILM = -1
    coeff%IMAXSH = -1

    call SHAPEG(LPOT, coeff%GSH, coeff%ILM, coeff%IMAXSH, WG, YRG, LMAX, NGSHD)

    coeff%NGSHD = NGSHD

    DEALLOCATECHECK(WG)
    DEALLOCATECHECK(YRG)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  subroutine destroyShapeGauntCoefficients(coeff)
    type(ShapeGauntCoefficients), intent(inout) :: coeff
    integer :: memory_stat

    DEALLOCATECHECK(coeff%GSH)
    DEALLOCATECHECK(coeff%ILM)
    DEALLOCATECHECK(coeff%IMAXSH)
  endsubroutine ! destroy

endmodule ShapeGauntCoefficients_mod
