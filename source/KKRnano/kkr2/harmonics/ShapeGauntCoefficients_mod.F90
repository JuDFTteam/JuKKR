!------------------------------------------------------------------------------
!> Module to calculate the Gaunt coefficients for Shape functions
!> This separated datastructure was necessary since it seems that a different
!> convention for the Gaunt coefficients was used in shape function related
!> code.
!> It wraps the previously used routines gaunt and shape
!> @author: Elias Rabel

module ShapeGauntCoefficients_mod
  implicit none

  type ShapeGauntCoefficients
    double precision, dimension(:), allocatable :: GSH
    integer, dimension(:,:), allocatable :: ILM
    integer, dimension(:),   allocatable :: IMAXSH
    integer :: NGSHD
    integer :: lmax
    integer :: lmpotd
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createShapeGauntCoefficients(coeff, lmax)
    implicit none
    type (ShapeGauntCoefficients), intent(inout) :: coeff
    integer, intent(in) :: lmax
    !---------------------------

    integer :: LASSLD
    integer :: LPOT
    integer :: LMPOTD
    integer :: NGSHD
    double precision, dimension(:),     allocatable :: WG
    double precision, dimension(:,:,:), allocatable :: YRG

    LPOT = 2*lmax
    LASSLD = 4*lmax

    LMPOTD = (LPOT+1) ** 2

    coeff%lmax = lmax
    coeff%lmpotd = LMPOTD

    allocate(WG(LASSLD))
    allocate(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG,YRG,lmax)
    call SHAPEG_count(LPOT,WG,YRG,LMAX, NGSHD)

    allocate(coeff%GSH(NGSHD))
    allocate(coeff%ILM(NGSHD,3))
    allocate(coeff%IMAXSH(0:LMPOTD))
    coeff%GSH = 0.0d0
    coeff%ILM = -1
    coeff%IMAXSH = -1

    call SHAPEG(LPOT,coeff%GSH,coeff%ILM,coeff%IMAXSH,WG,YRG,LMAX, NGSHD)

    coeff%NGSHD = NGSHD

    deallocate(WG)
    deallocate(YRG)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyShapeGauntCoefficients(coeff)
    implicit none
    type (ShapeGauntCoefficients), intent(inout) :: coeff

    deallocate(coeff%GSH)
    deallocate(coeff%ILM)
    deallocate(coeff%IMAXSH)
  end subroutine

end module ShapeGauntCoefficients_mod
