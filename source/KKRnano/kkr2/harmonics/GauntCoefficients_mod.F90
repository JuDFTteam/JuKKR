!------------------------------------------------------------------------------
!> Module to calculate the Gaunt coefficients
!> Declares a data structure GauntCoefficients that contains
!> Gaunt related data needed for a KKR calculation
!> It wraps the previously used routines gaunt and gaunt2
!> @author: Elias Rabel

module GauntCoefficients_mod
  implicit none

  type GauntCoefficients
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

    allocate(coeff%CLEB(NCLEB,2))
    allocate(coeff%ICLEB(NCLEB,3))
    allocate(coeff%JEND(LMPOTD,0:LMAX,0:LMAX))
    allocate(coeff%LOFLM(LM2D))
    coeff%CLEB = 0.0d0
    coeff%ICLEB = -1
    coeff%JEND = -1
    coeff%LOFLM = -1

    allocate(WG(LASSLD))
    allocate(YRG(LASSLD,0:LASSLD,0:LASSLD))

    call GAUNT2(WG,YRG,lmax)
    call GAUNT(lmax,LPOT,WG,YRG,coeff%CLEB,coeff%LOFLM, &
               coeff%ICLEB,coeff%IEND,coeff%JEND,coeff%NCLEB)

    deallocate(WG)
    deallocate(YRG)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyGauntCoefficients(coeff)
    implicit none
    type (GauntCoefficients), intent(inout) :: coeff

    deallocate(coeff%CLEB)
    deallocate(coeff%ICLEB)
    deallocate(coeff%JEND)
    deallocate(coeff%LOFLM)
  end subroutine

end module GauntCoefficients_mod
