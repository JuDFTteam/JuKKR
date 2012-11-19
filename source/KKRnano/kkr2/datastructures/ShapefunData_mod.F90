!> This module defines a datatype that contains shapefunction related data
!> @author Elias Rabel
module ShapefunData_mod
  implicit none

  type ShapefunData
    ! dimension params
    integer :: irid
    integer :: nfund
    integer :: lmmax_shape !< former name LMXSPD

    double precision, dimension(:,:), allocatable :: THETA
    integer, dimension(:), allocatable :: LLMSP
    integer, dimension(:), allocatable :: IFUNM
    integer, dimension(:), allocatable :: LMSP !< =0 if shape-function component is zero, otherwise =1
    integer :: NFU

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createShapefunData(shdata, irid, lmmax_shape, nfund)
    implicit none
    type (ShapeFunData), intent(inout) :: shdata
    integer, intent(in) :: irid
    integer, intent(in) :: lmmax_shape
    integer, intent(in) :: nfund

    shdata%irid = irid
    shdata%nfund = nfund
    shdata%lmmax_shape = lmmax_shape

    allocate(shdata%THETA(irid, nfund))
    allocate(shdata%LLMSP(nfund))
    allocate(shdata%IFUNM(lmmax_shape))
    allocate(shdata%LMSP(lmmax_shape))

    shdata%THETA = 0.0d0
    shdata%LLMSP = 0
    shdata%IFUNM = 0
    shdata%LMSP =  0
    shdata%NFU = 0

    ! TODO: check lmmax_shape <= nfund
  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyShapefunData(shdata)
    implicit none
    type (ShapeFunData), intent(inout) :: shdata

    deallocate(shdata%THETA)
    deallocate(shdata%LLMSP)
    deallocate(shdata%IFUNM)
    deallocate(shdata%LMSP)
  end subroutine

end module ShapefunData_mod
