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

    double precision :: max_muffin_tin !< maximum muffin tin radius in units of ALAT!!!
    integer :: num_faces !< number of faces of voronoi-cell

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
    shdata%max_muffin_tin = 0.0d0
    shdata%num_faces = 0

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

  !----------------------------------------------------------------------------
  !> Returns a string representation of ShapefunData.
  subroutine repr_ShapefunData(shdata, str)
    implicit none
    class (ShapefunData), intent(in) :: shdata
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(80) :: buffer
    integer :: ind, ifun

    nl = new_line(' ')

    str = ''
    write(buffer, *) "irid           = ", shdata%irid
    str = str // trim(buffer) // nl
    write(buffer, *) "nfund          = ", shdata%nfund
    str = str // trim(buffer) // nl
    write(buffer, *) "NFU            = ", shdata%NFU
    str = str // trim(buffer) // nl
    write(buffer, *) "lmmax_shape    = ", shdata%lmmax_shape
    str = str // trim(buffer) // nl
    write(buffer, *) "max_muffin_tin = ", shdata%max_muffin_tin !< maximum muffin tin radius in units of ALAT!!!
    str = str // trim(buffer) // nl
    write(buffer, *) "num_faces = ", shdata%num_faces
    str = str // trim(buffer) // nl

    write(buffer, *) "LLMSP = "
    str = str // trim(buffer) // nl
    do ind = 1, size(shdata%llmsp)
      write(buffer, *) shdata%llmsp(ind)
      str = str // trim(buffer) // nl
    end do

    write(buffer, '(A5, 2X, A5, 2X, A5)') "LM", "LMSP", "IFUNM"
    str = str // trim(buffer) // nl
    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl
    do ind = 1, size(shdata%lmsp)
      write(buffer, '(I5, 2X, I5, 2X, I5)') ind, shdata%lmsp(ind), shdata%ifunm(ind)
      str = str // trim(buffer) // nl
    end do

    write(buffer, '(A5, 2X, A5, 2X, A5, 2X, A5)') "ind", "ifun", "LM", "THETA"
    str = str // trim(buffer) // nl
    do ifun = 1, size(shdata%theta, 2)
      if (sum(abs(shdata%theta(:, ifun))) > 0.0) then
        do ind = 1, size(shdata%theta, 1)
          write(buffer, '(I5, 2X, I5, 2X, I5, 2X, E23.16)') ind, ifun, shdata%llmsp(ifun), shdata%THETA(ind, ifun)
          str = str // trim(buffer) // nl
        end do
      end if
    end do

  end subroutine

end module ShapefunData_mod
