!------------------------------------------------------------------------------
!> Data structure that contains cell specific data.
!>
!> Not very useful anymore, kept only for historical reasons.

module CellData_mod
  use ShapefunData_mod, only: ShapefunData
  implicit none
  private
  public :: CellData, create, destroy
  public :: createCellData, destroyCellData ! deprecated

  type CellData
    ! cell index?
    integer :: cell_index
    type (ShapefunData)   :: shdata
  end type

  interface create
    module procedure createCellData
  endinterface
  
  interface destroy
    module procedure destroyCellData
  endinterface
  
  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createCellData(cell, irid, lmmax_shape, nfund)
    use ShapefunData_mod, only: create
    type (CellData), intent(inout) :: cell
    integer, intent(in) :: irid
    integer, intent(in) :: lmmax_shape
    integer, intent(in) :: nfund

    cell%cell_index = -1
    call create(cell%shdata, irid, lmmax_shape, nfund)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyCellData(cell)
    use ShapefunData_mod, only: destroy
    type (CellData), intent(inout) :: cell

    call destroy(cell%shdata)

  end subroutine

end module CellData_mod
