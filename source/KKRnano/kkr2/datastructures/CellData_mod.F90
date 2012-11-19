! TODO: allow creation of separate files 'meshes' and 'shapes' ???
! TODO: how to initialise cell_index?

module CellData_mod
  use ShapefunData_mod
  implicit none

  type CellData
    ! cell index?
    integer :: cell_index
    type (ShapefunData)   :: shdata
  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createCellData(cell, irid, lmmax_shape, nfund)
    use ShapefunData_mod
    implicit none
    type (CellData), intent(inout) :: cell
    integer, intent(in) :: irid
    integer, intent(in) :: lmmax_shape
    integer, intent(in) :: nfund

    cell%cell_index = -1
    call createShapefunData(cell%shdata, irid, lmmax_shape, nfund)

  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyCellData(cell)
    use ShapefunData_mod
    implicit none
    type (CellData), intent(inout) :: cell

    call destroyShapefunData(cell%shdata)

  end subroutine

  !----------------------------------------------------------------------------
  !> Write cell data to direct access file 'fileunit' at record 'recnr'
  subroutine writeCellDataDA(cell, fileunit, recnr)

    implicit none
    type (CellData), intent(in) :: cell
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -889271554

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
                                cell%cell_index, &
                                cell%shdata%THETA, &
                                cell%shdata%LLMSP, &
                                cell%shdata%IFUNM, &
                                cell%shdata%LMSP, &
                                cell%shdata%NFU, &
                                MAGIC_NUMBER

  end subroutine

  !> Read cell data from direct access file 'fileunit' at record 'recnr'
  subroutine readCellDataDA(cell, fileunit, recnr)
    implicit none

    type (CellData), intent(inout) :: cell
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = -889271554
    integer :: magic, magic2

    read  (fileunit, rec=recnr) magic, &
                                cell%cell_index, &
                                cell%shdata%THETA, &
                                cell%shdata%LLMSP, &
                                cell%shdata%IFUNM, &
                                cell%shdata%LMSP, &
                                cell%shdata%NFU, &
                                magic2

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) then
      write (*,*) "ERROR: Invalid cell data read. ", __FILE__, __LINE__
      STOP
    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Opens CellData direct access file.
  subroutine openCellDataDAFile(cell, fileunit, filename)
    implicit none

    type (CellData), intent(in) :: cell
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    !------
    integer :: reclen
    integer, parameter :: MAGIC_NUMBER = -889271554

    inquire (iolength = reclen) MAGIC_NUMBER, &
                                cell%cell_index, &
                                cell%shdata%THETA, &
                                cell%shdata%LLMSP, &
                                cell%shdata%IFUNM, &
                                cell%shdata%LMSP, &
                                cell%shdata%NFU, &
                                MAGIC_NUMBER

    !write (*,*) reclen

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes CellData direct access file.
  subroutine closeCellDataDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine

end module CellData_mod
