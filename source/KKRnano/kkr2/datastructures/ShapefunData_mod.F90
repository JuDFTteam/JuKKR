!> This module defines a datatype that contains shapefunction related data
!> @author Elias Rabel
module ShapefunData_mod
  implicit none
  private
  public :: ShapefunData, create, destroy, represent

  type ShapefunData
    ! dimension params
    integer :: irid
    integer :: nfund
    integer :: lmmax_shape !< former name LMXSPD

    double precision, allocatable :: theta(:,:)
    integer, allocatable :: llmsp(:)
    integer, allocatable :: ifunm(:)
    integer, allocatable :: lmsp(:) !< =0 if shape-function component is zero, otherwise =1
    integer :: nfu

    double precision :: max_muffin_tin !< maximum muffin tin radius in units of ALAT!!!
    integer :: num_faces !< number of faces of voronoi-cell

    integer :: cell_index ! remainder from former type(CellData)
    
  endtype
  
  interface create
    module procedure createShapefunData
  endinterface
  
  interface destroy
    module procedure destroyShapefunData
  endinterface

  interface represent
    module procedure representShapefunData
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine createShapefunData(self, irid, lmmax_shape, nfund)
    type(ShapefunData), intent(inout) :: self
    integer, intent(in) :: irid
    integer, intent(in) :: lmmax_shape
    integer, intent(in) :: nfund

    self%irid = irid
    self%nfund = nfund
    self%lmmax_shape = lmmax_shape

    allocate(self%theta(irid,nfund), self%llmsp(nfund), self%ifunm(lmmax_shape), self%lmsp(lmmax_shape))

    self%theta = 0.d0
    self%llmsp = 0
    self%ifunm = 0
    self%lmsp =  0
    self%nfu = 0
    self%max_muffin_tin = 0.d0
    self%num_faces = 0

    self%cell_index = -1 ! ??
    
    ! TODO: check lmmax_shape <= nfund
  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyShapefunData(self)
    type(ShapefunData), intent(inout) :: self
    integer :: ist
    deallocate(self%theta, self%llmsp, self%ifunm, self%lmsp, stat=ist)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Returns a string representation of ShapefunData.
  subroutine representShapefunData(self, str)
    type(ShapefunData), intent(in) :: self
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(len=80) :: buffer
    integer :: ind, ifun

    nl = new_line(' ')

    str = ''
    write(buffer, *) "irid           = ", self%irid
    str = str // trim(buffer) // nl
    write(buffer, *) "nfund          = ", self%nfund
    str = str // trim(buffer) // nl
    write(buffer, *) "nfu            = ", self%nfu
    str = str // trim(buffer) // nl
    write(buffer, *) "lmmax_shape    = ", self%lmmax_shape
    str = str // trim(buffer) // nl
    write(buffer, *) "max_muffin_tin = ", self%max_muffin_tin !< maximum muffin tin radius in units of ALAT!!!
    str = str // trim(buffer) // nl
    write(buffer, *) "num_faces = ", self%num_faces
    str = str // trim(buffer) // nl

    write(buffer, *) "llmsp = "
    str = str // trim(buffer) // nl
    do ind = 1, size(self%llmsp)
      write(buffer, *) self%llmsp(ind)
      str = str // trim(buffer) // nl
    enddo ! ind

    write(buffer, '(A5, 2X, A5, 2X, A5)') "LM", "lmsp", "ifunm"
    str = str // trim(buffer) // nl
    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl
    do ind = 1, size(self%lmsp)
      write(buffer, '(I5, 2X, I5, 2X, I5)') ind, self%lmsp(ind), self%ifunm(ind)
      str = str // trim(buffer) // nl
    enddo ! ind

    write(buffer, '(A5, 2X, A5, 2X, A5, 2X, A5)') "ind", "ifun", "LM", "theta"
    str = str // trim(buffer) // nl
    do ifun = 1, size(self%theta, 2)
      if (sum(abs(self%theta(:, ifun))) > 0.0) then
        do ind = 1, size(self%theta, 1)
          write(buffer, '(I5, 2X, I5, 2X, I5, 2X, E23.16)') ind, ifun, self%llmsp(ifun), self%theta(ind, ifun)
          str = str // trim(buffer) // nl
        enddo ! ind
      endif
    enddo ! ifun

  endsubroutine ! represent

endmodule ! ShapefunData_mod
