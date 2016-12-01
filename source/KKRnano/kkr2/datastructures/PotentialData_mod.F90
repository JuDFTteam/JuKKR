!> Structure that stores the input and output potential for one cell.
!> @author Elias Rabel

module PotentialData_mod
  implicit none
  private
  public :: PotentialData, create, destroy, represent
  public :: getNumPotentialValues
  
  type PotentialData
    double precision, allocatable :: vins(:,:,:)        !< input potential - nonspherical parts only
    !> spherically averaged input potential
    !> ATTENTION: this differs by a factor of 1/sqrt(4pi) from the L=(0,0) component
    !> of the L-expansion
    !> Note: the L=(0,0) (index 1) component of the output potential vons does not
    !> contain this factor!
    double precision, allocatable :: visp(:,:)
    double precision, allocatable :: vons(:,:,:)        !< output potential - contains all L-components
    integer :: nspin
    integer :: lpot
    integer :: irmind
    integer :: irmd
    integer :: irnsd
    integer :: lmpot
  endtype ! PotentialData

  interface create
    module procedure createPotentialData
  endinterface
  
  interface destroy
    module procedure destroyPotentialData
  endinterface

  interface represent
    module procedure repr_PotentialData
  endinterface
  
  contains

  !----------------------------------------------------------------------------
  subroutine createPotentialData(potential, lpot, nspin, irmind, irmd)
    type(PotentialData), intent(inout) :: potential
    integer, intent(in) :: lpot
    integer, intent(in) :: nspin
    integer, intent(in) :: irmind  !< start of non-spherical region
    integer, intent(in) :: irmd  !< number of mesh points

    ! locals
    integer :: lmpot

    potential%nspin = nspin
    potential%lpot = lpot
    potential%irmind = irmind
    potential%irmd = irmd
    potential%irnsd = irmd - irmind

    lmpot = (lpot + 1)**2
    potential%lmpot = lmpot

    ! Unfortunately, for compatibility reasons,
    ! memory has to be allocated for both spin directions
    ! independent of nspin
    allocate(potential%vins(irmind:irmd,lmpot,2))
    allocate(potential%visp(irmd,2))
    allocate(potential%vons(irmd,lmpot,2))
    potential%vins = 0.d0
    potential%visp = 0.d0
    potential%vons = 0.d0

  endsubroutine ! create

  !----------------------------------------------------------------------------
  elemental subroutine destroyPotentialData(potential)
    type(PotentialData), intent(inout) :: potential

    deallocate(potential%vins)
    deallocate(potential%visp)
    deallocate(potential%vons)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Return the actual number of potential values.
  integer elemental function getNumPotentialValues(self)
    type(PotentialData), intent(in) :: self
    getNumPotentialValues = (self%irmd+(self%irnsd+1) * (self%lmpot-1)) * self%nspin
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns a string representation of PotentialData.
  subroutine repr_PotentialData(potential, str)
    type(PotentialData), intent(in) :: potential
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(len=80) :: buffer
    integer :: ind, lm, ispin

    nl = new_line(' ')

    str = ''
    write(buffer, *) "nspin  = ", potential%nspin
    str = str // trim(buffer) // nl
    write(buffer, *) "lpot   = ", potential%lpot
    str = str // trim(buffer) // nl
    write(buffer, *) "irmind = ", potential%irmind
    str = str // trim(buffer) // nl
    write(buffer, *) "irmd   = ", potential%irmd
    str = str // trim(buffer) // nl
    write(buffer, *) "irnsd  = ", potential%irnsd
    str = str // trim(buffer) // nl
    write(buffer, *) "lmpot  = ", potential%lmpot
    str = str // trim(buffer) // nl // nl

    write(buffer, *) "nr.    spin       visp"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%visp, 2)
      do ind = 1, size(potential%visp, 1)
        write(buffer, '(2(i5,2x),e23.16)') ind, ispin, potential%visp(ind, ispin)
        str = str // trim(buffer) // nl
      enddo ! ind
    enddo ! ispin

    str = str // nl

    write(buffer, *) "nr.      LM    spin  vins"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%vins, 3)
      do lm = 1, size(potential%vins, 2)
        ! print only non-zero potential components
        if (sum(abs(potential%vins(:, lm, ispin))) > 0.0) then
          do ind = lbound(potential%vins, 1), ubound(potential%vins, 1)
            write(buffer, '(3(i5,2x),e23.16)') ind, lm, ispin, potential%vins(ind, lm, ispin)
            str = str // trim(buffer) // nl
          enddo ! ind
        endif
      enddo ! lm
    enddo ! ispin

    str = str // nl

    write(buffer, *) "nr.      LM    spin  vons"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%vons, 3)
      do lm = 1, size(potential%vons, 2)
        ! print only non-zero potential components
        if (sum(abs(potential%vons(:, lm, ispin))) > 0.0) then
          do ind = lbound(potential%vons, 1), ubound(potential%vons, 1)
            write(buffer, '(3(i5,2x),e23.16)') ind, lm, ispin, potential%vons(ind, lm, ispin)
            str = str // trim(buffer) // nl
          enddo ! ind
        endif
      enddo ! lm
    enddo ! ispin
    
  endsubroutine ! represent

endmodule ! PotentialData_mod
