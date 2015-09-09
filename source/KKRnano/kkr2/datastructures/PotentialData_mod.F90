!> Structure that stores the input and output potential for one cell.
!> @author Elias Rabel

module PotentialData_mod
  implicit none
  private
  public :: PotentialData, create, destroy, represent
  public :: createPotentialData, destroyPotentialData, repr_PotentialData ! deprecated
  public :: getNumPotentialValues
  
  type PotentialData
    double precision, allocatable :: VINS(:,:,:)        !< input potential - nonspherical parts only
    !> spherically averaged input potential
    !> ATTENTION: this differs by a factor of 1/sqrt(4pi) from the L=(0,0) component
    !> of the L-expansion
    !> Note: the L=(0,0) (index 1) component of the output potential VONS does not
    !> contain this factor!
    double precision, allocatable :: VISP(:,:)
    double precision, allocatable :: VONS(:,:,:)        !< output potential - contains all L-components
    integer :: nspin
    integer :: lpot
    integer :: irmind
    integer :: irmd
    integer :: irnsd
    integer :: lmpot
  endtype

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
    type (PotentialData), intent(inout) :: potential
    integer, intent(in) :: lpot
    integer, intent(in) :: nspin
    integer, intent(in) :: irmind  !< start of non-spherical region
    integer, intent(in) :: irmd  !< number of mesh points

    ! ---- local ----------
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
    allocate(potential%VINS(IRMIND:IRMD,LMPOT,2))
    allocate(potential%VISP(IRMD,2))
    allocate(potential%VONS(IRMD,LMPOT,2))
    potential%VINS = 0.d0
    potential%VISP = 0.d0
    potential%VONS = 0.d0

  endsubroutine ! create


  !----------------------------------------------------------------------------
  subroutine destroyPotentialData(potential)
    type (PotentialData), intent(inout) :: potential

    deallocate(potential%VINS)
    deallocate(potential%VISP)
    deallocate(potential%VONS)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Return the actual number of potential values.
  integer function getNumPotentialValues(potential)
    type (PotentialData), intent(in) :: potential
    getNumPotentialValues = (potential%irmd+(potential%irnsd+1) * (potential%lmpot-1)) * potential%nspin
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Returns a string representation of PotentialData.
  subroutine repr_PotentialData(potential, str)
    class(PotentialData), intent(in) :: potential
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

    write(buffer, *) "nr.    spin       VISP"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%VISP, 2)
      do ind = 1, size(potential%VISP, 1)
        write(buffer, '(2(i5,2x),e23.16)') ind, ispin, potential%VISP(ind, ispin)
        str = str // trim(buffer) // nl
      enddo ! ind
    enddo ! ispin

    str = str // nl

    write(buffer, *) "nr.      LM    spin  VINS"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%VINS, 3)
      do lm = 1, size(potential%VINS, 2)
        ! print only non-zero potential components
        if (sum(abs(potential%VINS(:, lm, ispin))) > 0.0) then
          do ind = lbound(potential%VINS, 1), ubound(potential%VINS, 1)
            write(buffer, '(3(i5,2x),e23.16)') ind, lm, ispin, potential%VINS(ind, lm, ispin)
            str = str // trim(buffer) // nl
          enddo ! ind
        endif
      enddo ! lm
    enddo ! ispin

    str = str // nl

    write(buffer, *) "nr.      LM    spin  VONS"
    str = str // trim(buffer) // nl

    write(buffer, '(79("="))')
    str = str // trim(buffer) // nl

    do ispin = 1, size(potential%VONS, 3)
      do lm = 1, size(potential%VONS, 2)
        ! print only non-zero potential components
        if (sum(abs(potential%VONS(:, lm, ispin))) > 0.0) then
          do ind = lbound(potential%VONS, 1), ubound(potential%VONS, 1)
            write(buffer, '(3(i5,2x),e23.16)') ind, lm, ispin, potential%VONS(ind, lm, ispin)
            str = str // trim(buffer) // nl
          enddo ! ind
        endif
      enddo ! lm
    enddo ! ispin
    
  endsubroutine ! represent

endmodule PotentialData_mod
