!> Note: VONS is wastefully dimensioned... it is allocated for total number of
!> mesh points, but only non-spherical part is non-zero
!> also there is always space for 2 spin-directions.

module PotentialData_mod
  implicit none

  type PotentialData
    double precision, dimension(:,:,:), allocatable :: VINS        ! .. input potential
    double precision, dimension(:,:), allocatable :: VISP
    double precision, dimension(:,:,:), allocatable :: VONS        !     .. output potential
    integer :: nspin
    integer :: lpot
    integer :: irmind
    integer :: irmd
    integer :: irnsd
    ! derived
    integer :: lmpot

  end type

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine createPotentialData(potential, lpot, nspin, irmind, irmd)
    implicit none
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
    potential%VINS = 0.0d0
    potential%VISP = 0.0d0
    potential%VONS = 0.0d0

  end subroutine


  !----------------------------------------------------------------------------
  subroutine destroyPotentialData(potential)
    implicit none
    type (PotentialData), intent(inout) :: potential

    deallocate(potential%VINS)
    deallocate(potential%VISP)
    deallocate(potential%VONS)
  end subroutine

  !----------------------------------------------------------------------------
  !> Return the actual number of potential values.
  integer function getNumPotentialValues(potential)
    implicit none
    type (PotentialData), intent(in) :: potential
    getNumPotentialValues = (potential%irmd+(potential%irnsd+1) * &
                            (potential%lmpot-1)) &
                            * potential%nspin
  end function

  !----------------------------------------------------------------------------
  !> Returns a string representation of PotentialData.
  subroutine repr_PotentialData(potential, str)
    implicit none
    class (PotentialData), intent(in) :: potential
    character(len=:), allocatable, intent(inout) :: str

    character :: nl
    character(80) :: buffer
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
        write(buffer, '(I5, 2X, I5, 2X, E23.16)') ind, ispin, potential%VISP(ind, ispin)
        str = str // trim(buffer) // nl
      end do
    end do

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
            write(buffer, '(I5, 2X, I5, 2X, I5, 2X E23.16)') ind, lm, ispin, potential%VINS(ind, lm, ispin)
            str = str // trim(buffer) // nl
          end do
        end if
      end do
    end do

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
            write(buffer, '(I5, 2X, I5, 2X, I5, 2X E23.16)') &
              ind, lm, ispin, potential%VONS(ind, lm, ispin)
            str = str // trim(buffer) // nl
          end do
        end if
      end do
    end do
  end subroutine

end module PotentialData_mod
