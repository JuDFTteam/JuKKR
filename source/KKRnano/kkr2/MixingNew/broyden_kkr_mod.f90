module broyden_kkr_mod
  use brydbm_new_com_mod, only: BRYSH1_new, BRYSH2_new, BRYSH3_new
  use CalculationData_mod
  use broyden_second_mod
  use BroydenData_mod

  use, intrinsic :: ieee_features
  use, intrinsic :: ieee_arithmetic

  implicit none

  private :: calc_metric
  private :: collapse_input_potentials
  private :: collapse_output_potentials

  CONTAINS

  !----------------------------------------------------------------------------
  subroutine mix_broyden2_com(calc_data, iter, communicator)
    implicit none
    type (CalculationData), intent(inout) :: calc_data
    integer, intent(in) :: iter
    integer, intent(in) :: communicator

    integer :: length
    type (BroydenData), pointer :: broyden
    double precision, dimension(:), allocatable :: g_metric_all
    double precision, dimension(:), allocatable :: sm_input
    double precision, dimension(:), allocatable :: fm_output
    double precision :: nan

    length = getBroydenDim(calc_data)
    broyden => getBroyden(calc_data)
    allocate(sm_input(length))
    allocate(fm_output(length))
    allocate(g_metric_all(length))
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
    ! debug
    sm_input = nan
    fm_output = nan
    g_metric_all = nan
    ! end debug
    ! ASSERT length=ntird

    call collapse_input_potentials(calc_data, sm_input)
    call collapse_output_potentials(calc_data, fm_output)
    call calc_all_metrics(calc_data, g_metric_all)

    call broyden_second(sm_input, fm_output, broyden%sm1s, broyden%fm1s, &
                        broyden%ui2, broyden%vi2, g_metric_all, broyden%mixing, &
                        communicator, broyden%itdbryd, length, iter)

    call extract_mixed_potentials(sm_input, calc_data)

    if (any(isnan(sm_input))) then
      write(*,*) "NaN detected!"
      STOP
    end if

  end subroutine

  !----------------------------------------------------------------------------
  !> Returns if x is NaN.
  elemental logical function isnan(x)
    implicit none
    double precision, intent(in) :: x
    isnan = .not. (x == x)
  end function

  !----------------------------------------------------------------------------
  !> Collapse all local input potentials into one array
  subroutine collapse_input_potentials(calc_data, array)
    implicit none
    type (CalculationData), intent(in) :: calc_data
    double precision, dimension(*), intent(inout) :: array

    integer ilocal
    integer num_local_atoms
    integer imap, imap_new
    type (BasisAtom), pointer :: atomdata

    num_local_atoms = getNumLocalAtoms(calc_data)

    imap = 0
    do ilocal = 1, num_local_atoms
      atomdata     => getAtomData(calc_data, ilocal)
      call BRYSH3_new(array(imap+1),atomdata%potential%VISP,atomdata%potential%VINS, &
                      atomdata%potential%irmind,atomdata%potential%irmd, &
                      atomdata%potential%nspin, &
                      imap_new,atomdata%potential%lmpot, &
                      atomdata%potential%irmd, atomdata%potential%irnsd)

      imap = imap + imap_new
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> Collapse all local output potentials into one array
  subroutine collapse_output_potentials(calc_data, array)
    implicit none
    type (CalculationData), intent(in) :: calc_data
    double precision, dimension(*), intent(inout) :: array

    integer ilocal
    integer num_local_atoms
    integer imap, imap_new
    type (BasisAtom), pointer :: atomdata

    num_local_atoms = getNumLocalAtoms(calc_data)

    imap = 0
    do ilocal = 1, num_local_atoms
      atomdata     => getAtomData(calc_data, ilocal)
      call BRYSH1_new(array(imap+1),atomdata%potential%vons, &
                      atomdata%potential%irmind, atomdata%potential%irmd, &
                      atomdata%potential%nspin, imap_new, atomdata%potential%lmpot, &
                      atomdata%potential%irmd)

      imap = imap + imap_new
    end do
  end subroutine

  !----------------------------------------------------------------------------
  !> Extract (mixed) output potentials for all local atoms from 'array' and
  !> overwrite potentials in their respective datastructures (VONS)
  subroutine extract_mixed_potentials(array, calc_data)
    double precision, dimension(:), intent(in) :: array
    type (CalculationData), intent(inout) :: calc_data

    integer ilocal
    integer num_local_atoms
    integer ind, num
    type (BasisAtom), pointer :: atomdata

    num_local_atoms = getNumLocalAtoms(calc_data)

    ind = 1
    num = 0
    do ilocal = 1, num_local_atoms
      atomdata     => getAtomData(calc_data, ilocal)

      num = getNumPotentialValues(atomdata%potential)
      call BRYSH2_new(array(ind:ind+num-1),atomdata%potential%VONS, &
                      atomdata%potential%irmind, &
                      atomdata%potential%irmd, atomdata%potential%nspin, &
                      num, atomdata%potential%LMPOT, atomdata%potential%IRMD)

      ind = ind+num

    end do
  end subroutine

  !----------------------------------------------------------------------------
  subroutine calc_all_metrics(calc_data, g_metric_all)
    implicit none
    type (CalculationData), intent(in) :: calc_data
    double precision, dimension(:), intent(inout) :: g_metric_all

    integer ilocal
    integer num_local_atoms
    integer ind, num
    type (BasisAtom), pointer :: atomdata

    num_local_atoms = getNumLocalAtoms(calc_data)

    ind = 1
    num = 0
    do ilocal = 1, num_local_atoms
      atomdata     => getAtomData(calc_data, ilocal)

      num = getNumPotentialValues(atomdata%potential)
      call calc_metric(g_metric_all(ind:ind+num-1), atomdata%potential%LMPOT, &
                       atomdata%mesh_ptr%r, atomdata%mesh_ptr%drdi, &
                       atomdata%potential%irmd, atomdata%potential%irmind, &
                       atomdata%potential%nspin, num)

      ind = ind + num

    end do
  end subroutine

  !----------------------------------------------------------------------------
  subroutine calc_metric(g_metric, lmpot,r,drdi,irc,irmin,nspin,imap)
    implicit none

    double precision, intent(out) :: g_metric(imap)
    integer :: lmpot, imap, irc, irmin, nspin
    double precision :: r(:), drdi(:)

    integer ij, ir, lm
    integer isp
    double precision volinv

    ij = 0
    do isp = 1,nspin

      volinv = 3.0d0/(r(irc)**3)
      do ir = 1,irc
        ij = ij + 1
        g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
      end do
      !
      if (lmpot.gt.1) then

        do lm = 2,lmpot
          do ir = irmin,irc
            ij = ij + 1
            g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
          end do
        end do
      end if

   end do
 end subroutine

end module broyden_kkr_mod
