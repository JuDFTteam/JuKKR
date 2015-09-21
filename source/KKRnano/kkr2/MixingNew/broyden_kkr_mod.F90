module broyden_kkr_mod
  use, intrinsic :: ieee_features
  use, intrinsic :: ieee_arithmetic
  implicit none
  private
  public :: mix_broyden2_com

  contains

  !----------------------------------------------------------------------------
  subroutine mix_broyden2_com(calc_data, iter, communicator)
    use CalculationData_mod, only: CalculationData, getBroydenDim
    use broyden_second_mod, only: broyden_second

    type(CalculationData), intent(inout) :: calc_data
    integer, intent(in) :: iter
    integer, intent(in) :: communicator

    integer :: length
    double precision, allocatable :: g_metric_all(:)
    double precision, allocatable :: sm_input(:)
    double precision, allocatable :: fm_output(:)
    double precision :: nan

    length = getBroydenDim(calc_data)

    allocate(sm_input(length))
    allocate(fm_output(length))
    allocate(g_metric_all(length))

    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
    ! debug
    sm_input = nan
    fm_output = nan
    g_metric_all = nan
    ! enddebug
    ! ASSERT length=ntird

    call collapse_input_potentials(calc_data, sm_input)
    call collapse_output_potentials(calc_data, fm_output)
    call calc_all_metrics(calc_data, g_metric_all)

#define broyden calc_data%broyden
    call broyden_second(sm_input, fm_output, broyden%sm1s, broyden%fm1s, &
                        broyden%ui2, broyden%vi2, g_metric_all, broyden%mixing, &
                        communicator, broyden%itdbryd, length, iter)
#undef broyden
    call extract_mixed_potentials(sm_input, calc_data)

    if (any(isnan(sm_input))) then
      write(*,*) "NaN detected!"
      STOP
    endif

    deallocate(g_metric_all)
    deallocate(fm_output)
    deallocate(sm_input)

  endsubroutine

  !----------------------------------------------------------------------------
  !> Returns if x is NaN.
  elemental logical function isnan(x)
    double precision, intent(in) :: x
    isnan = .not. (x == x)
  endfunction

  !----------------------------------------------------------------------------
  !> Collapse all local input potentials into one array
  subroutine collapse_input_potentials(calc_data, array)
    use brydbm_new_com_mod, only: BRYSH3_new
    use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getAtomData
    use BasisAtom_mod, only: BasisAtom
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: array(*)

    integer ilocal
    integer num_local_atoms
    integer imap, imap_new
    type(BasisAtom), pointer :: atomdata

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
    enddo
  endsubroutine

  !----------------------------------------------------------------------------
  !> Collapse all local output potentials into one array
  subroutine collapse_output_potentials(calc_data, array)
    use brydbm_new_com_mod, only: BRYSH1_new
    use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getAtomData
    use BasisAtom_mod, only: BasisAtom
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: array(*)

    integer ilocal
    integer num_local_atoms
    integer imap, imap_new
    type(BasisAtom), pointer :: atomdata

    num_local_atoms = getNumLocalAtoms(calc_data)

    imap = 0
    do ilocal = 1, num_local_atoms
      atomdata     => getAtomData(calc_data, ilocal)
      call BRYSH1_new(array(imap+1),atomdata%potential%vons, &
                      atomdata%potential%irmind, atomdata%potential%irmd, &
                      atomdata%potential%nspin, imap_new, atomdata%potential%lmpot, &
                      atomdata%potential%irmd)

      imap = imap + imap_new
    enddo
  endsubroutine

  !----------------------------------------------------------------------------
  !> Extract (mixed) output potentials for all local atoms from 'array' and
  !> overwrite potentials in their respective datastructures (VONS)
  subroutine extract_mixed_potentials(array, calc_data)
    use brydbm_new_com_mod, only: BRYSH2_new
    use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getAtomData
    use BasisAtom_mod, only: BasisAtom
    use PotentialData_mod, only: getNumPotentialValues
  
    double precision, intent(in) :: array(:)
    type(CalculationData), intent(inout) :: calc_data

    integer :: ilocal
    integer :: num_local_atoms
    integer :: ind, num
    type(BasisAtom), pointer :: atomdata

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

    enddo
  endsubroutine

  !----------------------------------------------------------------------------
  subroutine calc_all_metrics(calc_data, g_metric_all)
    use CalculationData_mod, only: CalculationData, getNumLocalAtoms, getAtomData
    use BasisAtom_mod, only: BasisAtom
    use PotentialData_mod, only: getNumPotentialValues
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: g_metric_all(:)

    integer ilocal
    integer num_local_atoms
    integer ind, num
    type(BasisAtom), pointer :: atomdata

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

    enddo
  endsubroutine

  !----------------------------------------------------------------------------
  subroutine calc_metric(g_metric, lmpot,r,drdi,irc,irmin,nspin,imap)
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
      enddo
      !
      if (lmpot > 1) then

        do lm = 2, lmpot
          do ir = irmin, irc
            ij = ij + 1
            g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
          enddo
        enddo
      endif

   enddo
 endsubroutine

endmodule broyden_kkr_mod
