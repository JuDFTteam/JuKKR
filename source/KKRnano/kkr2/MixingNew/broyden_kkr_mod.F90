module broyden_kkr_mod
#ifndef __GFORTRAN__
  use, intrinsic :: ieee_features
  use, intrinsic :: ieee_arithmetic
#endif
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

    integer :: length, ist
    double precision, allocatable :: g_metric_all(:)
    double precision, allocatable :: sm_input(:)
    double precision, allocatable :: fm_output(:)
    double precision :: nan
#ifdef __GFORTRAN__
    double precision :: zero
#endif

    length = getBroydenDim(calc_data)

    allocate(sm_input(length))
    allocate(fm_output(length))
    allocate(g_metric_all(length))

#ifdef __GFORTRAN__
    zero = 0.d0
    nan = 1.d0/zero
#else
    nan = ieee_value(nan, IEEE_SIGNALING_NAN)
#endif

    ! debug
    sm_input = nan
    fm_output = nan
    g_metric_all = nan
    ! enddebug
    ! ASSERT length=ntird

    call collapse_input_potentials(calc_data, sm_input)
    call collapse_output_potentials(calc_data, fm_output)
    call calc_all_metrics(calc_data, g_metric_all)

#define broyden calc_data%Broyden
    call broyden_second(sm_input, fm_output, broyden%sm1s, broyden%fm1s, &
                        broyden%ui2, broyden%vi2, g_metric_all, broyden%mixing, &
                        communicator, broyden%itdbryd, length, iter)
#undef  broyden
    call extract_mixed_potentials(sm_input, calc_data)

    if (any(isnan(sm_input))) then
      write(*,*) "NaN detected!"
      stop
    endif

    deallocate(g_metric_all, fm_output, sm_input, stat=ist)
  endsubroutine ! mix_broyden2_com

  !----------------------------------------------------------------------------
  !> Returns if x is NaN.
  elemental logical function isnan(x)
    double precision, intent(in) :: x
    isnan = .not. (x == x)
  endfunction ! isnan

  !----------------------------------------------------------------------------
  !> Collapse all local input potentials into one array
  subroutine collapse_input_potentials(calc_data, array)
    use brydbm_new_com_mod, only: BRYSH3_new
    use CalculationData_mod, only: CalculationData, getAtomData
    use BasisAtom_mod, only: BasisAtom
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: array(*)

    integer :: ilocal, num_local_atoms, imap, imap_new
    type(BasisAtom), pointer :: atomdata

    num_local_atoms = calc_data%num_local_atoms

    imap = 0
    do ilocal = 1, num_local_atoms
      atomdata =>getAtomData(calc_data, ilocal)
      call BRYSH3_new(array(imap+1),atomdata%potential%VISP,atomdata%potential%VINS, &
                      atomdata%potential%irmind,atomdata%potential%irmd, &
                      atomdata%potential%nspin, &
                      imap_new,atomdata%potential%lmpot, &
                      atomdata%potential%irmd, atomdata%potential%irnsd)

      imap = imap + imap_new
    enddo ! ilocal
  endsubroutine ! collapse_input_potentials

  !----------------------------------------------------------------------------
  !> Collapse all local output potentials into one array
  subroutine collapse_output_potentials(calc_data, array)
    use brydbm_new_com_mod, only: BRYSH1_new
    use CalculationData_mod, only: CalculationData, getAtomData
    use BasisAtom_mod, only: BasisAtom
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: array(*)

    integer :: ilocal, num_local_atoms, imap, imap_new
    type(BasisAtom), pointer :: atomdata

    num_local_atoms = calc_data%num_local_atoms

    imap = 0
    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)
      call BRYSH1_new(array(imap+1),atomdata%potential%vons, &
                      atomdata%potential%irmind, atomdata%potential%irmd, &
                      atomdata%potential%nspin, imap_new, atomdata%potential%lmpot, &
                      atomdata%potential%irmd)

      imap = imap + imap_new
    enddo ! ilocal
  endsubroutine ! collapse_output_potentials

  !----------------------------------------------------------------------------
  !> Extract (mixed) output potentials for all local atoms from 'array' and
  !> overwrite potentials in their respective datastructures (VONS)
  subroutine extract_mixed_potentials(array, calc_data)
    use brydbm_new_com_mod, only: BRYSH2_new
    use CalculationData_mod, only: CalculationData, getAtomData
    use BasisAtom_mod, only: BasisAtom
    use PotentialData_mod, only: getNumPotentialValues
  
    double precision, intent(in) :: array(:)
    type(CalculationData), intent(inout) :: calc_data

    integer :: ilocal, num_local_atoms, ind, num
    type(BasisAtom), pointer :: atomdata

    num_local_atoms = calc_data%num_local_atoms

    ind = 1
    num = 0
    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)

      num = getNumPotentialValues(atomdata%potential)
      call BRYSH2_new(array(ind:ind+num-1),atomdata%potential%VONS, &
                      atomdata%potential%irmind, &
                      atomdata%potential%irmd, atomdata%potential%nspin, &
                      num, atomdata%potential%LMPOT, atomdata%potential%IRMD)

      ind = ind+num

    enddo ! ilocal
  endsubroutine ! extract_mixed_potentials

  !----------------------------------------------------------------------------
  subroutine calc_all_metrics(calc_data, g_metric_all)
    use CalculationData_mod, only: CalculationData, getAtomData
    use BasisAtom_mod, only: BasisAtom
    use PotentialData_mod, only: getNumPotentialValues
  
    type(CalculationData), intent(in) :: calc_data
    double precision, intent(inout) :: g_metric_all(:)

    integer :: ilocal, num_local_atoms, ind, num
    type(BasisAtom), pointer :: atomdata

    num_local_atoms = calc_data%num_local_atoms

    ind = 1
    num = 0
    do ilocal = 1, num_local_atoms
      atomdata => getAtomData(calc_data, ilocal)

      num = getNumPotentialValues(atomdata%potential)
      call calc_metric(g_metric_all(ind:ind+num-1), atomdata%potential%LMPOT, &
                       atomdata%mesh_ptr%r, atomdata%mesh_ptr%drdi, &
                       atomdata%potential%irmd, atomdata%potential%irmind, &
                       atomdata%potential%nspin, num)

      ind = ind + num
    enddo ! ilocal
  endsubroutine ! calc_all_metrics

  !----------------------------------------------------------------------------
  subroutine calc_metric(g_metric, lmpot, r, drdi, irc, irmin, nspin, nmap)
    double precision, intent(out) :: g_metric(nmap)
    integer, intent(in) :: lmpot, nmap, irc, irmin, nspin
    double precision :: r(:), drdi(:)

    integer :: ij, ir, lm, isp
    double precision :: volinv

    ij = 0
    do isp = 1,nspin

      volinv = 3.0d0/(r(irc)**3)
      do ir = 1,irc
        ij = ij + 1
        g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
      enddo ! ir
      !
      if (lmpot > 1) then

        do lm = 2, lmpot
          do ir = irmin, irc
            ij = ij + 1
            g_metric(ij) = volinv*r(ir)*r(ir)*drdi(ir)
          enddo ! ir
        enddo ! lm
      endif

    enddo ! isp
  endsubroutine ! calc_metric

endmodule ! broyden_kkr_mod
