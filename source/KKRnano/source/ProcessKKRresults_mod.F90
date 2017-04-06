#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"
#include "DebugHelpers/test_macros.h"

module ProcessKKRresults_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  use Logging_mod, only:    ! import no name, just mention it for the module dependency 
  use arraytest2_mod, only: ! import no name, just mention it for the module dependency 
  implicit none
  private

  public :: processKKRresults, output_forces

  integer, private, parameter :: MAX_MADELUNG_RADIUS_INDEX=101

  contains

  !------------------------------------------------------------------------------
  !> Returns 1 when target rms error has been reached,
  !> master rank adds 2 if STOP-file present.
  integer function processKKRresults(iter, calc, mp, emesh, dims, params, arrays, program_timer)

    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD

    use KKRnanoParallel_mod, only: KKRnanoParallel
    use EnergyMesh_mod, only: EnergyMesh
    use CalculationData_mod, only: CalculationData
    use TimerMpi_mod, only: TimerMpi, getTime, outTime
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use BasisAtom_mod, only: BasisAtom, openBasisAtomPotentialDAFile, resetPotentials, writeBasisAtomPotentialDA, closeBasisAtomPotentialDAFile
    use RadialMeshData_mod, only: RadialMeshData
    use DensityResults_mod, only: DensityResults
    use EnergyResults_mod, only: EnergyResults
    use CalculationData_mod, only: getDensities, getEnergies, getAtomData
    use NonCollinearMagnetismData_mod, only: store
    include 'mpif.h'

    integer, intent(in)                                 :: iter
    type(CalculationData), intent(inout)                :: calc
    type(KKRnanoParallel), intent(in)                   :: mp
    type(EnergyMesh), intent(inout)                     :: emesh
    type(Main2Arrays), intent(in)                       :: arrays
    type(DimParams), intent(in)                         :: dims
    type(InputParams), intent(in)                       :: params
    type(TimerMpi), intent(inout)                       :: program_timer

    ! locals
    type(BasisAtom) , pointer                           :: atomdata
    type(DensityResults), pointer                       :: densities
    type(EnergyResults), pointer                        :: energies

    type(RadialMeshData), pointer :: mesh
    integer :: atom_id, ila, num_local_atoms, ierr

    processKKRresults = 0
    num_local_atoms = calc%num_local_atoms

    densities => getDensities(calc, 1)
    energies  => getEnergies(calc, 1)

    ! kkr
    !  |
    !  v

    call calculateDensities(iter, calc, mp, dims, params, program_timer, arrays, emesh)

    ! |
    ! v
    ! modified: densities, emesh, energies (ESPV only)
    ! |
    ! v
    call calculatePotentials(iter, calc, mp, dims, params, program_timer, arrays)
    ! |
    ! v
    ! modified: atomdata, energies
    !
    ! |
    ! v

    ! mix_potential returns 1 when target_rms reached, otherwise 0

    ! ATTENTION: the spherical part of the potential is divided by sqrt(4*pi) here
    ! this is definitely not the optimal place to do this
    processKKRresults = mix_potential(calc, iter, params, dims, mp)

    ! |
    ! v
    ! modified: atomdata


    ! use any atomdata to open file - use reclen stored in calc
    atomdata => getAtomData(calc, 1)
    call openBasisAtomPotentialDAFile(atomdata, 37, "bin.vpotnew", calc%max_reclen_potential, action='write')

    !if (dims%korbit == 1 .and. mp%isMasterRank) then ! NOCO
    !  call store(calc%noco_data, 'bin.noco')
    !endif

    do ila = 1, num_local_atoms ! no OpenMP
      atomdata => getAtomData(calc, ila)
      energies => getEnergies(calc, ila)
      mesh => atomdata%mesh_ptr
      atom_id = calc%atom_ids(ila) ! get global atom_id from local index

      TESTARRAYLOG(3, atomdata%potential%VINS)
      TESTARRAYLOG(3, atomdata%potential%VISP)
      TESTARRAYLOG(3, atomdata%potential%VONS)

      !----------------------------------------------------------------------
      ! -->    reset to start new iteration
      ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      call resetPotentials(atomdata)

      ! ----------------------------------------------------- output_potential
      call writeBasisAtomPotentialDA(atomdata, 37, atom_id)
      ! ----------------------------------------------------- output_potential

    enddo ! ila

    call closeBasisAtomPotentialDAFile(37)

    call outTime(mp%isMasterRank, 'potential written ..............', getTime(program_timer), iter)

  ! Wait here in order to guarantee regular and non-errorneous output
  ! in RESULTS
    call outTime(mp%isMasterRank, 'barrier begin ..................', getTime(program_timer), iter)

    call MPI_Barrier(mp%mySEComm, IERR)

    call outTime(mp%isMasterRank, 'barrier end ....................', getTime(program_timer), iter)

    ! output of local density of states (task-local files)
    if (params%NPOL == 0) then
      do ila = 1, num_local_atoms
        densities => getDensities(calc, ila)
        atomdata  => getAtomData(calc, ila)
        atom_id = calc%atom_ids(ila) ! get global atom_id from local index

        ! TODO: Fermi energy written to DOS files is just a dummy value (99.99)
        call write_LDOS(densities%DEN, emesh%EZ, densities%lmaxd+1, emesh%IELAST, atomdata%core%ITITLE(:,1:dims%NSPIND), 99.99d0, &
                        emesh%E1, emesh%E2, params%ALAT, emesh%TK, dims%NSPIND, atom_id)
      enddo ! ila
    endif

  ! -----------------------------------------------------------------
  ! BEGIN: only MASTERRANK is working here
  ! -----------------------------------------------------------------
    if (mp%isMasterRank) then

      ! DOS was written to file 'results1' and read out here just
      ! to be written in routine wrldos (new: file complex.dos only)
      ! also other stuff is read from results1 (and results2)

      ! TODO: note: title written to complex.dos is not correct
      ! - taken from 1st local atom only
      ! TODO: Fermi energy written to complex.dos is not correct
      call RESULTS(dims%LRECRES2, densities%IEMXD, ITER, dims%LMAXD, &
      arrays%NAEZ, emesh%NPOL, &
      dims%NSPIND, params%KPRE, params%KTE, atomdata%potential%LPOT, &
      emesh%E1, emesh%E2, emesh%TK, emesh%EFERMI, &
      params%ALAT, atomdata%core%ITITLE(:,1:dims%NSPIND), &
      densities%total_charge_neutrality, &
      arrays%ZAT, emesh%EZ,&
!       emesh%WEZ,&
      params%LDAU, dims%iemxd, &
      dims%korbit)

      call outTime(mp%isMasterRank, 'results ........................', getTime(program_timer), iter)

      ! manual exit possible by creation of file 'STOP' in home directory
      processKKRresults = processKKRresults + 2*stopfile_flag()

    endif
  ! -----------------------------------------------------------------
  ! END: only MASTERRANK is working here
  ! -----------------------------------------------------------------

  endfunction ! processKKRresults


  !------------------------------------------------------------------------------
  !> Performs potential mixing.
  !>
  !> Returns 1 when params%target_rms is reached, otherwise 0
  !>
  !> Selection of algorithm according to value of params%imix
  !> imix = 0     straight mixing
  !> imix = 1, 2  not used - defaults to straight mixing
  !> imix = 3     Broyden's 1st method
  !> imix = 4     Broyden's 2nd method
  !> imix = 5     generalised Anderson
  !> imix = 6     Broyden's 2nd method with support for num_local_atoms > 1
  integer function mix_potential(calc, iter, params, dims, mp)

    use KKRnanoParallel_mod, only: KKRnanoParallel
    use CalculationData_mod, only: CalculationData, getAtomData
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData
    use broyden_kkr_mod, only: mix_broyden2_com
    use BRYDBM_new_com_mod, only: BRYDBM_new_com
    use wrappers_mod, only: MIXSTR_wrapper

    integer, intent(in)                   :: iter
    type(CalculationData), intent(inout)  :: calc
    type(KKRnanoParallel), intent(in)     :: mp
    type(DimParams), intent(in)           :: dims
    type(InputParams), intent(in)         :: params

    type(BasisAtom) , pointer             :: atomdata
    type(RadialMeshData), pointer         :: mesh
    double precision :: RMSAVQ, RMSAVM ! rms error charge dens. and mag. density (contribution of all local sites)
    double precision :: RMSAVQ_single, RMSAVM_single
    integer :: ila, num_local_atoms

    mix_potential = 0

    num_local_atoms = calc%num_local_atoms

    ! -->   calculation of RMS and final construction of the potentials (straight mixing)
    RMSAVQ = 0.d0
    RMSAVM = 0.d0

    ! straight/simple mixing
    !$omp parallel do reduction(+: RMSAVQ, RMSAVM) private(ila, atomdata, RMSAVQ_single, RMSAVM_single)
    do ila = 1, num_local_atoms
      atomdata => getAtomData(calc, ila)

      ! ATTENTION: the spherical part of the potential is divided by sqrt(4*pi) here - this is needed for the
      ! single site solver in the next iteration
      call MIXSTR_wrapper(atomdata, RMSAVQ_single, RMSAVM_single, params%MIXING, params%FCM)

      RMSAVQ = RMSAVQ + RMSAVQ_single
      RMSAVM = RMSAVM + RMSAVM_single
    enddo ! ila
    !$omp endparallel do

    ! summation and output of RMS error
    call RMSOUT_com(RMSAVQ, RMSAVM, ITER, dims%NSPIND, dims%NAEZ, mp%myAtomRank, mp%mySEComm)

    ! check if target rms error has been reached and set abort flag
    if (rmsavq <= params%target_rms) then
      mix_potential = 1
      if (mp%isMasterRank) write(*,*) "TARGET RMS ERROR REACHED..."
    endif

    ! straight mixing is called in any case before
    ! Broyden mixing - it is undone in Broyden routines

  ! -->   potential mixing procedures: Broyden or Andersen updating schemes
    if (params%IMIX >= 3 .and. params%IMIX < 6) then

      if (num_local_atoms > 1) die_here("Broyden mixing and num_local_atoms > 1 not supported.")

      ! Take data from 1st local atom, since only one local atom is supported
      atomdata     => getAtomData(calc, 1)
      mesh => atomdata%mesh_ptr

#define broyden calc%Broyden
      call BRYDBM_new_com(atomdata%potential%VISP,atomdata%potential%VONS, &
      atomdata%potential%VINS, &
      atomdata%potential%LMPOT,mesh%R,mesh%DRDI,broyden%MIXING, &
      mesh%IRC,mesh%IRMIN,atomdata%potential%NSPIN, &
      broyden%IMIX,ITER, &
      broyden%UI2,broyden%VI2,broyden%WIT,broyden%SM1S,broyden%FM1S, &
      mp%myAtomRank, &
      mp%mySEComm, &
      broyden%itdbryd, mesh%irmd, atomdata%potential%irnsd, &
      atomdata%potential%nspin)
#undef broyden

    ! this method supports num_local_atoms > 1
    else if (params%imix == 6) then
      ! use Broyden mixing that supports num_local_atoms > 1
      call mix_broyden2_com(calc, iter, mp%mySEComm)
    endif

  endfunction ! mix_potential

  !------------------------------------------------------------------------------
  !> Write forces to file 'forces'.
  !>
  !> Gather all forces at rank 'master', only this rank writes the file.
  !> Since the amount of data for forces is low this is a reasonable approach.
  subroutine output_forces(calc, master, rank, comm)
    use CalculationData_mod, only: CalculationData
    include 'mpif.h'
    type(CalculationData), intent(in) :: calc
    integer, intent(in) :: master, rank, comm

    integer :: ila, atom_id, num_local_atoms, nranks, ierr, max_local_atoms, ffu
    double precision, allocatable :: force_buffer(:,:), local_buffer(:,:)

    num_local_atoms = calc%num_local_atoms ! must be the same for all ranks, therefore get the maximum
    call MPI_Allreduce(num_local_atoms, max_local_atoms, 1, MPI_INTEGER, MPI_MAX, comm, ierr)

    if (rank == master) then
      call MPI_Comm_size(comm, nranks, ierr)
      allocate(force_buffer(3,max_local_atoms*nranks)) ! this could be large
    else
      allocate(force_buffer(1,1)) ! dummy allocation
    endif

    allocate(local_buffer(3,max_local_atoms))
    local_buffer = 0.d0
    do ila = 1, num_local_atoms
      local_buffer(:,ila) = calc%densities_a(ila)%force_flm(-1:1)
    enddo ! ila

    call MPI_Gather(local_buffer, 3*max_local_atoms, MPI_DOUBLE_PRECISION, &
                    force_buffer, 3*max_local_atoms, MPI_DOUBLE_PRECISION, &
                    master, comm, ierr)

    if (rank == master) then
      ffu = openForceFile()
      do atom_id = 1, max_local_atoms*nranks ! todo: limit to naez
        call writeForces(ffu, force_buffer(:,atom_id), atom_id)
      enddo ! atom_id
      close(ffu)
    endif

    deallocate(force_buffer, local_buffer, stat=ierr)

  endsubroutine ! output_forces

  !==============================================================================

  !------------------------------------------------------------------------------
  !> Calculate densities.
  !>
  !> output: emesh (Fermi-energy updated, renormalized weights), densities, ldau_data?
  !> files written: 'results1'
  subroutine calculateDensities(iter, calc, mp, dims, params, program_timer, arrays, emesh)

    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD

    use KKRnanoParallel_mod, only: KKRnanoParallel
    use EnergyMesh_mod, only: EnergyMesh
    use CalculationData_mod, only: CalculationData, getDensities, getAtomData
    use CalculationData_mod, only: getLDAUData, getEnergies
    use InputParams_mod, only: InputParams
    use Main2Arrays_mod, only: Main2Arrays
    use DimParams_mod, only: DimParams
    use TimerMpi_mod, only: TimerMpi, getTime, outTime, stopTimer, startTimer, outTimeStats
    
    use RadialMeshData_mod, only: RadialMeshData
    use ShapefunData_mod, only: ShapefunData
    use BasisAtom_mod, only: BasisAtom
    use LDAUData_mod, only: LDAUData
    use KKRresults_mod, only: KKRresults
    use DensityResults_mod, only: DensityResults
    use EnergyResults_mod, only: EnergyResults
    
    use lloyd0_new_mod, only: lloyd0_wrapper_com
    use wrappers_mod, only: rhoval_wrapper, RHOTOTB_wrapper, RHOMOM_NEW_wrapper
    use lloyds_formula_mod, only: renormalizeDOS

    integer, intent(in)                       :: iter
    type(CalculationData), intent(inout)      :: calc
    type(KKRnanoParallel), intent(in)         :: mp
    type(EnergyMesh), intent(inout)           :: emesh
    type(Main2Arrays), intent(in)             :: arrays
    type(DimParams), intent(in)               :: dims
    type(InputParams), intent(in)             :: params
    type(TimerMpi), intent(inout)             :: program_timer

    ! locals
    type(BasisAtom), pointer                  :: atomdata  ! not const
    type(LDAUData), pointer                   :: ldau_data ! not const
    type(DensityResults), pointer             :: densities ! not const
    type(EnergyResults), pointer              :: energies  ! not const
    type(RadialMeshData), pointer             :: mesh
    type(ShapefunData), pointer               :: cell

    double complex, parameter :: CZERO = (0.d0, 0.d0)
    logical :: LdoRhoEF
    integer :: atom_id
    double precision :: denEf !< charge density at Fermi level
    double precision :: chrgNt !< charge neutrality
    double precision :: denEf_local
    double precision :: chrgNt_local
    double precision :: new_fermi
    double precision :: CHRGSEMICORE !< total semicore charge over all atoms
    integer :: ila, r1fu
    integer :: num_local_atoms

    double complex, allocatable :: prefactors(:)  ! for Morgan charge test only

    num_local_atoms = calc%num_local_atoms

    atomdata  => getAtomData(calc, 1)
    ldau_data => getLDAUData(calc, 1)
    densities => getDensities(calc, 1)
    energies  => null()
    mesh      => null()
    cell      => null()

    atom_id = 0
    CHRGSEMICORE = 0.d0 !< Initialize semicore charge to zero

    if (dims%LLY /= 0 .and. num_local_atoms > 1) then
      if (mp%isMasterRank) write(*,*) "Lloyd's formula and num_local_atoms > 1 not supported."
      STOP
    endif

    ! out: emesh, RNORM
    call lloyd0_wrapper_com(atomdata, mp%mySEComm, calc%kkr_a(1)%LLY_GRDT, &
                            emesh, densities%RNORM, &
                            dims%LLY, params%ICST, params%NSRA, &
                            calc%kkr_a(1)%GMATN, calc%gaunts, ldau_data, params%Volterra)

    if (dims%LLY == 1) then
      TESTARRAYLOG(3, emesh%WEZRN)
      TESTARRAYLOG(3, densities%RNORM)
      call outTime(mp%isMasterRank, 'Lloyd processed ................', getTime(program_timer), ITER)
    endif

    ! now WEZRN stores the weights for E-integration

    DENEF = 0.d0
    CHRGNT = 0.d0

    LDORHOEF = (emesh%NPOL /= 0) ! needed in RHOVAL, 'L'ogical 'DO' RHO at 'E'-'F'ermi

    if (mp%isMasterRank .and. dims%korbit == 1) write(*,*) "Entering RHOVAL_wrapper! This might take some time, because NOCO (Bauer) solver is used..."

    call outTime(mp%isMasterRank, 'Lloyd done .....................', getTime(program_timer), ITER)

  !------------------------------------------------------------------------------
    !$omp parallel do reduction(+: chrgnt,denef, chrgsemicore) private(ila,atomdata,densities,energies,ldau_data,denef_local,chrgnt_local,atom_id)
    do ila = 1, num_local_atoms
      atomdata  => getAtomData(calc, ila)
      densities => getDensities(calc, ila)
      energies  => getEnergies(calc, ila)
      ldau_data => getLDAUData(calc, ila)
      atom_id = calc%atom_ids(ila) ! get global atom_id from local index
  !------------------------------------------------------------------------------

      ! has to be done after Lloyd
      ! output: RHO2NS, R2NEF, DEN, ESPV
      densities%DEN = CZERO
      densities%muorb = 0.0d0

      ! calculate valence charge density and band energies
      call RHOVAL_wrapper(atomdata, LdoRhoEF, params%ICST, params%NSRA, &
                          densities%RHO2NS, densities%R2NEF, &
                          densities%DEN, energies%ESPV, calc%kkr_a(ila)%GMATN, &
                          calc%gaunts, emesh, ldau_data, params%Volterra, &
                          dims%korbit, calc%noco_data%theta_noco(atom_id), calc%noco_data%phi_noco(atom_id), &
                          calc%noco_data%theta_noco_old(atom_id), calc%noco_data%phi_noco_old(atom_id), &
                          calc%noco_data%angle_fixed(atom_id), & 
                          calc%noco_data%moment_x(atom_id),calc%noco_data%moment_y(atom_id), calc%noco_data%moment_z(atom_id), &
                          densities%muorb, densities%iemxd, params)

      ! LDAU
      if (ldau_data%LDAU .and. ldau_data%NLDAU >= 1) then
          call LDAUWMAT(atom_id,ldau_data%NSPIND,ITER,params%MIXING,ldau_data%DMATLDAU, &
                        ldau_data%NLDAU,ldau_data%LLDAU, &
                        ldau_data%ULDAU,ldau_data%JLDAU, &
                        ldau_data%UMLDAU,ldau_data%WMLDAU, &
                        ldau_data%EULDAU,ldau_data%EDCLDAU, &
                        ldau_data%lmaxd)
      endif
      ! LDAU

    ! ----------------------------------------------------------------------
    ! -->   determine total charge expanded in spherical harmonics
    ! -------------------------------------------------------------- density
      ! output: CATOM, CATOM(1) = n_up + n_down, CATOM(2) = n_up - n_down
      !         core charge density added to densities%RHO2NS
      call RHOTOTB_wrapper(densities%CATOM, densities%RHO2NS, atomdata)

      CHRGNT_local = densities%CATOM(1) - atomdata%Z_nuclear

      if (dims%LLY == 1) then
        call renormalizeDOS(densities%DEN,densities%RNORM, &
                            densities%LMAXD+1,densities%IEMXD, &
                            densities%NSPIND,densities%IEMXD)
      endif

      ! calculate DOS at Fermi level
      DENEF_local = calcDOSatFermi(densities%DEN, densities%IEMXD, &
                                    densities%IEMXD, densities%LMAXD+1, &
                                    densities%NSPIND)

      ! ---> l/m_s/atom-resolved charges, output -> CHARGE
      ! Use WEZ or WEZRN ? - renormalisation already in DEN! (see renormalizeDOS)
      ! CHARGE -> written to result file

      ! Semicore contour included
      if (any(params%npntsemi > 0)) then

        call calcChargesLresSemi(densities%CHARGE, densities%CHRGSEMICORE_per_atom, densities%DEN, emesh%ielast, &
                            emesh%IESEMICORE, densities%LMAXD+1, densities%NSPIND, emesh%WEZ, densities%IEMXD)

      else
        ! Only valence contour
        call calcChargesLres(densities%CHARGE, densities%DEN, emesh%ielast, &
                            densities%LMAXD+1, densities%NSPIND, emesh%WEZ, densities%IEMXD)

      endif

      DENEF  = DENEF  + DENEF_local
      CHRGNT = CHRGNT + CHRGNT_local

      ! Add semicore charge for current atom
      CHRGSEMICORE = CHRGSEMICORE + densities%CHRGSEMICORE_per_atom

  !------------------------------------------------------------------------------
    enddo ! ila
    !$omp endparallel do
  !------------------------------------------------------------------------------

    call outTime(mp%isMasterRank, 'valence charge density .........', getTime(program_timer), ITER)
    call sumNeutralityDOSFermi_com(CHRGNT, DENEF, mp%mySEComm)

    ! write to 'results1' - only to be read in in results.f
    ! necessary for density of states calculation, otherwise
    ! only for informative reasons
    if (params%KTE >= 0) then
      r1fu = openResults1File(dims%IEMXD, dims%LMAXD, emesh%NPOL)

      do ila = 1, num_local_atoms
        atomdata  => getAtomData(calc, ila)
        densities => getDensities(calc, ila)
        atom_id = calc%atom_ids(ila) ! get global atom_id from local index

        call writeResults1File(r1fu, densities%CATOM, densities%CHARGE, densities%DEN, &
                              atomdata%core%ECORE, atom_id, emesh%NPOL, &
                              atomdata%core%QC_corecharge, densities%MUORB, &
                              calc%noco_data%phi_noco(atom_id), calc%noco_data%theta_noco(atom_id), &
                              calc%noco_data%phi_noco_old(atom_id), calc%noco_data%theta_noco_old(atom_id), &
                              calc%noco_data%angle_fixed(atom_id), &
                              calc%noco_data%moment_x(atom_id),calc%noco_data%moment_y(atom_id), &
                              calc%noco_data%moment_z(atom_id))
      enddo

      close(r1fu)
    endif

    call outTime(mp%isMasterRank, 'density calculated .............', getTime(program_timer), ITER)

  !------------------------------------------------------------------------------
    ! OMP had problems here, should be corrected now but not tested - therefore commented out
    !!!$omp parallel do private(ila, atomdata, densities, energies, mesh, cell) lastprivate(new_fermi)
    do ila = 1, num_local_atoms
      atomdata  => getAtomData(calc, ila)
      densities => getDensities(calc, ila)
      energies  => getEnergies(calc, ila)
      mesh      => atomdata%mesh_ptr
      cell      => atomdata%cell_ptr
      atom_id = calc%atom_ids(ila) ! get global atom_id from local index

  !=============== DEBUG: Morgan charge distribution test =======================
      if (params%DEBUG_morgan_electrostatics == 1) then
        if (mp%isMasterRank) call print_morgan_message()
        allocate(prefactors(size(arrays%rbasis,2)))
        call read_morgan_prefactors(prefactors)
        call overwrite_densities_gen_morgan(densities%RHO2NS, mesh%R, 2*dims%LMAXD, &
                      arrays%rbasis, arrays%rbasis(:,atom_id), arrays%bravais, prefactors)
        deallocate(prefactors)
        CHRGNT = 0.d0  ! don't do the Fermi energy correction
      endif
  !==============================================================================

      densities%total_charge_neutrality = CHRGNT

      new_fermi = emesh%E2

      ! allow only a maximal Fermi Energy shift of 0.03 Ry
      call doFermiEnergyCorrection(atomdata, mp%isMasterRank .and. (ila == 1), & ! show only once
                                  arrays%naez, 0.03d0, CHRGNT, DENEF, densities%R2NEF, &
                                  energies%ESPV, densities%RHO2NS, new_fermi)

      ! calculate multipole moments
      !output: CMOM, CMINST  ! only RHO2NS(:,:,1) passed (charge density)
      densities%CMOM   = 0.d0
      densities%CMINST = 0.d0

      call RHOMOM_NEW_wrapper(densities%CMOM,densities%CMINST, &
                              densities%RHO2NS(:,:,1), cell, mesh, calc%shgaunts)

  !------------------------------------------------------------------------------
    enddo ! ila
    !!!$omp endparallel do
  !------------------------------------------------------------------------------

    if (any(params%npntsemi > 0)) then
      ! --> Sum up semicore charges from different MPI ranks
      call sumChargeSemi_com(CHRGSEMICORE, mp%mySEComm)
      ! --> Recalculate the semicore contour factor FSEMICORE
      if (mp%isMasterRank) call calcFactorSemi(CHRGSEMICORE, emesh%FSEMICORE, params%fsemicore)
    endif

    emesh%E2 = new_fermi  ! Assumes that for every atom the same Fermi correction
                          ! was calculated !!!

    call outTime(mp%isMasterRank, 'rhomom .........................', getTime(program_timer), ITER)

#ifndef NOLOGGING
    ! log some results
    do ila = 1, num_local_atoms
      densities => getDensities(calc, ila)
      TESTARRAYLOG(3, calc%kkr_a(ila)%GMATN)
      TESTARRAYLOG(3, densities%CMOM)
      TESTARRAYLOG(3, densities%CMINST)
      TESTARRAYLOG(3, densities%RHO2NS)
    enddo ! ila
#endif

  endsubroutine ! calculateDensities


  !------------------------------------------------------------------------------
  !> Calculate potentials.
  !>
  !> Output: atomdata, ldau_data
  !> Files written: 'results2'
  subroutine calculatePotentials(iter, calc, mp, dims, params, program_timer, arrays)

    USE_LOGGING_MOD
    USE_ARRAYLOG_MOD

    use CalculationData_mod, only: CalculationData
    use KKRnanoParallel_mod, only: KKRnanoParallel
    use Main2Arrays_mod, only: Main2Arrays
    use DimParams_mod, only: DimParams
    use InputParams_mod, only: InputParams
    use TimerMpi_mod, only: TimerMpi, getTime, outTime
    
    use BasisAtom_mod, only: BasisAtom
    use LDAUData_mod, only: LDAUData
    use DensityResults_mod, only: DensityResults
    use EnergyResults_mod, only: EnergyResults
    use RadialMeshData_mod, only: RadialMeshData
    
    use CalculationData_mod, only: getDensities, getAtomData, getLDAUData, getEnergies, getDensities
    
    use NearField_calc_mod, only: add_near_field_corr
    use total_energy_mod, only: madelung_energy, madelung_ref_radius_correction, energy_electrostatic_wrapper
    use MadelungPotential_mod, only: addMadelungPotentialnew_com
    use muffin_tin_zero_mod, only: allreduceMuffinTinShift_com, printMuffinTinShift
    use AtomicForce_mod, only: force, forceh, forcxc
    use wrappers_mod, only: VINTRAS_wrapper, VXCDRV_wrapper, MTZERO_wrapper, CONVOL_wrapper
    
    integer, intent(in)                       :: iter
    type(CalculationData), intent(inout)      :: calc
    type(KKRnanoParallel), intent(in)         :: mp
    type(Main2Arrays), intent(in)             :: arrays
    type(DimParams), intent(in)               :: dims
    type(InputParams), intent(in)             :: params
    type(TimerMpi), intent(inout)             :: program_timer

    ! locals
    type(BasisAtom), pointer                  :: atomdata     ! not const
    type(LDAUData), pointer                   :: ldau_data    ! not const
    type(EnergyResults), pointer              :: energies     ! not const
    type(DensityResults), pointer             :: densities    ! not const
    type(RadialMeshData), pointer             :: mesh

    integer :: atom_id
    integer :: lcoremax
    double precision :: VAV0, VOL0
    double precision :: VAV0_local, VOL0_local
    double precision :: VBC_new(2)
    integer :: ila, ist, r2fu
    integer :: num_local_atoms
    logical :: calc_force
    double precision :: force_flmc(-1:1)

    double precision :: new_total_energy(2), new_total_energy_all(2)
    double precision, allocatable :: vons_temp(:,:,:)

    double complex, allocatable :: prefactors(:) ! for Morgan charge test only
    double precision :: direction(3)             ! for Morgan charge test only

    num_local_atoms = calc%num_local_atoms

    atomdata     => getAtomData(calc, 1)
    ldau_data    => getLDAUData(calc, 1)
    densities    => null()
    energies     => getEnergies(calc, 1)

    atom_id = 0

    new_total_energy = 0.d0
    
!#ifndef __GFORTRAN__
!    allocate(vons_temp, source=atomdata%potential%vons)
!#else
#define v atomdata%potential%vons
    allocate(vons_temp(size(v,1),size(v,2),size(v,3)))
!    vons_temp = v ! copy
!#undef v
!#endif
    
    calc_force = (params%KFORCE == 1) ! calculate force at each iteration

  !------------------------------------------------------------------------------
    !!!$omp parallel do private(ila, atomdata, densities)
    ! OpenMP problems?
    do ila = 1, num_local_atoms
      atomdata     => getAtomData(calc, ila)
      densities    => getDensities(calc, ila)
  !------------------------------------------------------------------------------

      !output: VONS
      call VINTRAS_wrapper(densities%RHO2NS(:,:,1), calc%shgaunts, atomdata)

  !------------------------------------------------------------------------------
    enddo ! ila
    !!!$omp endparallel do
  !------------------------------------------------------------------------------

    call outTime(mp%isMasterRank, 'vintras ........................', getTime(program_timer), ITER)

    ! perform near field correction for ALL local atoms
    if (params%near_field > 0) then
      call add_near_field_corr(calc, arrays, params%alat, mp%mySEComm)
      call outTime(mp%isMasterRank, 'near field .....................', getTime(program_timer),ITER)
    endif

#ifndef DEBUG_NO_ELECTROSTATICS
    ! output: VONS (changed), VMAD
    ! operation on all atoms! O(N**2)
    call addMadelungPotentialnew_com(calc, arrays%ZAT, mp%mySEComm)

    call outTime(mp%isMasterRank, 'vmadelblk ......................', getTime(program_timer), ITER)
#endif

  !=============== DEBUG: Morgan charge distribution test =======================
      if (params%DEBUG_morgan_electrostatics > 0 .and. mp%isMasterRank) then
        atomdata  => getAtomData(calc, 1)
        mesh => atomdata%mesh_ptr
        atom_id = calc%atom_ids(1) ! get global atom_id from local index
        allocate(prefactors(size(arrays%rbasis, 2)))
        call read_morgan_prefactors(prefactors)
        direction = read_direction()
!        call write_morgan_potential_dir(atomdata%potential%vons(:,:,1), mesh%R, direction)
!        call write_gen_morgan_potential_dir_analytical(mesh%R, arrays%rbasis, arrays%rbasis(:,atom_id), &
!                                                      arrays%bravais, prefactors, direction)
        deallocate(prefactors)
      endif
  !==============================================================================

    VAV0 = 0.d0
    VOL0 = 0.d0

    if (params%KTE >= 0) r2fu = openResults2File(dims%LRECRES2)

  !------------------------------------------------------------------------------
    !!!$nomp parallel do reduction(+: VAV0, VOL0) &
    !!!$nomp private(ila, atomdata, densities, energies, ldau_data, atom_id, &
    !!!$nomp         lcoremax, VAV0_local, VOL0_local, mesh, force_flmc)
    do ila = 1, num_local_atoms
      atomdata     => getAtomData(calc, ila)
      densities    => getDensities(calc, ila)
      energies     => getEnergies(calc, ila)
      ldau_data    => getLDAUData(calc, ila)
      mesh => atomdata%mesh_ptr
      atom_id = calc%atom_ids(ila) ! get global atom_id from local index
  !------------------------------------------------------------------------------

  ! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF  FORCES

      if (calc_force) then

        call forceh(densities%force_FLM,atomdata%potential%LPOT, &
                    densities%RHO2NS,atomdata%potential%VONS, &
                    mesh%R,mesh%DRDI, min(mesh%imt, MAX_MADELUNG_RADIUS_INDEX), atomdata%Z_nuclear,mesh%irmd)

        force_FLMC = 0.d0 ! temporary needed later in forcxc

        call force(densities%force_FLM,force_FLMC,atomdata%potential%LPOT, &
                  atomdata%potential%NSPIN, atomdata%core%RHOCAT, &
                  atomdata%potential%VONS, mesh%R, mesh%DRDI, &
                  mesh%IMT, mesh%irmd)

      endif

  ! Force Calculation stops here look after VXCDRV

  ! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF

    ! =====================================================================
      ! TODO: OpenMP critical !!! VXCDRV is most likely not threadsafe!
      ! output: vons_temp, EXC (exchange energy) (l-resolved)

      deallocate(vons_temp, stat=ist)
      allocate(vons_temp(size(v,1),size(v,2),size(v,3)))
      vons_temp = 0.d0 ! V_XC stored in temporary, must not add before energy calc.
      call VXCDRV_wrapper(vons_temp, energies%EXC, params%KXC, densities%RHO2NS, calc%shgaunts, atomdata)
  ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      call calculatePotentials_energies()
  ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE

      atomdata%potential%vons = atomdata%potential%vons + vons_temp

      deallocate(vons_temp, stat=ist)
      allocate(vons_temp(size(v,1),size(v,2),size(v,3)))
    ! Force calculation continues here

      if (calc_force) then
        call forcxc(densities%force_FLM,force_FLMC,atomdata%potential%LPOT, &
                    atomdata%potential%NSPIN, atomdata%core%RHOCAT, &
                    atomdata%potential%VONS, mesh%R, &
                    mesh%DRDI, mesh%IMT, mesh%irmd)
      endif

      ! writes some results (mostly energies) to direct access file 'results2',
      ! those are read in again in routine 'RESULTS'
      if (params%KTE >= 0) then
        call writeResults2File(r2fu, densities%CATOM, energies%ECOU, ldau_data%EDCLDAU, &
                              energies%EPOTIN, energies%ESPC, energies%ESPV, ldau_data%EULDAU, &
                              energies%EXC, atom_id, LCOREMAX, energies%VMAD)
      endif

      ! calculate new muffin-tin zero. output: VAV0, VOL0
      VAV0_local = 0.d0
      VOL0_local = 0.d0
      call MTZERO_wrapper(VAV0_local, VOL0_local, atomdata)
      VAV0 = VAV0 + VAV0_local
      VOL0 = VOL0 + VOL0_local

  !------------------------------------------------------------------------------
    enddo ! ila
    !!!$nomp endparallel do
  !------------------------------------------------------------------------------

    if (params%KTE >= 0) close(r2fu)

    call outTime(mp%isMasterRank, 'write results ..................', getTime(program_timer), ITER)

    call allreduceMuffinTinShift_com(mp%mySEComm, VAV0, VBC_new, VOL0)

    if (mp%isMasterRank) call printMuffinTinShift(VAV0, VBC_new, VOL0)

    !------------------------------------------------------------------------------
    !$omp parallel do private(ila, atomdata, energies, densities) reduction(+:new_total_energy)
    do ila = 1, num_local_atoms
      atomdata     => getAtomData(calc, ila)
      energies     => getEnergies(calc, ila)
      densities    => getDensities(calc, ila)
      !------------------------------------------------------------------------------

  ! -->   shift potential to muffin tin zero (average of interstitial potentials) and
  !       convolute potential with shape function for next iteration

  ! -->   shift potential by VBC and multiply with shape functions - output: VONS
  !       add also an optional muffin-tin-zero shift 'mt_zero_shift'

  !       Note: also the effect of the nuclear potential on the interstitial region
  !       is added in 'convol'
      energies%VBC = VBC_new + params%mt_zero_shift
      call CONVOL_wrapper(energies%VBC, calc%shgaunts, atomdata)

      ! MT-shift energy for Weinert energy only
      energies%e_shift = - energies%VBC(1) * densities%CATOM(1)
      energies%e_total(2) = energies%e_total(2) + energies%e_shift

      ! sum up over all atoms
      new_total_energy = new_total_energy + energies%e_total

      !write(*,*) energies%e_total
      !write(*,*) energies%e_vxc, energies%e_shift, energies%e_madelung

      !------------------------------------------------------------------------------
    enddo
    !$omp endparallel do
    !------------------------------------------------------------------------------

    call outTime(mp%isMasterRank, 'calculated pot .................', getTime(program_timer), ITER)

    call sum_total_energy_com(new_total_energy_all, new_total_energy, 0, mp%mySEComm)

    if (mp%isMasterRank) call printTotalEnergies(new_total_energy_all)

    deallocate(vons_temp, stat=ist)

    contains

    !----------------------------------------------------------------------------
    !> Inner subroutine: energy calculation.
    subroutine calculatePotentials_energies()
      use wrappers_mod, only: ESPCB_wrapper, EPOTINB_wrapper
      use total_energy_mod, only: energy_electrostatic_L_resolved_wrapper
      ! These energies have to be calculated BEFORE the XC-potential is added!
      ! calculate total energy and individual contributions if requested

      ! energies%e_total(1) stores result of total energy calculation using Harris functional
      !
      ! energies%e_total(2) stores the result of the method described in:
      ! Weinert, PRB 26, 4571 (1982)., which needs self-consistency to be accurate
      ! it does not use the input potential (EPOTIN) for the calculation

      double precision :: energy_temp

      energies%e_total = 0.d0

      ! core electron contribution
      call ESPCB_wrapper(energies%ESPC, LCOREMAX, atomdata)

      ! output: EPOTIN
      call EPOTINB_wrapper(energies%EPOTIN,densities%RHO2NS,atomdata)

      ! output: ECOU - l resolved Coulomb energy
      call energy_electrostatic_L_resolved_wrapper(energies%ECOU, atomdata%potential%vons, atomdata%Z_nuclear, densities%RHO2NS, calc%shgaunts, atomdata)

      ! coulomb energy and part of double counting
      energies%e_total(1) = energies%e_total(1) + sum(energies%ECOU) + energies%EPOTIN  ! Harris
      energies%e_total(2) = energies%e_total(2) - sum(energies%ECOU)                    ! Weinert

      ! madelung energy
      energy_temp = madelung_energy(atomdata%potential%vons(:,1,1), densities%rho2ns(:,1,1), &
                                    mesh%r, mesh%drdi, mesh%irmd, atomdata%Z_nuclear, min(mesh%imt, MAX_MADELUNG_RADIUS_INDEX))

      ! correction to madelung energy when reference radius is not muffin-tin radius
      energy_temp = energy_temp + madelung_ref_radius_correction(densities%rho2ns(:,1,1), &
                                    mesh%r, mesh%drdi, mesh%irmd, atomdata%Z_nuclear, min(mesh%imt, MAX_MADELUNG_RADIUS_INDEX), mesh%imt)

      energies%ECOU(0) = energies%ECOU(0) + energy_temp

      energies%e_total = energies%e_total + energy_temp

      energies%e_madelung = energy_temp

      ! core energies
      energies%e_total = energies%e_total + sum(sum(energies%ESPC, 2))

      ! single particle energies (band)
      energies%e_total = energies%e_total + sum(sum(energies%ESPV, 2))

      ! part of double counting energy that stems from V_XC (Weinert only)
      energies%e_vxc = -2.d0 * energy_electrostatic_wrapper(vons_temp, 0.d0, densities%RHO2NS, calc%shgaunts, atomdata)

      energies%e_total(2) = energies%e_total(2) + energies%e_vxc

      ! XC energy
      energies%e_total = energies%e_total + sum(energies%EXC)

      ! LDA+U contribution
      energies%e_total = energies%e_total + ldau_data%EULDAU - ldau_data%EDCLDAU

      ! missing contributions to energies%e_total(2), Weinert only: muffin-tin shift energy
    endsubroutine ! calculatePotentials_energies

  endsubroutine ! calculatePotentials

  !------------------------------------------------------------------------------
  integer function openForceFile() result(fu)
    integer :: reclen
    double precision :: dummy(-1:1)

    inquire(iolength=reclen) dummy ! get reclen for 3 doubles
    fu = 91
    open(unit=fu, access='direct', file='bin.forces', recl=reclen, form='unformatted', action='write')
  endfunction ! open

  !------------------------------------------------------------------------------
  subroutine writeForces(fu, forces_flm, recnr)
    integer, intent(in) :: fu !< file unit
    double precision, intent(in) :: forces_flm(-1:1)
    integer, intent(in) :: recnr

    write(unit=fu, rec=recnr) forces_flm
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Open Results2 File
  integer function openResults2File(lrecres2) result(fu)
    integer, intent(in) :: lrecres2
    fu = 72
    open(unit=fu, access='direct', recl=lrecres2, file='bin.results2', form='unformatted', action='write')
  endfunction ! open

  !----------------------------------------------------------------------------
  !> Write calculated stuff into historical 'results2' file
  subroutine writeResults2File(fu, catom, ecou, edcldau, epotin, espc, espv, euldau, exc, i1, lcoremax, vmad)
    integer, intent(in) :: fu !< file unit
    double precision, intent(in) :: catom(:), ecou(:), edcldau, epotin, espc(:,:), espv(:,:), euldau, exc(:)
    integer, intent(in) :: i1, lcoremax
    double precision, intent(in) :: vmad

    write(unit=fu, rec=i1) catom,vmad,ecou,epotin,espc,espv,exc,lcoremax,euldau,edcldau
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Open the file 'results1'
  integer function openResults1File(iemxd, lmaxd, npol) result(fu)
    integer, intent(in) :: iemxd, lmaxd, npol

    integer :: lrecres1

    !lrecres1 = 8*43 + 16*(lmaxd+2)
    !lrecres1 = 8*43 + 16*(lmaxd+2) + 8*(lmaxd+3)*3 + 8*2 + 1*1 ! NOCO with nonco angles and angle_fixed
    lrecres1 = 8*43 + 16*(lmaxd+2) + 8*(lmaxd+3)*3 + 8*4 + 1*1 + 8*3 ! NOCO with old and new nonco angles, angle_fixed and moments
    if (npol == 0) lrecres1 = lrecres1 + 32*(lmaxd+2)*iemxd

    fu = 71
    open(unit=fu, access='direct', recl=lrecres1, file='bin.results1', form='unformatted', action='write')
  endfunction ! open

  !----------------------------------------------------------------------------
  !> Write some stuff to the 'results1' file
  subroutine writeResults1File(fu, catom, charge, den, ecore, i1, npol, qc, &
                               muorb, phi_soc, theta_soc, phi_soc_old, &
                               theta_soc_old, angle_fixed, &
                               moment_x, moment_y, moment_z)
                           
    integer, intent(in) :: fu !< file unit
    double precision, intent(in) :: catom(:), charge(:,:)
    double complex, intent(in) :: den(:,:,:)
    double precision, intent(in) :: ecore(20,2)
    integer, intent(in) :: i1, npol
    double precision, intent(in) :: qc
    double precision, intent(in) :: muorb(:,:)    ! NOCO
    double precision, intent(in) :: phi_soc       ! NOCO
    double precision, intent(in) :: theta_soc     ! NOCO
    double precision, intent(in) :: phi_soc_old   ! NOCO
    double precision, intent(in) :: theta_soc_old ! NOCO
    integer (kind=1), intent(in) :: angle_fixed   ! NOCO
    double precision, intent(in) :: moment_x      ! NOCO
    double precision, intent(in) :: moment_y      ! NOCO
    double precision, intent(in) :: moment_z      ! NOCO

    if (npol == 0) then
      write(unit=fu, rec=i1) qc,catom,charge,ecore,muorb,phi_soc,theta_soc,phi_soc_old,theta_soc_old,angle_fixed, &
              moment_x,moment_y,moment_z,den  ! write density of states (den) only when certain options set
    else
      write(unit=fu, rec=i1) qc,catom,charge,ecore,muorb,phi_soc,theta_soc,phi_soc_old,theta_soc_old,angle_fixed, &
              moment_x,moment_y,moment_z
    endif
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Communicate and sum up contributions for charge neutrality and
  !> density of states at Fermi level.
  subroutine sumNeutralityDOSFermi_com(chrgnt, denef, communicator)
    include 'mpif.h'
    double precision, intent(inout) :: chrgnt, denef
    integer, intent(in) :: communicator

    double precision :: work1(2), work2(2)
    integer :: ierr

    work1(1:2) = [chrgnt, denef]

    call MPI_Allreduce(work1,work2,2,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

    chrgnt = work2(1); denef  = work2(2)

  endsubroutine ! sum

  !----------------------------------------------------------------------------
  !> Communicate and sum up contributions to semicore charge

  subroutine sumChargeSemi_com(chrgsemicore, communicator)
    include 'mpif.h'
    double precision, intent(inout) :: chrgsemicore
    integer, intent(in) :: communicator

    double precision :: buff
    integer :: ierr

    call MPI_Allreduce(chrgsemicore,buff,1,MPI_DOUBLE_PRECISION,MPI_SUM,communicator,ierr)

    chrgsemicore = buff

  endsubroutine ! sum

  !----------------------------------------------------------------------------
  !> Extracts the DOS at the Fermi level and stores result in DENEF
  !> @param IELAST energy index of Fermi level
  !> @param LMAXD1 lmax+1
  double precision function calcDOSatFermi(den, ielast, iemxd, lmaxd1, nspin) result(denef)
    use Constants_mod, only: pi
    double complex, intent(in) :: den(0:lmaxd1,iemxd,nspin)
    integer, intent(in) :: ielast
    integer, intent(in) :: iemxd
    integer, intent(in) :: lmaxd1
    integer, intent(in) :: nspin

    integer :: ispin, l

    ! get density of states at fermi-level
    denef = 0.d0
    do ispin = 1, nspin
      do l = 0, lmaxd1
        denef = denef - 2.d0*dimag(den(l,ielast,ispin))/(pi*dble(nspin))
      enddo ! l
    enddo ! ispin

  endfunction ! calcDOSatFermi

  !----------------------------------------------------------------------------
  !> Calculates the l-resolved charges.
  !> This is done by energy integration in the complex plane over the imaginary
  !> part of the diagonal of the structural Green's function (DEN)

  subroutine calcChargesLres(charge, den, ielast, lmaxd1, nspin, wez, iemxd)
    double precision, intent(out) :: charge(0:lmaxd1,2)
    double complex, intent(in) :: den(0:lmaxd1,iemxd,nspin)
    doublecomplex, intent(in) :: wez(iemxd)

    integer, intent(in) :: iemxd, nspin, lmaxd1, ielast

    integer :: ispin, l, ie

    ! ---> l/m_s/atom-resolved charges
    do ispin = 1, nspin
      do l = 0, lmaxd1
        charge(l,ispin) = 0.d0
        do ie = 1, ielast
          charge(l,ispin) = charge(l,ispin) + dimag(wez(ie)*den(l,ie,ispin))/dble(nspin)
        enddo ! ie
      enddo ! l
    enddo ! ispin
    
  endsubroutine ! calcChargesLres

   !----------------------------------------------------------------------------
  !> Calculates the l-resolved charges including the semicore charge.
  !> This is done by energy integration in the complex plane over the imaginary
  !> part of the diagonal of the structural Green's function (DEN)
  !> TO DO: Combine calcChargesLres and calcChargesLresSemi once semicore feature is stable!

  subroutine calcChargesLresSemi(charge, chrgsemicore, den, ielast, iesemicore, lmaxd1, nspin, wez, iemxd)
    double precision, intent(out) :: charge(0:lmaxd1,2)
    double precision, intent(out) :: chrgsemicore
    double complex, intent(in) :: den(0:lmaxd1,iemxd,nspin)
    doublecomplex, intent(in) :: wez(iemxd)

    integer, intent(in) :: iemxd
    integer, intent(in) :: nspin
    integer, intent(in) :: lmaxd1
    integer, intent(in) :: ielast
    integer, intent(in) :: iesemicore

    integer :: ispin, l, ie

    ! ---> l/m_s/atom-resolved charges
    chrgsemicore = 0.d0
    do ispin = 1, nspin
      do l = 0, lmaxd1
        charge(l,ispin) = 0.d0
        do ie = 1, ielast
          charge(l,ispin) = charge(l,ispin) + dimag(wez(ie)*den(l,ie,ispin))/dble(nspin)
          if (ie == iesemicore) chrgsemicore = chrgsemicore + charge(l,ispin)
        enddo ! ie
      enddo ! l
    enddo ! ispin
    
  endsubroutine ! calcChargesLresSemi

  !----------------------------------------------------------------------------
  !> Calculates the normalization factor for the semicore contour (FSEMICORE) analogously to the JM-Code

  subroutine calcFactorSemi(chrgsemicore, fsemicore, fsemicore_old)
    double precision, intent(inout) :: chrgsemicore ! semicore charge
    double precision, intent(inout) :: fsemicore    ! semicore factor to be updated
    double precision, intent(in) :: fsemicore_old    ! initial semicore factor (default = 1.0d0)
    
    integer :: nsb !> number of semicore bands

    if (chrgsemicore < 1d-10) chrgsemicore = 1d-10

    nsb = nint(chrgsemicore)
    fsemicore = dble(nsb)/chrgsemicore*fsemicore_old

    write(6, '(6X,"< SEMICORE > : ",/,21X,"charge found in semicore :",F10.6,/,21X,"new normalisation factor :",F20.16,/)') chrgsemicore, fsemicore
  endsubroutine ! calcFactorSemi

  !----------------------------------------------------------------------------
  !> correct Fermi-energy (charge neutrality).
  !>
  !> Modifies charge density, Fermi energy and valence band energy!
  subroutine doFermiEnergyCorrection(atomdata, output, naez, max_shift, chrgnt, denef, r2nef, espv, rho2ns, e2)
    use BasisAtom_mod, only: BasisAtom
    use RadialMeshData_mod, only: RadialMeshData
    use Constants_mod, only: pi

    type(BasisAtom), intent(in) :: atomdata
    logical, intent(in) :: output !< output to stdout - yes/no
    integer, intent(in) :: naez
    double precision, intent(in) :: max_shift !< maximally allowed fermi-energy shift (good: 0.03d0)
    double precision, intent(in) :: chrgnt, denef
    double precision, intent(in) :: r2nef(:,:,:)
    double precision, intent(inout) :: espv(0:,:), rho2ns(:,:,:)
    double precision, intent(inout) :: e2

!   type(RadialMeshData), pointer :: mesh
    double precision :: e2shift, efold, df
    integer :: nspind, lmpotd!, ispin, lm

    nspind = atomdata%nspin
    lmpotd = atomdata%potential%lmpot

!   mesh => atomdata%mesh_ptr
#define mesh atomdata%mesh_ptr
    
    ! --> determine new fermi level due to valence charge up to old fermi level e2 and density of states denef
    e2shift = min(max(-max_shift, chrgnt/denef), max_shift)
    efold = e2
    e2 = e2 - e2shift

    if (output) call printFermiEnergy(denef/dble(naez), e2, e2shift, efold)

    espv(0,1:nspind) = espv(0,1:nspind) - efold*chrgnt/dble(nspind*naez) ! get correct density and valence band energies
    df = 2.d0*e2shift/(pi*dble(nspind))
    rho2ns(1:mesh%irc,1:lmpotd,1:nspind) = rho2ns(1:mesh%irc,1:lmpotd,1:nspind) + df * r2nef(1:mesh%irc,1:lmpotd,1:nspind)

!    ! old code
!    do ispin = 1, nspind
! -->     get correct density and valence band energies
!       espv(0,ispin) = espv(0,ispin) - efold*chrgnt/dble(nspind*naez)
!       do lm = 1, lmpotd
!         call daxpy(mesh%irc, df, r2nef(1,lm,ispin), 1, rho2ns(1,lm,ispin), 1)
!       enddo ! lm
!    enddo ! ispin
#undef  mesh
  endsubroutine ! doFermiEnergyCorrection

  !---------------------------------------------------------------------------
  !> Checks for file 'STOP' in current working directory, returns 1 if it
  !> exists, otherwise 0.
  !>
  !> Should be called only by master rank!!!
  integer function stopfile_flag()
    logical :: stopfile_exists
    stopfile_exists = .false.
    
    inquire(file='STOP', exist=stopfile_exists)
    stopfile_flag = 0; if (stopfile_exists) stopfile_flag = 1
  endfunction ! stopfile_flag

  !----------------------------------------------------------------------------
  !> Print Fermi-Energy information to screen.
  subroutine printFermiEnergy(denef, e2, e2shift, efold)
    double precision, intent(in) :: denef, e2, e2shift, efold

    write(6,fmt="('                old E FERMI ',F12.6,' Delta E_F = ',f12.6)") efold, e2shift
    ! --> divided by naez because the weight of each atom has been already taken into account in 1c
    write(6,fmt="('                new E FERMI ',F12.6,'  DOS(E_F) = ',f12.6)") e2, denef
    write(6,'(79(1h+),/)')
  endsubroutine ! print

  !----------------------------------------------------------------------------
  !> Sums up the total energies of all processes - returns result in 'total' on master, 0.d0 on
  !> other ranks.
  subroutine sum_total_energy_com(total, total_energies, master, communicator)
    include 'mpif.h'
    integer, parameter :: ND=2
    double precision, intent(out) :: total(ND)
    double precision, intent(in) :: total_energies(ND)
    integer, intent(in) :: master
    integer, intent(in) :: communicator
    
    integer :: ierr

    total = 0.d0
    call MPI_Reduce(total_energies, total, ND, MPI_DOUBLE_PRECISION, MPI_SUM, master, communicator, ierr)

  endsubroutine ! sum

  !----------------------------------------------------------------------------
  !> Print total energy to screen (both methods: Harris and Weinert).
  subroutine printTotalEnergies(total_energies)
    double precision, intent(in) :: total_energies(2)

    write(*,"(/,3x,70('+'),/,15x,'TOTAL ENERGY in ryd. : ',        f21.8,/15x,'                 eV  : ',f21.8,/,3x,70('+'))") &
                  total_energies(1), total_energies(1)*13.6058D0
    write(*,"(/,3x,70('+'),/,15x,'Weinert total energy in ryd. : ',f21.8,/15x,'                 eV  : ',f21.8,/,3x,70('+'))") &
                  total_energies(2), total_energies(2)*13.6058D0
  endsubroutine ! print

!=============== DEBUG: Morgan charge distribution test =======================

  !----------------------------------------------------------------------------
  !> Print warning message.
  subroutine print_morgan_message()
    write(*,*) "==============================================================="
    write(*,*) "= DEBUG: MORGAN charge distribution test activated.           ="
    write(*,*) "= Set option DEBUG_morgan_electrostatics = 0 to deactivate    ="
    write(*,*) "= Results of calculation are wrong.                           ="
    write(*,*) "==============================================================="
    warn(6, "MORGAN charge distribution test activated")
  endsubroutine ! print

  !----------------------------------------------------------------------------
  !> Write results of potential in direction 'dir' to a file.
!  subroutine write_morgan_potential_dir(vons, mesh_points, dir)
!    use debug_morgan_mod, only: eval_expansion
!    
!    double precision, intent(in) :: vons(:,:), mesh_points(:), dir(3)
!    integer, parameter :: fu = 99
!
!    double precision :: vec(3), norm_dir(3), val
!    integer :: ii
!
!    norm_dir = dir/norm2(dir)
!
!    open(unit=fu, form='formatted', file='morgan_potential_dir.txt', action='write')
!    do ii = 1, size(mesh_points)
!      vec = norm_dir * mesh_points(ii)
!      if (norm2(vec) == 0.d0) vec(1) = 1d-6
!      val = eval_expansion(vons(ii,:), vec)
!      write(unit=fu, fmt=*) mesh_points(ii), val
!    enddo ! ii
!    close(unit=fu)
!
!  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Stores generalised Morgan test charge distribution into rho2ns_density.
  !>
  !> One can specify basis atoms.
  subroutine overwrite_densities_gen_morgan(rho2ns_density, mesh_points, lpot, rbasis, center, bravais, prefactors)
    use debug_morgan_mod, only: calc_gen_morgan_rho_expansion, calc_reciprocal_basis, calc_reciprocal_first_shell

    double precision, intent(inout) :: rho2ns_density(:,:,:)
    double precision, intent(in) :: mesh_points(:)
    integer, intent(in) :: lpot
    double precision, intent(in) :: rbasis(:,:)
    double precision, intent(in) :: center(3)
    double precision, intent(in) :: bravais(3,3)
    double complex, intent(in)   :: prefactors(:)

    double precision, allocatable :: reciprocals(:,:)
    double precision :: rec_basis(3,3)

    integer :: ii
    double complex, allocatable :: coeffs(:)

    allocate(coeffs(size(rho2ns_density, 2)))

    call calc_reciprocal_basis(rec_basis, bravais)
    call calc_reciprocal_first_shell(reciprocals, rec_basis) ! reciprocals gets allocated here

    do ii = 1, size(mesh_points)
      call calc_gen_morgan_rho_expansion(coeffs, reciprocals, prefactors, rbasis, center, mesh_points(ii), lpot)

      ! since we are using real spherical harmonics, one can just take the real part of the expansion coeffs
      rho2ns_density(ii,:,1) = dreal(coeffs)
    enddo ! ii

    ! multiply with r**2 (mesh_points = r)
    do ii = 1, size(rho2ns_density, 2)
      rho2ns_density(:,ii,1) = rho2ns_density(:,ii,1) * mesh_points * mesh_points
    enddo ! ii

    if (size(rho2ns_density, 3) > 1) rho2ns_density(:,:,2:) = 0.d0

    deallocate(reciprocals, coeffs, stat=ii)
  endsubroutine ! overwrite

  !----------------------------------------------------------------------------
  !> Read file 'morgan_prefactors.dat'.
  !>
  !> If the file does not exist, prefactors are set to 1.0
  !> for each basis atom there is a line with 1 complex prefactor
  subroutine read_morgan_prefactors(prefactors)
    double complex, intent(out) :: prefactors(:)

    double precision :: re, im
    integer :: ii, ios
    integer, parameter :: fu=99

    prefactors = dcmplx(1.d0, 0.d0)

    open(unit=fu, form='formatted', file='morgan_prefactors.dat', action='read', status='old', iostat=ios)
    if (ios /= 0) return ! (1.0, 0.0) for all prefactors
    
    do ii = 1, size(prefactors)
      read(unit=fu, fmt=*) re, im
      prefactors(ii) = dcmplx(re, im)
    enddo ! ii
    close(unit=fu)
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Read file a direction vector from file 'directions.dat'.
  !>
  !> If the file does not exist, prefactors are set to 1.0
  function read_direction() result(dir)
    double precision :: dir(3)
    integer, parameter :: fu=99
    open(unit=fu, form='formatted', file='directions.dat', action='read', status='old')
    read(unit=fu, fmt=*) dir(1:3)
    close(unit=fu)
  endfunction ! read

  !----------------------------------------------------------------------------
  !> Write analytical values of generalised morgan potential in direction 'dir' to a file.
!  subroutine write_gen_morgan_potential_dir_analytical(mesh_points, rbasis, center, bravais, prefactors, dir)
!    use debug_morgan_mod, only: eval_gen_morgan_potential, calc_reciprocal_basis, calc_reciprocal_first_shell
!    
!    double precision, intent(in) :: mesh_points(:)
!    double precision, intent(in) :: rbasis(:,:)
!    double precision, intent(in) :: center(3)
!    double precision, intent(in) :: bravais(3, 3)
!    double complex, intent(in) :: prefactors(:)
!    double precision, intent(in) :: dir(3)
!
!    integer, parameter :: fu=99
!    integer :: ii
!    double precision :: vec(3), norm_dir(3), val, rec_basis(3,3)
!    double precision, allocatable :: reciprocals(:,:)
!
!    call calc_reciprocal_basis(rec_basis, bravais)
!    call calc_reciprocal_first_shell(reciprocals, rec_basis) ! reciprocals gets allocated here
!
!    norm_dir = dir/norm2(dir)
!
!    ! also write the analytical solution
!    open(unit=fu, form='formatted', file='gen_morgan_potential_dir_analytical.txt', action='write')
!    do ii = 1, size(mesh_points)
!      vec = norm_dir * mesh_points(ii)
!      val = real(eval_gen_morgan_potential(reciprocals, vec, prefactors, rbasis, center))
!      write(unit=fu, fmt=*) mesh_points(ii), val
!    enddo ! ii
!    close(unit=fu)
!
!  endsubroutine ! write


! process with MYLRANK(LMPIC) == 0 and LMPIC == 1 writes results

  subroutine results(lrecres2, ielast, itscf, lmax, natoms, npol, nspin, kpre, compute_total_energy, lpot, e1, e2, tk, efermi, alat, ititle, chrgnt, zat, ez, &
!     wez, &
    ldau, iemxd, &
    korbit)
  use Constants_mod, only: pi
    integer, intent(in) :: iemxd
    integer, intent(in) :: ielast, itscf, lmax, natoms, npol, nspin
    integer, intent(in) :: kpre
    integer, intent(in) :: compute_total_energy ! former kte
    double precision, intent(in) :: e1, e2, tk, efermi
    double precision, intent(in) :: chrgnt, alat
    logical, intent(in) :: ldau
    double complex, intent(in) :: ez(iemxd)!, wez(iemxd)
    double precision, intent(in) :: zat(natoms)
    integer, intent(in) :: ititle(20,*)
    integer, intent(in) :: korbit ! NOCO
    
!   logical, external :: TEST
#define TEST(STRING) .false.

    !     .. locals ..
    double complex :: den(0:lmax+1,iemxd,nspin)
    double precision :: dostot(0:lmax+1,2)
    double precision :: ecou(0:lpot), epotin, euldau, edcldau, espc(0:3,nspin), espv(0:lmax+1,nspin), exc(0:lpot)
    double precision :: ecore(20,2)
    double precision :: charge(0:lmax+1,2)
    double precision :: catom(nspin), qc
    double precision :: vmad, totsmom 
    double precision muorb(0:LMAX+2,3) !NOCO
    double precision phi_noco !NOCO
    double precision theta_noco !NOCO
    double precision phi_noco_old !NOCO
    double precision theta_noco_old !NOCO
    integer (kind=1) angle_fixed !NOCO
    double precision moment_x !NOCO
    double precision moment_y !NOCO
    double precision moment_z !NOCO
    double precision max_delta_theta !NOCO
    double precision max_delta_phi !NOCO
    double precision max_delta_angle !NOCO
    double precision delta_angle !NOCO
    integer :: max_delta_atom !NOCO
    integer :: lrecres1, lrecres2
    integer :: lcoremax, i1, ispin, lpot
    character(len=*), parameter :: &
    F90="('  Atom ',I4,' charge in Wigner Seitz cell =',f10.6)", &
    F91="(7X,'spin moment in Wigner Seitz cell =',f10.6)", &
    F92="('      ITERATION',I4,' charge neutrality in unit cell = ',f12.6)", &
    F93="('                    TOTAL mag. moment in unit cell = ',f12.6)", &
    F86="('                    Largest spin moment direction change for atom   = ',i5.1)", &
    F87="('                    Angle between old and new moment (deg)  = ',f12.6)", &
    F88="('                    Change of angle theta (deg)  = ',f12.6)", &
    F89="('                    Change of angle phi (deg)  = ',f12.6)", &
    F94="(4X,'nuclear charge  ',F10.6,9X,'core charge =   ',F10.6)"
    
    integer :: npotd
    npotd = nspin*natoms

    delta_angle = 0.0d0    !NOCO
    max_delta_theta = 0.d0 !NOCO
    max_delta_phi = 0.d0   !NOCO
    max_delta_angle = 0.d0 !NOCO
    max_delta_atom = 1     !NOCO

    !lrecres1 = 8*43 + 16*(lmax+2) ! w/o NOCO
    !lrecres1 = 8*43 + 16*(lmax+2) + 8*(lmax+3)*3 + 8*2 + 1*1 ! NOCO with noco angles and angle_fixed option
    lrecres1 = 8*43 + 16*(lmax+2) + 8*(lmax+3)*3 + 8*4 + 1*1 + 8*3 ! NOCO with old and new noco angles, angle_fixed option and moments
    
    if (npol == 0) lrecres1 = lrecres1 + 32*(lmax+2)*iemxd ! dos calc.

    if (compute_total_energy >= 0) then
      open(71, access='direct', recl=lrecres1, file='bin.results1', form='unformatted', action='read', status='old')
      if (korbit == 1) open(13,file='nonco_angle_out.dat',form='formatted') ! NOCO
      if (korbit == 1) open(14,file='nonco_moment_out.txt',form='formatted') ! NOCO
    
      ! moments output
      do i1 = 1, natoms
        if (npol == 0) then 
          read(71, rec=i1) qc,catom,charge,ecore,muorb,phi_noco,theta_noco,phi_noco_old,theta_noco_old,angle_fixed, &
                  moment_x,moment_y,moment_z,den
        else
          read(71, rec=i1) qc,catom,charge,ecore,muorb,phi_noco,theta_noco,phi_noco_old,theta_noco_old,angle_fixed, &
                  moment_x,moment_y,moment_z
        endif
       
        call wrmoms(nspin, charge, muorb, i1, lmax, lmax+1, i1 == 1, i1 == natoms)! first=(i1 == 1), last=(i1 == natoms))
        if (korbit == 1) then ! NOCO

           delta_angle = acos(sin(theta_noco)*sin(theta_noco_old)*cos(phi_noco-phi_noco_old)+ &
                         cos(theta_noco)*cos(theta_noco_old))
           if (abs(delta_angle) >= max_delta_angle) then
             max_delta_atom = i1
!             write(*,*) 'max_delta_atom= ', i1
             max_delta_angle = abs(delta_angle)
             max_delta_theta = abs(theta_noco_old-theta_noco)
             max_delta_phi = abs(phi_noco_old-phi_noco)
           endif
!           max_delta_angle = max(max_delta_angle,abs(delta_angle)) 
!           max_delta_theta = max(max_delta_theta,abs(theta_noco-theta_noco_old)) 
!           max_delta_phi   = max(max_delta_phi,abs(phi_noco-phi_noco_old))

          ! save to 'nonco_angle_out.dat' in converted units (degrees)
          write(13,*) theta_noco/(2.0D0*PI)*360.0D0, &
                      phi_noco/(2.0D0*PI)*360.0D0, &
                      angle_fixed

          ! save extended information to 'nonco_moment_out.txt', e.g. for
          ! visualization
          write(14,"(6f12.5,1i5)") moment_x, &
                      moment_y, &
                      moment_z, &
                      sqrt(moment_x**2+moment_y**2+moment_z**2), &
                      theta_noco/(2.0D0*PI)*360.0D0, &
                      phi_noco/(2.0D0*PI)*360.0D0, &
                      angle_fixed
        endif
      enddo ! i1

      if (korbit == 1)  close(13)
      if (korbit == 1)  close(14)

      ! density of states output
      if (npol == 0) then
        do i1 = 1, natoms
          read(71, rec=i1) qc, catom, charge, ecore, den
!         call wrldos(den, ez, wez, lmax+1, iemxd, npotd, ititle, efermi, e1, e2, alat, tk, nspin, natoms, ielast, i1, dostot)
          call wrldos(den, ez, lmax+1, iemxd, ititle, efermi, e1, e2, alat, tk, nspin, natoms, ielast, i1, dostot)
        enddo ! i1
      endif


      totsmom = 0.d0
      do i1 = 1, natoms
        if (npol == 0) then
          read(71, rec=i1) qc, catom, charge, ecore, den
        else
          read(71, rec=i1) qc, catom, charge, ecore
        endif
        do ispin = 1, nspin
          if (ispin /= 1) then
            write(6, fmt=F91) catom(ispin)                  ! spin moments
          else
            write(6, fmt=F90) i1, catom(ispin)              ! atom charge
          endif
        enddo ! ispin
        write(6, fmt=F94) zat(i1), qc                        ! nuclear charge, total charge
        if (nspin == 2) totsmom = totsmom + catom(nspin)
      enddo ! i1
      write(6, '(79(1h+))')
      write(6, fmt=F92) itscf,chrgnt                        ! charge neutrality
      if (nspin == 2) write(6, fmt=F93) totsmom             ! total mag. moment
      if (korbit == 1) write(6, fmt=F86) max_delta_atom               ! atom with largest spin moment direction change, NOCO
      if (korbit == 1) write(6, fmt=F87) 180.0/PI*max_delta_angle     ! largest spin moment direction change, NOCO
      if (korbit == 1) write(6, fmt=F88) 180.0/PI*max_delta_theta     ! Corresponding theta angle change, NOCO
      if (korbit == 1) write(6, fmt=F89) 180.0/PI*max_delta_phi       ! Corresponding phi angle change, NOCO
      write(6, '(79(1h+))')

      close(71)
    endif

    if (compute_total_energy == 1) then
      open(72, access='direct', recl=lrecres2, file='bin.results2', form='unformatted', action='read', status='old')
      do i1 = 1, natoms
        read(72, rec=i1) catom, vmad, ecou, epotin, espc, espv, exc, lcoremax, euldau, edcldau
        ! output unfortunately integrated into etotb1
        ! etotb1 depends on 'saved' variables !!!
        call etotb1(ecou, epotin, espc, espv, exc, euldau, edcldau, ldau, kpre, lmax, lpot, lcoremax, nspin, i1, natoms)
      enddo ! i1
      close(72)
    endif

  endsubroutine ! results
 

endmodule ! ProcessKKRresults_mod

