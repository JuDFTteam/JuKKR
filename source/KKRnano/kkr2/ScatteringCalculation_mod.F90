#include "DebugHelpers/logging_macros.h"
#include "DebugHelpers/test_array_log.h"
#include "DebugHelpers/test_macros.h"

!> @author Modularisation: Elias Rabel
module ScatteringCalculation_mod

public  :: energyLoop
private :: printEnergyPoint
private :: calcDeltaTupTdown
private :: substractReferenceTmatrix

CONTAINS

!> Input:  *) ebalance_handler (properly initialised)
!>         *) jij_data, ldau_data (properly initialised !!!)
!>
!> Output: ebalance_handler (changed, updated timings and process distribution)
!>         KKRResults
!>         jij_data results of jij-calculation
!>         ldau_data LDA+U results
!>
!>         FILES WRITTEN:
!>         *) Logfiles (if requested)
!>         *) JIJ-Files (if requested)
!>         *) matrix dump (if requested)
subroutine energyLoop(iter, calc_data, emesh, params, dims, &
                      ebalance_handler, my_mpi, arrays)

  USE_LOGGING_MOD
  USE_ARRAYLOG_MOD

  use TimerMpi_mod
  use EBalanceHandler_mod

  use TEST_lcutoff_mod !TODO: remove

  use DimParams_mod
  use InputParams_mod
  use Main2Arrays_mod

  use CalculationData_mod
  use KKRresults_mod
  use GauntCoefficients_mod
  use BasisAtom_mod
  use JijData_mod
  use LDAUData_mod
  use EnergyMesh_mod
  use KKRnanoParallel_mod
  use EBalanceHandler_mod

  use KKRnano_Comm_mod
  use kloopz1_mod
  use TruncationZone_mod

  use wrappers_mod,     only: calctmat_wrapper, calcdtmat_wrapper

  implicit none

  integer, intent(in) :: iter
  type (CalculationData), intent(inout) :: calc_data
  type (KKRnanoParallel), intent(in)    :: my_mpi
  type (EBalanceHandler), intent(inout) :: ebalance_handler
  type (EnergyMesh), intent(in)         :: emesh
  type (Main2Arrays), intent(in)        :: arrays
  type (DimParams), intent(in)          :: dims
  type (InputParams), intent(in)        :: params

  !---------- locals ----------------------------

  type (BasisAtom), pointer             :: atomdata  ! referenced data does not change
  type (KKRresults), pointer            :: kkr       ! changes
  type (GauntCoefficients), pointer     :: gaunts    ! never changes
  type (LDAUData), pointer              :: ldau_data ! changes
  type (JijData), pointer               :: jij_data  ! changes
  type (TruncationZone), pointer        :: trunc_zone ! never changes

  double complex, parameter :: CZERO = (0.0d0, 0.0d0)
  type (TimerMpi) :: mult_scattering_timer
  type (TimerMpi) :: single_site_timer
  integer :: ie
  integer :: rf
  integer :: ispin
  integer :: prspin
  integer :: nmesh
  integer :: ekm
  logical :: xccpl
  double complex :: JSCAL ! scaling factor for Jij calculation
  integer :: I1
  integer, allocatable :: atom_indices(:)
  integer :: ilocal
  integer :: num_local_atoms
  integer :: lmmaxd
  double complex, allocatable, dimension(:,:,:) :: TMATLL !< all t-matrices
  double complex, allocatable, dimension(:,:,:) :: GmatN_buffer !< GmatN for all local atoms

  lmmaxd = (dims%lmaxd + 1) ** 2

  trunc_zone => getTruncationZone(calc_data)
  gaunts    => getGaunts(calc_data)
  atomdata  => getAtomData(calc_data, 1)
  I1 = atomdata%atom_index
  kkr       => null()
  ldau_data => getLDAUData(calc_data, 1)
  jij_data  => getJijData(calc_data, 1)

  num_local_atoms = getNumLocalAtoms(calc_data)

  ! allocate buffer for t-matrices
  allocate(TMATLL(lmmaxd, lmmaxd, trunc_zone%naez_trc))

  allocate(GmatN_buffer(lmmaxd,lmmaxd,num_local_atoms))
  allocate(atom_indices(num_local_atoms))

  if (params%jij .and. num_local_atoms > 1) then
    if (isMasterRank(my_mpi)) write(*,*) "Jij and num_local_atoms > 1 not supported."
    STOP
  endif

  if (params%ldau .and. num_local_atoms > 1) then
    if (isMasterRank(my_mpi)) write(*,*) "LDA+U and num_local_atoms > 1 not supported."
    STOP
  endif

  xccpl = .false.

  call resetTimer(mult_scattering_timer)
  call stopTimer(mult_scattering_timer)

  call resetTimer(single_site_timer)

  EKM = 0
  prspin = 1

  ! calculate exchange couplings only at last self-consistency step and when Jij=true
  if ((ITER==params%SCFSTEPS).and.params%JIJ) XCCPL = .true.

  if (XCCPL) then

    call CLSJIJ(I1,dims%NAEZ,arrays%RR,arrays%NR,arrays%RBASIS, &
                jij_data%RCUTJIJ,arrays%NSYMAT,arrays%ISYMINDEX, &
                jij_data%IXCP,jij_data%NXCP,jij_data%NXIJ,jij_data%RXIJ, &
                jij_data%RXCCLS,jij_data%ZKRXIJ, &
                arrays%nrd, jij_data%nxijd)

    jij_data%JXCIJINT = CZERO
    jij_data%GMATXIJ = CZERO

  endif

! IE ====================================================================
!     BEGIN do loop over energies (EMPID-parallel)
! IE ====================================================================
  do IE = 1, emesh%ielast
! IE ====================================================================
    if (getMyEnergyId(my_mpi)==ebalance_handler%EPROC(IE)) then
! IE ====================================================================
      call startEBalanceTiming(ebalance_handler, IE)

      WRITELOG(2, *) "Working on energy point ", IE

!------------------------------------------------------------------------------
      !$omp parallel do private(ilocal, kkr, RF)
      do ilocal = 1, num_local_atoms  ! not so smart, redundant calculations
        kkr => getKKR(calc_data, ilocal)
!------------------------------------------------------------------------------
        kkr%noiter = 0

        do RF = 1,arrays%NREF
          call TREF(emesh%EZ(IE),arrays%VREF(RF),arrays%LMAXD,arrays%RMTREF(RF), &
                    kkr%TREFLL(1,1,RF),kkr%DTREFLL(1,1,RF), dims%LLY)
        end do

        !TESTARRAYLOG(3, kkr%TREFLL)
        !TESTARRAYLOG(3, kkr%DTREFLL)


        call GREF(emesh%EZ(IE),params%ALAT,gaunts%IEND,arrays%NCLS,arrays%NAEZ, &
                      gaunts%CLEB,arrays%RCLS,arrays%ATOM,arrays%CLS,gaunts%ICLEB, &
                      gaunts%LOFLM,arrays%NACLS, arrays%REFPOT, &
                      kkr%TREFLL,kkr%DTREFLL,kkr%GREFN,kkr%DGREFN, &
                      kkr%LLY_G0TR(:,IE), &
                      arrays%lmaxd, arrays%naclsd, gaunts%ncleb, kkr%nrefd, kkr%nclsd, &
                      dims%LLY)

        !TESTARRAYLOG(3, kkr%GREFN)
        !TESTARRAYLOG(3, kkr%DGREFN)

!------------------------------------------------------------------------------
      end do  ! ilocal
      !$omp end parallel do
!------------------------------------------------------------------------------

! SPIN ==================================================================
!     BEGIN do loop over spins
! SPIN===================================================================
!------------------------------------------------------------------------------
!     beginning of SMPID-parallel section
!------------------------------------------------------------------------------
      spinloop: do ISPIN = 1,dims%NSPIND
        if(isWorkingSpinRank(my_mpi, ispin)) then

          if (dims%SMPID==1) then
            PRSPIN   = ISPIN
          else
            PRSPIN   = 1
          endif

!------------------------------------------------------------------------------
          !$omp parallel do private(ilocal, kkr, atomdata, ldau_data, jij_data, I1, RF)
          do ilocal = 1, num_local_atoms
            kkr => getKKR(calc_data, ilocal)
            atomdata => getAtomData(calc_data, ilocal)
            ldau_data => getLDAUData(calc_data, ilocal)
            jij_data  => getJijData(calc_data, ilocal)
            I1 = getAtomIndexOfLocal(calc_data, ilocal)
!------------------------------------------------------------------------------

            call CALCTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, &
                            params%NSRA, gaunts, kkr%TMATN, kkr%TR_ALPH, ldau_data)

            jij_data%DTIXIJ(:,:,ISPIN) = kkr%TMATN(:,:,ISPIN)  ! save t-matrix for Jij-calc.

            if(dims%LLY==1) then  ! calculate derivative of t-matrix for Lloyd's formula
              call CALCDTMAT_wrapper(atomdata, emesh, ie, ispin, params%ICST, &
                            params%NSRA, gaunts, kkr%DTDE, kkr%TR_ALPH, ldau_data)
            end if

            RF = arrays%REFPOT(I1)
            call substractReferenceTmatrix(kkr%TMATN(:,:,ISPIN), &
                                           kkr%TREFLL(:,:,RF), kkr%LMMAXD)
            ! do the same for derivative of T-matrix
            call substractReferenceTmatrix(kkr%DTDE(:,:,ISPIN), &
                                           kkr%DTREFLL(:,:,RF), kkr%LMMAXD)

            ! TMATN now contains Delta t = t - t_ref !!!
            ! DTDE now contains Delta dt !!!

            ! renormalize TR_ALPH
            kkr%TR_ALPH(ISPIN) = kkr%TR_ALPH(ISPIN) - kkr%LLY_G0TR(arrays%CLS(I1), IE)

            call rescaleTmatrix(kkr%TMATN(:,:,ISPIN), kkr%lmmaxd, params%alat)

!------------------------------------------------------------------------------
          end do ! ilocal
          !$omp end parallel do
!------------------------------------------------------------------------------

          NMESH = arrays%KMESH(IE)

          if( getMyAtomRank(my_mpi)==0 ) then
            if (params%KTE >= 0) call printEnergyPoint(emesh%EZ(IE), IE, ISPIN, NMESH)
          end if

          call stopTimer(single_site_timer)
          call resumeTimer(mult_scattering_timer)

! <<>> Multiple scattering part

!                 if (atom_indices(1) == 1 .and. IE == params%IELAST .and. ISPIN == 1) then
!                    DEBUG_dump_matrix = .true.
!                 else
!                    DEBUG_dump_matrix = .false.
!                 endif

          ! gather t-matrices from own truncation zone
          call gatherTmatrices_com(calc_data, TMATLL, ispin, &
                                   getMySEcommunicator(my_mpi))

          TESTARRAYLOG(3, TMATLL)

!------------------------------------------------------------------------------

          do ilocal = 1, num_local_atoms
            atom_indices(ilocal) = getAtomIndexOfLocal(calc_data, ilocal)
            atom_indices(ilocal) = trunc_zone%index_map(atom_indices(ilocal))
            CHECKASSERT(atom_indices(ilocal) > 0)
          end do

          ! problem: reference Green's functions
          ! here: known by all atoms - therefore pick any atom (nr.1)
          kkr => getKKR(calc_data, 1)
          jij_data => getJijData(calc_data, 1)

          call KLOOPZ1_new(GmatN_buffer, params%ALAT, &
          trunc_zone%NAEZ_trc,arrays%NOFKS(NMESH),arrays%VOLBZ(NMESH), &
          arrays%BZKP(:,:,NMESH),arrays%VOLCUB(:,NMESH), trunc_zone%CLS_trc, &
          arrays%NACLS,arrays%RR,trunc_zone%EZOA_trc,trunc_zone%ATOM_trc, &
          kkr%GREFN, &
          arrays%NSYMAT,arrays%DSYMLL, &
          TMATLL, &
          trunc_zone%NUMN0_trc,trunc_zone%INDN0_trc,atom_indices, &
          params%QMRBOUND, &
          arrays%lmmaxd, arrays%naclsd,  &
          arrays%nrd)

!------------------------------------------------------------------------------

          ! copy results from buffer: G_LL'^NN (E, spin) =
          !                           GmatN_buffer_LL'^N(ilocal) N(ilocal)
          do ilocal = 1, num_local_atoms
            kkr => getKKR(calc_data, ilocal)
            kkr%GMATN(:,:,ie,ispin) = GmatN_buffer(:,:,ilocal)
          end do

          call stopTimer(mult_scattering_timer)
          call resumeTimer(single_site_timer)

        endif
      end do spinloop                          ! ISPIN = 1,NSPIN
!------------------------------------------------------------------------------
!        End of SMPID-parallel section
!------------------------------------------------------------------------------
! SPIN ==================================================================
!     END do loop over spins
! SPIN===================================================================

! =====================================================================
! Calculate Jij for the in CLSJIJ predefined atom pairs i,j
! xccpl

      if (XCCPL) then
        call jijSpinCommunication_com(my_mpi, jij_data%GMATXIJ, jij_data%DTIXIJ)

        ! calculate DTIXIJ = T_down - T_up
        call calcDeltaTupTdown(jij_data%DTIXIJ)

        JSCAL = emesh%WEZ(IE)/DBLE(jij_data%NSPIND)

        call jijLocalEnergyIntegration(my_mpi, JSCAL, jij_data%GMATXIJ, &
                                        jij_data%DTIXIJ(:,:,1), jij_data%RXIJ,&
                                        jij_data%NXIJ, jij_data%IXCP, &
                                        jij_data%RXCCLS, jij_data%JXCIJINT)
      end if

! xccpl
! End of Jij calculation
! =====================================================================

      call stopEBalanceTiming(ebalance_handler, ie)

! IE ====================================================================
    endif
! IE ====================================================================

! for initial guess calculate sparse indices combining IE.KPT
    EKM = EKM + arrays%NOFKS(arrays%KMESH(IE))

  end do                   ! IE = 1,IELAST

! IE ====================================================================
!     END do loop over energies (EMPID-parallel)
! IE ====================================================================

  call stopTimer(single_site_timer)

!=======================================================================
!communicate information of 1..EMPID and 1..SMPID processors to MASTERGROUP
  do ilocal = 1, num_local_atoms
    kkr => getKKR(calc_data, ilocal)
    call collectMSResults_com(my_mpi, kkr%GMATN, kkr%LLY_GRDT, &
                              ebalance_handler%EPROC)
  end do
!=======================================================================

! TIME
  call OUTTIME(isMasterRank(my_mpi),'Single Site took.....', &
               getElapsedTime(single_site_timer),ITER)
  call OUTTIME(isMasterRank(my_mpi),'Mult. Scat. took.....', &
               getElapsedTime(mult_scattering_timer),ITER)

!=======================================================================
!     output of Jij's
!=======================================================================
  if (XCCPL) then

    call jijReduceIntResults_com(my_mpi, jij_data%JXCIJINT)

    if (isInMasterGroup(my_mpi)) then
      call writeJiJs(I1,jij_data%RXIJ,jij_data%NXIJ,jij_data%IXCP, &
                     jij_data%RXCCLS,jij_data%JXCIJINT, jij_data%nxijd)
    end if
  endif

!=======================================================================
!     on the basis of new timings determine now new distribution of
!     work to 1 .. EMPID processors - all processes SYNCED
!=======================================================================
  call updateEBalance_com(ebalance_handler, my_mpi)

!=======================================================================
!     in case of IGUESS and EMPID > 1 initial guess arrays might
!     have to be adjusted to new distributions
!=======================================================================
  if ((dims%IGUESSD==1) .and. (dims%EMPID>1)) then

    do ISPIN = 1,dims%NSPIND
      if(isWorkingSpinRank(my_mpi, ispin)) then

        if (dims%SMPID==1) then
          PRSPIN   = ISPIN
        else
          PRSPIN   = 1
        endif

        WRITELOG(3, *) "EPROC:     ", ebalance_handler%EPROC
        WRITELOG(3, *) "EPROC_old: ", ebalance_handler%EPROC_old

      do ilocal = 1, num_local_atoms
        kkr => getKKR(calc_data, ilocal)
        ! Note: MPI overhead does not increase because number of procs goes down
        call redistributeInitialGuess_com(my_mpi, kkr%PRSC(:,:,PRSPIN), &
             ebalance_handler%EPROC, ebalance_handler%EPROC_old, &
             arrays%KMESH, arrays%NofKs)
      end do

      endif
    enddo

  endif  ! IGUESS == 1 .and. EMPID > 1
!=======================================================================

  deallocate(atom_indices)
  deallocate(GmatN_buffer)
  deallocate(TMATLL)

end subroutine

! =============================================================================
! Helper routines
! =============================================================================

!----------------------------------------------------------------------------
!> Print info about Energy-Point currently treated.
!>
subroutine printEnergyPoint(EZ_point, IE, ISPIN, NMESH)
  implicit none
  double complex :: EZ_point
  integer :: IE
  integer :: ISPIN
  integer :: NMESH
  write (6,'(A,I3,A,2(1X,F10.6),A,I3,A,I3)')  &
  ' ** IE = ',IE,' ENERGY =',EZ_point, &
  ' KMESH = ', NMESH,' ISPIN = ',ISPIN
end subroutine

!----------------------------------------------------------------------------
!> Calculate \Delta T_up - T_down for exchange couplings calculation.
!> The result is stored in DTIXIJ(:,:,1)
subroutine calcDeltaTupTdown(DTIXIJ)
  implicit none
  double complex, intent(inout) :: DTIXIJ(:,:,:)
  integer :: LMMAXD

  integer :: LM1
  integer :: LM2

  lmmaxd = size(DTIXIJ,1)

  do LM2 = 1,LMMAXD
    do LM1 = 1,LMMAXD
      DTIXIJ(LM1,LM2,1) = DTIXIJ(LM1,LM2,2) - DTIXIJ(LM1,LM2,1)
    enddo
  enddo

end subroutine

!----------------------------------------------------------------------------
!> Substract diagonal reference T matrix of certain spin channel
!> from real system's T matrix.
subroutine substractReferenceTmatrix(TMATN, TREFLL, LMMAXD)
  implicit none
  integer :: LM1
  integer :: LMMAXD
  double complex :: TMATN(:,:)
  double complex :: TREFLL(:,:)

  ! Note: TREFLL is diagonal! - spherical reference potential
  do LM1 = 1,LMMAXD
    TMATN(LM1,LM1) =  TMATN(LM1,LM1) - TREFLL(LM1,LM1)
  end do

end subroutine

!------------------------------------------------------------------------------
!> Rescale and symmetrise T-matrix.
subroutine rescaleTmatrix(tsst_local, lmmaxd, alat)
  implicit none

  double complex, intent(inout), dimension(lmmaxd, lmmaxd) :: tsst_local
  integer, intent(in) :: lmmaxd
  double precision, intent(in) :: alat

  integer :: lm1, lm2
  double precision :: RFCTOR

  !     RFCTOR=A/(2*PI) conversion factor to p.u.
    RFCTOR = ALAT/(8.D0*ATAN(1.0D0))           ! = ALAT/(2*PI)

! --> convert inverted delta_t-matrices to p.u.
!     Also a symmetrisation of the matrix is performed

    do LM2 = 1,LMMAXD
        do LM1 = 1,LM2
            TSST_LOCAL(LM1,LM2) = 0.5D0/RFCTOR * &
            ( TSST_LOCAL(LM1,LM2) + TSST_LOCAL(LM2,LM1) )
            TSST_LOCAL(LM2,LM1) = TSST_LOCAL(LM1,LM2)
        end do
    end do
end subroutine

!------------------------------------------------------------------------------
!> Gather all t-matrices for 'ispin'-channel (from truncation zone only).
!>
!> Uses MPI-RMA
subroutine gatherTmatrices_com(calc_data, TMATLL, ispin, communicator)
  use CalculationData_mod
  use KKRresults_mod
  use TruncationZone_mod
  use one_sided_commZ_mod
  implicit none
  include 'mpif.h'

  type (CalculationData), intent(in) :: calc_data
  double complex, dimension(:,:,:), intent(inout) :: TMATLL
  integer, intent(in) :: ispin
  integer, intent(in) :: communicator

  type (KKRresults), pointer :: kkr
  type (TruncationZone), pointer :: trunc_zone

  type (ChunkIndex), dimension(:), allocatable :: chunk_inds

  integer :: ii
  integer :: ilocal
  integer :: num_local_atoms
  integer :: ierr
  integer :: lmmaxd
  integer :: naez_trc ! number of atoms in trunc. zone
  integer :: naez
  integer :: nranks
  integer :: atom_requested
  integer :: win
  integer :: chunk_size
  double complex, allocatable, dimension(:,:,:) :: TSST_LOCAL

  num_local_atoms = getNumLocalAtoms(calc_data)
  trunc_zone => getTruncationZone(calc_data)
  lmmaxd = size(TMATLL, 1)

  call MPI_Comm_size(communicator, nranks, ierr)

  allocate(TSST_LOCAL(lmmaxd, lmmaxd, num_local_atoms))

  chunk_size = size(TSST_LOCAL, 1) * size(TSST_LOCAL, 2)

  do ilocal = 1, num_local_atoms
    kkr => getKKR(calc_data, ilocal)
    TSST_LOCAL(:,:,ilocal) = kkr%TMATN(:,:,ispin)
  end do

  naez_trc = trunc_zone%naez_trc

  naez = num_local_atoms * nranks
  CHECKASSERT( naez == size(trunc_zone%index_map) )

  allocate(chunk_inds(naez_trc))

  do ii = 1, naez_trc
    atom_requested = trunc_zone%trunc2atom_index(ii)
    chunk_inds(ii)%owner = getOwner(atom_requested, naez, nranks)
    chunk_inds(ii)%local_ind = getLocalInd(atom_requested, naez, nranks)
  end do

  call exposeBufferZ(win, TSST_LOCAL, size(TSST_LOCAL), chunk_size, communicator)
  call copyChunksZ(TMATLL, win, chunk_inds, chunk_size)
  call hideBufferZ(win)

  deallocate(chunk_inds)
  deallocate(TSST_LOCAL)

end subroutine

end module
