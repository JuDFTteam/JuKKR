module KKRnano_Comm_mod
  implicit none

  CONTAINS

  !> Collect the results from the multiple scattering part at
  !> the corresponding atom process of the master group.
  !> Master Group: (Spin, Energy)-id = 1
  subroutine collectMultScatteringResults(my_mpi, GMATN_ALL, LLY_GRDT_ALL, EPROC)
    use KKRnanoParallel_mod
    use comm_patternsZ_mod
    implicit none

    type (KKRnanoParallel) :: my_mpi

    double complex, dimension(:,:,:,:), intent(inout) ::  GMATN_ALL
    double complex, dimension(:,:), intent(inout) ::  LLY_GRDT_ALL
    integer, dimension(:), intent(in) :: EPROC

    integer :: iemxd
    integer :: nspind
    integer :: lmmaxd

    integer :: ie
    integer :: ispin

    integer :: my_atom_rank

    integer, dimension(:), allocatable :: owning_ranks
    integer :: receiver

    lmmaxd = size(GMATN_ALL,1)
    iemxd = size(GMATN_ALL,3)
    nspind = size(GMATN_ALL,4)

    ! TODO: check allocate
    allocate(owning_ranks(iemxd * nspind))

    my_atom_rank = getMyAtomRank(my_mpi)

    do ispin = 1, nspind
      do ie = 1, iemxd

        owning_ranks( (ispin-1)*iemxd + ie ) = mapToWorldRank(my_mpi, my_atom_rank, ispin, EPROC(ie))

      end do
    end do

    receiver = mapToWorldRankSE(my_mpi, my_atom_rank, 1)

    call comm_gatherZ(my_world_rank, GMATN_ALL, lmmaxd*lmmaxd, owning_ranks, receiver)
    call comm_gatherZ(my_world_rank, LLY_GRDT_ALL, 1, owning_ranks, receiver)


    ! TODO: check dimensions

    deallocate(owning_ranks)

  end subroutine

end module KKRnano_Comm_mod
