!------------------------------------------------------------------------------
!> Module for two-sided MPI communication
!>
!> @author Paul Baumeister
!>
!> Only for equally sized chunks of data which have basic numeric datatypes.
!>
!> Data layout:
!>
!> \verbatim
!>
!> A distributed array is created (very similar to Co-array Fortran),
!> which consists of equally sized chunks.
!>
!> Each rank then copies a slice of this array into its local memory.
!>
!> Note: this corresponds to the "partitioned global address space" (PGAS)
!>       programming model - however in our case only READ access is allowed.
!>
!> Each chunk has an "owner"-index (MPI-rank starting at 0)
!> and a local index (starting at 1 !)
!>
!> chunk index = (rank, local index)
!>
!>  ___________|_____ _____|
!> |     |     |     |     |
!> |(0,1)|(0,2)|(1,1)|(1,2)| usw...
!> |_____|_____|_____|_____|
!>    rank 0   |  rank 1   |
!>
!> |atm 1|atm 2|atm 3|atm 4| ...
!>
!> To access the chunks by a continuous index (called 'atom index')
!> one can use the routines getOwner and getLocalInd
!> to convert to a chunk index
!>
!> \endverbatim

#define NUMBERZ double complex
#define NUMBERMPIZ MPI_DOUBLE_COMPLEX
#define NUMBERC complex
#define NUMBERMPIC MPI_COMPLEX
#define NUMBERD double precision
#define NUMBERMPID MPI_DOUBLE_PRECISION
#define NUMBERI integer
#define NUMBERMPII MPI_INTEGER

#ifndef NUMBER_TYPE
! defaults
#define NUMBER_TYPE integer
#define NUMBERMPI_TYPE MPI_INTEGER
#endif

module two_sided_comm_TYPE_mod
  use ExchangeTable_mod, only: ExchangeTable

#define assert(condition) if(.not. (condition)) stop __LINE__

  implicit none
  private
  
  public :: reference_sys_com_TYPE

  contains

  subroutine reference_sys_com_TYPE(self, ncount, Gout, Ginp)
    type(ExchangeTable), intent(in) :: self
    integer, intent(in) :: ncount
    NUMBER_TYPE, intent(out) :: Gout(ncount,*) !> dim(ncount,num_trunc_atoms)
    NUMBER_TYPE, intent(in)  :: Ginp(ncount,*) !< dim(ncount,num_local_atoms)
    
    ! In order to implement point-to-point communication also with max_local_atoms > 1,
    ! we need to set up a sparse matrix that describes 
    !   which references system (col) we need to receive from which rank (pair)
    ! then, we need to communicate these matrices
    !   and create the tables about which local reference system we need to send to which rank (pair)

    integer, allocatable :: rreq(:), sreq(:), rstats(:,:), sstats(:,:)
    integer :: ipair, rank, tag, ierr, ist, inz, iinp, iout, myrank
    include 'mpif.h' ! only: MPI_STATUS_SIZE, MPI_INTEGER, MPI_REQUEST_NULL

    assert( self%comm /= 0 )
    
    call MPI_Comm_rank(self%comm, myrank, ierr)

    allocate(sreq(self%send_n), sstats(MPI_STATUS_SIZE,self%send_n), &
             rreq(self%recv_n), rstats(MPI_STATUS_SIZE,self%recv_n), stat=ist) ! ToDo: catch status
    sreq(:) = MPI_REQUEST_NULL
    rreq(:) = MPI_REQUEST_NULL

    ! now communication works like this
    do ipair = 1, self%npairs
      rank = self%pair2rank(ipair)

      if (rank /= myrank) then
      
        do inz = self%send_start(ipair), self%send_start(ipair + 1) - 1
          tag = inz - self%send_start(ipair)
          iinp = self%send_index(inz)
          call MPI_Isend(Ginp(:,iinp), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, sreq(inz), ierr)
        enddo ! inz

        do inz = self%recv_start(ipair), self%recv_start(ipair + 1) - 1
          tag = inz - self%recv_start(ipair)
          iout = self%recv_index(inz)
          call MPI_Irecv(Gout(:,iout), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, rreq(inz), ierr)
        enddo ! inz
        
      else  ! rank /= myrank

        assert( self%send_start(ipair) - self%send_start(ipair + 1) == self%recv_start(ipair) - self%recv_start(ipair + 1) )
        do inz = self%recv_start(ipair), self%recv_start(ipair + 1) - 1
          iout = self%recv_index(inz)
          ist = inz - self%recv_start(ipair) + self%send_start(ipair)
          iinp = self%send_index(ist)
          Gout(:,iout) = Ginp(:,iinp) ! local copy
        enddo ! inz
        
      endif ! rank /= myrank

    enddo ! ipair
    call MPI_Waitall(self%send_n, sreq, sstats, ierr) ! wait until all sends have finished
    call MPI_Waitall(self%recv_n, rreq, rstats, ierr) ! wait until all receives have finished

    deallocate(rreq, rstats, sreq, sstats, stat=ist) ! ignore status
  endsubroutine ! reference_sys_com

  subroutine reference_sys_com_TYPE_var_sizes(self, Gout, Ginp)
    ! if we have additional information about the sizes stored in self, we can extend this method to variable message sizes
    type(ExchangeTable), intent(in) :: self
    NUMBER_TYPE, intent(out) :: Gout(:,:,:,:,:) !> dim(lmmaxd,lmmaxd,0:LLy,naclsd,num_trunc_atoms)
    NUMBER_TYPE, intent(in)  :: Ginp(:,:,:,:,:) !< dim(lmmaxd,lmmaxd,0:LLy,naclsd,num_local_atoms)
    
    ! In order to implement point-to-point communication also with max_local_atoms > 1,
    ! we need to set up a sparse matrix that describes 
    !   which references system (col) we need to receive from which rank (pair)
    ! then, we need to communicate these matrices
    !   and create the tables about which local reference system we need to send to which rank (pair)

    integer, allocatable :: rreq(:), sreq(:), rstats(:,:), sstats(:,:)
    integer :: ncount, ipair, rank, tag, ierr, ist, inz, iinp, iout
    include 'mpif.h' ! only: MPI_STATUS_SIZE, MPI_INTEGER, MPI_REQUEST_NULL

    assert( self%comm /= 0 )
    
    ncount = size(Ginp(:,:,:,:,1))
    do ist = 1, 4
      assert( size(Ginp, ist) == size(Gout, ist) )
    enddo ! ist

    allocate(sreq(self%send_n), sstats(MPI_STATUS_SIZE,self%send_n), &
             rreq(self%recv_n), rstats(MPI_STATUS_SIZE,self%recv_n), stat=ist) ! ToDo: catch status
    sreq(:) = MPI_REQUEST_NULL
    rreq(:) = MPI_REQUEST_NULL

    ! now communication works like this
    do ipair = 1, self%npairs
      rank = self%pair2rank(ipair)
      
      do inz = self%send_start(ipair), self%send_start(ipair + 1) - 1
        tag = inz - self%send_start(ipair)
        iinp = self%send_index(inz)
        ! ncount = size(Ginp, 1)*size(Ginp, 2)*size(Ginp, 3)*self%send_size(inz) ! ToDo
        call MPI_Isend(Ginp(:,:,:,:,iinp), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, sreq(inz), ierr)
      enddo ! inz

      do inz = self%recv_start(ipair), self%recv_start(ipair + 1) - 1
        tag = inz - self%recv_start(ipair)
        iout = self%recv_index(inz)
        ! ncount = size(Gout, 1)*size(Gout, 2)*size(Gout, 3)*self%recv_size(inz) ! ToDo
        call MPI_Irecv(Gout(:,:,:,:,iout), ncount, NUMBERMPI_TYPE, rank, tag, self%comm, rreq(inz), ierr)
      enddo ! inz

    enddo ! ipair
    call MPI_Waitall(self%send_n, sreq, sstats, ierr) ! wait until all sends have finished
    call MPI_Waitall(self%recv_n, rreq, rstats, ierr) ! wait until all receives have finished

    deallocate(rreq, rstats, sreq, sstats, stat=ist) ! ignore status
  endsubroutine ! reference_sys_com

endmodule ! two_sided_comm_TYPE_mod
