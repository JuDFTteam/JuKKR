module mod_types

      Use mod_datatypes, Only: dp
implicit none


   type :: type_opt_test

      integer :: Nparams = 2
      CHARACTER*8 OPTC(32),TESTC(32)

   end type type_opt_test
   


   type (type_opt_test), save :: t_optt


contains

#ifdef CPP_MPI
   subroutine bcast_t_opt_test(t_optt)

    use mpi
    use mod_mympi,   only: myrank, nranks, master
    implicit none

    type(type_opt_test), intent(inout) :: t_optt

    integer :: blocklen1(t_inc%Nparams), etype1(t_inc%Nparams), myMPItype1 ! for parameter from t_inc
    integer :: ierr
    integer(kind=MPI_ADDRESS_KIND) :: disp1(t_inc%Nparams),  base

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !broadcast parameters from t_optt
    call MPI_Get_address(t_inc%Nparams,       disp1(1), ierr)
    call MPI_Get_address(t_inc%LMMAXD,        disp1(2), ierr)
    base  = disp1(1)
    disp1 = disp1 - base

    blocklen1(1:15)=1

    etype1(1:14) = MPI_INTEGER
    etype1(15) = MPI_LOGICAL

    call MPI_Type_create_struct(t_inc%Nparams, blocklen1, disp1, etype1, myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'Problem in create_mpimask_t_inc'

    call MPI_Type_commit(myMPItype1, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error commiting create_mpimask_t_inc'

    call MPI_Bcast(t_inc%Nparams, 1, myMPItype1, master, MPI_COMM_WORLD, ierr)
    if(ierr/=MPI_SUCCESS) stop 'error brodcasting t_inc'

    call MPI_Type_free(myMPItype1, ierr)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   end subroutine bcast_t_opt_test
#endif


end module
