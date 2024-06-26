!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the handling of the wavefunctions
!> Author:
!> Wrapper module for the handling of the wavefunctions
!------------------------------------------------------------------------------------
module mod_save_wavefun

  use :: mod_datatypes

  implicit none

  !-------------------------------------------------------------------------------
  !> Summary: Definition of the variables used to store and handle the wavefunctions
  !> Author:
  !> Category: single-site, KKRhost
  !> Deprecated: False
  !> Definition of the variables used to store and handle the wavefunctions
  !-------------------------------------------------------------------------------
  type :: type_wavefunctions

    integer :: nwfsavemax          !! number of wavefunctions (na*ne) that this rank has to store
    integer :: maxmem_units        !! units of maxmem_number (1024=kbyte, 1024**2=Mbyte, 1024**3=Gbyte)
    integer :: maxmem_number       !! maximal amount of memory available to store wavefunctions (maxmem_number*maxmem_units byte)
    integer :: nth                 ! number of OpenMP tasks used without MPI (not the number of tasks in hybrid mode) in OpenMP mode only

    ! allocatable arrays
    integer, dimension(:,:), allocatable :: isave_wavefun !! 0 if (iat_myrank, ie_myrank) pair is not saved, 1...Nwfsavemax otherwise, for first, second, ... Nwfsavemax-th saved wavefunction on this rank; (Nat_myrank, Ne_myrank)
    logical :: save_rll, save_sll, save_rllleft, save_sllleft !! logicals that say which wavefunctions are stored (used to reduce amount of memory that is used
    complex (kind=dp), dimension(:, :, :, :, :),  allocatable :: rll !! regular right wavefunction; (Nwfsavemax, NSRA*lmmaxd, lmmaxd, IRMDNEW, 0:nth-1)
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: rllleft !! regular left wavefunction; (Nwfsavemax, NSRA*lmmaxd, lmmaxd, IRMDNEW, 0:nth-1)
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: sll !! iregular right wavefunction; (Nwfsavemax, NSRA*lmmaxd, lmmaxd, IRMDNEW, 0:nth-1)
    complex (kind=dp), dimension(:, :, :, :, :), allocatable :: sllleft !! iregular left wavefunction; (Nwfsavemax, NSRA*lmmaxd, lmmaxd, IRMDNEW, 0:nth-1)

  end type type_wavefunctions

  type (type_wavefunctions), save :: t_wavefunctions

  ! switch to enable some additional output
  integer, parameter :: iprint_debug = 0

contains

  !-------------------------------------------------------------------------------
  !> Summary: Find how many wavefunctions can be stored and set `isave_wavefun` array to determine which wavefunctions are stored
  !> Author:
  !> Category: input-output, profiling, initialization, KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> Find how many wavefunctions can be stored and set `isave_wavefun` array to
  !> determine which wavefunctions are stored
  !-------------------------------------------------------------------------------
  subroutine find_isave_wavefun(t_wavefunctions)

#ifdef CPP_OMP
    use :: omp_lib
#endif
    use :: mod_mympi, only: myrank, nranks, master
    use :: mod_types, only: t_inc
#ifdef CPP_MPI
    use :: mod_types, only: t_mpi_c_grid
#endif
    use global_variables, only: natypd, lmmaxd

    implicit none

    type (type_wavefunctions), intent (inout) :: t_wavefunctions

    integer :: nth, nat, ne, ierr, i, ie, iat, maxpos(1), nwfsave_count, nat_myrank, ne_myrank, ioff_at, ioff_ie, nsave
    real (kind=dp) :: delta_mem

    integer, allocatable :: kmesh_priority(:), my_kmesh(:)
    logical, allocatable :: mask(:)

    logical :: firstrun

    character (len=200) :: message

    ! >>>>>>>  set some parameters  >>>>>>>!
    nth = 1
#ifdef CPP_OMP
    nth = omp_get_num_threads()
#endif
    t_wavefunctions%nth = nth

    nat = natypd
    ne = t_inc%ielast
#ifdef CPP_MPI
    ! local parameter for myrank
    nat_myrank = t_mpi_c_grid%ntot1
    ne_myrank = t_mpi_c_grid%ntot2
    ioff_at = t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie)
    ioff_ie = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
#else
    nat_myrank = nat
    ne_myrank = ne
    ioff_at = 0
    ioff_ie = 0
#endif

    ! check if this is the first run or not
    firstrun = .true.
    if (allocated(t_wavefunctions%isave_wavefun)) firstrun = .false.

    if (firstrun) then

      allocate (t_wavefunctions%isave_wavefun(nat,ne), stat=ierr)
      if (ierr/=0) stop '[find_isave_wavefun] Error allocating isave_wavefun'
      ! <<<<<<<  set some parameters  <<<<<<<!
      ! >>>>>>>  find numer of wavefunctions that can be stored and allocate store arrays  >>>>>>>!
      ! memory demand for one atom and one energy point in Mbyte
      delta_mem = real(t_inc%nsra*lmmaxd*lmmaxd*t_inc%irmdnew*nth*16, kind=dp)/(1024.0_dp**2)

      ! number of wavefunction (rll, sll, rllleft, sllleft) that are needed to be stored
      nsave = 0
      if (t_wavefunctions%save_rll) nsave = nsave + 1
      if (t_wavefunctions%save_sll) nsave = nsave + 1
      if (t_wavefunctions%save_rllleft) nsave = nsave + 1
      if (t_wavefunctions%save_sllleft) nsave = nsave + 1

      delta_mem = delta_mem*nsave
      ! avoid division by zero
      if (delta_mem<1.0_dp) delta_mem = 1.0_dp

      ! find numer of wavefunctions that can be stored
      t_wavefunctions%nwfsavemax = t_wavefunctions%maxmem_number*1024**(t_wavefunctions%maxmem_units-2)/int(delta_mem)
      if (nsave==0) then
        write (1337, '(A)') '[find_isave_wavefun] Warning: Nsave is zero!!!'
        write (1337, '(A)') 'reset Nwfsavemax to 0'
        t_wavefunctions%nwfsavemax = 0
      end if

      ! avoid unnessesary large allocations of arrays
      if (t_wavefunctions%nwfsavemax>nat_myrank*ne_myrank) then
        t_wavefunctions%nwfsavemax = nat_myrank*ne_myrank
        if (t_inc%i_write>0) write (1337, '(A,I5,A,I9)') '  rank', myrank, ' reset Nwfsavemax to maximal needed number for this thread:', t_wavefunctions%nwfsavemax
      end if

      ! write out info about allocation for wf storage
      if (myrank==master) then

        write (message, '(A)') '   ==> find_isave_wavefun '
        write (1337, '(A)') trim(message)
        if (t_wavefunctions%nwfsavemax>0) write (*, '(A)') trim(message)
        write (message, '(A,F15.2,A)') '   (maxmem given per rank for storage:', real((1024**(t_wavefunctions%maxmem_units-2))*t_wavefunctions%maxmem_number, kind=dp), 'MB'
        write (1337, '(A)') trim(message)
        if (t_wavefunctions%nwfsavemax>0) write (*, '(A)') trim(message)
        write (message, '(A,I11)') '    number of wavefunctions that fit in:', t_wavefunctions%nwfsavemax
        write (1337, '(A)') trim(message)
        if (t_wavefunctions%nwfsavemax>0) write (*, '(A)') trim(message)
        write (message, '(A,F15.2,A)') '    total memory needed for storage:', nranks*t_wavefunctions%nwfsavemax*delta_mem, 'MB'
        write (1337, '(A)') trim(message)
        if (t_wavefunctions%nwfsavemax>0) write (*, '(A)') trim(message)
        write (message, '(A,I11,A,I11,A)') '    Ne=', ne, ', Nat=', nat, ''
        write (1337, '(A)') trim(message)
        if (t_wavefunctions%nwfsavemax>0) write (*, '(A)') trim(message)

        write (message, '(A,F15.2,A)') '    memory demand per atom and energy point for rll, rllleft, sll and sllleft respectively:', delta_mem, 'MB)  <=='
        write (1337, '(A)') trim(message)
        ! write(*, '(A)') trim(message)

        write (message, '(A,I3,4(A,L1))') '    Number of saved wavefunctions per atom and energy:', nsave, '; save rll:', t_wavefunctions%save_rll, '; save sll:', &
          t_wavefunctions%save_sll, '; save rllleft:', t_wavefunctions%save_rllleft, '; save sllleft:', t_wavefunctions%save_sllleft
        write (1337, '(A)') trim(message)
        ! write(*, '(A)') trim(message)
      end if
      ! allocate store arrays in t_wavefunctions
      if (t_wavefunctions%nwfsavemax>0) then
        if (t_wavefunctions%save_rll) then
          allocate (t_wavefunctions%rll(t_wavefunctions%nwfsavemax,t_inc%nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:nth-1), stat=ierr)
          if (ierr/=0) stop '[find_isave_wavefun] Error allocating rll'
        end if
        if (t_wavefunctions%save_sll) then
          allocate (t_wavefunctions%rllleft(t_wavefunctions%nwfsavemax,t_inc%nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:nth-1), stat=ierr)
          if (ierr/=0) stop '[find_isave_wavefun] Error allocating rllleft'
        end if
        if (t_wavefunctions%save_rllleft) then
          allocate (t_wavefunctions%sll(t_wavefunctions%nwfsavemax,t_inc%nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:nth-1), stat=ierr)
          if (ierr/=0) stop '[find_isave_wavefun] Error allocating sll'
        end if
        if (t_wavefunctions%save_sllleft) then
          allocate (t_wavefunctions%sllleft(t_wavefunctions%nwfsavemax,t_inc%nsra*lmmaxd,lmmaxd,t_inc%irmdnew,0:nth-1), stat=ierr)
          if (ierr/=0) stop '[find_isave_wavefun] Error allocating sllleft'
        end if
      end if
      ! <<<<<<<  find numer of wavefunctions that can be stored and allocate store arrays  <<<<<<<!
      if (t_wavefunctions%nwfsavemax>0) then
        ! >>>>>>>  find priority of saving wavefunctions according to kmesh  >>>>>>>!
        ! kmesh( 1 : ne )
        allocate (my_kmesh(ne_myrank), stat=ierr)
        if (ierr/=0) stop '[find_isave_wavefun] Error allocating my_kmesh'
        my_kmesh(:) = t_inc%kmesh_ie(1+ioff_ie:ne_myrank+ioff_ie)
        maxpos(1) = 0
        do i = 1, ne_myrank
          if (my_kmesh(i)>maxpos(1)) maxpos(1) = my_kmesh(i)
        end do
        my_kmesh(:) = maxpos(1) - my_kmesh(:)

        allocate (kmesh_priority(ne_myrank), stat=ierr)
        if (ierr/=0) stop '[find_isave_wavefun] Error allocating kmesh_priority'
        allocate (mask(ne_myrank), stat=ierr)
        if (ierr/=0) stop '[find_isave_wavefun] Error allocating mask'

        mask(:) = .true.
        do i = 1, ne_myrank
          ! find position of maximal value in my_kmesh and fill array kmesh_priority
          maxpos = maxloc(my_kmesh, mask=mask)
          kmesh_priority(maxpos(1)) = i
          ! update mask to exclude previously found position
          mask(maxpos(1)) = .false.
        end do                     ! i=1,ne_myrank

        if (any(kmesh_priority==-1)) then
          write (1337, *) '[find_isave_wavefun] kmesh_priority not found correctly; kmesh_priority:', kmesh_priority, '; mask:', mask
          stop
        end if
        ! <<<<<<<  find priority of saving wavefunctions according to kmesh  <<<<<<<!
        ! >>>>>>>  set isave_wavefun array  >>>>>>>!
        t_wavefunctions%isave_wavefun(:, :) = 0
        nwfsave_count = t_wavefunctions%nwfsavemax
        do iat = 1, nat_myrank
          t_wavefunctions%isave_wavefun(iat+ioff_at, 1+ioff_ie:ioff_ie+ne_myrank) = kmesh_priority(:) + ne_myrank*(iat-1)
          do ie = 1, ne_myrank
            if ((kmesh_priority(ie)>nwfsave_count)) then
              t_wavefunctions%isave_wavefun(iat+ioff_at, ie+ioff_ie) = 0
            end if
          end do
          nwfsave_count = nwfsave_count - ne_myrank
        end do
        ! <<<<<<<  set isave_wavefun array  <<<<<<<!
        ! deallocate work arrays
        deallocate (mask, kmesh_priority, my_kmesh, stat=ierr)
        if (ierr/=0) stop '[find_isave_wavefun] Error dellocating arrays'
      end if                       ! (t_wavefunctions%Nwfsavemax>0)

    end if                         ! firstrun

  end subroutine find_isave_wavefun

  !-------------------------------------------------------------------------------
  !> Summary: Saves wavefunction of atom iat and energypoint ie if enough memory is given
  !> Author:
  !> Category: single-site, KKRhost
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !> Saves wavefunction of atom iat and energypoint ie if enough memory is given
  !-------------------------------------------------------------------------------
  subroutine save_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,iat,ie,nsra,     &
    lmmaxd,irmdnew,ith)

    implicit none
    type (type_wavefunctions), intent (inout) :: t_wavefunctions
    integer, intent (in) :: iat, ie, nsra, lmmaxd, irmdnew, ith
    complex (kind=dp), intent (in) :: rll(:, :, :, 0:), rllleft(:, :, :, 0:)
    complex (kind=dp), intent (in) :: sll(:, :, :, 0:), sllleft(:, :, :, 0:)

    integer :: isave

    isave = t_wavefunctions%isave_wavefun(iat, ie)

    if (isave>0) then

      if (t_wavefunctions%save_rll) then

        if (iprint_debug>0) write (*, *) 'save rll', ie, iat
        t_wavefunctions%rll(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = rll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
      end if

      if (t_wavefunctions%save_rllleft) then
        if (iprint_debug>0) write (*, *) 'write out rllleft', ie, iat
        t_wavefunctions%rllleft(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = rllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
      end if

      if (t_wavefunctions%save_sll) then
        if (iprint_debug>0) write (*, *) 'write out sll', ie, iat
        t_wavefunctions%sll(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = sll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
      end if

      if (t_wavefunctions%save_sllleft) then
        if (iprint_debug>0) write (*, *) 'write out sllleft', ie, iat
        t_wavefunctions%sllleft(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = sllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
      end if

    end if

  end subroutine save_wavefunc

  !-------------------------------------------------------------------------------
  !> Summary: Reads wavefunction of atom iat and energypoint ie if it was stored
  !> Author:
  !> Category: single-site, KKRhost
  !> Deprecated: False
  !> Reads wavefunction of atom iat and energypoint ie if it was stored
  !-------------------------------------------------------------------------------
  subroutine read_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,iat,ie,nsra,     &
    lmmaxd,irmdnew,ith,nth,read_in_rll,read_in_sll,read_in_rllleft,read_in_sllleft)

    implicit none
    type (type_wavefunctions), intent (inout) :: t_wavefunctions
    integer, intent (in) :: iat, ie, nsra, lmmaxd, irmdnew, ith, nth
    complex (kind=dp), intent (out) :: rll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, 0:(nth-1)), rllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, 0:(nth-1)), &
      sll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, 0:(nth-1)), sllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, 0:(nth-1))
    logical, intent (out) :: read_in_rll, read_in_sll, read_in_rllleft, read_in_sllleft ! true or false, if wavefunction of this (iat,ie) pair was read in or not, respectively

    integer :: isave

    isave = t_wavefunctions%isave_wavefun(iat, ie)

    read_in_rll = .false.
    read_in_sll = .false.
    read_in_rllleft = .false.
    read_in_sllleft = .false.

    if (isave>0) then

      if (t_wavefunctions%save_rll) then
        if (iprint_debug>0) write (*, *) 'read in rll', ie, iat
        rll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = t_wavefunctions%rll(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
        read_in_rll = .true.
      end if

      if (t_wavefunctions%save_rllleft) then
        if (iprint_debug>0) write (*, *) 'read in rllleft', ie, iat
        rllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = t_wavefunctions%rllleft(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
        read_in_rllleft = .true.
      end if

      if (t_wavefunctions%save_sll) then
        if (iprint_debug>0) write (*, *) 'read in sll', ie, iat
        sll(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = t_wavefunctions%sll(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
        read_in_sll = .true.
      end if

      if (t_wavefunctions%save_sllleft) then
        if (iprint_debug>0) write (*, *) 'read in sllleft', ie, iat
        sllleft(1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith) = t_wavefunctions%sllleft(isave, 1:nsra*lmmaxd, 1:lmmaxd, 1:irmdnew, ith)
        read_in_sllleft = .true.
      end if

    end if


  end subroutine read_wavefunc

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Broadcast parameters for memory demand from master to all other ranks
  !> Author:
  !> Category: communication, KKRhost
  !> Deprecated: False
  !> Broadcast parameters for memory demand from master to all other ranks
  !-------------------------------------------------------------------------------
  subroutine bcast_params_savewf(t_wavefunctions)
    ! broadcast parameters for memory demand from master to all other ranks
    use :: mpi
    implicit none

    type (type_wavefunctions), intent (inout) :: t_wavefunctions
    integer :: ierr

    call mpi_bcast(t_wavefunctions%maxmem_number, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting maxmem_number'
    call mpi_bcast(t_wavefunctions%maxmem_units, 1, mpi_integer, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting maxmem_units'
    call mpi_bcast(t_wavefunctions%save_rll, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_rll'
    call mpi_bcast(t_wavefunctions%save_sll, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_sll'
    call mpi_bcast(t_wavefunctions%save_rllleft, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_rllleft'
    call mpi_bcast(t_wavefunctions%save_sllleft, 1, mpi_logical, 0, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop '[bcast_params_savewf] Error broadcasting save_sllleft'

  end subroutine bcast_params_savewf
#endif

end module mod_save_wavefun
