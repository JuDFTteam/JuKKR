module mod_save_wavefun

implicit none

   type :: type_wavefunctions
 
      integer :: Nwfsavemax    ! number of wavefunctions (na*ne) that this rank has to store
      integer :: maxmem_units  ! units of maxmem_number (1024=kbyte, 1024**2=Mbyte, 1024**3=Gbyte)
      integer :: maxmem_number ! maximal amount of memory available to store wavefunctions (maxmem_number*maxmem_units byte)
      integer :: nth           ! number of OpenMP tasks used without MPI (not the number of tasks in hybrid mode) in OpenMP mode only
      
      ! allocatable arrays
      integer, allocatable :: isave_wavefun(:,:)    ! 0 if (iat_myrank, ie_myrank) pair is not saved, 1...Nwfsavemax otherwise, for first, second, ... Nwfsavemax-th saved wavefunction on this rank; (Nat_myrank, Ne_myrank)
      double complex, allocatable :: rll(:,:,:,:,:)     ! regular right wavefunction; (Nwfsavemax, NSRA*LMMAXSO, LMMAXSO, IRMDNEW, 0:nth-1)
      double complex, allocatable :: rllleft(:,:,:,:,:) ! regular left wavefunction; (Nwfsavemax, NSRA*LMMAXSO, LMMAXSO, IRMDNEW, 0:nth-1)
      double complex, allocatable :: sll(:,:,:,:,:)     ! iregular right wavefunction; (Nwfsavemax, NSRA*LMMAXSO, LMMAXSO, IRMDNEW, 0:nth-1)
      double complex, allocatable :: sllleft(:,:,:,:,:) ! iregular left wavefunction; (Nwfsavemax, NSRA*LMMAXSO, LMMAXSO, IRMDNEW, 0:nth-1)

   end type type_wavefunctions

   type (type_wavefunctions), save :: t_wavefunctions

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine find_isave_wavefun(t_wavefunctions)
   ! find how many wavefunctions can be stored and set isave_wavefun array to determine which wavefunctions are stored
   
#ifdef CPP_OMP
      use omp_lib
#endif
      use mod_mympi, only: myrank, master, nranks
      use mod_types, only: t_inc, t_mpi_c_grid
      
      implicit none
   
      type(type_wavefunctions), intent(inout) :: t_wavefunctions
      
      integer :: nth, nat, ne, ierr, i, ie, iat, maxpos(1), Nwfsave_count, nat_myrank, ne_myrank, ioff_at, ioff_ie
      double precision :: delta_mem
      
      integer, allocatable :: kmesh_priority(:), my_kmesh(:)
      logical, allocatable :: mask(:)
      
      
      !>>>>>>>  set some parameters  >>>>>>>!
      nth = 1
#ifdef CPP_OMP
      nth = omp_get_num_threads()
#endif
      t_wavefunctions%nth = nth

      nat = t_inc%natyp
      ne =  t_inc%ielast
#ifdef CPP_MPI
      ! local parameter for myrank
      nat_myrank = t_mpi_c_grid%ntot1
      ne_myrank  = t_mpi_c_grid%ntot2
      ioff_at = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)
      ioff_ie = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
#else
      nat_myrank = nat
      ne_myrank  = ne
      ioff_at = 0
      ioff_ie = 0
#endif
      
      allocate(t_wavefunctions%isave_wavefun(nat, ne), stat=ierr)
      if(ierr/=0) stop '[find_isave_wavefun] Error allocating isave_wavefun'
      !<<<<<<<  set some parameters  <<<<<<<!
      
      
      !>>>>>>>  find numer of wavefunctions that can be stored and allocate store arrays  >>>>>>>!
      ! memory demand for one atom and one energy point in Mbyte
      delta_mem = dfloat(t_inc%NSRA*t_inc%LMMAXSO * t_inc%LMMAXSO * t_inc%IRMDNEW * nth * 16) / (1024.0d0**2)
      ! find numer of wavefunctions that can be stored
      t_wavefunctions%Nwfsavemax = t_wavefunctions%maxmem_number * 1024**(t_wavefunctions%maxmem_units-2) / dint( 4*delta_mem )
      if(myrank==master) then
         write(1337,"(A,/A,F15.2,A,/A,I11,/A,I11,A,I11,A)") '   ==> find_isave_wavefun ', &
         &'   (maxmem given for storage:',dfloat((1024**(t_wavefunctions%maxmem_units-2))*t_wavefunctions%maxmem_number), 'MB', &
         &'    number of wavefunctions (factor 4 for rll, sll left and right already in) that fit in:',t_wavefunctions%Nwfsavemax, &
         &'    Ne=',ne,', Nat=',nat,''
         write(1337,'(A,F15.2,A)') '    memory demand per atom and energy point for rll, rllleft, sll and sllleft respectively:', delta_mem, 'MB)  <=='
      end if
      
      ! avoid unnessesary large allocations of arrays
      if( t_wavefunctions%Nwfsavemax > nat_myrank*ne_myrank ) then
         t_wavefunctions%Nwfsavemax = nat_myrank*ne_myrank
         if(t_inc%i_write>0) write(1337,'(A,I5,A,I9)') '  rank',myrank,' reset Nwfsavemax to maximal needed number for this thread:',t_wavefunctions%Nwfsavemax
      end if

      
      if(t_wavefunctions%Nwfsavemax>0) then
         ! allocate store arrays int t_wavefunctions
         allocate(t_wavefunctions%rll(t_wavefunctions%Nwfsavemax, t_inc%NSRA*t_inc%LMMAXSO, t_inc%LMMAXSO, t_inc%IRMDNEW, 0:nth-1), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating rll'
         allocate(t_wavefunctions%rllleft(t_wavefunctions%Nwfsavemax, t_inc%NSRA*t_inc%LMMAXSO, t_inc%LMMAXSO, t_inc%IRMDNEW, 0:nth-1), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating rllleft'
         allocate(t_wavefunctions%sll(t_wavefunctions%Nwfsavemax, t_inc%NSRA*t_inc%LMMAXSO, t_inc%LMMAXSO, t_inc%IRMDNEW, 0:nth-1), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating sll'
         allocate(t_wavefunctions%sllleft(t_wavefunctions%Nwfsavemax, t_inc%NSRA*t_inc%LMMAXSO, t_inc%LMMAXSO, t_inc%IRMDNEW, 0:nth-1), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating sllleft'
      end if
      !<<<<<<<  find numer of wavefunctions that can be stored and allocate store arrays  <<<<<<<!


      if(t_wavefunctions%Nwfsavemax>0) then
         !>>>>>>>  find priority of saving wavefunctions according to kmesh  >>>>>>>!
         ! kmesh( 1 : ne )
         allocate(my_kmesh(ne_myrank), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating my_kmesh'
         my_kmesh(:) = t_inc%kmesh_ie( 1+ioff_ie : ne_myrank+ioff_ie)
         maxpos(1) = 0
         do i=1,ne_myrank
            if(my_kmesh(i)>maxpos(1)) maxpos(1) = my_kmesh(i)
         end do
         my_kmesh(:) = maxpos(1) - my_kmesh(:)
         
         allocate(kmesh_priority(ne_myrank), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating kmesh_priority'
         allocate(mask(ne_myrank), stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error allocating mask'
         
         mask(:) = .true.
         do i=1,ne_myrank
            ! find position of maximal value in my_kmesh and fill array kmesh_priority
            maxpos = maxloc(my_kmesh, mask=mask)
            kmesh_priority(maxpos(1)) = i
            ! update mask to exclude previously found position
            mask(maxpos(1)) = .false.
         end do ! i=1,ne_myrank

         if(any(kmesh_priority==-1)) then
            write(1337,*) '[find_isave_wavefun] kmesh_priority not found correctly; kmesh_priority:', kmesh_priority, '; mask:', mask
            stop 
         end if
         !<<<<<<<  find priority of saving wavefunctions according to kmesh  <<<<<<<!
         
         
         !>>>>>>>  set isave_wavefun array  >>>>>>>!
         t_wavefunctions%isave_wavefun(:,:) = 0
         Nwfsave_count = t_wavefunctions%Nwfsavemax
         do iat=1,nat_myrank
            t_wavefunctions%isave_wavefun(iat+ioff_at,1+ioff_ie:ioff_ie+ne_myrank) = kmesh_priority(:) + ne_myrank*(iat-1)
            do ie=1,ne_myrank
               if( (kmesh_priority(ie)>Nwfsave_count) ) then
                  t_wavefunctions%isave_wavefun(iat+ioff_at, ie+ioff_ie) = 0
               end if
            end do
            Nwfsave_count = Nwfsave_count - ne_myrank
         end do
         
!          do ie=1,ne_myrank
!             do iat=1,nat_myrank
!                write(*,*) 'isave_wavefun',myrank, iat, ie, t_wavefunctions%isave_wavefun(iat+ioff_at, ie+ioff_ie), my_kmesh(ie)
!             end do
!          end do
         !<<<<<<<  set isave_wavefun array  <<<<<<<!
              
              
         ! deallocate work arrays
         deallocate(mask, kmesh_priority, my_kmesh, stat=ierr)
         if(ierr/=0) stop '[find_isave_wavefun] Error dellocating arrays'
      end if
      
   end subroutine find_isave_wavefun
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   subroutine save_wavefunc(t_wavefunctions, rll, rllleft, sll, sllleft, iat, ie, NSRA, LMMAXSO, IRMDNEW, ith, nth)
      !saves wavefunction of atom iat and energypoint ie if enough memory is given

      use mod_mympi, only: myrank
      implicit none
      type(type_wavefunctions), intent(inout) :: t_wavefunctions
      integer, intent(in) :: iat, ie, NSRA, LMMAXSO, IRMDNEW, ith, nth
      double complex, intent(in) :: rll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), rllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), sll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), sllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1))
      
      integer :: isave
      
      isave = t_wavefunctions%isave_wavefun(iat, ie)
      
      if(isave>0) then
         t_wavefunctions%rll(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = rll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         t_wavefunctions%rllleft(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = rllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         t_wavefunctions%sll(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = sll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         t_wavefunctions%sllleft(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = sllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
      end if

   end subroutine save_wavefunc
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   subroutine read_wavefunc(t_wavefunctions, rll, rllleft, sll, sllleft, iat, ie, NSRA, LMMAXSO, IRMDNEW, ith, nth, read_in)
      !reads wavefunction of atom iat and energypoint ie if it was stored
   
      use mod_mympi, only: myrank
      implicit none
      type(type_wavefunctions), intent(inout) :: t_wavefunctions
      integer, intent(in) :: iat, ie, NSRA, LMMAXSO, IRMDNEW, ith, nth
      double complex, intent(out) :: rll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), rllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), sll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1)), sllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,0:(nth-1))
      logical, intent(out) :: read_in ! true or false, if wavefunction of this (iat,ie) pair was read in or not, respectively
      
      integer :: isave
      
      isave = t_wavefunctions%isave_wavefun(iat, ie)
      
      if(isave>0) then
         rll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = t_wavefunctions%rll(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         rllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = t_wavefunctions%rllleft(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         sll(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = t_wavefunctions%sll(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         sllleft(1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith) = t_wavefunctions%sllleft(isave,1:NSRA*LMMAXSO,1:LMMAXSO,1:IRMDNEW,ith)
         read_in = .true.
      else
         read_in = .false.
      end if


   end subroutine read_wavefunc
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
#ifdef CPP_MPI
   subroutine bcast_params_savewf(t_wavefunctions)
      !broadcast parameters for memory demand from master to all other ranks
      use mpi
      use mod_mympi, only: myrank
      implicit none
      
      type(type_wavefunctions), intent(inout) :: t_wavefunctions
      integer :: ierr
      
      call MPI_Bcast(t_wavefunctions%maxmem_number, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop '[bcast_params_savewf] Error broadcasting maxmem_number'
      call MPI_Bcast(t_wavefunctions%maxmem_units, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if(ierr/=MPI_SUCCESS) stop '[bcast_params_savewf] Error broadcasting maxmem_units'
            
   end subroutine bcast_params_savewf
#endif
   
         

end module mod_save_wavefun
