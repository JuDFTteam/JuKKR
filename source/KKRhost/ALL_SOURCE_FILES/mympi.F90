module mod_mympi
  ! ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann

  use :: mod_datatypes, only: dp

  implicit none

  private
  public :: myrank, nranks, master, mympi_init, mpiatom, mpiadapt, distribute_work_atoms, distribute_work_energies
#ifdef CPP_MPI
  public :: distribute_linear_on_tasks, find_dims_2d, create_newcomms_group_ie, mympi_main1c_comm, mympi_main1c_comm_newsosol, mympi_main1c_comm_newsosol2, &
    check_communication_pattern
#endif

  integer, save :: myrank = -1
  integer, save :: nranks = -1
  integer, save :: master = -1
  logical, save :: mpiatom = .false.
  integer, save :: mpiadapt = -1

contains

  subroutine mympi_init()

#ifdef CPP_MPI
    use :: mpi

    integer :: ierr
#endif

    master = 0

#ifdef CPP_MPI
    call mpi_comm_rank(mpi_comm_world, myrank, ierr)
    call mpi_comm_size(mpi_comm_world, nranks, ierr)
#else
    myrank = master
    nranks = 1
#endif

  end subroutine mympi_init

#ifdef CPP_MPI
  subroutine distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pt, ioff_pt, output, fill_rest)
    ! distributes 'ntot' points on 'nranks' tasks and returns the number of points per task, ntot_pT, and the offsets of the tasks, ioff_pT.

    implicit none

    integer, intent (in) :: nranks, myrank, master, ntot
    logical, intent (in) :: output
    logical, intent (in), optional :: fill_rest
    integer, intent (out) :: ntot_pt(0:nranks-1), ioff_pt(0:nranks-1)

    integer :: irest, irank

    ntot_pt = int(ntot/nranks)
    ioff_pt = int(ntot/nranks)* [ (irank,irank=0,nranks-1) ]
    irest = ntot - int(ntot/nranks)*nranks

    if (irest>0) then

      do irank = 0, irest - 1
        ntot_pt(irank) = ntot_pt(irank) + 1
        ioff_pt(irank) = ioff_pt(irank) + irank
      end do                       ! irank

      do irank = irest, nranks - 1
        ioff_pt(irank) = ioff_pt(irank) + irest
      end do                       ! irank

    end if                         ! irest>0

    if (present(fill_rest)) then
      if (fill_rest) then
        if (ntot<nranks) then
          ! write(*,*) 'set rest',myrank,ntot,nranks
          do irank = ntot, nranks - 1
            ioff_pt(irank) = ioff_pt(irank-1)
            ntot_pt(irank) = ntot_pt(irank-1)
          end do                   ! irank
        end if                     ! ntot<nranks
      end if                       ! fill_rest
    end if                         ! present(fill_rest)

    if (myrank==master .and. output) then
      write (1337, *) '==== DISTRIBUTION OF POINTS ON TASKS: ===='
      do irank = 0, nranks - 1
        write (1337, '("Task ",I0," treats points ",I0," to ",I0,",#of points= ",I0)') irank, ioff_pt(irank) + 1, ioff_pt(irank) + ntot_pt(irank), ntot_pt(irank)
      end do                       ! irank
      write (1337, *) '=========================================='
    end if                         ! myrank==master

  end subroutine distribute_linear_on_tasks
#endif


#ifdef CPP_MPI
  subroutine find_dims_2d(nranks, ntot1, ntot2, dims, mpiatom)
    ! find dimensions to create cartesian communicator
    ! input:  nranks, ntot1 is N_atom, ntot2 is N_E
    ! output: dims(2), dims(1) is N_atomranks, dims(2) is N_Eranks
    use :: mpi

    implicit none
    integer, intent (in) :: nranks, ntot1, ntot2
    logical, intent (in) :: mpiatom
    integer, intent (out) :: dims(2)
    integer :: ierr

    if (.not. mpiatom) then
      if (nranks<=ntot2) then
        dims(1) = 1
        dims(2) = nranks
      else
        dims(1) = nranks/ntot2
        dims(2) = ntot2
      end if
    else
      if (nranks<=ntot1) then
        dims(2) = 1
        dims(1) = nranks
      else
        dims(2) = nranks/ntot1
        dims(1) = ntot1
      end if
    end if

    if (nranks>(ntot1*(ntot2+1)-1)) then
      if (myrank==master) write (*, '(A,I3,A,I3,A,I5)') 'Error for', ntot1, ' atoms and', ntot2, ' energy points you use too many processors. Nranks=', nranks
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize(ierr)
      stop 'Error: too many ranks'
    end if

  end subroutine find_dims_2d
#endif


#ifdef CPP_MPI
  subroutine create_newcomms_group_ie(master, nranks, myrank, nat, ne, nkmesh, kmesh, mympi_comm_ie, myrank_ie, nranks_ie, mympi_comm_at, myrank_at, nranks_at, myrank_atcomm, &
    nranks_atcomm)
    ! takes vector kmesh with mesh/timing information and finds number of rest procs that are devided in fractions given in ktake for optimal division of work

    use :: mpi
    implicit none

    integer, intent (in) :: master, nranks, myrank, ne, nat, nkmesh
    integer, intent (in) :: kmesh(nkmesh)
    integer, intent (out) :: mympi_comm_ie, mympi_comm_at, nranks_atcomm, nranks_at, nranks_ie, myrank_ie, myrank_at, myrank_atcomm

    integer :: rest, k(nkmesh-1), ktake(nkmesh), myg, ie, iat, ierr, ik, i1, i2, i3
    real (kind=dp) :: f(nkmesh-1), q, qmin
    integer, allocatable :: groups(:, :), mygroup(:)
    integer :: mympi_group_ie, mympi_group_world, mympi_comm_grid

    rest = 0
    qmin = -1
    if (myrank==master) write (1337, *) 'create_newcomms_group_ie input:', nranks, ne, nat

    ktake(:) = 0

    if ((ne*nat)<nranks .and. (ne>1)) then ! .and. nat>1)) then

      if (nkmesh<=1) then
        if (myrank==master) write (*, '(A,2I7)') &
          'no load imbalance found (all energy points have the same k-mesh), please use regular grid to not waste any resources.#E,#atoms = ', ne, nat
      end if

      rest = nranks - int(nranks/(ne*nat))*ne*nat
      if (myrank==master) write (1337, *) 'rest:', rest, ne, nat, nranks
      if (myrank==master) write (1337, *) 'kmesh:', kmesh

      if (rest>0 .and. nkmesh>1) then
        ! find fraction of k:l:m
        do ik = 1, nkmesh - 1
          k(ik) = int(real(kmesh(1))/real(kmesh(nkmesh-ik+1))-1.)
          if (k(ik)==0) k(ik) = 1
        end do

        do ik = 1, nkmesh - 2
          f(ik) = real(k(ik+1))/real(k(ik))
        end do
        f(nkmesh-1) = real(k(nkmesh-1))/real(k(1))

        if (myrank==master) write (1337, *) 'set k,i:', k, 'f', f

        ! brute force look for optimal division of rest ranks after N_E*N_at are already
        ! assigned to rectangular part of processor matrix:
        ! N_E=8
        ! --------------->
        ! ^ ( | | | | | | | )    example for 49 processors,
        ! | ( | | | | | | | )    devided according to:
        ! N_at=5  | ( | | | | | | | )         N_E = 8, N_at = 5
        ! | ( | | | | | | | )        rest = 9 = 5+3+1
        ! v ( | | | | | | | )                   k+l+m
        ! ^    ( | | ) m=1  ^
        ! l=3 |      ( | )      |
        ! v      ( | )      | k=5
        ! ( )      |
        ! ( )      v

        if (nkmesh==4) then
          ktake(:) = 0
          do i1 = 1, rest
            do i2 = 0, rest - i1
              i3 = rest - i1 - i2
              if (i1>=i2 .and. i2>=i3) then
                if (i3==0 .and. i2==0) then
                  q = sqrt((f(1)-real(i2)/real(i1))**2+(f(2)-1.)**2+(f(3)-real(i3)/real(i1))**2)
                else
                  q = sqrt((f(1)-real(i2)/real(i1))**2+(f(2)-real(i3)/real(i2))**2+(f(3)-real(i3)/real(i1))**2)
                end if
                if (q<qmin .or. qmin==-1) then
                  ktake = [ i1, i2, i3 ]
                  qmin = q
                end if
              end if
            end do
          end do
        else if (nkmesh==3) then
          ktake(:) = 0
          do i1 = 1, rest
            i2 = rest - i1
            if (i1>=i2) then
              q = sqrt((f(1)-real(i2)/real(i1))**2)
              if (q<qmin .or. qmin==-1) then
                ktake = [ i1, i2 ]
                qmin = q
              end if
            end if
          end do
        else if (nkmesh==2) then
          ktake(1) = rest
        else
          stop 'ERROR: nkmesh>4 not implemented yet'
        end if

        ! special case when only one additional rank
        if (rest==1) ktake(1) = rest

        if (myrank==master) write (1337, *) 'found ktake', ktake, 'with', qmin

      else if (rest>0) then
        ktake(1) = rest
      end if                       ! if(rest>0 .and. nkmesh>1)

      ! find processor groups according to non-uniform division
      allocate (groups(ne+1,2))
      groups(:, :) = -1

      do ie = 1, ne + 1
        if (ie==ne-2 .and. nkmesh>3) then
          groups(ie, 1) = nat + ktake(3)
        else if (ie==ne-1 .and. nkmesh>2) then
          groups(ie, 1) = nat + ktake(2)
        else if (ie==ne-0 .and. nkmesh>=1) then
          groups(ie, 1) = nat + ktake(1)
        else
          groups(ie, 1) = nat
        end if
        if (ie==1) then
          groups(ie, 2) = 0
        else
          groups(ie, 2) = groups(ie-1, 1) + groups(ie-1, 2)
        end if
      end do

      if (myrank==master) write (1337, *) 'groups(1:Ne), number of ranks:', groups(1:ne, 1)
      if (myrank==master) write (1337, *) 'groups(1:Ne), ie offset:', groups(1:ne, 2)

      ! find my group
      myg = -1
      do ie = 1, ne
        do iat = groups(ie, 2), groups(ie+1, 2) - 1
          if (myrank==iat) myg = ie
        end do
      end do

      if (myg==-1) then
        write (1337, *) 'no group found for rank', myrank
        stop
      end if

      ! get group of processors in my group
      allocate (mygroup(groups(myg,1)))


      ie = 0
      do iat = groups(myg, 2), groups(myg+1, 2) - 1
        ie = ie + 1
        mygroup(ie) = iat
      end do


      ! create new communicator from group
      call mpi_comm_group(mpi_comm_world, mympi_group_world, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_GROUP with code ', ierr
        write (*, *) 'Error in MPI_COMM_GROUP with code ', ierr
        stop 'Error in MPI_COMM_GROUP'
      end if

      call mpi_group_incl(mympi_group_world, groups(myg,1), mygroup, mympi_group_ie, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_GROUP_INCL with code ', ierr
        write (*, *) 'Error in MPI_GROUP_INCL with code ', ierr
        stop 'Error in MPI_GROUP_INCL'
      end if
      call mpi_comm_create(mpi_comm_world, mympi_group_ie, mympi_comm_ie, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_CREATE with code ', ierr
        write (*, *) 'Error in MPI_COMM_CREATE with code ', ierr
        stop 'Error in MPI_COMM_CREATE'
      end if
      call mpi_group_free(mympi_group_ie, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_GROUP_FREE with code ', ierr
        write (*, *) 'Error in MPI_GROUP_FREE with code ', ierr
        stop 'Error in MPI_GROUP_FREE'
      end if

      ! get rank and size in new communicator
      call mpi_comm_rank(mympi_comm_ie, myrank_ie, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_RANK(ie) with code ', ierr
        write (*, *) 'Error in MPI_COMM_RANK(ie) with code ', ierr
        stop 'Error in MPI_COMM_RANK(ie)'
      end if
      call mpi_comm_size(mympi_comm_ie, nranks_ie, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_SIZE(ie) with code ', ierr
        write (*, *) 'Error in MPI_COMM_SIZE(ie) with code ', ierr
        stop 'Error in MPI_COMM_SIZE(ie)'
      end if

      ! create communicator to communicate between differen energies (i.e. different groups)
      call mpi_comm_split(mpi_comm_world, myrank_ie, myg, mympi_comm_at, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_SPLIT', ierr
        write (*, *) 'Error in MPI_COMM_SPLIT ', ierr
        stop 'Error in MPI_COMM_SPLIT'
      end if
      call mpi_comm_rank(mympi_comm_at, myrank_atcomm, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_RANK(at) with code ', ierr
        write (*, *) 'Error in MPI_COMM_RANK(at) with code ', ierr
        stop 'Error in MPI_COMM_RANK(at)'
      end if
      call mpi_comm_size(mympi_comm_at, nranks_atcomm, ierr)
      if (ierr/=mpi_success) then
        write (*, *) 'Error in MPI_COMM_SIZE(at) with code ', ierr
        write (*, *) 'Error in MPI_COMM_SIZE(at) with code ', ierr
        stop 'Error in MPI_COMM_SIZE(at)'
      end if

      nranks_at = ne
      myrank_at = myg - 1

    else if ((ne*nat)==nranks .and. (ne>1 .and. nat>1)) then

      rest = 0

      if (myrank==master) write (1337, *) 'create cartesian grid:', ne, nat, nranks
      call mpi_cart_create(mpi_comm_world, 2, [ne,nat], [.false.,.false.], [.true.,.true.], mympi_comm_grid, ierr)

      if (myrank==master) write (1337, *) 'MPI_Cart_sub'
      call mpi_cart_sub(mympi_comm_grid, [.true.,.false.], mympi_comm_at, ierr) ! row communicator
      call mpi_cart_sub(mympi_comm_grid, [.false.,.true.], mympi_comm_ie, ierr) ! col communicator

      if (myrank==master) write (1337, *) 'MPI_Comm_rank'
      call mpi_comm_rank(mympi_comm_ie, myrank_ie, ierr)
      call mpi_comm_rank(mympi_comm_at, myrank_at, ierr)

      if (myrank==master) write (1337, *) 'MPI_Comm_size'
      call mpi_comm_size(mympi_comm_ie, nranks_ie, ierr)
      call mpi_comm_size(mympi_comm_at, nranks_at, ierr)

      myrank_atcomm = myrank_at
      nranks_atcomm = nranks_at

    else if (ne==1 .and. nat==nranks) then

      rest = 0

      mympi_comm_ie = mpi_comm_world
      myrank_ie = myrank
      nranks_ie = nranks
      mympi_comm_at = mpi_comm_self
      myrank_at = 0
      nranks_at = 1

      myrank_atcomm = myrank_at
      nranks_atcomm = nranks_at

    else

      rest = 0

      mympi_comm_at = mpi_comm_world
      myrank_at = myrank
      nranks_at = nranks
      mympi_comm_ie = mpi_comm_self
      myrank_ie = 0
      nranks_ie = 1

      myrank_atcomm = myrank_at
      nranks_atcomm = nranks_at


    end if


    if (myrank==master) then
      write (1337, '(A)') '=================================================='
      write (1337, '(A,I5,A)') '    MPI parallelization: use', nranks, ' ranks'
      write (1337, '(A,I3,A,I4)') '    create processor array of size (nat x ne) ', nat, ' x', ne
      write (1337, '(A,I5,A,I5)') '    nranks_at: ', nranks_at, ', nranks_ie:', nranks_ie
      if (rest>0) write (1337, '(A,I3)') '                                   with rest', rest
      if (rest>0) write (1337, '(A,10I3)') '    divide rest onto last energy points (k,l,m):', ktake
      write (1337, '(A)') '                N_E'
      write (1337, '(A)') '         <--------------->'
      write (1337, '(A)') '       ^ ( | | | | | | | )'
      write (1337, '(A)') '       | ( | | | | | | | )'
      write (1337, '(A)') '  N_at | ( | | | | | | | )'
      write (1337, '(A)') '       | ( | | | | | | | )'
      if (rest==0) write (1337, '(A)') '       v ( | | | | | | | )'
      if (rest>0) write (1337, '(A)') '       v ( | | | | | | | )....'
      if (rest>0) write (1337, '(A)') '              ^    ( | | ) m  ^'
      if (rest>0) write (1337, '(A)') '            l |      ( | )    |'
      if (rest>0) write (1337, '(A)') '              v......( | )    | k'
      if (rest>0) write (1337, '(A)') '                       ( )    |'
      if (rest>0) write (1337, '(A)') '                       ( )....v'
    end if

  end subroutine create_newcomms_group_ie
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm(irmd, lmpotd, natypd, lmaxd, lmaxd1, lmmaxd, npotd, ielast, mmaxd, idoldau, natyp, krel, lmomvec, nmvecmax, nqdos, rho2ns, r2nef, espv, den, denlm, &
    denmatc, denef, denefat, rhoorb, muorb, mvevi, mvevil, mvevief, mympi_comm)

    use :: mpi
    implicit none
    integer, intent (in) :: irmd, lmpotd, natypd, lmaxd, lmmaxd, ielast, mmaxd, idoldau, natyp, krel, nmvecmax, npotd, lmaxd1, nqdos
    logical, intent (in) :: lmomvec
    integer, intent (in) :: mympi_comm
    real (kind=dp), intent (inout) :: rho2ns(irmd, lmpotd, natypd, 2), r2nef(irmd, lmpotd, natypd, 2), espv(0:lmaxd1, npotd), denef, denefat(natypd), &
      rhoorb(irmd*krel+(1-krel), natypd), muorb(0:lmaxd1+1, 3, natypd)
    complex (kind=dp), intent (inout) :: den(0:lmaxd1, ielast, npotd, nqdos), denlm(lmmaxd, ielast, npotd, nqdos), denmatc(mmaxd, mmaxd, npotd), mvevi(natypd, 3, nmvecmax), &
      mvevil(0:lmaxd, natypd, 3, nmvecmax), mvevief(natypd, 3, nmvecmax)

    integer :: idim, ierr          ! , myrank_comm
    integer, parameter :: master = 0
    real (kind=dp), allocatable :: work1(:), work2(:, :), work3(:, :, :), work4(:, :, :, :)
    complex (kind=dp), allocatable :: work3c(:, :, :), work4c(:, :, :, :)


    allocate (work4(irmd,lmpotd,natypd,2), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4 = 0.d0
    idim = irmd*lmpotd*natypd*2
    call mpi_allreduce(rho2ns, work4, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work4, 1, rho2ns, 1)
    deallocate (work4)


    allocate (work4(irmd,lmpotd,natypd,2), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4 = 0.d0
    idim = irmd*lmpotd*natypd*2
    call mpi_allreduce(r2nef, work4, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work4, 1, r2nef, 1)
    deallocate (work4)

    ! ESPV needs integration over atoms and energies -> MPI_COMM_WORLD
    allocate (work2(0:lmaxd+1,npotd), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work2 = 0.d0
    idim = (lmaxd+2)*npotd
    call mpi_allreduce(espv, work2, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work2, 1, espv, 1)
    deallocate (work2)

    allocate (work4c(0:lmaxd+1,ielast,npotd,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4c = (0.d0, 0.d0)
    idim = ielast*(lmaxd+2)*npotd*nqdos
    call mpi_allreduce(den, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(idim, work4c, 1, den, 1)
    deallocate (work4c)

    allocate (work4c(ielast,lmmaxd,npotd,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4c = (0.d0, 0.d0)
    idim = ielast*(lmmaxd)*npotd*nqdos
    call mpi_allreduce(denlm, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(idim, work4c, 1, denlm, 1)
    deallocate (work4c)

    if (idoldau==1) then
      allocate (work3c(mmaxd,mmaxd,npotd), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work3c = (0.d0, 0.d0)
      idim = mmaxd*mmaxd*npotd
      call mpi_allreduce(denmatc, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
      call zcopy(idim, work3c, 1, denmatc, 1)
      deallocate (work3c)
    end if

    allocate (work1(1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work1 = 0.d0
    idim = 1
    call mpi_allreduce(denef, work1, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work1, 1, denef, 1)
    deallocate (work1)

    allocate (work1(natyp), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work1 = 0.d0
    idim = natyp
    call mpi_allreduce(denefat, work1, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work1, 1, denefat, 1)
    deallocate (work1)

    if (krel==1) then
      allocate (work2(irmd,natypd), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work2 = 0.d0
      idim = irmd*natypd
      call mpi_allreduce(rhoorb, work2, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
      call dcopy(idim, work2, 1, rhoorb, 1)
      deallocate (work2)

      allocate (work3(0:lmaxd+2,natypd,3), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work3 = 0.d0
      idim = (lmaxd+3)*natypd*3
      call mpi_allreduce(muorb, work3, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
      call dcopy(idim, work3, 1, muorb, 1)
      deallocate (work3)

      if (lmomvec) then
        allocate (work3c(natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work3c = (0.d0, 0.d0)
        idim = natypd*3*nmvecmax
        call mpi_allreduce(mvevi, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work3c, 1, mvevi, 1)
        deallocate (work3c)

        allocate (work4c(lmaxd+1,natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work4c = (0.d0, 0.d0)
        idim = (lmaxd+1)*natypd*3*nmvecmax
        call mpi_allreduce(mvevil, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work4c, 1, mvevil, 1)
        deallocate (work4c)

        allocate (work3c(natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work3c = (0.d0, 0.d0)
        idim = natypd*3*nmvecmax
        call mpi_allreduce(mvevief, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work3c, 1, mvevief, 1)
        deallocate (work3c)
      end if                       ! LMOMVEC
    end if                         ! KREL.EQ.1

  end subroutine mympi_main1c_comm
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm_newsosol(nspin, korbit, irmdnew, lmpotd, lmaxd, lmaxd1, lmmaxd, lmmaxso, ielast, nqdos, den, denlm, gflle, rho2nsc, r2nefc, rho2int, espv, muorb, denorbmom, &
    denorbmomsp, denorbmomlm, denorbmomns, mympi_comm)

    use :: mpi
    implicit none
    integer, intent (in) :: nspin, korbit, irmdnew, lmpotd, lmaxd, lmaxd1, lmmaxd, lmmaxso, ielast, nqdos
    integer, intent (in) :: mympi_comm
    complex (kind=dp), intent (inout) :: r2nefc(irmdnew, lmpotd, nspin*(1+korbit)), rho2nsc(irmdnew, lmpotd, nspin*(1+korbit)), den(0:lmaxd1, ielast, nqdos, nspin), denlm(lmmaxd, ielast, nqdos, nspin), rho2int(nspin*(1+korbit)), &
      gflle(lmmaxso, lmmaxso, ielast, nqdos)
    real (kind=dp), intent (inout) :: espv(0:lmaxd1, 2), muorb(0:lmaxd1+1, 3), denorbmom(3), denorbmomsp(2, 3), denorbmomlm(0:lmaxd, 3), denorbmomns(3)

    integer :: ierr, idim
    real (kind=dp), allocatable :: work(:, :, :, :)
    complex (kind=dp), allocatable :: workc(:, :, :, :)


    ! all with reduce instead of allreduce:
    ! complex (kind=dp) arrays
    idim = irmdnew*lmpotd*4
    allocate (workc(irmdnew,lmpotd,4,1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, r2nefc'
    workc = (0.d0, 0.d0)
    call mpi_reduce(r2nefc, workc(:,:,:,1), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for r2nefc'
    call zcopy(idim, workc, 1, r2nefc, 1)
    deallocate (workc)

    idim = irmdnew*lmpotd*4
    allocate (workc(irmdnew,lmpotd,4,1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, rho2nsc'
    workc = (0.d0, 0.d0)
    call mpi_reduce(rho2nsc, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for rho2nsc'
    call zcopy(idim, workc, 1, rho2nsc, 1)
    deallocate (workc)

    idim = (lmaxd1+1)*ielast*2*nqdos
    allocate (workc(0:lmaxd1,ielast,2,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, den'
    workc = (0.d0, 0.d0)
    call mpi_reduce(den, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for den'
    call zcopy(idim, workc, 1, den, 1)
    deallocate (workc)

    idim = lmmaxd*ielast*2*nqdos
    allocate (workc(lmmaxd,ielast,2,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, denlm'
    workc = (0.d0, 0.d0)
    call mpi_reduce(denlm, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denlm'
    call zcopy(idim, workc, 1, denlm, 1)
    deallocate (workc)

    idim = 4
    allocate (workc(4,1,1,1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, rho2int'
    workc = (0.d0, 0.d0)
    call mpi_reduce(rho2int, workc(:,1,1,1), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for rho2int'
    call zcopy(idim, workc, 1, rho2int, 1)
    deallocate (workc)

    idim = lmmaxso*lmmaxso*ielast*nqdos
    allocate (workc(lmmaxso,lmmaxso,ielast,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, gflle'
    workc = (0.d0, 0.d0)
    call mpi_reduce(gflle, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for gflle'
    call zcopy(idim, workc, 1, gflle, 1)
    deallocate (workc)

    ! real (kind=dp) arrays
    idim = (lmaxd1+1)*2
    allocate (work(0:lmaxd1,2,1,1))
    work = 0.d0
    call mpi_reduce(espv, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for espv'
    call dcopy(idim, work, 1, espv, 1)
    deallocate (work)

    idim = (lmaxd1+2)*3
    allocate (work(0:lmaxd1+1,3,1,1))
    work = 0.d0
    call mpi_reduce(muorb, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for muorb'
    call dcopy(idim, work, 1, muorb, 1)
    deallocate (work)

    idim = 3
    allocate (work(3,1,1,1))
    work = 0.d0
    call mpi_reduce(denorbmom, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denobrmom'
    call dcopy(idim, work, 1, denorbmom, 1)
    deallocate (work)

    idim = 2*4
    allocate (work(2,4,1,1))
    work = 0.d0
    call mpi_reduce(denorbmomsp, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomsp'
    call dcopy(idim, work, 1, denorbmomsp, 1)
    deallocate (work)

    idim = 3
    allocate (work(3,1,1,1))
    work = 0.d0
    call mpi_reduce(denorbmomns, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomns'
    call dcopy(idim, work, 1, denorbmomns, 1)
    deallocate (work)

    idim = (lmaxd+1)*3
    allocate (work(0:lmaxd1,3,1,1))
    work = 0.d0
    call mpi_reduce(denorbmomlm, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomlm'
    call dcopy(idim, work, 1, denorbmomlm, 1)
    deallocate (work)

  end subroutine mympi_main1c_comm_newsosol
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm_newsosol2(lmaxd1, lmmaxd, ielast, nqdos, npotd, natypd, lmpotd, irmd, mmaxd, den, denlm, muorb, espv, r2nef, rho2ns, denefat, denef, denmatn, &
    angles_new, mympi_comm)

    use :: mpi
    implicit none
    integer, intent (in) :: lmaxd1, ielast, nqdos, npotd, natypd, lmpotd, irmd, lmmaxd, mmaxd
    integer, intent (in) :: mympi_comm
    complex (kind=dp), intent (inout) :: den(0:lmaxd1, ielast, nqdos, npotd), denlm(lmmaxd, ielast, nqdos, npotd), denmatn(mmaxd, mmaxd, 2, 2, natypd)
    real (kind=dp), intent (inout) :: muorb(0:lmaxd1+1, 3, natypd), espv(0:lmaxd1, npotd), r2nef(irmd, lmpotd, natypd, 2), rho2ns(irmd, lmpotd, natypd, 2), denefat(natypd), denef, &
      angles_new(natypd, 2)

    integer :: ierr, idim
    real (kind=dp), allocatable :: work(:, :, :, :)
    complex (kind=dp), allocatable :: workc(:, :, :, :)
    complex (kind=dp), allocatable :: workc1(:, :, :, :, :)

    ! complex (kind=dp) arrays
    idim = (1+lmaxd1)*ielast*nqdos*npotd
    allocate (workc(0:lmaxd1,ielast,nqdos,npotd))
    workc = (0.d0, 0.d0)
    call mpi_reduce(den, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(den,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for den'
    call zcopy(idim, workc, 1, den, 1)
    deallocate (workc)

    idim = lmmaxd*ielast*nqdos*npotd
    allocate (workc(lmmaxd,ielast,nqdos,npotd))
    workc = (0.d0, 0.d0)
    call mpi_reduce(denlm, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denlm,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denlm'
    call zcopy(idim, workc, 1, denlm, 1)
    deallocate (workc)

    idim = mmaxd*mmaxd*npotd
    allocate (workc1(mmaxd,mmaxd,2,2,natypd))
    workc1 = (0.d0, 0.d0)
    call mpi_reduce(denmatn, workc1(:,:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denmatn,workc1(:,:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denmatn'
    call zcopy(idim, workc1, 1, denmatn, 1)
    deallocate (workc1)

    ! real (kind=dp) arrays
    idim = (lmaxd1+2)*3*natypd
    allocate (work(0:lmaxd1+1,3,natypd,1))
    work = 0.d0
    call mpi_reduce(muorb, work(:,:,:,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(muorb,work(:,:,:,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for muorb'
    call dcopy(idim, work, 1, muorb, 1)
    deallocate (work)

    idim = (lmaxd1+1)*npotd
    allocate (work(0:lmaxd1,npotd,1,1))
    work = 0.d0
    call mpi_reduce(espv, work(:,:,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(ESPV,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for espv'
    call dcopy(idim, work, 1, espv, 1)
    deallocate (work)

    idim = irmd*lmpotd*natypd*2
    allocate (work(irmd,lmpotd,natypd,2))
    work = 0.d0
    call mpi_reduce(r2nef, work(:,:,:,:), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(r2nef,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for r2nef'
    call dcopy(idim, work, 1, r2nef, 1)
    deallocate (work)

    idim = irmd*lmpotd*natypd*2
    allocate (work(irmd,lmpotd,natypd,2))
    work = 0.d0
    call mpi_reduce(rho2ns, work(:,:,:,:), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(RHO2NS,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for rho2ns'
    call dcopy(idim, work, 1, rho2ns, 1)
    deallocate (work)

    idim = natypd
    allocate (work(natypd,1,1,1))
    work = 0.d0
    call mpi_reduce(denefat, work(:,1,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denefat,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denefat'
    call dcopy(idim, work, 1, denefat, 1)
    deallocate (work)

    idim = 1
    allocate (work(1,1,1,1))
    work = 0.d0
    call mpi_reduce(denef, work(:,1,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denef,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denef'
    call dcopy(idim, work, 1, denef, 1)
    deallocate (work)

    idim = 2*natypd
    allocate (work(2,natypd,1,1))
    work = 0.d0
    call mpi_reduce(angles_new, work(:,:,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(angles_new,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for angles_new'
    call dcopy(idim, work, 1, angles_new, 1)
    deallocate (work)

  end subroutine mympi_main1c_comm_newsosol2
#endif


#ifdef CPP_MPI
  subroutine check_communication_pattern(mpiatom, mpiadapt, timings_1a, timings_1b, load_imbalance, nkmesh, kmesh_ie)

    use :: mpi

    implicit none

    integer, intent (inout) :: mpiadapt
    logical, intent (inout) :: mpiatom
    real (kind=dp), intent (inout) :: timings_1a(:, :), timings_1b(:)
    integer, intent (out) :: load_imbalance(:)
    integer, intent (in) :: nkmesh, kmesh_ie(:)

    real (kind=dp), allocatable :: t_average(:), work(:, :)
    integer, allocatable :: kmesh_n(:)
    integer :: nat, ne, ne_1b, ierr, ie, iat, ik, iwork


    ! find some dimensions
    nat = size(timings_1a(1,:))
    ne = size(timings_1a(:,1))
    ne_1b = size(timings_1b)
    if (ne/=ne_1b) stop '[check_communication_pattern] Error in shapes of timing arrays'


    ! communicate timing arrays
    ! timings_1a
    iwork = ne*nat
    allocate (work(ne,nat), stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] error allocating work'
    call mpi_allreduce(timings_1a, work, iwork, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    call dcopy(iwork, work, 1, timings_1a, 1)
    deallocate (work, stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] error deallocating work'
    ! timings_1b
    iwork = ne
    allocate (work(ne,1), stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] error allocating work'
    call mpi_allreduce(timings_1b, work(:,1), iwork, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
    call dcopy(iwork, work(:,1), 1, timings_1b, 1)
    deallocate (work, stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] error deallocating work'


    ! first find average over atoms of timings of 1a
    allocate (t_average(ne), stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] Error allocating t_average'

    do ie = 1, ne
      t_average(ie) = 0.0d0
      do iat = 1, nat
        t_average(ie) = t_average(ie) + timings_1a(ie, iat)/dfloat(nat)
      end do
      ! if(myrank==master) write(1337,'(A,i9,2ES23.16)') '[check_communication_pattern]: ie, time for 1a and 1b', ie, timings_1b(ie),t_average(ie)
      ! if(myrank==master .and. t_inc%i_write) write(1337,'(A,i9,2ES23.16)') '[check_communication_pattern]: ie, time for 1a and 1b', ie, timings_1b(ie),t_average(ie)
    end do


    ! find how many different energy points have the same kmesh
    allocate (kmesh_n(nkmesh), stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] Error allocating kmesh_n'
    kmesh_n(:) = 0
    ik = 1
    do ie = 1, ne
      ik = kmesh_ie(ie)
      kmesh_n(ik) = kmesh_n(ik) + 1
    end do

    ! average timings of energypoints in different kmesh
    load_imbalance(:) = 0
    do ie = 1, ne
      ik = kmesh_ie(ie)            ! multiply with large number to have nice distiguishable integers
      load_imbalance(ik) = load_imbalance(ik) + (int(10000.0d0*((t_average(ie)+timings_1b(ie))/t_average(ie))))/kmesh_n(ik)
    end do
    if (myrank==master) write (*, *) 'load_imbalance', load_imbalance
    ! if(myrank==master .and. t_inc%i_write) write(1337,'(A,i9,1000I9)') '[check_communication_pattern] load imbalance:', load_imbalance


    ! set MPIatom and MPIadapt accordingly
    ! ...rest_at, rest_e => which fits actual load imbalance better => set MPIatom and MPIadapt
    if (myrank==master) write (*, *) mpiatom, mpiadapt



    ! finally deallocate work arrays
    deallocate (t_average, kmesh_n, stat=ierr)
    if (ierr/=0) stop '[check_communication_pattern] Error deallocating work arrays'

  end subroutine check_communication_pattern
#endif


  !! wrapper for work distribution over energies, also saves parallelization
  !! information for later use in t_mpi_c_grid
  subroutine distribute_work_energies(n_work, distribute_rest)

    use :: mod_types, only: t_mpi_c_grid
    use :: mod_profiling, only: memocc

    implicit none

    !! number of energies to be distributed
    integer, intent (in) :: n_work
    !! decide wether or not to distribute the rest rank or leave them idle
    logical, optional, intent (in) :: distribute_rest
    ! locals
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt
    integer :: i_stat

#ifdef CPP_MPI
    if (present(distribute_rest)) then
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, n_work, ntot_pt, ioff_pt, .true., distribute_rest)
    else
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_at, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, n_work, ntot_pt, ioff_pt, .true., .true.)
    end if

    ! store in t_mpi_c_grid for later use
    if (.not. (allocated(t_mpi_c_grid%ntot_pt2) .or. allocated(t_mpi_c_grid%ioff_pt2))) then
      allocate (t_mpi_c_grid%ntot_pt2(0:t_mpi_c_grid%nranks_at-1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt2))*kind(t_mpi_c_grid%ntot_pt2), 't_mpi_c_grid%ntot_pT2', 'distribute_work_energies')
      allocate (t_mpi_c_grid%ioff_pt2(0:t_mpi_c_grid%nranks_at-1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt2))*kind(t_mpi_c_grid%ioff_pt2), 't_mpi_c_grid%ioff_pT2', 'distribute_work_energies')
    end if
    t_mpi_c_grid%ntot_pt2 = ntot_pt
    t_mpi_c_grid%ioff_pt2 = ioff_pt
    t_mpi_c_grid%ntot2 = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    if (.not. (allocated(t_mpi_c_grid%ntot_pt2) .or. allocated(t_mpi_c_grid%ioff_pt2))) then
      allocate (t_mpi_c_grid%ntot_pt2(1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt2))*kind(t_mpi_c_grid%ntot_pt2), 't_mpi_c_grid%ntot_pT2', 'distribute_work_energies')
      allocate (t_mpi_c_grid%ioff_pt2(1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt2))*kind(t_mpi_c_grid%ioff_pt2), 't_mpi_c_grid%ioff_pT2', 'distribute_work_energies')
    end if
    t_mpi_c_grid%ntot2 = n_work
    t_mpi_c_grid%ntot_pt2 = n_work
    t_mpi_c_grid%ioff_pt2 = 0
#endif

  end subroutine distribute_work_energies


  !! wrapper for work distribution over atoms, also saves parallelization
  !! information for later use in t_mpi_c_grid
  subroutine distribute_work_atoms(n_work, i1_start, i1_end, distribute_rest)

    use :: mod_types, only: t_mpi_c_grid
    use :: mod_profiling, only: memocc

    implicit none

    !! number of energies to be distributed
    integer, intent (in) :: n_work
    !! decide wether or not to distribute the rest rank or leave them idle
    logical, optional, intent (in) :: distribute_rest
    !! number of energies to be distributed
    integer, intent (out) :: i1_start
    !! number of energies to be distributed
    integer, intent (out) :: i1_end
    ! locals
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt
    integer :: i_stat

#ifdef CPP_MPI
    if (present(distribute_rest)) then
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, n_work, ntot_pt, ioff_pt, .true., distribute_rest)
    else
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, n_work, ntot_pt, ioff_pt, .true., .true.)
    end if

    i1_start = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
    i1_end = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)

    ! store in t_mpi_c_grid for later use
    t_mpi_c_grid%ntot1 = ntot_pt(t_mpi_c_grid%myrank_ie)

    if (.not. (allocated(t_mpi_c_grid%ntot_pt1) .and. allocated(t_mpi_c_grid%ioff_pt1))) then
      allocate (t_mpi_c_grid%ntot_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt1))*kind(t_mpi_c_grid%ntot_pt1), 't_mpi_c_grid%ntot_pT1', 'distribute_work_atoms')
      allocate (t_mpi_c_grid%ioff_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt1))*kind(t_mpi_c_grid%ioff_pt1), 't_mpi_c_grid%ioff_pT1', 'distribute_work_atoms')
    end if
    t_mpi_c_grid%ntot_pt1 = ntot_pt
    t_mpi_c_grid%ioff_pt1 = ioff_pt
#else
    if (.not. (allocated(t_mpi_c_grid%ntot_pt1) .or. allocated(t_mpi_c_grid%ioff_pt1))) then
      allocate (t_mpi_c_grid%ntot_pt1(1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt1))*kind(t_mpi_c_grid%ntot_pt1), 't_mpi_c_grid%ntot_pT1', 'distribute_work_atoms')
      allocate (t_mpi_c_grid%ioff_pt1(1), stat=i_stat)
      call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt1))*kind(t_mpi_c_grid%ioff_pt1), 't_mpi_c_grid%ioff_pT1', 'distribute_work_atoms')
    end if
    t_mpi_c_grid%ntot1 = n_work
    t_mpi_c_grid%ntot_pt1 = n_work
    t_mpi_c_grid%ioff_pt1 = 0

    i1_start = 1
    i1_end = n_work
#endif

  end subroutine distribute_work_atoms


end module mod_mympi
