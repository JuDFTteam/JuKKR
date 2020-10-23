!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper for the definition of the MPI helper functions
!> Author: 
!> Wrapper for the definition of the MPI helper functions
!------------------------------------------------------------------------------------
!> @note ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann
!> @endnote
!------------------------------------------------------------------------------------
module mod_mympi

  implicit none

  private
  public :: myrank, nranks, master, mympi_init, mpiatom, mpiadapt, distribute_work_atoms, distribute_work_energies
#ifdef CPP_MPI
  public :: distribute_linear_on_tasks, find_dims_2d, create_newcomms_group_ie, mympi_main1c_comm, mympi_main1c_comm_newsosol, mympi_main1c_comm_newsosol2, &
    check_communication_pattern, bcast_global_variables
#endif

  integer, save :: myrank   = -1
  integer, save :: nranks   = -1
  integer, save :: master   = -1
  logical, save :: mpiatom  = .false.
  integer, save :: mpiadapt = -1

contains

  !-------------------------------------------------------------------------------
  !> Summary: Initialization of the MPI process
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Initialization of the MPI process
  !-------------------------------------------------------------------------------
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
  !-------------------------------------------------------------------------------
  !> Summary: Distributes 'ntot' points on 'nranks' tasks and returns the number of points per task, `ntot_pT`, and the offsets of the tasks, `ioff_pT`.
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Distributes 'ntot' points on 'nranks' tasks and returns the number of points 
  !> per task, `ntot_pT`, and the offsets of the tasks, `ioff_pT`.
  !-------------------------------------------------------------------------------
  subroutine distribute_linear_on_tasks(nranks,myrank,master,ntot,ntot_pt,ioff_pt,  &
    output,fill_rest)

    implicit none

    integer, intent (in) :: nranks, myrank, master, ntot
    logical, intent (in) :: output
    logical, intent (in), optional :: fill_rest
    integer, dimension(0:nranks-1), intent (out) :: ntot_pt, ioff_pt

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
  !-------------------------------------------------------------------------------
  !> Summary: Find dimensions to create cartesian communicator
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Find dimensions to create cartesian communicator
  !-------------------------------------------------------------------------------
  subroutine find_dims_2d(nranks, ntot1, ntot2, dims, mpiatom)
    ! find dimensions to create cartesian communicator
    ! input:  nranks, ntot1 is N_atom, ntot2 is N_E
    ! output: dims(2), dims(1) is N_atomranks, dims(2) is N_Eranks
    use :: mpi

    implicit none
    integer, intent (in) :: nranks, ntot1, ntot2
    logical, intent (in) :: mpiatom
    integer, dimension(2), intent (out) :: dims
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
  !-------------------------------------------------------------------------------
  !> Summary: Takes vector `kmesh` with mesh/timing information and finds number of rest procs that are devided in fractions given in ktake for optimal division of work
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Takes vector `kmesh` with mesh/timing information and finds number of rest 
  !> procs that are devided in fractions given in ktake for optimal division of work
  !-------------------------------------------------------------------------------
  subroutine create_newcomms_group_ie(master,nranks,myrank,nat,ne,nkmesh,kmesh,     &
    mympi_comm_ie,myrank_ie,nranks_ie,mympi_comm_at,myrank_at,nranks_at,            &
    myrank_atcomm,nranks_atcomm)

    use :: mpi
    use :: mod_datatypes, only: dp
    implicit none

    integer, intent (in) :: master, nranks, myrank, ne, nat, nkmesh
    integer, dimension(nkmesh), intent (in) :: kmesh
    integer, intent (out) :: mympi_comm_ie, mympi_comm_at, nranks_atcomm, nranks_at
    integer, intent (out) :: nranks_ie, myrank_ie, myrank_at, myrank_atcomm

    integer :: rest,myg, ie, iat, ierr, ik, i1, i2, i3
    integer, dimension(nkmesh-1)  :: k
    integer, dimension(nkmesh)    :: ktake
    real (kind=dp) :: q, qmin
    real (kind=dp), dimension(nkmesh-1) :: f
    integer, dimension(:), allocatable :: mygroup
    integer, dimension(:,:), allocatable :: groups
    integer :: mympi_group_ie, mympi_group_world, mympi_comm_grid

    rest = 0
    qmin = -1
    if (myrank==master) write (1337, *) 'create_newcomms_group_ie input:', nranks, ne, nat

    ktake(:) = 0

    ! now check different cases:
    ! 1: mor than one energy group with rest ranks (tackle load imbalance)
    ! 2: no rest ranks (use cartesian grid)
    ! 3: only atom parallelization without rest ranks
    ! 4: only atom paralelization with rest ranks (discard rest ranks and use cartesian grid)
    ! 5: only energy parallelization without rest ranks

    if ((ne*nat)<nranks .and. (ne>1)) then ! .and. nat>1)) then
      ! find parallelization with rest ranks divided over energy group with
      ! denset k-mesh

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
                  ktake(1:3) = [ i1, i2, i3 ]
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
                ktake(1:2) = [ i1, i2 ]
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
      ! no rest ranks, use cartesian grid

      rest = 0

      if (myrank==master) write (1337, *) 'create cartesian grid:', ne, nat, nranks
      call mpi_cart_create(mpi_comm_world, 2, [ne,nat], [.false.,.false.], .true., mympi_comm_grid, ierr)

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
      ! no rest ranks, can use mpi_comm world entirely for atom parallelization

      rest = 0

      mympi_comm_ie = mpi_comm_world
      myrank_ie = myrank
      nranks_ie = nranks
      mympi_comm_at = mpi_comm_self
      myrank_at = 0
      nranks_at = 1

      myrank_atcomm = myrank_at
      nranks_atcomm = nranks_at

    else if (nat<nranks .and. (ne==1 .and. nat>1)) then
      ! technically have rest ranks but do not use them since ne==1 (wasts some resources!)
      if (myrank==master) then
        write (*, '(A,I5,A)') 'WARNING: nat<nranks but ne==1: ', nranks-nat, ' ranks are unused!'
        write (*, '(A)') 'Please consider modifying your jobsript to not waste any resources!'
        write (1337, '(A,I5,A)') 'WARNING: nat<nranks but ne==1: ', nranks-nat, ' ranks are unused!'
        write (1337, '(A)') 'Please consider modifying your jobsript to not waste any resources!'
      end if

      ! set rest artificially to zero (affects printout only)
      rest = 0
      if (myrank<nat) then
        call mpi_comm_split(mpi_comm_world, 1, myrank, mympi_comm_ie, ierr)
      else
        call mpi_comm_split(mpi_comm_world, myrank+1000, 0, mympi_comm_ie, ierr)
      end if

      call mpi_comm_rank(mympi_comm_ie, myrank_ie, ierr)
      call mpi_comm_size(mympi_comm_ie, nranks_ie, ierr)
      if (myrank>=nat) then
        myrank_ie = -1 ! overwrite this so that myrank==master checks return false for unused ranks
      end if

      ! set atom communicator manually to mpi_comm_self
      mympi_comm_at = mpi_comm_self
      nranks_at = 1
      myrank_at = 0
      ! fill atcomm helpes (only used for rest>0)
      nranks_atcomm = nranks_at
      myrank_atcomm = myrank_at

    else
      ! fallback: use all ranks for energy parallelization without rest ranks

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

    ! print output
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
  !-------------------------------------------------------------------------------
  !> Summary: MPI communication for the `main1c` subroutine
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> MPI communication for the `main1c` subroutine
  !-------------------------------------------------------------------------------
  subroutine mympi_main1c_comm(irmd,lmpotd,natypd,lmaxd,lmaxd1,lmmaxd,npotd,ielast, &
    mmaxd,idoldau,natyp,krel,lmomvec,nmvecmax,nqdos,rho2ns,r2nef,espv,den,denlm,    &
    denmatc,denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,mympi_comm)

    use :: mpi
    use :: mod_datatypes, only: dp
    implicit none

    integer, intent (in) :: krel !! Switch for non- (or scalar-) relativistic/relativistic (Dirac) program (0/1). Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
    integer, intent (in) :: irmd    !! Maximum number of radial points
    integer, intent (in) :: npotd   !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: nqdos
    integer, intent (in) :: natyp   !! Number of kinds of atoms in unit cell
    integer, intent (in) :: mmaxd
    integer, intent (in) :: lmaxd   !! Maximum l component in wave function expansion
    integer, intent (in) :: lmpotd  !! (lpot+1)**2
    integer, intent (in) :: natypd  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: lmmaxd  !! (KREL+KORBIT+1)*(LMAX+1)**2
    integer, intent (in) :: ielast
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: idoldau !! flag to perform LDA+U
    logical, intent (in) :: lmomvec
    integer, intent (in) :: nmvecmax
    integer, intent (in) :: mympi_comm
    ! .. In/Out variables
    real (kind=dp), intent (inout) :: denef
    real (kind=dp), dimension(natypd), intent (inout) :: denefat
    real (kind=dp), dimension(0:lmaxd1, npotd), intent (inout)            :: espv
    real (kind=dp), dimension(irmd*krel+(1-krel), natypd), intent (inout) :: rhoorb !! orbital density
    real (kind=dp), dimension(0:lmaxd1+1, 3, natypd), intent (inout)  :: muorb    !! orbital magnetic moment
    real (kind=dp), dimension(irmd, lmpotd, natypd, 2), intent (inout) :: r2nef   !! rho at FERMI energy
    real (kind=dp), dimension(irmd, lmpotd, natypd, 2), intent (inout) :: rho2ns  !! radial density
    complex (kind=dp), dimension(natypd, 3, nmvecmax), intent (inout) :: mvevi
    complex (kind=dp), dimension(natypd, 3, nmvecmax), intent (inout) :: mvevief
    complex (kind=dp), dimension(mmaxd, mmaxd, npotd), intent (inout) :: denmatc
    complex (kind=dp), dimension(0:lmaxd1, ielast, npotd, nqdos), intent (inout)  :: den
    complex (kind=dp), dimension(lmmaxd, ielast, npotd, nqdos), intent (inout)    :: denlm
    complex (kind=dp), dimension(0:lmaxd, natypd, 3, nmvecmax), intent (inout)    :: mvevil
    ! .. Local variables
    integer :: idim, ierr          ! , myrank_comm
    integer, parameter :: master = 0
    real (kind=dp), dimension(:), allocatable :: work1
    real (kind=dp), dimension(:, :), allocatable :: work2
    real (kind=dp), dimension(:, :, :), allocatable :: work3
    real (kind=dp), dimension(:, :, :, :), allocatable :: work4
    complex (kind=dp), dimension(:,:,:), allocatable :: work3c
    complex (kind=dp), dimension(:,:,:,:), allocatable :: work4c

    allocate (work4(irmd,lmpotd,natypd,2), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4 = 0.0_dp
    idim = irmd*lmpotd*natypd*2
    call mpi_allreduce(rho2ns, work4, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work4, 1, rho2ns, 1)
    deallocate (work4)

    allocate (work4(irmd,lmpotd,natypd,2), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4 = 0.0_dp
    idim = irmd*lmpotd*natypd*2
    call mpi_allreduce(r2nef, work4, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work4, 1, r2nef, 1)
    deallocate (work4)

    ! ESPV needs integration over atoms and energies -> MPI_COMM_WORLD
    allocate (work2(0:lmaxd+1,npotd), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work2 = 0.0_dp
    idim = (lmaxd+2)*npotd
    call mpi_allreduce(espv, work2, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work2, 1, espv, 1)
    deallocate (work2)

    allocate (work4c(0:lmaxd+1,ielast,npotd,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4c = (0.0_dp, 0.0_dp)
    idim = ielast*(lmaxd+2)*npotd*nqdos
    call mpi_allreduce(den, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(idim, work4c, 1, den, 1)
    deallocate (work4c)

    allocate (work4c(ielast,lmmaxd,nqdos,npotd), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work4c = (0.0_dp, 0.0_dp)
    idim = ielast*(lmmaxd)*npotd*nqdos
    call mpi_allreduce(denlm, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
    call zcopy(idim, work4c, 1, denlm, 1)
    deallocate (work4c)

    if (idoldau==1) then
      allocate (work3c(mmaxd,mmaxd,npotd), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work3c = (0.0_dp, 0.0_dp)
      idim = mmaxd*mmaxd*npotd
      call mpi_allreduce(denmatc, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
      call zcopy(idim, work3c, 1, denmatc, 1)
      deallocate (work3c)
    end if

    allocate (work1(1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work1 = 0.0_dp
    idim = 1
    call mpi_allreduce(denef, work1, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work1, 1, denef, 1)
    deallocate (work1)

    allocate (work1(natyp), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
    work1 = 0.0_dp
    idim = natyp
    call mpi_allreduce(denefat, work1, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
    call dcopy(idim, work1, 1, denefat, 1)
    deallocate (work1)

    if (krel==1) then
      allocate (work2(irmd,natypd), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work2 = 0.0_dp
      idim = irmd*natypd
      call mpi_allreduce(rhoorb, work2, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
      call dcopy(idim, work2, 1, rhoorb, 1)
      deallocate (work2)

      allocate (work3(0:lmaxd+2,natypd,3), stat=ierr)
      if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
      work3 = 0.0_dp
      idim = (lmaxd+3)*natypd*3
      call mpi_allreduce(muorb, work3, idim, mpi_double_precision, mpi_sum, mympi_comm, ierr)
      call dcopy(idim, work3, 1, muorb, 1)
      deallocate (work3)

      if (lmomvec) then
        allocate (work3c(natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work3c = (0.0_dp, 0.0_dp)
        idim = natypd*3*nmvecmax
        call mpi_allreduce(mvevi, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work3c, 1, mvevi, 1)
        deallocate (work3c)

        allocate (work4c(lmaxd+1,natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work4c = (0.0_dp, 0.0_dp)
        idim = (lmaxd+1)*natypd*3*nmvecmax
        call mpi_allreduce(mvevil, work4c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work4c, 1, mvevil, 1)
        deallocate (work4c)

        allocate (work3c(natypd,3,nmvecmax), stat=ierr)
        if (ierr/=0) stop '[mympi_main1c_comm] error allocating work array'
        work3c = (0.0_dp, 0.0_dp)
        idim = natypd*3*nmvecmax
        call mpi_allreduce(mvevief, work3c, idim, mpi_double_complex, mpi_sum, mympi_comm, ierr)
        call zcopy(idim, work3c, 1, mvevief, 1)
        deallocate (work3c)
      end if                       ! LMOMVEC
    end if                         ! KREL.EQ.1

  end subroutine mympi_main1c_comm
#endif

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: MPI communication for the new solver in the `main1c` subroutine 
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> MPI communication for the new solver in the `main1c` subroutine 
  !-------------------------------------------------------------------------------
  subroutine mympi_main1c_comm_newsosol(nspin, korbit, irmdnew, lmpotd, lmaxd, lmaxd1, lmmax0d, lmmaxd, ielast, nqdos, den, denlm, gflle, rho2nsc, r2nefc, rho2int, espv, muorb, denorbmom, &
    denorbmomsp, denorbmomlm, denorbmomns, mympi_comm)

    use :: mpi
    use :: mod_datatypes, only: dp
    implicit none
    ! .. Input variables
    integer, intent (in) :: nqdos
    integer, intent (in) :: lmaxd   !! Maximum l component in wave function expansion
    integer, intent (in) :: lmpotd  !! (lpot+1)**2
    integer, intent (in) :: lmmax0d  !! (LMAX+1)**2
    integer, intent (in) :: ielast
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: irmdnew
    integer, intent (in) :: lmmaxd !! (krel+korbit+1)*(LMAX+1)**2
    integer, intent (in) :: nspin
    integer, intent (in) :: korbit
    integer, intent (in) :: mympi_comm
    ! .. In/Out variables
    complex (kind=dp), dimension(nspin*(1+korbit)), intent (inout) :: rho2int
    complex (kind=dp), dimension(irmdnew, lmpotd, nspin*(1+korbit)), intent (inout) :: r2nefc
    complex (kind=dp), dimension(irmdnew, lmpotd, nspin*(1+korbit)), intent (inout) :: rho2nsc
    complex (kind=dp), dimension(0:lmaxd1, ielast, nqdos, nspin), intent (inout) :: den
    complex (kind=dp), dimension(lmmax0d, ielast, nqdos, nspin), intent (inout) :: denlm
    complex (kind=dp), dimension(lmmaxd, lmmaxd, ielast, nqdos), intent (inout) :: gflle
    real (kind=dp), dimension(3), intent (inout) :: denorbmom
    real (kind=dp), dimension(3), intent (inout) :: denorbmomns
    real (kind=dp), dimension(2, 3), intent (inout) :: denorbmomsp
    real (kind=dp), dimension(0:lmaxd1, 2), intent (inout) :: espv
    real (kind=dp), dimension(0:lmaxd1+1, 3), intent (inout) :: muorb !! orbital magnetic moment
    real (kind=dp), dimension(0:lmaxd, 3), intent (inout) :: denorbmomlm
    
    ! .. Local variables
    integer :: ierr, idim
    real (kind=dp), dimension(:, :, :, :), allocatable :: work
    complex (kind=dp), dimension(:, :, :, :), allocatable :: workc

    ! all with reduce instead of allreduce:
    ! complex (kind=dp) arrays
    idim = irmdnew*lmpotd*nspin*(1+korbit)
    allocate (workc(irmdnew,lmpotd,nspin*(1+korbit),1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, r2nefc'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(r2nefc, workc(:,:,:,1), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for r2nefc'
    call zcopy(idim, workc, 1, r2nefc, 1)
    deallocate (workc)

    idim = irmdnew*lmpotd*nspin*(1+korbit)
    allocate (workc(irmdnew,lmpotd,nspin*(1+korbit),1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, rho2nsc'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(rho2nsc, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for rho2nsc'
    call zcopy(idim, workc, 1, rho2nsc, 1)
    deallocate (workc)

    idim = (lmaxd1+1)*ielast*nspin*nqdos
    allocate (workc(0:lmaxd1,ielast,nspin,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, den'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(den, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for den'
    call zcopy(idim, workc, 1, den, 1)
    deallocate (workc)

    idim = lmmax0d*ielast*nspin*nqdos
    allocate (workc(lmmax0d,ielast,nspin,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, denlm'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(denlm, workc, idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denlm'
    call zcopy(idim, workc, 1, denlm, 1)
    deallocate (workc)

    idim = nspin*(1+korbit)
    allocate (workc(nspin*(1+korbit),1,1,1), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, rho2int'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(rho2int, workc(:,1,1,1), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for rho2int'
    call zcopy(idim, workc, 1, rho2int, 1)
    deallocate (workc)

    idim = lmmaxd*lmmaxd*ielast*nqdos
    allocate (workc(lmmaxd,lmmaxd,ielast,nqdos), stat=ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error allocating workc, gflle'
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(gflle, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for gflle'
    call zcopy(idim, workc, 1, gflle, 1)
    deallocate (workc)

    ! real (kind=dp) arrays
    idim = (lmaxd1+1)*2
    allocate (work(0:lmaxd1,2,1,1))
    work = 0.0_dp
    call mpi_reduce(espv, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for espv'
    call dcopy(idim, work, 1, espv, 1)
    deallocate (work)

    idim = (lmaxd1+2)*3
    allocate (work(0:lmaxd1+1,3,1,1))
    work = 0.0_dp
    call mpi_reduce(muorb, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for muorb'
    call dcopy(idim, work, 1, muorb, 1)
    deallocate (work)

    idim = 3
    allocate (work(3,1,1,1))
    work = 0.0_dp
    call mpi_reduce(denorbmom, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denobrmom'
    call dcopy(idim, work, 1, denorbmom, 1)
    deallocate (work)

    idim = 2*3
    allocate (work(2,3,1,1))
    work = 0.0_dp
    call mpi_reduce(denorbmomsp, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomsp'
    call dcopy(idim, work, 1, denorbmomsp, 1)
    deallocate (work)

    idim = 3
    allocate (work(3,1,1,1))
    work = 0.0_dp
    call mpi_reduce(denorbmomns, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomns'
    call dcopy(idim, work, 1, denorbmomns, 1)
    deallocate (work)

    idim = (lmaxd+1)*3
    allocate (work(0:lmaxd1,3,1,1))
    work = 0.0_dp
    call mpi_reduce(denorbmomlm, work, idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol] Error in MPI_REDUCE for denorbmomlm'
    call dcopy(idim, work, 1, denorbmomlm, 1)
    deallocate (work)

  end subroutine mympi_main1c_comm_newsosol
#endif

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: MPI communication for the new solver in the `main1c` subroutine 
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> MPI communication for the new solver in the `main1c` subroutine 
  !-------------------------------------------------------------------------------
  subroutine mympi_main1c_comm_newsosol2(lmaxd1,lmmaxd,ielast,nqdos,npotd,natypd,   &
    lmpotd,irmd,mmaxd,den,denlm,muorb,espv,r2nef,rho2ns,denefat,denef,denmatn,      &
    angles_new,totmoment,mympi_comm)

    use :: mpi
    use :: mod_datatypes, only: dp
    implicit none

    integer, intent (in) :: irmd    !! Maximum number of radial points
    integer, intent (in) :: npotd   !! (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
    integer, intent (in) :: nqdos
    integer, intent (in) :: mmaxd
    integer, intent (in) :: lmpotd  !! (lpot+1)**2
    integer, intent (in) :: natypd  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: lmmaxd  !! (KREL+KORBIT+1)*(LMAX+1)**2
    integer, intent (in) :: ielast
    integer, intent (in) :: lmaxd1
    integer, intent (in) :: mympi_comm
    ! .. In/Out variables
    real (kind=dp), intent (inout) :: denef
    complex (kind=dp), dimension(0:lmaxd1, ielast, nqdos, npotd), intent (inout) :: den
    complex (kind=dp), dimension(lmmaxd, ielast, nqdos, npotd), intent (inout) :: denlm
    complex (kind=dp), dimension(mmaxd, mmaxd, 2, 2, natypd), intent (inout) :: denmatn
    real (kind=dp), dimension(natypd), intent (inout)     :: denefat
    real (kind=dp), dimension(natypd, 2), intent (inout)  :: angles_new
    real (kind=dp), dimension(natypd), intent (inout)  :: totmoment
    real (kind=dp), dimension(0:lmaxd1, npotd), intent (inout) :: espv
    real (kind=dp), dimension(0:lmaxd1+1, 3, natypd), intent (inout) :: muorb !! orbital magnetic moment
    real (kind=dp), dimension(irmd, lmpotd, natypd, 2), intent (inout) :: r2nef
    real (kind=dp), dimension(irmd, lmpotd, natypd, 2), intent (inout) :: rho2ns !! radial density

    integer :: ierr, idim
    real (kind=dp), dimension(:,:,:,:), allocatable :: work
    complex (kind=dp), dimension(:,:,:,:), allocatable :: workc
    complex (kind=dp), dimension(:,:,:,:,:), allocatable :: workc1

    ! complex (kind=dp) arrays
    idim = (1+lmaxd1)*ielast*nqdos*npotd
    allocate (workc(0:lmaxd1,ielast,nqdos,npotd))
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(den, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(den,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for den'
    call zcopy(idim, workc, 1, den, 1)
    deallocate (workc)

    idim = lmmaxd*ielast*nqdos*npotd
    allocate (workc(lmmaxd,ielast,nqdos,npotd))
    workc = (0.0_dp, 0.0_dp)
    call mpi_reduce(denlm, workc(:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denlm,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denlm'
    call zcopy(idim, workc, 1, denlm, 1)
    deallocate (workc)

    idim = mmaxd*mmaxd*npotd
    allocate (workc1(mmaxd,mmaxd,2,2,natypd))
    workc1 = (0.0_dp, 0.0_dp)
    call mpi_reduce(denmatn, workc1(:,:,:,:,:), idim, mpi_double_complex, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denmatn,workc1(:,:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denmatn'
    call zcopy(idim, workc1, 1, denmatn, 1)
    deallocate (workc1)

    ! real (kind=dp) arrays
    idim = (lmaxd1+2)*3*natypd
    allocate (work(0:lmaxd1+1,3,natypd,1))
    work = 0.0_dp
    call mpi_reduce(muorb, work(:,:,:,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(muorb,work(:,:,:,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for muorb'
    call dcopy(idim, work, 1, muorb, 1)
    deallocate (work)

    idim = (lmaxd1+1)*npotd
    allocate (work(0:lmaxd1,npotd,1,1))
    work = 0.0_dp
    call mpi_reduce(espv, work(:,:,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(ESPV,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for espv'
    call dcopy(idim, work, 1, espv, 1)
    deallocate (work)

    idim = irmd*lmpotd*natypd*2
    allocate (work(irmd,lmpotd,natypd,2))
    work = 0.0_dp
    call mpi_reduce(r2nef, work(:,:,:,:), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(r2nef,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for r2nef'
    call dcopy(idim, work, 1, r2nef, 1)
    deallocate (work)

    idim = irmd*lmpotd*natypd*2
    allocate (work(irmd,lmpotd,natypd,2))
    work = 0.0_dp
    call mpi_reduce(rho2ns, work(:,:,:,:), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(RHO2NS,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for rho2ns'
    call dcopy(idim, work, 1, rho2ns, 1)
    deallocate (work)

    idim = natypd
    allocate (work(natypd,1,1,1))
    work = 0.0_dp
    call mpi_reduce(denefat, work(:,1,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denefat,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denefat'
    call dcopy(idim, work, 1, denefat, 1)
    deallocate (work)

    idim = 1
    allocate (work(1,1,1,1))
    work = 0.0_dp
    call mpi_reduce(denef, work(:,1,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(denef,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for denef'
    call dcopy(idim, work, 1, denef, 1)
    deallocate (work)

    idim = 2*natypd
    allocate (work(2,natypd,1,1))
    work = 0.0_dp
    call mpi_reduce(angles_new, work(:,:,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    ! CALL MPI_ALLREDUCE(angles_new,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for angles_new'
    call dcopy(idim, work, 1, angles_new, 1)
    deallocate (work)

    idim = natypd
    allocate (work(natypd,1,1,1))
    work = 0.0_dp
    call mpi_reduce(totmoment, work(:,1,1,1), idim, mpi_double_precision, mpi_sum, master, mympi_comm, ierr)
    if (ierr/=0) stop '[mympi_main1c_comm_newsosol2] Error in MPI_REDUCE for totmoment'
    call dcopy(idim, work, 1, totmoment, 1)
    deallocate (work)

  end subroutine mympi_main1c_comm_newsosol2
#endif


#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: Subroutine to check if the MPI communication is working properly
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Subroutine to check if the MPI communication is working properly
  !-------------------------------------------------------------------------------
  subroutine check_communication_pattern(mpiatom,mpiadapt,timings_1a,timings_1b,    &
    load_imbalance,nkmesh,kmesh_ie)

    use :: mpi
    use :: mod_datatypes, only: dp

    implicit none

    integer, intent (inout) :: mpiadapt
    logical, intent (inout) :: mpiatom
    real (kind=dp), dimension(:), intent (inout) :: timings_1b
    real (kind=dp), dimension(:,:), intent (inout) :: timings_1a
    integer, dimension(:), intent (out) :: load_imbalance
    integer, intent (in) :: nkmesh
    integer, dimension(:), intent (in) :: kmesh_ie

    real (kind=dp), dimension(:), allocatable :: t_average
    real (kind=dp), dimension(:,:), allocatable :: work
    integer, dimension(:), allocatable :: kmesh_n
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
      t_average(ie) = 0.0_dp
      do iat = 1, nat
        t_average(ie) = t_average(ie) + timings_1a(ie, iat)/real(nat, kind=dp)
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
      load_imbalance(ik) = load_imbalance(ik) + (int(10000.0_dp*((t_average(ie)+timings_1b(ie))/t_average(ie))))/kmesh_n(ik)
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

  !-------------------------------------------------------------------------------
  !> Summary: Wrapper for work distribution over energies, also saves parallelization information for later use in `t_mpi_c_grid`
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> Wrapper for work distribution over energies, also saves parallelization
  !> information for later use in `t_mpi_c_grid`
  !-------------------------------------------------------------------------------
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


  !-------------------------------------------------------------------------------
  !> Summary: wrapper for work distribution over atoms, also saves parallelization information for later use in `t_mpi_c_grid`
  !> Author: 
  !> Category: communication, KKRhost 
  !> Deprecated: False
  !> wrapper for work distribution over atoms, also saves parallelization
  !> information for later use in `t_mpi_c_grid`
  !-------------------------------------------------------------------------------
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

#ifdef CPP_MPI
  !-------------------------------------------------------------------------------
  !> Summary: MPI Briadcast of global variables
  !> Author: Jonathan Chico, Philipp Rüßmann
  !> Category: KKRhost, communication, initialization
  !> Deprecated: False ! This needs to be set to True for deprecated subroutines
  !>
  !> MPI broadcast routine for global variables (i.e. array dimensions etc.)
  !-------------------------------------------------------------------------------
  subroutine bcast_global_variables()
    use :: mpi
    use :: global_variables ! use all parameters
    implicit none
  
    integer :: n !! number of paramters that are broadcasted
    integer, allocatable :: blocklen1(:), etype1(:) !! blocklength of variuables in derived data type and list of MPI datatypes
    integer :: ierr !! error status
    integer :: mympitype1 !! derived data type for collective communication
    integer (kind=mpi_address_kind), allocatable :: disp1(:) !! MPI addresses
    integer (kind=mpi_address_kind) :: base !! base address of first entry
  
  
    n = 65
    allocate (blocklen1(n), etype1(n), disp1(n), stat=ierr)
    if (ierr/=0) stop 'error allocating arrays in bcast_global_variables'
  
    call mpi_get_address(n, disp1(1), ierr)
    call mpi_get_address(irid, disp1(2), ierr)
    call mpi_get_address(krel, disp1(3), ierr)
    call mpi_get_address(nfund, disp1(4), ierr)
    call mpi_get_address(ipand, disp1(5), ierr)
    call mpi_get_address(ngshd, disp1(6), ierr)
    call mpi_get_address(ncleb, disp1(7), ierr)
    call mpi_get_address(knoco, disp1(8), ierr)
    call mpi_get_address(iemxd, disp1(9), ierr)
    call mpi_get_address(irnsd, disp1(10), ierr)
    call mpi_get_address(nmaxd, disp1(11), ierr)
    call mpi_get_address(ishld, disp1(12), ierr)
    call mpi_get_address(naclsd, disp1(13), ierr)
    call mpi_get_address(nspotd, disp1(14), ierr)
    call mpi_get_address(ntperd, disp1(15), ierr)
    call mpi_get_address(ntrefd, disp1(16), ierr)
    call mpi_get_address(nsheld, disp1(17), ierr)
    call mpi_get_address(ncelld, disp1(18), ierr)
    call mpi_get_address(nspind, disp1(19), ierr)
    call mpi_get_address(knosph, disp1(20), ierr)
    call mpi_get_address(korbit, disp1(21), ierr)
    call mpi_get_address(kpoibz, disp1(22), ierr)
    call mpi_get_address(wlength, disp1(23), ierr)
    call mpi_get_address(nprincd, disp1(24), ierr)
    call mpi_get_address(nlayerd, disp1(25), ierr)
    call mpi_get_address(natomimpd, disp1(26), ierr)
    call mpi_get_address(lmaxd, disp1(27), ierr)
    call mpi_get_address(lmmaxd, disp1(28), ierr)
    call mpi_get_address(lmgf0d, disp1(29), ierr)
    call mpi_get_address(alm, disp1(30), ierr)
    call mpi_get_address(almgf0, disp1(31), ierr)
    call mpi_get_address(ndim_slabinv, disp1(32), ierr)
    call mpi_get_address(nembd, disp1(33), ierr)
    call mpi_get_address(nembd1, disp1(34), ierr)
    call mpi_get_address(nembd2, disp1(35), ierr)
    call mpi_get_address(nrd, disp1(36), ierr)
    call mpi_get_address(lm2d, disp1(37), ierr)
    call mpi_get_address(nclsd, disp1(38), ierr)
    call mpi_get_address(mmaxd, disp1(39), ierr)
    call mpi_get_address(npotd, disp1(40), ierr)
    call mpi_get_address(lmxspd, disp1(41), ierr)
    call mpi_get_address(lassld, disp1(42), ierr)
    call mpi_get_address(irmind, disp1(43), ierr)
    call mpi_get_address(nofgij, disp1(44), ierr)
    call mpi_get_address(nspindd, disp1(45), ierr)
    call mpi_get_address(nsatypd, disp1(46), ierr)
    call mpi_get_address(nrefd, disp1(47), ierr)
    call mpi_get_address(irmd, disp1(48), ierr)
    call mpi_get_address(naezd, disp1(49), ierr)
    call mpi_get_address(natypd, disp1(50), ierr)
    call mpi_get_address(lmpotd, disp1(51), ierr)
    call mpi_get_address(ntotd, disp1(52), ierr)
    call mpi_get_address(nrmaxd, disp1(53), ierr)
    call mpi_get_address(lpotd, disp1(54), ierr)
    call mpi_get_address(nchebd, disp1(55), ierr)
    call mpi_get_address(maxmshd, disp1(56), ierr)
    call mpi_get_address(kBdG, disp1(57), ierr)
    call mpi_get_address(ninit_broydenspin, disp1(58), ierr)
    call mpi_get_address(memlen_broydenspin, disp1(59), ierr)
    call mpi_get_address(nsimplemixfirst, disp1(60), ierr)
    call mpi_get_address(linterface, disp1(61), ierr)
    call mpi_get_address(lnc, disp1(62), ierr)
    call mpi_get_address(pot_ns_cutoff, disp1(63), ierr)
    call mpi_get_address(mixfac_broydenspin, disp1(64), ierr)
    call mpi_get_address(qbound_spin, disp1(65), ierr)
  
    ! find displacements of variables
    base = disp1(1)
    disp1 = disp1 - base
  
    ! set length of variables in derived data type
    blocklen1(1:n) = 1
  
    ! set datatype of variables
    etype1(1:n-5) = mpi_integer
    etype1(n-4:n-3) = mpi_logical
    etype1(n-3:n) = mpi_double_precision
  
    ! create new Type structure for derived data type
    call mpi_type_create_struct(n, blocklen1, disp1, etype1, mympitype1, ierr)
    if (ierr/=mpi_success) stop 'Problem in create_mpimask_t_inc'
  
    ! commit new type
    call mpi_type_commit(mympitype1, ierr)
    if (ierr/=mpi_success) stop 'error commiting create_mpimask_t_inc'
  
    ! broadcast derived data type
  
    call mpi_bcast(n, 1, mympitype1, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error brodcasting t_inc'
  
    ! finally free auxiliary type and deallocate working arrays
    call mpi_type_free(mympitype1, ierr)
    deallocate (blocklen1, etype1, disp1, stat=ierr)
    if (ierr/=0) stop 'error deallocating arrays in bcast_global_variables'

  end subroutine bcast_global_variables
#endif


end module mod_mympi
