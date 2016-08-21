module mod_mympi
!ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann

implicit none

  private
  public :: myrank, nranks, master, mympi_init, MPIatom
#ifdef CPP_MPI
  public :: distribute_linear_on_tasks, find_dims_2d, create_newcomms_group_ie, mympi_main1c_comm, mympi_main1c_comm_newsosol, mympi_main1c_comm_newsosol2
#endif

  integer, save :: myrank = -1
  integer, save :: nranks = -1
  integer, save :: master = -1
  logical, save :: MPIatom = .false.

contains

  subroutine mympi_init()

#ifdef CPP_MPI
    use mpi
#endif

    integer :: ierr

    master = 0

#ifdef CPP_MPI
    call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size ( MPI_COMM_WORLD, nranks, ierr )
#else
    myrank = master
    nranks = 1
#endif

  end subroutine mympi_init

#ifdef CPP_MPI
  subroutine distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pT, ioff_pT, output, fill_rest)
    !distributes 'ntot' points on 'nranks' tasks and returns the number of points per task, ntot_pT, and the offsets of the tasks, ioff_pT.

    implicit none

    integer, intent(in)  :: nranks, myrank, master, ntot
    logical, intent(in)  :: output
    logical, intent(in), optional :: fill_rest
    integer, intent(out) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)

    integer :: irest, irank

    ntot_pT = int(ntot/nranks)
    ioff_pT = int(ntot/nranks)*(/ (irank, irank=0,nranks-1) /)
    irest   = ntot-int(ntot/nranks)*nranks

    if(irest>0) then

      do irank=0,irest-1
        ntot_pT(irank) = ntot_pT(irank) + 1
        ioff_pT(irank) = ioff_pT(irank) + irank
      end do!irank

      do irank=irest,nranks-1
        ioff_pT(irank) = ioff_pT(irank) + irest
      end do!irank

    end if!irest>0
    
    if(present(fill_rest)) then
    if(fill_rest) then
    if(ntot<nranks) then
!        write(*,*) 'set rest',myrank,ntot,nranks
       do irank=ntot,nranks-1
          ioff_pT(irank) = ioff_pT(irank-1)
          ntot_pT(irank) = ntot_pT(irank-1)
       end do ! irank
    endif!ntot<nranks
    endif!fill_rest
    endif!present(fill_rest)

    if(myrank==master .and. output) then
      write(1337,*) '==== DISTRIBUTION OF POINTS ON TASKS: ===='
      do irank=0,nranks-1
        write(1337,'("Task ",I0," treats points ",I0," to ",I0,", #of points= ",I0)') irank, ioff_pT(irank)+1, ioff_pT(irank)+ntot_pT(irank), ntot_pT(irank)
      end do!irank
      write(1337,*) '=========================================='
    end if!myrank==master

  end subroutine distribute_linear_on_tasks
#endif


#ifdef CPP_MPI
  subroutine find_dims_2d(nranks,ntot1,ntot2,dims,MPIatom)
    !find dimensions to create cartesian communicator
    !input:  nranks, ntot1 is N_atom, ntot2 is N_E
    !output: dims(2), dims(1) is N_atomranks, dims(2) is N_Eranks
    use mpi

    implicit none
    integer, intent(in)  :: nranks,ntot1,ntot2
    logical, intent(in)  :: MPIatom
    integer, intent(out) :: dims(2)
    integer ierr
    
    if(.not. MPIatom) then
       if(nranks.le.ntot2) then
          dims(1) = 1
          dims(2) = nranks
       else
          dims(1) = nranks/ntot2
          dims(2) = ntot2
       end if
    else
       if(nranks.le.ntot1) then
          dims(2) = 1
          dims(1) = nranks
       else
          dims(2) = nranks/ntot1
          dims(1) = ntot1
       end if
    end if
    
    if(nranks>(ntot1*(ntot2+1)-1)) then
       if(myrank==master) write(*,'(A,I3,A,I3,A,I5)') 'Error for',ntot1,'atoms and',ntot2,'energy points you use too many processors. Nranks=',nranks
       call MPI_BARRIER(MPI_COMM_WORLD,ierr)
       call MPI_Finalize(ierr)
       stop 'Error: too many ranks'
    endif
        
  end subroutine find_dims_2d
#endif


#ifdef CPP_MPI
  subroutine create_newcomms_group_ie(nranks,myrank,nat,ne,nkmesh,kmesh,mympi_comm_ie,  &
  &                                  myrank_ie,nranks_ie,mympi_comm_at,myrank_at, nranks_at, myrank_atcomm, nranks_atcomm)
  !takes vector kmesh with mesh/timing information and finds number of rest procs that are devided in fractions given in ktake for optimal division of work
  
  use mpi
  implicit none
  
  integer, intent(in) :: nranks,myrank,ne,nat,nkmesh
  integer, intent(in) :: kmesh(nkmesh)
  integer, intent(out) :: mympi_comm_ie, mympi_comm_at, nranks_atcomm, nranks_at, nranks_ie, myrank_ie, myrank_at, myrank_atcomm
  
  integer :: rest, k(nkmesh-1), ktake(nkmesh-1), myg, ie, iat, ierr, ik, i1,i2,i3
  double precision :: f(nkmesh-1), q, qmin
  integer, allocatable :: groups(:,:), mygroup(:)
  integer :: mympi_group_ie, mympi_group_world, myMPI_comm_grid

  rest = 0
  qmin = -1
  if(myrank==0) write(1337,*) 'create_newcomms_group_ie input:',nranks,ne,nat
  
  ktake(:) = 0

  if((ne*nat)<nranks .and. (ne>1)) then ! .and. nat>1)) then
  
    if(nkmesh<=1) stop 'no load imbalance, please use regular grid'
  
    rest = nranks-int(nranks/(ne*nat))*ne*nat
    if(myrank==0) write(1337,*) 'rest:',rest,ne,nat,nranks
   
    !find fraction of k:l:m
    do ik=1,nkmesh-1
      k(ik) = (real(kmesh(1))/real(kmesh(nkmesh-ik+1))-1.)
    end do
    
    do ik=1,nkmesh-2
      f(ik) = real(k(ik+1))/real(k(ik))
    end do
    f(nkmesh-1) = real(k(nkmesh-1))/real(k(1))
    
    if(myrank==0) write(1337,*) 'set k,i:',k,'f',f

    !brute force look for optimal division of rest ranks after N_E*N_at are already
    !assigned to rectangular part of processor matrix:
    !                  N_E=8
    !            <--------------->
    !          ^ ( | | | | | | | )    example for 49 processors,
    !          | ( | | | | | | | )    devided according to:
    !  N_at=5  | ( | | | | | | | )         N_E = 8, N_at = 5 
    !          | ( | | | | | | | )        rest = 9 = 5+3+1
    !          v ( | | | | | | | )                   k+l+m  
    !                 ^    ( | | ) m=1  ^
    !             l=3 |      ( | )      |
    !                 v      ( | )      | k=5
    !                          ( )      |
    !                          ( )      v
    if(rest>0) then

    if(nkmesh==4) then
     ktake(:) = 0
     do i1=1,rest
      do i2=0,rest-i1
        i3 = rest-i1-i2
        if (i1>=i2 .and. i2>=i3) then
           if(i3==0 .and. i2==0) then
             q = sqrt((f(1)-real(i2)/real(i1))**2+(f(2)-1.)**2+(f(3)-real(i3)/real(i1))**2)
           else
             q = sqrt((f(1)-real(i2)/real(i1))**2+(f(2)-real(i3)/real(i2))**2+(f(3)-real(i3)/real(i1))**2)
           endif
           if(q<qmin .or. qmin==-1) then
               ktake = (/ i1,i2,i3 /)
               qmin = q
           end if
        end if
      enddo
     enddo
    elseif(nkmesh==3) then
     ktake(:) = 0
     do i1=1,rest
        i2 = rest-i1
        if (i1>=i2) then
           q = sqrt((f(1)-real(i2)/real(i1))**2)
           if(q<qmin .or. qmin==-1) then
               ktake = (/ i1,i2 /)
               qmin = q
           end if
        end if
     enddo
    elseif(nkmesh==2) then
     ktake(1) = rest
    else
     stop 'ERROR: nkmesh>4 not implemented yet'
    end if
    
    ! special case when only one additional rank
    if(rest==1) ktake(1) = rest

    if(myrank==0) write(1337,*) 'found ktake',ktake,'with',qmin
    end if
   
    !find processor groups according to non-uniform division
    allocate(groups(ne+1,2))
    groups(:,:) = -1
   
    do ie=1,ne+1
     if(ie==ne-2 .and. nkmesh>3) then
      groups(ie,1) = nat+ktake(3)
     elseif(ie==ne-1 .and. nkmesh>2) then
      groups(ie,1) = nat+ktake(2)
     elseif(ie==ne-0 .and. nkmesh>1) then
      groups(ie,1) = nat+ktake(1)
     else
      groups(ie,1) = nat
     endif
     if(ie==1) then
      groups(ie,2) = 0
     else
      groups(ie,2) = groups(ie-1,1)+groups(ie-1,2)
     endif
    enddo

    if(myrank==0) write(1337,*) 'groups:',groups(1:ne,1)
    if(myrank==0) write(1337,*) 'groups:',groups(1:ne,2)
   
    !find my group
    myg = -1
    do ie=1,ne
      do iat=groups(ie,2),groups(ie+1,2)-1
        if (myrank==iat) myg = ie
      end do
    end do
    
    if(myg==-1) then
      write(1337,*) 'no group found for rank', myrank
      stop
    end if

    !get group of processors in my group
    allocate(mygroup(groups(myg,1)))
   
   
    ie = 0
    do iat=groups(myg,2),groups(myg+1,2)-1
      ie = ie + 1
      mygroup(ie) = iat
    end do


    !create new communicator from group
    call MPI_COMM_GROUP(MPI_COMM_WORLD,mympi_group_world,ierr)
    call MPI_GROUP_INCL(mympi_group_world,groups(myg,1),mygroup,mympi_group_ie,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,mympi_group_ie,mympi_comm_ie,ierr)
    call MPI_GROUP_FREE(mympi_group_ie,ierr)
   
    !get rank and size in new communicator
    call MPI_COMM_RANK(mympi_comm_ie,myrank_ie,ierr)
    call MPI_COMM_SIZE(mympi_comm_ie,nranks_ie,ierr)
   
    !create communicator to communicate between differen energies (i.e. different groups)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,myrank_ie,myg,mympi_comm_at,ierr)
    call MPI_COMM_RANK(mympi_comm_at,myrank_atcomm,ierr)
    call MPI_COMM_SIZE(mympi_comm_at,nranks_atcomm,ierr)
    
    nranks_at = ne    
    myrank_at = myg-1

  elseif((ne*nat)==nranks) then ! .and. (ne>1 .and. nat>1)) then

    rest = 0

    call MPI_Cart_create( MPI_COMM_WORLD, 2, (/ ne, nat /), (/ .false., .false. /), (/ .true., .true. /), myMPI_comm_grid, ierr )

    call MPI_Cart_sub( myMPI_comm_grid, (/ .true., .false. /), myMPI_comm_at, ierr ) ! row communicator
    call MPI_Cart_sub( myMPI_comm_grid, (/ .false., .true. /), myMPI_comm_ie, ierr ) ! col communicator

    call MPI_Comm_rank( myMPI_comm_ie, myrank_ie, ierr )
    call MPI_Comm_rank( myMPI_comm_at, myrank_at, ierr )
    
    call MPI_Comm_size ( myMPI_comm_ie, nranks_ie, ierr )
    call MPI_Comm_size ( myMPI_comm_at, nranks_at, ierr )
    
    myrank_atcomm = myrank_at
    nranks_atcomm = nranks_at

  else
  
    rest = 0

    mympi_comm_at = MPI_COMM_WORLD
    myrank_at     = myrank
    nranks_at     = nranks
    mympi_comm_ie = MPI_COMM_SELF
    myrank_ie     = 0
    nranks_ie     = 1
    
    myrank_atcomm = myrank_at
    nranks_atcomm = nranks_at


  end if
  

  if(myrank==0) then
    write(1337,'(A)') '=================================================='  
    write(1337,'(A,I5,A)') '    MPI parallelization: use',nranks,' ranks'
    write(1337,'(A,I3,A,I4)') '    create processor array of size',nat,' x',ne
    if(rest>0) write(1337,'(A,I3)') '                                   with rest',rest
    if(rest>0) write(1337,'(A,10I3)') '    divide rest onto last energy points (k,l,m):',ktake
    write(1337,'(A)') '                N_E'
    write(1337,'(A)') '         <--------------->'
    write(1337,'(A)') '       ^ ( | | | | | | | )'
    write(1337,'(A)') '       | ( | | | | | | | )'
    write(1337,'(A)') '  N_at | ( | | | | | | | )'
    write(1337,'(A)') '       | ( | | | | | | | )'
    if(rest==0) write(1337,'(A)')'       v ( | | | | | | | )'
    if(rest>0) write(1337,'(A)') '       v ( | | | | | | | )....'
    if(rest>0) write(1337,'(A)') '              ^    ( | | ) m  ^'
    if(rest>0) write(1337,'(A)') '            l |      ( | )    |'
    if(rest>0) write(1337,'(A)') '              v......( | )    | k'
    if(rest>0) write(1337,'(A)') '                       ( )    |'
    if(rest>0) write(1337,'(A)') '                       ( )....v'
  end if

  end subroutine create_newcomms_group_ie
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm(IRMD,LMPOTD,NATYPD,LMAXD,LMAXD1,LMMAXD,NPOTD,IEMXD,MMAXD,IDOLDAU,NATYP,KREL,  &
                             & LMOMVEC,NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,denef,denefat,  &
                             & rhoorb,muorb,mvevi,mvevil,mvevief,mympi_comm)

    use mpi
    implicit none
    integer, intent(in) :: irmd,lmpotd,natypd,lmaxd,lmmaxd,iemxd,mmaxd,idoldau,natyp,krel,nmvecmax,npotd,lmaxd1,nqdos
    logical, intent(in) :: lmomvec
    integer, intent(in) :: mympi_comm
    double precision, intent(inout) :: RHO2NS(IRMD,LMPOTD,NATYPD,2), R2NEF(IRMD,LMPOTD,NATYPD,2),    &
                                     & ESPV(0:LMAXD1,NPOTD), DENEF, DENEFAT(NATYPD), RHOORB(IRMD*KREL + (1-KREL),NATYPD),   &
                                     & MUORB(0:LMAXD1+1,3,NATYPD) 
    double complex, intent(inout)   :: DEN(0:LMAXD1,IEMXD,NPOTD,NQDOS), DENLM(LMMAXD,IEMXD,NPOTD,NQDOS), DENMATC(MMAXD,MMAXD,NPOTD), MVEVI(NATYPD,3,NMVECMAX),   &
                                     & MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX), MVEVIEF(NATYPD,3,NMVECMAX)
    
    integer :: idim, ierr
    double precision, allocatable :: work1(:), work2(:,:), work3(:,:,:), work4(:,:,:,:)
    double complex,   allocatable :: work1c(:), work2c(:,:), work3c(:,:,:), work4c(:,:,:,:)
    
        allocate(work4(IRMD,LMPOTD,NATYPD,2))
        work4 = 0.d0
        IDIM = IRMD*LMPOTD*NATYPD*2
        CALL MPI_ALLREDUCE(RHO2NS,work4,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
        CALL DCOPY(IDIM,WORK4,1,RHO2NS,1)
        deallocate(work4)

        allocate(work4(IRMD,LMPOTD,NATYPD,2))
        work4 = 0.d0
        IDIM = IRMD*LMPOTD*NATYPD*2
        CALL MPI_ALLREDUCE(R2NEF,work4,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
        CALL DCOPY(IDIM,WORK4,1,R2NEF,1)
        deallocate(work4)

        !ESPV needs integration over atoms and energies -> MPI_COMM_WORLD
        allocate(work2(0:LMAXD+1,NPOTD))
        work2 = 0.d0
        IDIM = (LMAXD+2)*NPOTD
        CALL MPI_ALLREDUCE(ESPV,work2,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
        CALL DCOPY(IDIM,WORK2,1,ESPV,1)
        deallocate(work2)

        allocate(work4c(0:LMAXD+1,IEMXD,NPOTD,NQDOS))
        work4c = (0.d0, 0.d0)
        IDIM = IEMXD*(LMAXD+2)*NPOTD*NQDOS
        CALL MPI_ALLREDUCE(DEN,work4c,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
        CALL zCOPY(IDIM,WORK4c,1,DEN,1)
        deallocate(work4c)

        allocate(work4c(IEMXD,LMMAXD,NPOTD,NQDOS))
        work4c = (0.d0, 0.d0)
        IDIM = IEMXD*(LMMAXD)*NPOTD*NQDOS
        CALL MPI_ALLREDUCE(DENLM,work4c,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
        CALL zCOPY(IDIM,WORK4c,1,DENLM,1)
        deallocate(work4c)

        IF (IDOLDAU.EQ.1) THEN 
           allocate(work3c(MMAXD,MMAXD,NPOTD))
           work3c = (0.d0, 0.d0)
           IDIM = MMAXD*MMAXD*NPOTD
           CALL MPI_ALLREDUCE(DENMATC,work3c,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
           CALL zCOPY(IDIM,WORK3c,1,DENMATC,1)
           deallocate(work3c)
        END IF
        

        allocate(work1(1))
        work1 = 0.d0
        IDIM = 1
        CALL MPI_ALLREDUCE(DENEF,work1,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
        CALL DCOPY(IDIM,WORK1,1,DENEF,1)
        deallocate(work1)

        allocate(work1(NATYP))
        work1 = 0.d0
        IDIM = NATYP
        CALL MPI_ALLREDUCE(DENEFAT,work1,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
        CALL DCOPY(IDIM,WORK1,1,DENEFAT,1)
        deallocate(work1)

        IF (KREL.EQ.1) THEN 
          allocate(work2(IRMD,NATYPD))
          work2 = 0.d0
          IDIM = IRMD*NATYPD
          CALL MPI_ALLREDUCE(RHOORB,work2,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
          CALL DCOPY(IDIM,WORK2,1,RHOORB,1)
          deallocate(work2)

          allocate(work3(0:LMAXD+2,NATYPD,3))
          work3 = 0.d0
          IDIM = (LMAXD+3)*NATYPD*3
          CALL MPI_ALLREDUCE(MUORB,work3,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
          CALL DCOPY(IDIM,WORK3,1,MUORB,1)
          deallocate(work3)

          IF (LMOMVEC) THEN
             allocate(work3c(NATYPD,3,NMVECMAX))
             work3c = (0.d0, 0.d0)
             IDIM = NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVI,work3c,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
             CALL zCOPY(IDIM,WORK3c,1,MVEVI,1)
             deallocate(work3c)

             allocate(work4c(LMAXD+1,NATYPD,3,NMVECMAX))
             work4c = (0.d0, 0.d0)
             IDIM = (LMAXD+1)*NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVIL,work4c,IDIM, MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
             CALL zCOPY(IDIM,WORK4c,1,MVEVIL,1)
             deallocate(work4c)

             allocate(work3c(NATYPD,3,NMVECMAX))
             work3c = (0.d0, 0.d0)
             IDIM = NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVIEF,work3c,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
             CALL zCOPY(IDIM,WORK3c,1,MVEVIEF,1)
             deallocate(work3c)
             
          END IF   ! LMOMVEC
        END IF     ! KREL.EQ.1
        
  end subroutine mympi_main1c_comm
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm_newsosol(IRMDNEW,LMPOTD,LMAXD,LMAXD1,LMMAXD,  &
     &                                LMMAXSO,IEMXD,IELAST,NQDOS,            &
     &                                den,denlm,gflle,rho2nsc,r2nefc,        &
     &                                rho2int,espv,muorb,denorbmom,          &
     &                                denorbmomsp,denorbmomlm,denorbmomns,   &
     &                                mympi_comm)
     
     use mpi
     implicit none
     integer, intent(in) :: IRMDNEW, LMPOTD, LMAXD, LMAXD1, LMMAXD, LMMAXSO, IEMXD, IELAST, NQDOS
     integer, intent(in) :: mympi_comm
     double complex, intent(inout)   :: R2NEFC(IRMDNEW,LMPOTD,4), RHO2NSC(IRMDNEW,LMPOTD,4), DEN(0:LMAXD1,IEMXD,NQDOS,2), DENLM(LMMAXD,IEMXD,NQDOS,2), RHO2INT(4), GFLLE(LMMAXSO,LMMAXSO,IELAST,NQDOS)
     double precision, intent(inout) :: ESPV(0:LMAXD1,2), MUORB(0:LMAXD1+1,3), DENORBMOM(3), DENORBMOMSP(2,4), DENORBMOMLM(0:LMAXD,3), DENORBMOMNS(3)
     
     integer :: ierr, idim
     double precision, allocatable :: work(:,:,:,:)
     double complex, allocatable :: workc(:,:,:,:)

      !double complex arrays
     IDIM = IRMDNEW*LMPOTD*4
     allocate(workc(IRMDNEW,LMPOTD,4,1), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, r2nefc'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(r2nefc,workc(:,:,:,1),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,R2NEFC,1)
     deallocate(workc)

     IDIM = IRMDNEW*LMPOTD*4
     allocate(workc(IRMDNEW,LMPOTD,4,1), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, rho2nsc'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(RHO2NSC,workc,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,RHO2NSC,1)
     deallocate(workc)

     IDIM = (LMAXD1+1)*IEMXD*2*NQDOS
     allocate(workc(0:LMAXD1,IEMXD,2,NQDOS), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, den'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(DEN,workc,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,DEN,1)
     deallocate(workc)

     IDIM = LMMAXD*IEMXD*2*NQDOS
     allocate(workc(LMMAXD,IEMXD,2,NQDOS), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, denlm'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(DENLM,workc,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,DENLM,1)
     deallocate(workc)

     IDIM = 4
     allocate(workc(4,1,1,1), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, rho2int'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(RHO2INT,workc(:,1,1,1),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,RHO2INT,1)
     deallocate(workc)

     IDIM = LMMAXSO*LMMAXSO*IELAST*NQDOS
     allocate(workc(LMMAXSO,LMMAXSO,IELAST,NQDOS), stat=ierr)
     if(ierr/=0) stop 'mympi_main1c_comm_newsosol] Error allocating workc, gflle'
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(GFLLE,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,GFLLE,1)
     deallocate(workc)
     
     !double precision arrays
     IDIM = (LMAXD1+1)*2
     allocate(work(0:LMAXD1,2,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(ESPV,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,ESPV,1)
     deallocate(work)

     IDIM = (LMAXD1+2)*3
     allocate(work(0:LMAXD1+1,3,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(MUORB,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,MUORB,1)
     deallocate(work)

     IDIM = 3
     allocate(work(3,1,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(DENORBMOM,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENORBMOM,1)
     deallocate(work)

     IDIM = 2*4
     allocate(work(2,4,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(DENORBMOMSP,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENORBMOMSP,1)
     deallocate(work)

     IDIM = 3
     allocate(work(3,1,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(DENORBMOMNS,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENORBMOMNS,1)
     deallocate(work)

     IDIM = (LMAXD+1)*3
     allocate(work(0:LMAXD1,3,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(DENORBMOMLM,work,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENORBMOMLM,1)
     deallocate(work)
        
  end subroutine mympi_main1c_comm_newsosol
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm_newsosol2(LMAXD1,LMMAXD,IEMXD,NQDOS,      &
     &                                   NPOTD,NATYPD,LMPOTD,IRMD,MMAXD, &
     &                                   den, denlm, muorb, espv, r2nef, &
     &                                   rho2ns, denefat, denef,denmatn, &
     &                                   angles_new, mympi_comm)
     
     use mpi
     implicit none
     integer, intent(in) :: LMAXD1,IEMXD,NQDOS,NPOTD,NATYPD,LMPOTD,IRMD,LMMAXD,MMAXD
     integer, intent(in) :: mympi_comm
     double complex, intent(inout)   :: den(0:LMAXD1,IEMXD,NQDOS,NPOTD), denlm(LMMAXD,IEMXD,NQDOS,NPOTD),denmatn(MMAXD,MMAXD,2,2,NATYPD)
     double precision, intent(inout) :: muorb(0:LMAXD1+1,3,NATYPD), espv(0:LMAXD1,NPOTD), r2nef(IRMD,LMPOTD,NATYPD,2), rho2ns(IRMD,LMPOTD,NATYPD,2), denefat(NATYPD), denef, angles_new(NATYPD,2)
     
     integer :: ierr, idim
     double precision, allocatable :: work(:,:,:,:)
     double complex, allocatable :: workc(:,:,:,:)
     double complex, allocatable :: workc1(:,:,:,:,:)
     
     !double complex arrays
     IDIM = (1+LMAXD1)*IEMXD*NQDOS*NPOTD
     allocate(workc(0:LMAXD1,IEMXD,NQDOS,NPOTD))
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(den,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,DEN,1)
     deallocate(workc)
     
     IDIM = LMMAXD*IEMXD*NQDOS*NPOTD
     allocate(workc(LMMAXD,IEMXD,NQDOS,NPOTD))
     workc = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(denlm,workc(:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC,1,DENLM,1)
     deallocate(workc)

     IDIM = MMAXD*MMAXD*NPOTD
     allocate(workc1(MMAXD,MMAXD,2,2,NATYPD))
     workc1 = (0.d0, 0.d0)
     CALL MPI_ALLREDUCE(denmatn,workc1(:,:,:,:,:),IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,mympi_comm,IERR)
     CALL ZCOPY(IDIM,WORKC1,1,DENMATN,1)
     deallocate(workc1)

     !double precision arrays
     IDIM = (LMAXD1+2)*3*NATYPD
     allocate(work(0:LMAXD1+1,3,NATYPD,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(muorb,work(:,:,:,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,MUORB,1)
     deallocate(work)

     IDIM = (LMAXD1+1)*NPOTD
     allocate(work(0:LMAXD1,NPOTD,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(ESPV,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,ESPV,1)
     deallocate(work)

     IDIM = IRMD*LMPOTD*NATYPD*2
     allocate(work(IRMD,LMPOTD,NATYPD,2))
     work = 0.d0
     CALL MPI_ALLREDUCE(r2nef,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,R2NEF,1)
     deallocate(work)

     IDIM = IRMD*LMPOTD*NATYPD*2
     allocate(work(IRMD,LMPOTD,NATYPD,2))
     work = 0.d0
     CALL MPI_ALLREDUCE(RHO2NS,work(:,:,:,:),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,RHO2NS,1)
     deallocate(work)

     IDIM = NATYPD
     allocate(work(NATYPD,1,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(denefat,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENEFAT,1)
     deallocate(work)

     IDIM = 1
     allocate(work(1,1,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(denef,work(:,1,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,DENEF,1)
     deallocate(work)

     IDIM = 2*NATYPD
     allocate(work(2,NATYPD,1,1))
     work = 0.d0
     CALL MPI_ALLREDUCE(angles_new,work(:,:,1,1),IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,mympi_comm,IERR)
     CALL DCOPY(IDIM,WORK,1,angles_new,1)
     deallocate(work)
         
  end subroutine mympi_main1c_comm_newsosol2
#endif

end module mod_mympi
