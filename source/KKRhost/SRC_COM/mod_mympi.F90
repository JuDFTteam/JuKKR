module mod_mympi
!ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann

implicit none

  private
  public :: myrank, nranks, master, mympi_init
#ifdef CPP_MPI
  public :: distribute_linear_on_tasks, find_dims_2d, create_newcomms_group_ie, mympi_main1c_comm, mympi_main1c_comm_newsosol
#endif

  integer, save :: myrank = -1
  integer, save :: nranks = -1
  integer, save :: master = -1

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
  subroutine distribute_linear_on_tasks(nranks, myrank, master, ntot, ntot_pT, ioff_pT, output)
    !distributes 'ntot' points on 'nranks' tasks and returns the number of points per task, ntot_pT, and the offsets of the tasks, ioff_pT.

    implicit none

    integer, intent(in)  :: nranks, myrank, master, ntot
    logical, intent(in)  :: output
    integer, intent(out) :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)

    integer :: irest, irank
    
!     write(*,*) 'distribute call:',myrank, nranks,ntot

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

    if(myrank==master .and. output) then
      write(*,*) '==== DISTRIBUTION OF POINTS ON TASKS: ===='
      do irank=0,nranks-1
        write(*,'("Task ",I0," treats points ",I0," to ",I0,", #of points= ",I0)') irank, ioff_pT(irank)+1, ioff_pT(irank)+ntot_pT(irank), ntot_pT(irank)
      end do!irank
      write(*,*) '=========================================='
    end if!myrank==master

  end subroutine distribute_linear_on_tasks
#endif


#ifdef CPP_MPI
  subroutine find_dims_2d(nranks,ntot1,ntot2,dims)
    !find dimensions to create cartesian communicator
    implicit none
    integer, intent(in)  :: nranks,ntot1,ntot2
    integer, intent(out) :: dims(2)
    

    if(nranks.le.ntot2) then
       dims(1) = 1
       dims(2) = nranks
    else
       ! rest not implemented!!!
!        stop 'ERROR: only parallelisation with maximally the number of energy points can be used!'
       dims(1) = nranks/ntot2
       dims(2) = ntot2
!        if(mod(float(nranks)/float(ntot1),1.).ne.0) stop 'ERROR in find_dims_2d'
    end if
    
!     write(*,*) 'find_dims',myrank,nranks,ntot1,ntot2,dims
    
  end subroutine find_dims_2d
#endif


#ifdef CPP_MPI
  subroutine create_newcomms_group_ie(nranks,myrank,nat,ne,nkmesh,kmesh,mympi_comm_ie,  &
  &                                  myrank_ie,nranks_ie,mympi_comm_at,myrank_at, nranks_at, myrank_atcomm, nranks_atcomm)
  !takes vector kmesh with mesh/timing information and finds number of rest procs that are devided in fractions given in ktake for optimal division of work
  
  use mpi
!   use mod_types, only: t_inc
  implicit none
  
  integer, intent(in) :: nranks,myrank,ne,nat,nkmesh
  integer, intent(in) :: kmesh(nkmesh)
  integer, intent(out) :: mympi_comm_ie, mympi_comm_at, nranks_atcomm, nranks_at, nranks_ie, myrank_ie, myrank_at, myrank_atcomm
  
  integer :: rest, k(nkmesh-1), ktake(nkmesh-1), myg, ie, iat, ierr, ik, i1,i2,i3
  double precision :: f(nkmesh-1), q, qmin
  integer, allocatable :: groups(:,:), mygroup(:)
  integer :: mympi_group_ie, mympi_group_world, myMPI_comm_grid
  
  if(myrank==0) write(*,*) 'create_newcomms_group_ie input:',nranks,ne,nat
  
  
!   if(nranks>(t_inc%natyp*t_inc%ielast)) stop 'you can only use less or equal number of processors that energypoints*atoms'
  
  if((ne*nat)<nranks .and. (ne>1 .and. nat>1)) then
  
    rest = nranks-int(nranks/(ne*nat))*ne*nat
    if(myrank==0) write(*,*) 'rest:',rest
   
    !find fraction of k:l:m
    do ik=1,nkmesh-1
      k(ik) = (real(kmesh(1))/real(kmesh(nkmesh-ik+1))-1.)
    end do
    
    do ik=1,nkmesh-2
      f(ik) = real(k(ik+1))/real(k(ik))
    end do
    f(nkmesh-1) = real(k(nkmesh-1))/real(k(1))
    
    if(myrank==0) write(*,*) 'set k,i:',k,'f',f

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
     qmin = -1
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
     qmin = -1
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

    if(myrank==0) write(*,*) 'found ktake',ktake,'with',qmin
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

    if(myrank==0) write(*,*) 'groups:',groups(1:ne,1)
    if(myrank==0) write(*,*) 'groups:',groups(1:ne,2)
   
    !find my group
    myg = -1
    do ie=1,ne
      do iat=groups(ie,2),groups(ie+1,2)-1
        if (myrank==iat) myg = ie
      end do
    end do
    
    if(myg==-1) then
      write(*,*) 'no group found for rank', myrank
      stop
    end if

    !get group of processors in my group
    allocate(mygroup(groups(myg,1)))
   
   
    ie = 0
    do iat=groups(myg,2),groups(myg+1,2)-1
      ie = ie + 1
      mygroup(ie) = iat
    end do

    write(*,'(A,I3,A,I3,A,I3,A,100I3)') 'rank ',myrank ,' found group: ',myg,' of size',groups(myg,1),' with group members:',mygroup
 
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
    
  elseif((ne*nat)==nranks .and. (ne>1 .and. nat>1)) then
  
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
    write(*,*) 'comm_test',myMPI_comm_at==MPI_COMM_WORLD
    myrank_at     = myrank
    nranks_at     = nranks
    mympi_comm_ie = MPI_COMM_SELF
    myrank_ie     = 0
    nranks_ie     = 1
    
    myrank_atcomm = myrank_at
    nranks_atcomm = nranks_at

  end if
  
      write(*,'(A,100I3)') 'my_newcomm_props:',myrank,nranks,nranks_ie,myrank_ie,myrank_at,nranks_at,nranks_atcomm,myrank_atcomm
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!       write(*,*) 'newcomm_at:',myMPI_comm_at==MPI_COMM_WORLD,myMPI_comm_at==MPI_COMM_NULL,myMPI_comm_at==MPI_COMM_SELF
!       write(*,*) 'newcomm_ie:',myMPI_comm_ie==MPI_COMM_WORLD,myMPI_comm_ie==MPI_COMM_NULL,myMPI_comm_ie==MPI_COMM_SELF
!       call MPI_BARRIER(MPI_COMM_WORLD, ierr)

  if(myrank==0) then
    write(*,'(A)') '=================================================='  
    write(*,'(A,I5,A)') '    MPI parallelization: use',nranks,' ranks'
    !write(*,'(AI3AI5A)') '    devide these onto Ne=',ne,' energy points and ',nat,' atoms'
    write(*,'(A,I3,A,I4,A,I3)') '    create processor array of size',nat,' x',ne,' with rest',rest
    write(*,'(A,10I3)') '    divide rest onto last energy points (k,l,m):',ktake
    write(*,'(A)') '                N_E'
    write(*,'(A)') '         <--------------->'
    write(*,'(A)') '       ^ ( | | | | | | | )'
    write(*,'(A)') '       | ( | | | | | | | )'
    write(*,'(A)') '  N_at | ( | | | | | | | )'
    write(*,'(A)') '       | ( | | | | | | | )'
    write(*,'(A)') '       v ( | | | | | | | )'
    write(*,'(A)') '              ^    ( | | ) m  ^'
    write(*,'(A)') '            l |      ( | )    |'
    write(*,'(A)') '              v      ( | )    | k'
    write(*,'(A)') '                       ( )    |'
    write(*,'(A)') '                       ( )    v'
  end if

  end subroutine create_newcomms_group_ie
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm(IRMD,LMPOTD,NATYPD,LMAXD,LMAXD1,LMMAXD,NPOTD,IEMXD,MMAXD,IDOLDAU,NATYP,KREL,  &
                             & LMOMVEC,NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,denef,denefat,  &
                             & rhoorb,muorb,mvevi,mvevil,mvevief)

    use mpi
    implicit none
    integer, intent(in) :: irmd,lmpotd,natypd,lmaxd,lmmaxd,iemxd,mmaxd,idoldau,natyp,krel,nmvecmax,npotd,lmaxd1,nqdos
    logical, intent(in) :: lmomvec
    double precision, intent(inout) :: RHO2NS(IRMD,LMPOTD,NATYPD,2), R2NEF(IRMD,LMPOTD,NATYPD,2),    &
                                     & ESPV(0:LMAXD1,NPOTD), DENEF, DENEFAT(NATYPD), RHOORB(IRMD*KREL + (1-KREL),NATYPD),   &
                                     & MUORB(0:LMAXD1+1,3,NATYPD) 
    double complex, intent(inout)   :: DEN(0:LMAXD1,IEMXD,NPOTD,NQDOS), DENLM(LMMAXD,IEMXD,NPOTD,NQDOS), DENMATC(MMAXD,MMAXD,NPOTD), MVEVI(NATYPD,3,NMVECMAX),   &
                                     & MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX), MVEVIEF(NATYPD,3,NMVECMAX)
    
    integer :: idim, ierr
    double complex, allocatable :: work(:,:,:,:)
    
        allocate(work(IRMD,LMPOTD,NATYPD,2),stat=ierr)
        if(ierr.ne.0) stop 'problem allocating work array in mympi_main1c_comm'
    
    
        IDIM = IRMD*LMPOTD*NATYPD*2
        CALL MPI_ALLREDUCE(RHO2NS,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL DCOPY(IDIM,WORK,1,RHO2NS,1)
        
        IDIM = IRMD*LMPOTD*NATYPD*2
        CALL MPI_ALLREDUCE(R2NEF,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL DCOPY(IDIM,WORK,1,R2NEF,1)
        
        
        IDIM = (LMAXD+2)*NPOTD
        CALL MPI_ALLREDUCE(ESPV,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL DCOPY(IDIM,WORK,1,ESPV,1)

        IDIM = IEMXD*(LMAXD+2)*NPOTD*NQDOS
        CALL MPI_ALLREDUCE(DEN,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL ZCOPY(IDIM,WORK,1,DEN,1)
        
        IDIM = IEMXD*(LMMAXD)*NPOTD*NQDOS
        CALL MPI_ALLREDUCE(DENLM,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL ZCOPY(IDIM,WORK,1,DENLM,1)

        IF (IDOLDAU.EQ.1) THEN 
           IDIM = MMAXD*MMAXD*NPOTD
           CALL MPI_ALLREDUCE(DENMATC,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
           CALL ZCOPY(IDIM,WORK,1,DENMATC,1)
        END IF
        
       
        IDIM = 1
        CALL MPI_ALLREDUCE(DENEF,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL DCOPY(IDIM,WORK,1,DENEF,1)

        IDIM = NATYP
        CALL MPI_ALLREDUCE(DENEFAT,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL DCOPY(IDIM,WORK,1,DENEFAT,1)
                

        IF (KREL.EQ.1) THEN 
          IDIM = IRMD*NATYPD
          CALL MPI_ALLREDUCE(RHOORB,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
          CALL DCOPY(IDIM,WORK,1,RHOORB,1)

          IDIM = (LMAXD+3)*NATYPD*3
          CALL MPI_ALLREDUCE(MUORB,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
          CALL DCOPY(IDIM,WORK,1,MUORB,1)

          IF (LMOMVEC) THEN
             IDIM = NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVI,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
             CALL ZCOPY(IDIM,WORK,1,MVEVI,1)

             IDIM = (LMAXD+1)*NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVIL,WORK,IDIM, MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
             CALL ZCOPY(IDIM,WORK,1,MVEVIL,1)

             IDIM = NATYPD*3*NMVECMAX
             CALL MPI_ALLREDUCE(MVEVIEF,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
             CALL ZCOPY(IDIM,WORK,1,MVEVIEF,1)
             
          END IF   ! LMOMVEC
        END IF     ! KREL.EQ.1
        
        deallocate(work)
        
  end subroutine mympi_main1c_comm
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm_newsosol(IRMDNEW,LMPOTD,LMAXD,LMAXD1,LMMAXD,  &
     &                                LMMAXSO,IEMXD,IELAST,NQDOS,            &
     &                                den,denlm,gflle,rho2nsc,r2nefc,   &
     &                                rho2int,espv,muorb,denorbmom,    &
     &                                denorbmomsp,denorbmomlm,denorbmomns)
     
     use mpi
     implicit none
     integer, intent(in) :: IRMDNEW, LMPOTD, LMAXD, LMAXD1, LMMAXD, LMMAXSO, IEMXD, IELAST, NQDOS
     double complex, intent(inout)   :: R2NEFC(IRMDNEW,LMPOTD,4), RHO2NSC(IRMDNEW,LMPOTD,4), DEN(0:LMAXD1,IEMXD,NQDOS,2), DENLM(LMMAXD,IEMXD,NQDOS,2), RHO2INT(4), GFLLE(LMMAXSO,LMMAXSO,IELAST,NQDOS)
     double precision, intent(inout) :: ESPV(0:LMAXD1,2), MUORB(0:LMAXD1+1,3), DENORBMOM(3), DENORBMOMSP(2,4), DENORBMOMLM(0:LMAXD,3), DENORBMOMNS(3)
     
     integer :: ierr, idim
     double complex, allocatable :: work(:,:,:,:)
     
     allocate(work(IRMDNEW,IEMXD,IEMXD,NQDOS),stat=ierr)
     if(ierr.ne.0) stop 'problem allocating work array in mympi_main1c_comm_newsosol'

     !double complex arrays
     IDIM = IRMDNEW*LMPOTD*4
     CALL MPI_ALLREDUCE(r2nefc,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,r2nefc,1)
!      if (myrank==master) write(*,*) 'mpireduce done for r2nefc'
     
     IDIM = IRMDNEW*LMPOTD*4
     CALL MPI_ALLREDUCE(RHO2NSC,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,RHO2NSC,1)
!      if (myrank==master) write(*,*) 'mpireduce done for rho2nsc'
     
     IDIM = (LMAXD1+1)*IEMXD*2*NQDOS
     CALL MPI_ALLREDUCE(DEN,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DEN,1)
!      if (myrank==master) write(*,*) 'mpireduce done for den'
     
     IDIM = LMMAXD*IEMXD*2*NQDOS
     CALL MPI_ALLREDUCE(DENLM,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENLM,1)
!      if (myrank==master) write(*,*) 'mpireduce done for denlm'
     
     IDIM = 4
     CALL MPI_ALLREDUCE(RHO2INT,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,RHO2INT,1)
!      if (myrank==master) write(*,*) 'mpireduce done for rho2int'
     
     IDIM = LMMAXSO*LMMAXSO*IELAST*NQDOS
     CALL MPI_ALLREDUCE(GFLLE,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,GFLLE,1)
!      if (myrank==master) write(*,*) 'mpireduce done for gflle'
     
     !double precision arrays
     IDIM = (LMAXD1+1)*2
     CALL MPI_ALLREDUCE(ESPV,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,ESPV,1)
!      if (myrank==master) write(*,*) 'mpireduce done for espv'
     
     IDIM = (LMAXD1+2)*3
     CALL MPI_ALLREDUCE(MUORB,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,MUORB,1)
!      if (myrank==master) write(*,*) 'mpireduce done for muorb'
     
     IDIM = 3
     CALL MPI_ALLREDUCE(DENORBMOM,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOM,1)     
!      if (myrank==master) write(*,*) 'mpireduce done for denorbmom'
     
     IDIM = 2*4
     CALL MPI_ALLREDUCE(DENORBMOMSP,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMSP,1)  
!      if (myrank==master) write(*,*) 'mpireduce done for denorbmomsp'
     
     IDIM = 3
     CALL MPI_ALLREDUCE(DENORBMOMNS,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMNS,1)     
!      if (myrank==master) write(*,*) 'mpireduce done for denorbmomns'
     
     IDIM = (LMAXD+1)*3
     CALL MPI_ALLREDUCE(DENORBMOMLM,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMLM,1)
!      if (myrank==master) write(*,*) 'mpireduce done for denorbmomlm'
     
     deallocate(work)
        
  end subroutine mympi_main1c_comm_newsosol
#endif


end module mod_mympi
