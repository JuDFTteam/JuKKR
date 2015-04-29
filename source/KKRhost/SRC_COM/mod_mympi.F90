module mod_mympi
!ruess: taken from Pkkr_sidebranch2D_2014_12_16, created by Bernd Zimmermann

implicit none

  private
  public :: myrank, nranks, master, mympi_init, distribute_linear_on_tasks, find_dims_2d, create_subarr_comm, mympi_main1c_comm, mympi_main1c_comm_newsosol

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
       dims(1) = 0
       dims(2) = nranks
    else
       dims(1) = ntot1
       dims(2) = nranks/ntot1
       if(mod(float(nranks)/float(ntot1),1.).ne.0) stop 'ERROR in find_dims_2d'
    end if
    
    write(*,*) 'find_dims',myrank,nranks,ntot1,ntot2,dims
    
  end subroutine find_dims_2d
#endif


#ifdef CPP_MPI
  subroutine create_subarr_comm( subarr_dim, myMPI_comm_grid, myMPI_comm_row,myMPI_comm_col, myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col)
    use mpi
    implicit none
    integer, intent(in) :: subarr_dim(2)
    integer, intent(out) :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col,myrank_grid, myrank_row, myrank_col, nranks_row, nranks_col

    integer :: ierr
    logical :: logic2(2)
    logical, parameter :: periodic(2) = .false., reorder(2) = .false.

 
    if (subarr_dim(2).le.1) then
      if(myrank==master) write(*,*) 'option 1 in create_subarr_comm'
      myMPI_comm_grid = MPI_COMM_WORLD
      myMPI_comm_row  = MPI_COMM_WORLD
      myMPI_comm_col  = MPI_COMM_WORLD
      myrank_grid     = myrank
      myrank_row      = myrank
      myrank_col      = 0
      nranks_row      = nranks
      nranks_col      = 1
    elseif (subarr_dim(1).le.1) then
      if(myrank==master) write(*,*) 'option 2 in create_subarr_comm'
      myMPI_comm_grid = MPI_COMM_WORLD
      myMPI_comm_row  = MPI_COMM_WORLD
      myMPI_comm_col  = MPI_COMM_WORLD
      myrank_grid     = myrank
      myrank_row      = 0
      myrank_col      = myrank
      nranks_row      = 1
      nranks_col      = nranks
    else
      if(myrank==master) write(*,*) 'option 3 in create_subarr_comm'

      call MPI_Cart_create( MPI_COMM_WORLD, 2, subarr_dim, periodic, reorder, myMPI_comm_grid, ierr )
      call MPI_Comm_rank( myMPI_comm_grid, myrank_grid, ierr )

      logic2 = (/ .true., .false. /)
      call MPI_Cart_sub( myMPI_comm_grid, logic2, myMPI_comm_row, ierr ) ! row communicator
      logic2 = (/ .false., .true. /)
      call MPI_Cart_sub( myMPI_comm_grid, logic2, myMPI_comm_col, ierr ) ! col communicator

      call MPI_Comm_rank( myMPI_comm_row, myrank_row, ierr )
      call MPI_Comm_rank( myMPI_comm_col, myrank_col, ierr )
      
      call MPI_Comm_size ( myMPI_comm_row, nranks_row, ierr )
      call MPI_Comm_size ( myMPI_comm_col, nranks_col, ierr )

    end if
    
          write(*,*) 'in create_subarr_comm',myrank,nranks,subarr_dim,myrank_row,myrank_col,nranks_row,nranks_col


  end subroutine create_subarr_comm
#endif

#ifdef CPP_MPI
  subroutine mympi_main1c_comm(IRMD,LMPOTD,NATYPD,LMAXD,LMAXD1,NPOTD,IEMXD,MMAXD,IDOLDAU,NATYP,KREL,  &
                             & LMOMVEC,NMVECMAX,rho2ns,r2nef,espv,den,denmatc,denef,denefat,  &
                             & rhoorb,muorb,mvevi,mvevil,mvevief)

    use mpi
    implicit none
    integer, intent(in) :: irmd,lmpotd,natypd,lmaxd,iemxd,mmaxd,idoldau,natyp,krel,nmvecmax,npotd,lmaxd1
    logical, intent(in) :: lmomvec
    double precision, intent(inout) :: RHO2NS(IRMD,LMPOTD,NATYPD,2), R2NEF(IRMD,LMPOTD,NATYPD,2),    &
                                     & ESPV(0:LMAXD1,NPOTD), DENEF, DENEFAT(NATYPD), RHOORB(IRMD*KREL + (1-KREL),NATYPD),   &
                                     & MUORB(0:LMAXD1+1,3,NATYPD) 
    double complex, intent(inout)   :: DEN(0:LMAXD1,IEMXD,NPOTD), DENMATC(MMAXD,MMAXD,NPOTD), MVEVI(NATYPD,3,NMVECMAX),   &
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

        IDIM = IEMXD*(LMAXD+2)*NPOTD
        CALL MPI_ALLREDUCE(DEN,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
        CALL ZCOPY(IDIM,WORK,1,DEN,1)

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
     &                                LMMAXSO,IEMXD,NQDOS,            &
     &                                den,denlm,gflle,rho2nsc,r2nefc,   &
     &                                rho2int,espv,muorb,denorbmom,    &
     &                                denorbmomsp,denorbmomlm,denorbmomns)
     
     use mpi
     implicit none
     integer, intent(in) :: IRMDNEW, LMPOTD, LMAXD, LMAXD1, LMMAXD, LMMAXSO, IEMXD, NQDOS
     double complex, intent(inout)   :: R2NEFC(IRMDNEW,LMPOTD,4), RHO2NSC(IRMDNEW,LMPOTD,4), DEN(0:LMAXD1,IEMXD,2,NQDOS), DENLM(LMMAXD,IEMXD,2,NQDOS), RHO2INT(4), GFLLE(LMMAXSO,LMMAXSO,IEMXD,NQDOS)
     double precision, intent(inout) :: ESPV(0:LMAXD1,2), MUORB(0:LMAXD1+1,3), DENORBMOM(3), DENORBMOMSP(2,4), DENORBMOMLM(0:LMAXD,3), DENORBMOMNS(3)
     
     integer :: ierr, idim
     double complex, allocatable :: work(:,:,:,:)
     
     allocate(work(IRMDNEW,IEMXD,IEMXD,NQDOS),stat=ierr)
     if(ierr.ne.0) stop 'problem allocating work array in mympi_main1c_comm_newsosol'

     !double complex arrays
     IDIM = IRMDNEW*LMPOTD*4
     CALL MPI_ALLREDUCE(r2nefc,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,r2nefc,1)
     if (myrank==master) write(*,*) 'mpireduce done for r2nefc'
     
     IDIM = IRMDNEW*LMPOTD*4
     CALL MPI_ALLREDUCE(RHO2NSC,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,RHO2NSC,1)
     if (myrank==master) write(*,*) 'mpireduce done for rho2nsc'
     
     IDIM = (LMAXD1+1)*IEMXD*2*NQDOS
     CALL MPI_ALLREDUCE(DEN,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DEN,1)
     if (myrank==master) write(*,*) 'mpireduce done for den'
     
     IDIM = LMMAXD*IEMXD*2*NQDOS
     CALL MPI_ALLREDUCE(DENLM,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENLM,1)
     if (myrank==master) write(*,*) 'mpireduce done for denlm'
     
     IDIM = 4
     CALL MPI_ALLREDUCE(RHO2INT,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,RHO2INT,1)
     if (myrank==master) write(*,*) 'mpireduce done for rho2int'
     
     IDIM = LMMAXSO*LMMAXSO*IEMXD*NQDOS
     CALL MPI_ALLREDUCE(GFLLE,WORK,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,GFLLE,1)
     if (myrank==master) write(*,*) 'mpireduce done for gflle'
     
     !double precision arrays
     IDIM = (LMAXD1+1)*2
     CALL MPI_ALLREDUCE(ESPV,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,ESPV,1)
     if (myrank==master) write(*,*) 'mpireduce done for espv'
     
     IDIM = (LMAXD1+2)*3
     CALL MPI_ALLREDUCE(MUORB,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,MUORB,1)
     if (myrank==master) write(*,*) 'mpireduce done for muorb'
     
     IDIM = 3
     CALL MPI_ALLREDUCE(DENORBMOM,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOM,1)     
     if (myrank==master) write(*,*) 'mpireduce done for denorbmom'
     
     IDIM = 2*4
     CALL MPI_ALLREDUCE(DENORBMOMSP,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMSP,1)  
     if (myrank==master) write(*,*) 'mpireduce done for denorbmomsp'
     
     IDIM = 3
     CALL MPI_ALLREDUCE(DENORBMOMNS,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMNS,1)     
     if (myrank==master) write(*,*) 'mpireduce done for denorbmomns'
     
     IDIM = (LMAXD+1)*3
     CALL MPI_ALLREDUCE(DENORBMOMLM,WORK,IDIM,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
     CALL ZCOPY(IDIM,WORK,1,DENORBMOMLM,1)
     if (myrank==master) write(*,*) 'mpireduce done for denorbmomlm'
     
     deallocate(work)
        
  end subroutine mympi_main1c_comm_newsosol
#endif


end module mod_mympi
