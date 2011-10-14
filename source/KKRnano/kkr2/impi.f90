subroutine IMPI( &
  NAEZ, &                                     ! > in
  MYRANK,NROFNODES, &                         ! < out
  LMPIC,MYLRANK,LGROUP,LCOMM,LSIZE, &         ! < out
  LSMPIB,LSMPIC,LSRANK,LSMYRANK, &            ! < out
  SMPIB,SMPIC,SRANK,SMYRANK, &                ! < out
  EMPIB,EMPIC,ERANK,EMYRANK, &                ! < out
  MYACTVRANK,ACTVGROUP,ACTVCOMM,ACTVSIZE, &   ! < out
  ! new input parameters after inc.p removal
  lmpid, smpid, empid, nthrds)                ! < in

  ! ======================================================================
  !                       build MPI_Groups
  ! ======================================================================

  ! the required information for l-parallelization are
  ! intitialized.

  ! all depends on the additional switch in inc.p:

  ! MPI_L = 1  .... still one process per atom
  ! MPI_L > 1  .... LMPI number of processes running per atom

  ! parallelization strategy:
  ! only the solution of the Dyson equation, resp. the inversion of
  ! the Lloyd-matrix are parallelized of LM. All tasks before are
  ! executed anyway, to save communication. After INVIT and LLYINVIT
  ! finished, the results for different L are communicated to the
  ! root process by MPI_REDUCE(...,MPI_SUM,...) and from here on only
  ! the rout processes are running.
  !                               Alexander Thiess, 19th of November 2009

  !                               f90 conversion: Elias Rabel, Oct 2011
  ! =======================================================================

  use mpi
  implicit none

  integer, intent(in) :: lmpid
  integer, intent(in) :: smpid
  integer, intent(in) :: empid
  integer, intent(in) :: nthrds

  integer::I1
  integer::NAEZ
  integer::ISPIN
  integer::IRANK
  !     .. MPI ..
  integer::WGROUP
  integer::IERR
  !     .. N-MPI
  integer::MYRANK
  integer::NROFNODES
  !     .. L-MPI

  integer::MYLRANK    (LMPID*SMPID*EMPID)
  integer::LRANKS(NAEZ,LMPID*SMPID*EMPID)
  integer::LCOMM      (LMPID*SMPID*EMPID)
  integer::LGROUP     (LMPID*SMPID*EMPID)
  integer::LSIZE      (LMPID*SMPID*EMPID)
  integer::LMPI
  integer::LMPIC

  !     .. LS-MPI
  integer::LSMYRANK(LMPID,NAEZ*SMPID*EMPID)
  integer::LSRANK  (LMPID,NAEZ*SMPID*EMPID)
  integer::LSMPIC
  integer::LSMPIB

  !     .. S-MPI
  integer::SMYRANK(SMPID,NAEZ*LMPID*EMPID)
  integer::SRANK  (SMPID,NAEZ*LMPID*EMPID)
  integer::SMPIC
  integer::SMPIB

  !     .. E-MPI
  integer::EMYRANK(EMPID,NAEZ*LMPID*SMPID)
  integer::ERANK  (EMPID,NAEZ*LMPID*SMPID)
  integer::EMPI
  integer::EMPIC
  integer::EMPIB

  !     .. ACTV-MPI
  integer::MYACTVRANK
  integer::ACTVRANKS(NAEZ*LMPID*SMPID*EMPID)
  integer::ACTVCOMM
  integer::ACTVGROUP
  integer::ACTVSIZE


  ! intialize MPI and standard groups and communicator

  call MPI_INIT(IERR)
  call MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)

  ! define ranks for all non-standard groups of processors

  LMPIC = 0

  LSMPIC = 0
  LSMPIB = 0

  SMPIC = 0
  SMPIB = 0

  EMPIC = 0
  EMPIB = 0


  do I1 = 1, NAEZ
    do ISPIN = 1, SMPID
      do LMPI = 1, LMPID
        do EMPI = 1, EMPID
                
          IRANK = (ISPIN-1)*(NAEZ*LMPID*EMPID) + &
          (LMPI-1)*NAEZ*EMPID + &
          (I1-1)*EMPID + &
          EMPI - 1
                
          LRANKS(I1,(ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID+EMPI) = IRANK
                
          LSMYRANK(LMPI,(ISPIN-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)= IRANK
          LSRANK(LMPI,(ISPIN-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)  = LMPI-1
                
          SMYRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI) = IRANK
          SRANK(ISPIN,(LMPI-1)*NAEZ*EMPID+(I1-1)*EMPID+EMPI)   = ISPIN-1
                
          EMYRANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)= IRANK
          ERANK(EMPI,(ISPIN-1)*NAEZ*LMPID+(I1-1)*LMPID+LMPI)  = LMPI-1
                
          ACTVRANKS(IRANK+1) = IRANK
                
          if (MYRANK == IRANK) then
                    
            LMPIC  = (ISPIN-1)*LMPID*EMPID+(LMPI-1)*EMPID + EMPI
                    
            LSMPIC = (ISPIN-1)*NAEZ*EMPID +(I1-1)  *EMPID + EMPI
            LSMPIB = LMPI
                    
            SMPIC  = (LMPI-1)*NAEZ*EMPID  +(I1-1)*EMPID   + EMPI
            SMPIB  = ISPIN
                    
            EMPIC  = (ISPIN-1)*NAEZ*LMPID +(I1-1)*LMPID   + LMPI
            EMPIB  = EMPI
                    
          endif
                
        enddo
      enddo
    enddo
  enddo

  ! build groups and communicators


  ! ACTIVE GROUP (ACTVGROUP) ...............................................

  call MPI_COMM_GROUP(MPI_COMM_WORLD,WGROUP,IERR) ! get default
  call MPI_GROUP_INCL(WGROUP,NAEZ*LMPID*SMPID*EMPID, &
                      ACTVRANKS(1), &
                      ACTVGROUP,IERR) !create a group ACTVGROUP
  call MPI_COMM_CREATE(MPI_COMM_WORLD,ACTVGROUP,ACTVCOMM,IERR)

  MYACTVRANK=-2
  call MPI_GROUP_RANK(ACTVGROUP,MYACTVRANK,IERR)
  call MPI_GROUP_SIZE(ACTVGROUP,ACTVSIZE,IERR)

  ! LGROUP ...............................................................

  do LMPI=1,LMPID*SMPID*EMPID
    call MPI_COMM_GROUP(MPI_COMM_WORLD, WGROUP, IERR) ! get default
    call MPI_GROUP_INCL(WGROUP, NAEZ, LRANKS(1,LMPI), &
                        LGROUP(LMPI), IERR) !create groups LGROUP
    call MPI_COMM_CREATE(MPI_COMM_WORLD, LGROUP(LMPI), &
                         LCOMM(LMPI), IERR) !create communicators for
  enddo

  do LMPI=1,LMPID*SMPID*EMPID
    MYLRANK(LMPI) = -2
    call MPI_GROUP_RANK(LGROUP(LMPI),MYLRANK(LMPI),IERR)
    call MPI_GROUP_SIZE(LGROUP(LMPI),LSIZE(LMPI),IERR)
  enddo

  !.......................................................................

  if (MYRANK == 0) then
    write(6,'(79(1H=))')
    write(6,*) '  groups of processes created'
    write(6,*) '  NMPI                    = ',NAEZ
    write(6,*) '  LMPI                    = ',LMPID
    write(6,*) '  SMPI                    = ',SMPID
    write(6,*) '  EMPI                    = ',EMPID
    write(6,*) '  NTHRDS                  = ',NTHRDS
    write(6,*) '  MPI-processes           = ',NAEZ*LMPID*SMPID*EMPID
    write(6,*) '  total no. of processors = ',NAEZ*LMPID*SMPID*EMPID*NTHRDS
    write(6,'(79(1H=))')
  endif


end subroutine IMPI
