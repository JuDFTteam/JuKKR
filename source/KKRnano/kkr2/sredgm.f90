!> Communication between the different spin processes.
!> The diagonal part of the Green's function (needed for charge density) is
!> thrown with Lloyd's formula results into one array and communicated
!> (if necessary)
subroutine SREDGM(  NSPIN,IELAST, &           ! >
                    MYRANK, &                 ! >
                    SMPIC,SMYRANK, &          ! >
                    EMPIC,EMYRANK,EPROC, &    ! >
                    GMATN,LLY_GRDT, &         ! >
                    GMATN_ALL,LLY_GRDT_ALL, & ! < output
                    ! new input parameters after inc.p removal
                    naez, lmax, lmpid, smpid, empid, iemxd)

! ======================================================================
! this routine preforms an MPI_ALLREDUCE command for
! a given set of ranks
!                               Alexander Thiess, 7th of December 2009
!
! F90 conversion: Elias Rabel, Oct. 2011
! =======================================================================

  !use mpi
  implicit none

  include 'mpif.h'

  integer, intent(in) :: naez
  integer, intent(in) :: lmax
  integer, intent(in) :: lmpid
  integer, intent(in) :: smpid
  integer, intent(in) :: empid
  integer, intent(in) :: iemxd

  !integer::        LMMAXD
  !parameter      (LMMAXD= (LMAXD+1)**2)

  !     .. scalars ..
  integer::ISPIN
  integer::NSPIN
  !     .. arrays ..
  !double complex :: GMATN(LMMAXD,LMMAXD,IEMXD,NSPIND)
  !double complex :: GMATN_ALL(LMMAXD,LMMAXD,IEMXD,NSPIND)
  double complex :: GMATN((LMAX+1)**2,(LMAX+1)**2,IEMXD,NSPIN)
  double complex :: GMATN_ALL((LMAX+1)**2,(LMAX+1)**2,IEMXD,NSPIN)

  double complex :: LLY_GRDT(IEMXD,NSPIN)
  double complex :: LLY_GRDT_ALL(IEMXD,NSPIN)

  !     .. local scalars ..
  integer::IE
  integer::IELAST
  integer::LM1
  integer::LM2
  integer::SEND
  integer::RECV

  !     .. local arrays ..
  !double complex :: GSEND(LMMAXD*LMMAXD*IEMXD+IEMXD)
  !double complex :: GRECV(LMMAXD*LMMAXD*IEMXD+IEMXD)
  !double complex :: GESEND(LMMAXD*LMMAXD+1,NSPIND)
  !double complex :: GERECV(LMMAXD*LMMAXD+1,NSPIND)

  double complex :: GSEND((LMAX+1)**2 * (LMAX+1)**2 * IEMXD + IEMXD)
  double complex :: GRECV((LMAX+1)**2 * (LMAX+1)**2 * IEMXD + IEMXD)
  double complex :: GESEND((LMAX+1)**2 * (LMAX+1)**2 + 1, NSPIN)
  double complex :: GERECV((LMAX+1)**2 * (LMAX+1)**2 + 1, NSPIN)

  !     .. MPI ..
  integer, dimension(MPI_STATUS_SIZE) :: STATUS
  integer::IERR
  integer::MAPBLOCK
  !     .. N-MPI
  integer::MYRANK
  !     .. LS-MPI
  integer::SMYRANK(SMPID,NAEZ*LMPID*EMPID)
  integer::SMPI
  integer::SMPIC

  !     .. E-MPI
  integer::EMYRANK(EMPID,NAEZ*LMPID*SMPID)
  integer::EMPI
  integer::EMPIC
  integer::EPROC(IEMXD)

  integer::LMMAXD

  LMMAXD= (LMAX+1)**2

  !     fist copy the complete array to destination array:

  call ZCOPY(LMMAXD*LMMAXD*IEMXD*NSPIN,GMATN,1,GMATN_ALL,1)

  do ISPIN=1,NSPIN
    do IE=1,IELAST
      LLY_GRDT_ALL(IE,ISPIN) = LLY_GRDT(IE,ISPIN)
    enddo
  enddo


!=======================================================================
!     now take care that all processors get the results of
!     all energy points
!=======================================================================

  if (EMPID > 1) then
    
! IE
    do IE = 1, IELAST
        
      !       identify sending proc:
      SEND  = EMYRANK(MAPBLOCK(EPROC(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
        
        
! EMPI     begin loop over E-parallel processors
      do EMPI = 1, EMPID
! EMPI

        RECV = EMYRANK(MAPBLOCK(EMPI,1,EMPID,1,0,EMPID-1)+1,EMPIC)


            
        ! check if communication is necessary
        if (RECV == SEND .and. RECV == MYRANK) then
          !             in this case: no communication necessary
          do ISPIN=1,NSPIN
            do LM1=1,LMMAXD
              do LM2=1,LMMAXD
                GMATN_ALL(LM1,LM2,IE,ISPIN)=GMATN(LM1,LM2,IE,ISPIN)
              enddo
            enddo
            LLY_GRDT_ALL(IE,ISPIN)=LLY_GRDT(IE,ISPIN)
          enddo
                
        else
                
          if (MYRANK == SEND) then
            ! i am the sender
            do ISPIN=1,NSPIN
              call ZCOPY(LMMAXD*LMMAXD,GMATN(1,1,IE,ISPIN),1, &
              GESEND(1,ISPIN),1)
              GESEND(LMMAXD*LMMAXD+1,ISPIN)= &
              LLY_GRDT(IE,ISPIN)
            enddo
                    
            call MPI_SEND(GESEND,(LMMAXD*LMMAXD+1)*NSPIN, &
            MPI_DOUBLE_COMPLEX, &
            RECV,96,MPI_COMM_WORLD,IERR)
                    
          endif
                
          if (MYRANK == RECV) then
            ! i am the receiver
            call MPI_RECV(GERECV,(LMMAXD*LMMAXD+1)*NSPIN, &
            MPI_DOUBLE_COMPLEX, &
            SEND,96,MPI_COMM_WORLD,STATUS,IERR)
                    
            do ISPIN=1,NSPIN
              call ZCOPY(LMMAXD*LMMAXD,GERECV(1,ISPIN),1, &
              GMATN_ALL(1,1,IE,ISPIN),1)
              LLY_GRDT_ALL(IE,ISPIN)=GERECV(LMMAXD*LMMAXD+1,ISPIN)
            enddo
                    
          endif
                
        endif
            
! EMPI
      enddo
! EMPI
        
        
      if (SMPID == 2) then
        do ISPIN=1,NSPIN
          do LM1=1,LMMAXD
            do LM2=1,LMMAXD
              GMATN(LM1,LM2,IE,ISPIN) = GMATN_ALL(LM1,LM2,IE,ISPIN)
            enddo
          enddo
          LLY_GRDT(IE,ISPIN) = LLY_GRDT_ALL(IE,ISPIN)
        enddo
      endif
! IE
    enddo
! IE
    
  endif ! EMPID > 1

!=======================================================================
!=======================================================================




!=======================================================================
!     as a third step exchange information of different spin channels
!=======================================================================

  if (SMPID == 2) then
    
    do SMPI=1, SMPID
        
      SEND = SMYRANK(MAPBLOCK(SMPI,1,SMPID,1,0,SMPID-1)+1,SMPIC)
      RECV = SMYRANK(MAPBLOCK(1+MOD(SMPI,2),1,SMPID,1,0,SMPID-1)+1,SMPIC)
        
      if (MYRANK == SEND) then
            
        call ZCOPY(LMMAXD*LMMAXD*IEMXD,GMATN_ALL(1,1,1,SMPI),1, &
        GSEND,1)
        call ZCOPY(IEMXD,LLY_GRDT_ALL(1,SMPI),1, &
        GSEND(LMMAXD*LMMAXD*IEMXD+1),1)
            
        call MPI_SEND(GSEND,LMMAXD*LMMAXD*IEMXD+IEMXD, &
        MPI_DOUBLE_COMPLEX, &
        RECV,97,MPI_COMM_WORLD,IERR)
            
      endif
        
      if (MYRANK == RECV) then
            
        call MPI_RECV(GRECV,LMMAXD*LMMAXD*IEMXD+IEMXD, &
        MPI_DOUBLE_COMPLEX, &
        SEND,97,MPI_COMM_WORLD,STATUS,IERR)
            
        call ZCOPY(LMMAXD*LMMAXD*IEMXD,GRECV,1, &
        GMATN_ALL(1,1,1,SMPI),1)
        call ZCOPY(IEMXD,GRECV(LMMAXD*LMMAXD*IEMXD+1),1, &
        LLY_GRDT_ALL(1,SMPI),1)
            
      endif
        
    enddo
    
  endif

end subroutine SREDGM
