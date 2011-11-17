subroutine SREDGX(ISPIN,NSPIN, &                 ! >
                  MYRANK, &                      ! >
                  SMPIC,SMYRANK, &               ! >
                  GMATXIJ, &                     ! >
                  GXIJ_ALL, &                    ! < output
                  ! new input parameters after inc.p removal
                  naez, lmax, lmpid, empid, smpid, nxijd)

! ======================================================================
! this routine preforms an MPI_ALLREDUCE command for
! a given set of ranks
!                               Alexander Thiess, 7th of December 2009
! =======================================================================

  !use mpi

  implicit none
  include 'mpif.h'

  integer, intent(in) :: naez
  integer, intent(in) :: lmax
  integer, intent(in) :: lmpid
  integer, intent(in) :: smpid
  integer, intent(in) :: empid
  integer, intent(in) :: nxijd

  !integer::        LMMAXD
  !parameter      (LMMAXD= (LMAXD+1)**2)

  !     .. scalar arguments ..
  integer::ISPIN
  integer::NSPIN
  !     .. array arguments ..
  !double complex :: GMATXIJ(LMMAXD,LMMAXD,NXIJD,NSPIND)
  !double complex :: GXIJ_ALL(LMMAXD,LMMAXD,NXIJD,NSPIND)
  double complex :: GMATXIJ ((LMAX+1)**2,(LMAX+1)**2,NXIJD,NSPIN)
  double complex :: GXIJ_ALL((LMAX+1)**2,(LMAX+1)**2,NXIJD,NSPIN)

  !     .. local scalars ..
  integer::IX
  integer::LM1
  integer::LM2
  integer::SEND
  integer::RECV

  !     .. local arrays ..
  !double complex :: GSEND(LMMAXD,LMMAXD,NXIJD)
  !double complex :: GRECV(LMMAXD,LMMAXD,NXIJD)
  double complex :: GSEND((LMAX+1)**2,(LMAX+1)**2,NXIJD)
  double complex :: GRECV((LMAX+1)**2,(LMAX+1)**2,NXIJD)

  !     .. MPI ..
  integer, dimension(MPI_STATUS_SIZE) :: STATUS
  integer::IERR
  integer::MAPBLOCK
  !     .. N-MPI
  integer::        MYRANK
  !     .. LS-MPI
  integer::SMYRANK(SMPID,NAEZ*LMPID*EMPID)
  integer::SMPI
  integer::SMPIC

  integer::LMMAXD

  LMMAXD= (LMAX+1)**2

  do ISPIN=1,NSPIN
    do IX=1,NXIJD
      do LM1=1,LMMAXD
        do LM2=1,LMMAXD
          GXIJ_ALL(LM1,LM2,IX,ISPIN)=GMATXIJ(LM1,LM2,IX,ISPIN)
        enddo
      enddo
    enddo
  enddo

  if (SMPID == 2) then
    
    do SMPI=1, SMPID
        
      SEND = SMYRANK(MAPBLOCK(SMPI,1,SMPID,1,0,SMPID-1)+1,SMPIC)
      RECV = SMYRANK(MAPBLOCK(1+MOD(SMPI,2),1,SMPID,1,0,SMPID-1)+1,SMPIC)
        
      if (MYRANK == SEND) then
            
        do IX=1,NXIJD
          do LM1=1,LMMAXD
            do LM2=1,LMMAXD
              GSEND(LM1,LM2,IX)=GMATXIJ(LM1,LM2,IX,SMPI)
            enddo
          enddo
        enddo
            
        call MPI_SEND(GSEND,LMMAXD*LMMAXD*NXIJD, &
        MPI_DOUBLE_COMPLEX, &
        RECV,97,MPI_COMM_WORLD,IERR)
            
      endif
        
      if (MYRANK == RECV) then
            
        call MPI_RECV(GRECV,LMMAXD*LMMAXD*NXIJD, &
        MPI_DOUBLE_COMPLEX, &
        SEND,97,MPI_COMM_WORLD,STATUS,IERR)
            
        do IX=1,NXIJD
          do LM1=1,LMMAXD
            do LM2=1,LMMAXD
              GXIJ_ALL(LM1,LM2,IX,SMPI)=GRECV(LM1,LM2,IX)
            enddo
          enddo
        enddo
            
      endif
        
    enddo
    
  endif

end subroutine SREDGX
