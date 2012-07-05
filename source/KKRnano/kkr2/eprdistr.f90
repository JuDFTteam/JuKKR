! TODO: This communication pattern could be wrong
    subroutine EPRDIST(IELAST,KMESH,NOFKS, &
                       PRSC,CNVFAC, &
                       MYRANK,EMPIC,EMYRANK, &
                       EPROC,EPROCO, &
                       ! new input parameters after inc.p removal
                       lmpid, smpid, empid, naez, lmax, nguessd, ekmd, iemxd)
! ======================================================================

! this routine redistributes on demand the initial guess
! information to the individual processors EMPI=1,...,EMPID

!                               Alexander Thiess, 14th of December 2009
! =======================================================================

    implicit none
    include 'mpif.h'

    integer, intent(in) :: lmpid
    integer, intent(in) :: smpid
    integer, intent(in) :: empid
    integer, intent(in) :: naez
    integer, intent(in) :: lmax
    integer, intent(in) :: nguessd
    integer, intent(in) :: ekmd
    integer, intent(in) :: iemxd

    !     .. Parameters ..

    integer::         LMMAXD
    !parameter        (LMMAXD= (LMAXD+1)**2)
    integer::         MAXMSHD
    parameter        (MAXMSHD=8)

    !     .. scalar arguments ..
    integer::          IELAST
    !     .. array arguments ..
   !complex::          PRSC(NGUESSD*LMMAXD,EKMD)
    complex::          PRSC(NGUESSD*(LMAX+1)**2,EKMD)
    double precision:: CNVFAC(EKMD)

    integer::          EPROC(IEMXD)
    integer::          EPROCO(IEMXD)
    integer::          KMESH(IEMXD)
    integer::          NOFKS(MAXMSHD)
    !     .. local scalars ..
    integer::IDIM
    integer::IE
    integer::JEKM  ! index for initial guess array
    integer::RECV
    integer::SEND
    !     .. local arrays ..


    !     .. MPI ..
    integer, dimension(MPI_STATUS_SIZE) :: STATUS
    integer::IERR
    integer::MAPBLOCK
    !     .. N-MPI
    integer::MYRANK
    !     .. E-MPI
    integer::EMYRANK(EMPID,NAEZ*LMPID*SMPID)
    integer::EMPIC

    LMMAXD = (LMAX+1)**2

    JEKM = 0

    do IE= 1, IELAST
    
    !       EPROC = E-process EPROCO = E-process old
    !       Check if the energy process which has worked on the energy
    !       point IE is not the same as in the previous self-
    !       consistency step

        if (EPROC(IE) /= EPROCO(IE)) then
        
            SEND  = &
            EMYRANK(MAPBLOCK(EPROCO(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
        
            RECV  = &
            EMYRANK(MAPBLOCK(EPROC(IE),1,EMPID,1,0,EMPID-1)+1,EMPIC)
        
        
        !         first exchange PRSC (initial guess)
        
            IDIM = NGUESSD*LMMAXD*NOFKS(KMESH(IE))
        
            if (MYRANK == SEND) then
                call MPI_SEND(PRSC(1,JEKM+1),IDIM, &
                MPI_COMPLEX, &
                RECV,92,MPI_COMM_WORLD,IERR)
            elseif (MYRANK == RECV) then
                call MPI_RECV(PRSC(1,JEKM+1),IDIM, &
                MPI_COMPLEX, &
                SEND,92,MPI_COMM_WORLD,STATUS,IERR)
            endif

        !         third step: exchange CNVFAC
        
!            IDIM = NOFKS(KMESH(IE))
        
!            if (MYRANK == SEND) then
!                call MPI_SEND(CNVFAC(JEKM+1),IDIM, &
!                MPI_DOUBLE_PRECISION, &
!                RECV,93,MPI_COMM_WORLD,IERR)
!            elseif (MYRANK == RECV) then
!                call MPI_RECV(CNVFAC(JEKM+1),IDIM, &
!                MPI_DOUBLE_PRECISION, &
!                SEND,93,MPI_COMM_WORLD,STATUS,IERR)
!            endif
        
        endif
    
        JEKM = JEKM + NOFKS(KMESH(IE))
    
    enddo

    end subroutine EPRDIST
