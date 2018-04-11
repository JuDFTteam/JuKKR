subroutine operators_for_FScode(KORBIT)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  ! interface routine to normcoeff routines that prepare operators for
  ! use in FScode (compuation of spin expectation value etc.)
  !
  ! first wavefuncitons are read in and converted to old mesh and then
  ! normcoeff_* routines are called
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# ifdef CPP_MPI
  use mod_types, only: t_inc, t_mpi_c_grid
  use mpi
  use mod_mympi, only: nranks, master, myrank, distribute_linear_on_tasks
# else
  use mod_types, only: t_inc
  use mod_mympi, only: nranks, master, myrank
# endif
  use mod_wunfiles, only: t_params
  use mod_save_wavefun, only: t_wavefunctions, read_wavefunc

  implicit none

  integer, intent(in) :: KORBIT

  ! read in wavefunctions 
  logical :: rll_was_read_in, sll_was_read_in,  rllleft_was_read_in, sllleft_was_read_in
  double complex, allocatable :: RLL(:,:,:,:),SLL(:,:,:,:), RLLLEFT(:,:,:,:),SLLLEFT(:,:,:,:), RLLTEMP(:,:), PNSTEMP(:,:), PNS_SO(:,:,:,:), PNS_SO_ALL(:,:,:,:,:) 

  ! loop counter etc.
  integer :: ie, ie_start, ie_end, ie_num, lm1, lm2, ir, i1, i1_start, i1_end, ierr

  ! array dimensions
  integer :: lmmaxd, irmd, natyp, nsra, ncheb, ntotd

  ! arrays for rmeshes (old and new), and nonco_angles
  integer, allocatable :: irws(:), npan_tot(:), ipan_intervall(:,:)
  double precision, allocatable :: rmesh(:,:), rpan_intervall(:,:), theta(:), phi(:)

  ! constants
  double complex, parameter :: czero=(0.0d0, 0.0d0)

# ifdef CPP_MPI
  ! communcate PNS_SO_ALL for OPERATOR option
  integer :: ihelp
  integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1)
  double complex, allocatable :: work(:,:,:,:,:)
# endif

  ! first fill scalar and array parameters that are used here
  !  call get_params_operators(lmmaxd, irmd, natyp, nsra, ncheb, ntot, irws, 
  ! scalars
  lmmaxd = t_inc%lmmaxd
  irmd = t_params%irmd
  natyp = t_params%natypd
  nsra = t_params%nsra 
  ncheb = t_params%ncheb 
  ntotd = t_params%ntotd
  ! arrays
  allocate(irws(natyp))
  irws = t_params%irws
  allocate(rmesh(irmd, natyp))
  rmesh = t_params%r
  allocate(npan_tot(natyp))
  npan_tot = t_params%npan_tot
  allocate(rpan_intervall(0:ntotd, natyp))
  rpan_intervall = t_params%rpan_intervall
  allocate(ipan_intervall(0:ntotd, natyp))
  ipan_intervall = t_params%ipan_intervall
  allocate(theta(natyp), phi(natyp))
  theta = t_params%theta
  phi = t_params%phi

  ALLOCATE(PNS_SO_ALL(LMMAXD,LMMAXD,IRMD,2,NATYP))
  PNS_SO_ALL = CZERO

# ifdef CPP_MPI
  call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master,natyp,ntot_pT,ioff_pT,.true.)
  i1_start = ioff_pT(t_mpi_c_grid%myrank_ie) + 1
  i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie) + ntot_pT(t_mpi_c_grid%myrank_ie)
  t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)
  t_mpi_c_grid%ntot_pT1 = ntot_pT
  t_mpi_c_grid%ioff_pT1 = ioff_pT
# else
  i1_start = 1
  i1_end   = NATYP
# endif

  DO I1=i1_start, i1_end

    ALLOCATE(RLL(NSRA*LMMAXD,LMMAXD,t_inc%IRMDNEW,0:0))
    ALLOCATE(SLL(NSRA*LMMAXD,LMMAXD,t_inc%IRMDNEW,0:0))
    ALLOCATE(RLLLEFT(NSRA*LMMAXD,LMMAXD,t_inc%IRMDNEW,0:0))
    ALLOCATE(SLLLEFT(NSRA*LMMAXD,LMMAXD,t_inc%IRMDNEW,0:0))
    ALLOCATE(PNS_SO(LMMAXD,LMMAXD,IRMD,2))
    RLL = CZERO
    SLL = CZERO
    RLLLEFT = CZERO
    SLLLEFT = CZERO
    PNS_SO = CZERO

#   ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
    ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#   else
    ie_start = 0 ! offset
    ie_end   = t_params%IELAST
#   endif

    DO ie_num=1,ie_end

      IE = ie_start+ie_num

      ! make sure only calculated at the Fermi level
      !if(ie_end==1 .or. ie==1) then
      if(ie==1) then

         if(t_wavefunctions%Nwfsavemax>0) then ! read wavefunctions?
           ! read in wavefunction from memory
           call read_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,i1,ie,NSRA,LMMAXD,t_inc%IRMDNEW,0,1,rll_was_read_in,sll_was_read_in,rllleft_was_read_in,sllleft_was_read_in)
         end if ! t_wavefunctions%Nwfsavemax

      end if ! ie==1

    END DO !ie_num=1,ie_end

    ! transform radial wavefunction back to old mesh
    ALLOCATE(RLLTEMP(t_inc%IRMDNEW,LMMAXD))
    ALLOCATE(PNSTEMP(IRWS(I1),LMMAXD))
    DO LM1=1,LMMAXD
     RLLTEMP=CZERO
     PNSTEMP=CZERO
     IR = 0
     DO IR=1,t_inc%IRMDNEW
      DO LM2=1,LMMAXD
       RLLTEMP(IR,LM2)=RLL(LM1,LM2,IR,0)
      ENDDO
     ENDDO
     CALL CHEB2OLDGRID(IRWS(I1),t_inc%IRMDNEW,LMMAXD,RMESH(:,I1), NCHEB, NPAN_TOT(I1),RPAN_INTERVALL(:,I1), IPAN_INTERVALL(:,I1),RLLTEMP,PNSTEMP,IRMD)
     DO IR=1,IRWS(I1)
      DO LM2=1,LMMAXD
       PNS_SO(LM1,LM2,IR,1)=PNSTEMP(IR,LM2)
      ENDDO
     ENDDO
    ENDDO ! LM1
    ! for small component
    IF (NSRA.EQ.2) THEN
     DO LM1=1,LMMAXD
      RLLTEMP=CZERO
      PNSTEMP=CZERO
      DO IR=1,t_inc%IRMDNEW
       DO LM2=1,LMMAXD
        RLLTEMP(IR,LM2)=RLL(LM1+LMMAXD,LM2,IR,0)
       ENDDO
      ENDDO
      CALL CHEB2OLDGRID(IRWS(I1),t_inc%IRMDNEW,LMMAXD,RMESH(:,I1), NCHEB, NPAN_TOT(I1),RPAN_INTERVALL(:,I1), IPAN_INTERVALL(:,I1),RLLTEMP,PNSTEMP,IRMD)
      DO IR=1,IRWS(I1)
       DO LM2=1,LMMAXD
        PNS_SO(LM1,LM2,IR,2)=PNSTEMP(IR,LM2)
       ENDDO
      ENDDO
     ENDDO ! LM1
    ENDIF ! NSRA.EQ.2

    ! rotate radial wavefunction to global frame
    DO IR=1,IRMD
      CALL ROTATEMATRIX(PNS_SO(1,1,IR,1), THETA(I1), PHI(I1), LMMAXD/2,0)
      CALL ROTATEMATRIX(PNS_SO(1,1,IR,2), THETA(I1), PHI(I1), LMMAXD/2,0)
    ENDDO

    ! finally collect wavefuncitons in global frame and old mesh for all atoms to be used in normcoeff-routines below
    PNS_SO_ALL(:,:,:,:,i1) = PNS_SO(:,:,:,:)

    DEALLOCATE(RLL,SLL,RLLLEFT,SLLLEFT,PNS_SO)
    DEALLOCATE(RLLTEMP)
    DEALLOCATE(PNSTEMP)

  END DO ! I1=i1_start, i1_end

# ifdef CPP_MPI
  ! finally gather PNS_SO_ALL on master in case of MPI run
  allocate(work(LMMAXD,LMMAXD,IRMD,2,NATYP), stat=ierr)
  if(ierr.ne.0) stop 'Error allocating work for MPI comm of PNS_SO_ALL in main1a'
  ihelp = LMMAXD*LMMAXD*IRMD*2*NATYP
  call MPI_ALLREDUCE(PNS_SO_ALL, work, ihelp, MPI_DOUBLE_COMPLEX, MPI_SUM, t_mpi_c_grid%myMPI_comm_ie,ierr)

  if(ierr.ne.MPI_SUCCESS) stop 'Error in MPI comm of PNS_SO_ALL in main1a'
  PNS_SO_ALL(:,:,:,:,:) = work(:,:,:,:,:)
  deallocate(work, stat=ierr)
  if(ierr.ne.0) stop 'Error deallocating work for MPI comm of PNS_SO_ALL in main1a'
# endif

  ! done with preparations, call normcoeff routines that construct operators
  if(myrank==master) WRITE(*,*) 'Computing spin operator'
  CALL NORMCOEFF_SO(t_params%IRCUT,t_params%LMMAXD/2,PNS_SO_ALL,t_params%THETAS,t_params%NTCELL,t_params%IFUNM,t_params%IPAN,t_params%LMSP,t_inc%KVREL,t_params%CLEB,t_params%ICLEB,t_params%IEND,t_params%DRDI,t_params%IRWS,1+KORBIT)
     
  if(myrank==master) WRITE(*,*) 'Computing torq operator'
  CALL NORMCOEFF_SO_TORQ(t_params%IRCUT,t_params%LMMAXD/2,PNS_SO_ALL,t_params%NTCELL,t_params%IFUNM,t_params%IPAN,t_params%LMSP,t_inc%KVREL,t_params%CLEB,t_params%ICLEB,t_params%IEND,t_params%DRDI,t_params%IRWS,t_params%VISP,t_inc%NSPIN,t_params%VINS,t_params%IRMIN)

  if(myrank==master) WRITE(*,*) 'Computing spinflux operator'
  CALL NORMCOEFF_SO_SPINFLUX(t_params%IRCUT,t_params%LMMAXD/2,PNS_SO_ALL,t_inc%KVREL,t_params%DRDI)

  ! deallocate temporary arrays
  deallocate(irws, rmesh, npan_tot, rpan_intervall, ipan_intervall, theta, phi)

end subroutine operators_for_FScode
