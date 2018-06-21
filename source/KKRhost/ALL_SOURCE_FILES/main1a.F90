!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN1A
!> @brief Wrapper module for the calculation of the T-matrix for the JM-KKR package
!> @details The code uses the information obtained in the main0 module, this is
!> mostly done via the get_params_1a() call, that obtains parameters of the type
!> t_params and passes them to local variables
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!> @note
!> - Jonathan Chico Jan. 2018: Removed inc.p dependencies and rewrote to Fortran90
!-------------------------------------------------------------------------------
module MOD_MAIN1A

   use Profiling
   use Constants
   use global_variables
      Use mod_datatypes, Only: dp

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main1a
   !> @brief Main subroutine regarding the calculation of the t-matrix
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !----------------------------------------------------------------------------
   subroutine main1a()

#ifdef CPP_MPI
      use mpi
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif

      use mod_types, only: t_tgmat, t_inc, t_lloyd, t_dtmatJij,init_t_dtmatJij,&
                           init_t_dtmatJij_at
#ifdef CPP_MPI
      use mod_types, only: gather_tmat, gather_lly_dtmat,t_mpi_c_grid,         &
                           save_t_mpi_c_grid,get_ntot_pT_ioff_pT_2D
#endif
      use mod_mympi, only: nranks, master, myrank
#ifdef CPP_MPI
      use mod_mympi, only: find_dims_2d,distribute_linear_on_tasks
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_wunfiles
      use mod_jijhelp, only: set_Jijcalc_flags

      use mod_main0

      ! .. Local variables
      integer :: I1
      integer :: IPOT
      integer :: ILTMP
      integer :: ISPIN
      integer :: ITMPDIR
      character(len=80) :: TMPDIR
      logical :: OPT
      logical :: TEST
      logical :: LREFSYS
      integer :: LRECTMT
      integer :: LRECTRA
      ! .. Local arrays
      real (kind=dp), dimension(NATYPD) :: PHI
      real (kind=dp), dimension(NATYPD) :: THETA
      real (kind=dp), dimension(:,:,:), allocatable :: VINSNEW

#ifdef CPP_MPI
      integer :: ntot1, mytot, ii
      integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT,ntot_all, ioff_all
      ! communication of dtmat in case of lloyd
      integer :: iwork
      complex (kind=dp), dimension(:,:,:,:), allocatable :: work_jij
#endif
      integer :: i1_start, i1_end, ierr,i_stat,i_all

      !-------------------------------------------------------------------------
      !     .. External Subroutines ..
      !-------------------------------------------------------------------------
      external :: TBREF,CALCTMAT,OPT
      !-------------------------------------------------------------------------
      !     .. Intrinsic Functions ..
      !-------------------------------------------------------------------------
      intrinsic :: ATAN,MOD
      !     ..
      !data TOLRDIF /1.5D0/ ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      !data LLY /0/
      !
      LLY=0
      TOLRDIF=1.5d0
      !LRECTMT=WLENGTH*kind(czero)*LMMAXD*LMMAXD
      !LRECTRA=WLENGTH*kind(czero)
      LRECTMT=WLENGTH*4*LMMAXD*LMMAXD
      LRECTRA=WLENGTH*4

      allocate(VINSNEW(NRMAXD,LMPOTD,NSPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VINSNEW))*kind(VINSNEW),'VINSNEW','main1a')
      VINSNEW=0.0D0

      ! Consistency check
      if ( (KREL.lt.0) .or. (KREL.gt.1) ) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      if ( (KREL.eq.1) .and. (NSPIND.eq.2) ) stop ' set NSPIN = 1 for KREL = 1 in the inputcard'
      !-------------------------------------------------------------------------
      ! This routine previously used to read from unformatted files created by
      ! the main0 module, now  instead of unformatted files take parameters from
      ! types defined in wunfiles.F90
      !-------------------------------------------------------------------------
      call get_params_1a(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,NCLSD,NREFD,&
         NCLEB,NEMB,NAEZD,LM2D,NSRA,INS,NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,    &
         NINEQ,IDOLDAU,LLY,KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,    &
         IEND,EZ,VINS,IRMIN,ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,  &
         VISP,RMTREF,VREF,CLEB,CSCL,SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,    &
         SOLVER,TMPDIR,DELTAE,TOLRDIF,NPAN_LOG_AT,NPAN_EQ_AT,NCHEB,NPAN_TOT,     &
         IPAN_INTERVALL,RPAN_INTERVALL,RNEW,NTOTD,NRMAXD,R_LOG,NTLDAU,     &
         ITLDAU,LOPT,VTREL,BTREL,DRDIREL,R2DRDIREL,RMREL,IRMIND,LMPOTD,     &
         NSPOTD,NPOTD,JWSREL,ZREL,ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
      !
      if ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0

      !-------------------------------------------------------------------------
      !                       End read in variables
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! LDA+U treatment
      !-------------------------------------------------------------------------
      if ( IDOLDAU.eq.1 ) then
         open (67,FILE='ldau.unformatted',FORM='unformatted')
         read (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         close(67)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate Coulomb matrix ULDAU it calculates U matrix only once.
         ! Remove the next IF statement to have U calculated for each iteration anew.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!! IF ( ITRUNLDAU.LE.0 ) THEN
         call INITLDAU(NSRA,NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU,VISP,NSPIN, &
            RMESH,DRDI,ZAT,IPAN,IRCUT,PHILDAU,ULDAU)
         !!!!!!!! END IF
      end if
      !-------------------------------------------------------------------------
      ! End of LDA+U setup
      !-------------------------------------------------------------------------
#ifdef CPP_MPI
      ! MPI:
      ntot1 = t_inc%NATYP
#endif
      !-------------------------------------------------------------------------
      ! No need to recalculate the reference system in SCF decimation case
      !-------------------------------------------------------------------------
      ! ITSCF is initialised to 0 in main0
      LREFSYS = .TRUE.
      if (  OPT('DECIMATE').and.(ITSCF.gt.0) ) LREFSYS = .FALSE.
      if (  OPT('rigid-ef').and.(ITSCF.gt.0) ) LREFSYS = .FALSE.
      if ( TEST('no-neutr').and.(ITSCF.gt.0) ) LREFSYS = .FALSE.
      if (  OPT('no-neutr').and.(ITSCF.gt.0) ) LREFSYS = .FALSE.
      if ( TEST('lrefsysf').or.OPT('lrefsysf') ) LREFSYS = .FALSE.
      !

      if (t_tgmat%tmat_to_file) then
         call OPENDAFILE(69,'tmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if

      if (LLY.ne.0) then
         if(t_lloyd%dtmat_to_file) then
            call OPENDAFILE(691,'dtmatde',7,LRECTMT,TMPDIR,ITMPDIR,ILTMP) ! LLY
         end if
         if(t_lloyd%tralpha_to_file) then
            call OPENDAFILE(692,'tralpha',7,LRECTRA,TMPDIR,ITMPDIR,ILTMP) ! LLY
         end if
      endif
      !

#ifdef CPP_MPI
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,  &
         t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master, &
         ntot1,ntot_pT,ioff_pT,.true.,.true.)

      i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
      i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)

      if (.not. (allocated(t_mpi_c_grid%ntot_pT1) .and.  &
         allocated(t_mpi_c_grid%ioff_pT1))) then
         allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
         call memocc(i_stat,product(shape(t_mpi_c_grid%ntot_pT1))*kind(t_mpi_c_grid%ntot_pT1),&
            't_mpi_c_grid%ntot_pT1','main1a')
         allocate(t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_ie-1),stat=i_stat)
         call memocc(i_stat,product(shape(t_mpi_c_grid%ioff_pT1))*kind(t_mpi_c_grid%ioff_pT1),&
            't_mpi_c_grid%ioff_pT1','main1a')
      endif
      t_mpi_c_grid%ntot_pT1 = ntot_pT
      t_mpi_c_grid%ioff_pT1 = ioff_pT
#else
      i1_start = 1
      i1_end   = NATYP
#endif
     i1 = 1
     ipot = 1 
     write(2222+myrank,*) 'line1',ICST,INS,IELAST,NSRA,NSPIN,I1,EZ,DRDI(1,I1)
     write(2222+myrank,*) 'line2',RMESH(1,I1),VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH))
     write(2222+myrank,*) 'line3',VISP(1,IPOT),ZAT(I1),IRMIN(I1),IPAN(I1),IRCUT(0,I1),CLEB
     write(2222+myrank,*) 'line4',LOFLM,ICLEB,IEND,SOLVER,SOCSCL(1,KREL*I1+(1-KREL))
     write(2222+myrank,*) 'line5',CSCL(1,KREL*I1+(1-KREL)),VTREL(1,I1),BTREL(1,I1),RMREL(1,I1)
     write(2222+myrank,*) 'line6',DRDIREL(1,I1),R2DRDIREL(1,I1),ZREL(I1),JWSREL(I1),IDOLDAU
     write(2222+myrank,*) 'line7',LOPT(I1),WLDAU(1,1,1,I1),LLY,DELTAE
#ifdef CPP_MPI
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

      !skip this part with GREENIMP option
      if(opt('GREENIMP') .or. TEST('IMP_ONLY')) then
         if(myrank==master) write(*,*) 'Skipping atom loop in main1a'
         i1_start = 1
         i1_end = 0
      end if

      if (.not.OPT('NEWSOSOL')) then
         do I1 = i1_start,i1_end
            do ISPIN = 1,NSPIN
               IPOT=NSPIN*(I1-1)+ispin
               !
               call CALCTMAT(ICST,INS,IELAST,NSRA,ISPIN,NSPIN,I1,EZ,DRDI(1,I1),  &
                  RMESH(1,I1),VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),             &
                  VISP(1,IPOT),ZAT(I1),IRMIN(I1),IPAN(I1),IRCUT(0,I1),CLEB,      &
                  LOFLM,ICLEB,IEND,SOLVER,SOCSCL(1,KREL*I1+(1-KREL)),            &
                  CSCL(1,KREL*I1+(1-KREL)),VTREL(1,I1),BTREL(1,I1),RMREL(1,I1),  &
                  DRDIREL(1,I1),R2DRDIREL(1,I1),ZREL(I1),JWSREL(I1),IDOLDAU,     &
                  LOPT(I1),WLDAU(1,1,1,I1),LLY,DELTAE) ! LLY
               !
            end do
         end do

      else
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! For calculation of Jij-tensor: create array for additional t-matrices and
         ! set atom-dependent flags which indicate if t-matrix is needed
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call init_t_dtmatJij(t_inc,t_dtmatJij)
         if(OPT('XCPL    '))then
            call set_Jijcalc_flags(t_dtmatJij,NATYP,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
         end if !OPT('XCPL')

         ! nonco angles: defined in mod_wunfiles
         call read_angles(t_params,NATYP,THETA,PHI)

         ! Interpolate potential
         call INTERPOLATE_POTEN(LPOT,IRMD,IRNSD,NATYP,IPAND,LMPOTD,NSPOTD,NTOTD,   &
            NCHEB,NTOTD*(NCHEB+1),NSPIN,RMESH,IRMIND,IRMIN,IRWS,IRCUT,VINS,VISP,   &
            NPAN_LOG_AT,NPAN_EQ_AT,NPAN_TOT,RNEW,IPAN_INTERVALL,VINSNEW)

         do I1=i1_start,i1_end

            IPOT=NSPIN*(I1-1)+1

            call TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),SOCSCALE(I1),EZ, &
               NSRA,CLEB(:,1),ICLEB,IEND,NCHEB,NPAN_TOT(I1),               &
               RPAN_INTERVALL(:,I1),IPAN_INTERVALL(:,I1),RNEW(:,I1),       &
               VINSNEW,THETA(I1),PHI(I1),I1,IPOT,LMPOTD,LLY,DELTAE,IDOLDAU, &
               LOPT(I1),WLDAU(:,:,:,I1),t_dtmatJij(I1))

         enddo !I1, atom loop

      endif !NEWSOSOL
      !
      if ( IDOLDAU.eq.1 ) then
         open (67,FILE='ldau.unformatted',FORM='unformatted')
         write (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         close(67)
      end if
      !
      close (69)

      if (LLY.ne.0) then
         if(t_lloyd%dtmat_to_file) close(691)
         if(t_lloyd%tralpha_to_file) close(692)
      endif


#ifdef CPP_MPI
      !skip this part with GREENIMP option
      if(.not.(opt('GREENIMP').or.TEST('IMP_ONLY'))) then

         if(.not.t_tgmat%tmat_to_file) then
            do ii=0,t_mpi_c_grid%nranks_ie-1
               ntot_all(ii) = t_mpi_c_grid%ntot_pT1(ii)
               ioff_all(ii) = t_mpi_c_grid%ioff_pT1(ii)
            end do
            mytot = t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
            call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,mytot,&
               t_mpi_c_grid%mympi_comm_ie, t_mpi_c_grid%nranks_ie)
         end if
         
         if(lly/=0 .and. .not.t_lloyd%dtmat_to_file) then
            call gather_lly_dtmat(t_mpi_c_grid, t_lloyd, lmmaxd,t_mpi_c_grid%mympi_comm_ie)
         endif
         
         !-------------------------------------------------------------------------
         !for calculation of Jij-tensor
         !-------------------------------------------------------------------------
         if(OPT('XCPL    ').and.OPT('NEWSOSOL'))then
            do I1=1,t_inc%NATYP
               !initialize t_dtmatJij on other tasks
               !   t_dtmatJij was already allocated for certain atoms within the atom loop
               !   (in tmat_newsolver). This initialization cannot be made before tmat_newsolver,
               !   because the division of the enegry loop (done in there) influences t_dtmatJij.
               call init_t_dtmatJij_at(t_inc,t_mpi_c_grid,t_dtmatJij(I1))
         
               !communicate
               if(t_dtmatJij(I1)%calculate)then
                  iwork = product(shape(t_dtmatJij(I1)%dtmat_xyz))
                  allocate(work_jij(iwork,1,1,1),stat=i_stat)
                  call memocc(i_stat,product(shape(work_jij))*kind(work_jij),'work_jij','main1a')
         
                  call MPI_ALLREDUCE(t_dtmatJij(I1)%dtmat_xyz,work_jij, &
                     iwork,MPI_DOUBLE_COMPLEX,MPI_SUM,                  &
                     t_mpi_c_grid%mympi_comm_ie,IERR)
                  if(ierr/=MPI_SUCCESS) stop 'error communicating t_dtmatJij'
                  call ZCOPY(iwork,work_jij,1,t_dtmatJij(I1)%dtmat_xyz,1)
                  i_all=-product(shape(work_jij))*kind(work_jij)
                  deallocate(work_jij, stat=i_stat)
                  call memocc(i_stat,i_all,'work_jij','main1a')
               end if !t_dtmatJij(I1)%calculate
         
            end do !I1=1,t_inc%NATYP
         
         end if !OPT('XCPL    ').and.OPT('NEWSOSOL')

      end if ! .not.opt('GREENIMP')
      ! end skip this part with GREENIMP option
#endif
      !-------------------------------------------------------------------------
      ! End of calculation of Jij-tensor
      !-------------------------------------------------------------------------


#ifdef CPP_TIMING
      call timing_start('main1a - tbref')
#endif
      if ( LREFSYS ) then
         call TBREF(EZ,IELAST,ALAT,VREF,IEND,LMAX,NCLS,NINEQ,NREF,CLEB,RCLS,  &
            ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,RMTREF,TOLRDIF,TMPDIR,ITMPDIR,  &
            ILTMP,NAEZ,LLY) ! LLY Lloyd
      endif
#ifdef CPP_TIMING
      call timing_stop('main1a - tbref')
#endif

      if(t_inc%i_write>0) write (1337,'(79("="),/,30X,"< KKR1a finished >",/,79("="),/)')

      ! Deallocate leftover arrays
      if (allocated(VINSNEW)) then
         i_all=-product(shape(VINSNEW))*kind(VINSNEW)
         deallocate(VINSNEW,stat=i_stat)
         call memocc(i_stat,i_all,'VINSNEW','main1a')
      endif

   end subroutine main1a

end module MOD_MAIN1A
