!-------------------------------------------------------------------------------
! MODULE: MOD_MAIN1A
!> @brief Wrapper module for the calculation of the T-matrix for the JM-KKR package
!> @details The code uses the information obtained in the main0 module, this is
!> mostly done via the get_params_1a() call, that obtains parameters of the type
!> t_params and passes them to local variables
!> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
!> and many others ...
!-------------------------------------------------------------------------------
module MOD_MAIN1A

   use Profiling

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

      ! *********************************************************************
      ! * For KREL = 1 (relativistic mode)                                  *
      ! *                                                                   *
      ! *  NPOTD = 2 * NATYPD                                               *
      ! *  LMMAXD = 2 * (LMAXD+1)^2                                         *
      ! *  NSPIND = 1                                                       *
      ! *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
      ! *          function, set up in the spin-independent non-relativstic *
      ! *          (l,m_l)-representation                                   *
      ! *                                                                   *
      ! *********************************************************************
      !     .. Parameters ..
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER MMAXD
      PARAMETER (MMAXD = 2*LMAXD+1)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) +(1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NRMAXD
      PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
      INTEGER LRECTMT
      PARAMETER (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
      INTEGER LRECTRA
      PARAMETER (LRECTRA=WLENGTH*4)
      !     ..
      !     .. Local Scalars ..
      integer :: I1
      integer :: INS
      integer :: LLY ! LLY <> 0: apply Lloyds formula
      integer :: ICST
      integer :: IEND
      integer :: IPOT
      integer :: NCLS
      integer :: LMAX
      integer :: NREF
      integer :: NSRA
      integer :: NAEZ
      integer :: NATYP
      integer :: ISPIN
      integer :: ITSCF
      integer :: NINEQ
      integer :: NSPIN
      integer :: NCHEB
      integer :: IELAST
      double precision :: PI
      double precision :: ALAT
      double precision :: R_LOG
      double precision :: TOLRDIF
      double complex :: DELTAE  ! Energy difference for numerical derivative
      integer, dimension(NATYP) :: NPAN_EQ
      integer, dimension(NATYP) :: NPAN_LOG
      integer, dimension(NATYP) :: NPAN_TOT
      integer, dimension(0:NTOTD,NATYP) :: IPAN_INTERVALL
      double precision, dimension(NATYP) :: PHI
      double precision, dimension(NATYP) :: THETA
      double precision, dimension(NRMAXD,NATYP) :: RNEW
      double precision, dimension(0:NTOTD,NATYP) :: RPAN_INTERVALL
      double precision, dimension(:,:,:), allocatable :: VINSNEW
      !     ..
      !     .. Local Arrays ..
      integer, dimension(NAEZ+NEMB) :: CLS
      integer, dimension(NATYP) :: IPAN
      integer, dimension(NATYP) :: IRWS
      integer, dimension(NATYP) :: IRMIN
      integer, dimension(LM2D) :: LOFLM
      integer, dimension(NCLSD) :: NACLS
      integer, dimension(NAEZ+NEMB) :: REFPOT
      integer, dimension(NACLSD,NAEZ+NEMB) :: ATOM
      integer, dimension(NCLEB,4) :: ICLEB
      integer, dimension(0:IPAND,NATYP) :: IRCUT
      double precision, dimension(NATYP) :: ZAT
      double precision, dimension(NREF) :: VREF
      double precision, dimension(NREF) :: RMTREF
      double precision, dimension(NCLEB,2) :: CLEB
      double precision, dimension(IRM,NPOTD) :: VISP
      double precision, dimension(IRM,NATYP) :: DRDI
      double precision, dimension(IRM,NATYP) :: RMESH
      double precision, dimension(3,NACLSD,NCLSD) :: RCLS
      double precision, dimension(:,:,:), allocatable :: VINS
      double complex, dimension(IEMXD) :: EZ
      !-------------------------------------------------------------------------
      !     RELATIVISTIC MODE
      !-------------------------------------------------------------------------
      integer :: ILTMP
      integer :: ITMPDIR
      character(len=10) :: SOLVER
      character(len=80) :: TMPDIR
      logical :: OPT
      logical :: TEST
      logical :: LREFSYS
      integer, dimension(NATYP) :: ZREL
      integer, dimension(NATYP) :: JWSREL
      double precision, dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)) :: SOCSCL
      double precision, dimension(NATYP) :: SOCSCALE
      double precision, dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)) :: CSCL
      double precision, dimension(IRM*KREL+(1-KREL),NATYP) :: VTREL
      double precision, dimension(IRM*KREL+(1-KREL),NATYP) :: BTREL
      double precision, dimension(IRM*KREL+(1-KREL),NATYP) :: RMREL
      double precision, dimension(IRM*KREL+(1-KREL),NATYP) :: DRDIREL
      double precision, dimension(IRM*KREL+(1-KREL),NATYP) :: R2DRDIREL
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      integer :: NTLDAU
      integer :: IDOLDAU
      integer :: ITRUNLDAU
      integer, dimension(NATYP) :: LOPT
      integer, dimension(NATYP) :: ITLDAU
      double precision, dimension(NATYP) :: UEFF
      double precision, dimension(NATYP) :: JEFF
      double precision, dimension(NATYP) :: EREFLDAU
      double precision, dimension(MMAXD,MMAXD,NSPIND,NATYP) :: WLDAU
      double precision, dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYP) :: ULDAU
      double complex, dimension(IRM,NATYP) :: PHILDAU
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      integer :: NATOMIMP
      integer, dimension(NATYP) :: IQAT
      integer, dimension(NATOMIMPD) :: ATOMIMP

      ! Assignment of values to parameters
      PARAMETER  (PI=4.d0*datan(1.d0))

#ifdef CPP_MPI
      integer :: ntot1, mytot, ii
      integer, dimension(0:nranks-1) :: ntot_pT, ioff_pT,ntot_all, ioff_all
      ! communication of dtmat in case of lloyd
      integer :: iwork
      double complex, dimension(:,:,:,:), allocatable :: work_jij
#endif
      integer :: i1_start, i1_end, ierr,i_stat,i_all

      !-------------------------------------------------------------------------
      !     .. External Subroutines ..
      !-------------------------------------------------------------------------
      EXTERNAL TBREF,CALCTMAT,OPT
      !-------------------------------------------------------------------------
      !     .. Intrinsic Functions ..
      !-------------------------------------------------------------------------
      INTRINSIC ATAN,MOD
      !     ..
      DATA TOLRDIF /1.5D0/ ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      DATA LLY /0/
      !

      allocate(VINSNEW(NRMAXD,LMPOTD,NSPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VINSNEW))*kind(VINSNEW),'VINSNEW','main1a')
      allocate(VINS(IRMIND:IRM,LMPOTD,NSPOTD),stat=i_stat)
      call memocc(i_stat,product(shape(VINS))*kind(VINS),'VINS','main1a')

      ! Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) ) STOP ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) STOP ' set NSPIN = 1 for KREL = 1 in the inputcard'
      !-------------------------------------------------------------------------
      ! This routine previously used to read from unformatted files created by
      ! the main0 module, now  instead of unformatted files take parameters from
      ! types defined in wunfiles.F90
      !-------------------------------------------------------------------------
      call get_params_1a(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,   &
         NCLSD,NREFD,NCLEB,NEMBD,NAEZD,LM2D,NSRA,INS,NAEZ,NATYP,     &
         NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY,     &
         KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,IEND,EZ,VINS,   &
         IRMIN,ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,VISP,    &
         RMTREF,VREF,CLEB,CSCL,SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,   &
         SOLVER,TMPDIR,DELTAE,tolrdif,NPAN_LOG,NPAN_EQ,              &
         NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,LMAXD,    &
         NTOTD,NRMAXD,R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,  &
         R2DRDIREL,RMREL,IRMIND,LMPOTD,NSPOTD,NPOTD,JWSREL,ZREL,     &
         ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
      !
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0

      !-------------------------------------------------------------------------
      !                       End read in variables
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      ! LDA+U treatment
      !-------------------------------------------------------------------------
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         READ (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! Calculate Coulomb matrix ULDAU it calculates U matrix only once.
         ! Remove the next IF statement to have U calculated for each iteration anew.
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         !!!!!!!! IF ( ITRUNLDAU.LE.0 ) THEN
         CALL INITLDAU(NSRA,NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU,&
            VISP,NSPIN,RMESH,DRDI,ZAT,IPAN,IRCUT,                 &
            PHILDAU,ULDAU)
         !!!!!!!! END IF
      END IF
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
      IF (  OPT('DECIMATE').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF (  OPT('rigid-ef').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF ( TEST('no-neutr').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF (  OPT('no-neutr').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF ( TEST('lrefsysf').OR.OPT('lrefsysf') ) LREFSYS = .FALSE.
      !

      if (t_tgmat%tmat_to_file) then
         CALL OPENDAFILE(69,'tmat',4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)
      end if

      IF (LLY.NE.0) THEN
         if(t_lloyd%dtmat_to_file) then
            CALL OPENDAFILE(691,'dtmatde',7,LRECTMT,TMPDIR,ITMPDIR,ILTMP) ! LLY
         end if
         if(t_lloyd%tralpha_to_file) then
            CALL OPENDAFILE(692,'tralpha',7,LRECTRA,TMPDIR,ITMPDIR,ILTMP) ! LLY
         end if
      ENDIF
      !

#ifdef CPP_MPI
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,  &
         t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master, &
         ntot1,ntot_pT,ioff_pT,.true.,.true.)

      i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
      i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+ntot_pT(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)

      if (.not. (allocated(t_mpi_c_grid%ntot_pT1) .and.   &
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

      !skip this part with GREENIMP option
      if(opt('GREENIMP')) then
         if(myrank==master) write(*,*) 'Skipping atom loop in main1a'
         i1_start = 1
         i1_end = 1
      end if

      IF (.NOT.OPT('NEWSOSOL')) THEN
         DO I1 = i1_start,i1_end
            DO ISPIN = 1,NSPIN
               IPOT=NSPIN*(I1-1)+ispin
               !
               CALL CALCTMAT(ICST,INS,IELAST,                  &
                     NSRA,ISPIN,NSPIN,                         &
                     I1,EZ,                                    &
                     DRDI(1,I1),RMESH(1,I1),                   &
                     VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),    &
                     VISP(1,IPOT),ZAT(I1),IRMIN(I1),IPAN(I1),  &   ! Added IRMIN 1.7.2014
                     IRCUT(0,I1),CLEB,LOFLM,ICLEB,IEND,SOLVER, &
                     SOCSCL(1,KREL*I1+(1-KREL)),               &
                     CSCL(1,KREL*I1+(1-KREL)),                 &
                     VTREL(1,I1),BTREL(1,I1),                  &
                     RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),&
                     ZREL(I1),JWSREL(I1),                      &
                     IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1),         &
                     LLY,DELTAE) ! LLY
               !
            END DO
         END DO

      ELSE
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! For calculation of Jij-tensor: create array for additional t-matrices and
         ! set atom-dependent flags which indicate if t-matrix is needed
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         call init_t_dtmatJij(t_inc,t_dtmatJij)
         IF(OPT('XCPL    '))THEN
            call set_Jijcalc_flags(t_dtmatJij,NATYPD,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
         END IF!OPT('XCPL')

         ! nonco angles: defined in mod_wunfiles
         call read_angles(t_params,NATYP,THETA,PHI)

         ! Interpolate potential
         CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYPD,IPAND, &
            NSPOTD,NTOTD,NCHEBD,NTOTD*(NCHEBD+1),              &
            NSPIN,RMESH,IRMIN,IRWS,IRCUT,VINS,                 &
            VISP,NPAN_LOG,NPAN_EQ,NPAN_TOT,                    &
            RNEW,IPAN_INTERVALL,                               &
            VINSNEW)

         DO I1=i1_start,i1_end

            IPOT=NSPIN*(I1-1)+1

            CALL TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),  &
               SOCSCALE(I1),EZ,NSRA,CLEB(1,1),ICLEB,IEND,   &
               NCHEB,NPAN_TOT(I1),                          &
               RPAN_INTERVALL(0,I1),IPAN_INTERVALL(0,I1),   &
               RNEW(1,I1),VINSNEW,THETA(I1),PHI(I1),I1,IPOT,&
               LLY,DELTAE,                                  &  !LLY
               IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1),            &  ! LDAU
               t_dtmatJij(I1))

         ENDDO !I1, atom loop

         !      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         !      TTTTTTTTTTTTT TESTOUTPUT   Dij-implementation TTTTTTTTTTTTTTTTTT
         !      DO I1=i1_start,i1_end
         !       ie_start=t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         !       ie_end  =t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
         !       DO ie_num=1,ie_end
         !        IE = ie_start+ie_num
         !        DO II=1,3
         !       write(filename,'(A,I3.3,A,I3.3,A,I1,A)') 'test_dtmat_iat=',I1,
         !    +                                                  '_ie=',IE,
         !    +                                                  'xyz=',II
         !       write(*,*) 'writing testfile:', trim(filename)
         !       if(allocated(t_dtmatJij(I1)%dtmat_xyz))then
         !       open(unit=13626,file=trim(filename), action='write',
         !    +       form='formatted')
         !      write(13626,'(2ES25.16)') t_dtmatJij(I1)%dtmat_xyz(:,:,II,ie_num)
         !       close(13626)
         !       end if
         !        ENDDO!II
         !       ENDDO!IE
         !      ENDDO !I1, atom loop
         !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      ENDIF !NEWSOSOL
      !
      IF ( IDOLDAU.EQ.1 ) THEN
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
      !
      CLOSE (69)

      IF (LLY.NE.0) THEN
         if(t_lloyd%dtmat_to_file) CLOSE(691)
         if(t_lloyd%tralpha_to_file) CLOSE(692)
      ENDIF


#ifdef CPP_MPI
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
            end if!t_dtmatJij(I1)%calculate

         end do!I1=1,t_inc%NATYP
         !      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
         !      TTTTTTTTTTTTT TESTOUTPUT   Dij-implementation TTTTTTTTTTTTTTTTTT
         !      DO I1=1,t_inc%NATYP
         !       ie_start=t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
         !       ie_end  =t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
         !       DO ie_num=1,ie_end
         !        IE = ie_start+ie_num
         !        DO II=1,3
         !       write(filename,'(A,I3.3,A,I3.3,A,I1,A,I2)') 'test_dtmat_iat=',I1
         !    +                                                  ,'_ie=',IE,
         !    +                                                  'xyz=',II,
         !    +                                               'myrank=',myrank
         !       write(*,*) 'writing testfile:', trim(filename)
         !       if(allocated(t_dtmatJij(I1)%dtmat_xyz))then
         !       open(unit=13626,file=trim(filename), action='write',
         !    +       form='formatted')
         !      write(13626,'(2ES25.16)') t_dtmatJij(I1)%dtmat_xyz(:,:,II,ie_num)
         !       close(13626)
         !       end if
         !        ENDDO!II
         !       ENDDO!IE
         !      ENDDO !I1, atom loop
         !      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
         !      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      end if!OPT('XCPL    ').and.OPT('NEWSOSOL')
      !-------------------------------------------------------------------------
      ! End of calculation of Jij-tensor
      !-------------------------------------------------------------------------
#endif


#ifdef CPP_TIMING
      call timing_start('main1a - tbref')
#endif
      if ( LREFSYS ) then
         call TBREF(EZ,IELAST,ALAT,VREF,IEND,LMAX,NCLS,  &
            NINEQ,NREF,CLEB,RCLS,ATOM,CLS,ICLEB,         &
            LOFLM,NACLS,REFPOT,RMTREF,TOLRDIF,           &
            TMPDIR,ITMPDIR,ILTMP,                        &
            NAEZ,LLY) ! LLY Lloyd
      endif
#ifdef CPP_TIMING
      call timing_stop('main1a - tbref')
#endif

      if(t_inc%i_write>0) write (1337,'(79(1H=),/,30X,"< KKR1a finished >",/,79(1H=),/)')

      ! Deallocate leftover arrays
      !> @note JC: These arrays are allocated here, but referenced in the main0
      !> and then in the other modules, why?
      if (allocated(VINS)) then
         i_all=-product(shape(VINS))*kind(VINS)
         deallocate(VINS,stat=i_stat)
         call memocc(i_stat,i_all,'VINS','main1a')
      endif
      if (allocated(VINSNEW)) then
         i_all=-product(shape(VINSNEW))*kind(VINSNEW)
         deallocate(VINSNEW,stat=i_stat)
         call memocc(i_stat,i_all,'VINSNEW','main1a')
      endif

   end subroutine main1a

end module MOD_MAIN1A
