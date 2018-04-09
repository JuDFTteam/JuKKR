      module MOD_MAIN1A
            
      implicit none
      
      contains
      
      
      subroutine main1a()
      
#ifdef CPP_MPI
      use mpi
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif
      
      use mod_types, only: t_tgmat, t_inc, t_lloyd, t_dtmatJij,
     &                     init_t_dtmatJij, init_t_dtmatJij_at
#ifdef CPP_MPI
     &                     ,gather_tmat, gather_lly_dtmat,
     &                      t_mpi_c_grid, save_t_mpi_c_grid,
     &                      get_ntot_pT_ioff_pT_2D 
#endif
      use mod_mympi, only: nranks, master, myrank
#ifdef CPP_MPI
     &                     ,find_dims_2d,
     &                     distribute_linear_on_tasks
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_wunfiles
      use mod_jijhelp, only: set_Jijcalc_flags
      use mod_save_wavefun, only: t_wavefunctions, read_wavefunc
     
      INCLUDE 'inc.p'
C
C *********************************************************************
C * For KREL = 1 (relativistic mode)                                  *
C *                                                                   *
C *  NPOTD = 2 * NATYPD                                               *
C *  LMMAXD = 2 * (LMAXD+1)^2                                         *
C *  NSPIND = 1                                                       *
C *  LMGF0D = (LMAXD+1)^2 dimension of the reference system Green     *
C *          function, set up in the spin-independent non-relativstic *
C *          (l,m_l)-representation                                   *
C *                                                                   *
C *********************************************************************
C
C     .. Parameters ..
C
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      INTEGER MMAXD
      PARAMETER (MMAXD = 2*LMAXD+1)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NRMAXD
      PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
      INTEGER LRECTMT
      PARAMETER (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
      INTEGER LRECTRA
      PARAMETER (LRECTRA=WLENGTH*4)
      DOUBLE COMPLEX, PARAMETER :: CZERO=(0.0d0, 0.0d0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,PI
      PARAMETER  (PI=4.d0*datan(1.d0))
      INTEGER I1,ICST,IELAST,IEND,INS,IPOT,ISPIN,ITSCF,LMAX,NATYP,NAEZ, ! LLY added NAEZ
     +        NCLS,NINEQ,NREF,NSPIN,NSRA!,IE
      INTEGER LLY ! LLY <> 0: apply Lloyds formula
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION R_LOG,THETA(NATYPD),PHI(NATYPD),TOLRDIF
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     +                 RNEW(NRMAXD,NATYPD)!,
!      +                 VINSNEW(NRMAXD,LMPOTD,NSPOTD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
      double precision, allocatable :: vinsnew(:,:,:)
C     ..
C     .. Local Arrays ..
      DOUBLE COMPLEX EZ(IEMXD)
      DOUBLE PRECISION DRDI(IRMD,NATYPD),RMESH(IRMD,NATYPD),ZAT(NATYPD),
     +                 RCLS(3,NACLSD,NCLSD),RMTREF(NREFD),VREF(NREFD),
     +                 VISP(IRMD,NPOTD),
     +                 CLEB(NCLEB,2)
      double precision, allocatable :: vins(:,:,:)
      INTEGER ATOM(NACLSD,NAEZD+NEMBD),CLS(NAEZD+NEMBD),ICLEB(NCLEB,4),
     +        IPAN(NATYPD),IRCUT(0:IPAND,NATYPD),
     +        LOFLM(LM2D),NACLS(NCLSD),REFPOT(NAEZD+NEMBD),
     +        IRWS(NATYPD),IRMIN(NATYPD)
C-----------------------------------------------------------------------
C     RELATIVISTIC MODE 
C
      CHARACTER*10 SOLVER
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION SOCSCALE(NATYPD)
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 RMREL(IRMD*KREL+(1-KREL),NATYPD)
      INTEGER JWSREL(NATYPD),ZREL(NATYPD)
      INTEGER ITMPDIR,ILTMP
      CHARACTER*80 TMPDIR!, filename
      LOGICAL LREFSYS,OPT,TEST
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
      INTEGER NATOMIMP, ATOMIMP(NATOMIMPD), IQAT(NATYPD)

      ! read in wavefunctions 
      logical :: rll_was_read_in, sll_was_read_in,  
     &           rllleft_was_read_in, sllleft_was_read_in
      ! for normcoeff (OPERATOR option) stuff
      DOUBLE COMPLEX, ALLOCATABLE :: RLL(:,:,:,:),SLL(:,:,:,:),
     &                               RLLLEFT(:,:,:,:),SLLLEFT(:,:,:,:),
     &                               RLLTEMP(:,:), PNSTEMP(:,:),
     &                               PNS_SO(:,:,:,:),
     &                               PNS_SO_ALL(:,:,:,:,:) 
      INTEGER :: ie, ie_start, ie_end, ie_num, lm1, lm2, ir

#ifdef CPP_MPI
      integer :: ntot1, mytot, ii
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),
     &           ntot_all(0:nranks-1), ioff_all(0:nranks-1)
      ! communication of dtmat in case of lloyd
      integer :: iwork
      DOUBLE COMPLEX, allocatable :: work_jij(:,:,:,:)
      ! communcate PNS_SO_ALL for OPERATOR option
      integer :: ihelp
      double complex, allocatable :: work(:,:,:,:,:)
#endif
      integer :: i1_start, i1_end, ierr

C-----------------------------------------------------------------------
C     ..
C     .. External Subroutines ..
      EXTERNAL TBREF,CALCTMAT,OPT      
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,MOD
C     ..
      DATA TOLRDIF /1.5D0/ ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      DATA LLY /0/
C


      allocate(VINSNEW(NRMAXD,LMPOTD,NSPOTD),
     &         VINS(IRMIND:IRMD,LMPOTD,NSPOTD))
     
Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C
C
C ======================================================================
C =             read in variables from unformatted files               =
C               > instead of unformatted files take parameters from 
C                 type defined in wunfiles.f
C ======================================================================
C
      call get_params_1a(t_params,IPAND,NATYPD,IRMD,NACLSD,IELAST,
     &       NCLSD,NREFD,NCLEB,NEMBD,NAEZD,LM2D,NSRA,INS,NAEZ,NATYP,
     &       NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY,
     &       KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,IEND,EZ,VINS,
     &       IRMIN,ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,VISP,
     &       RMTREF,VREF,CLEB,CSCL,SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,
     &       SOLVER,TMPDIR,DELTAE,tolrdif,NPAN_LOG,NPAN_EQ,
     &       NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,LMAXD,
     &       NTOTD,NRMAXD,R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,
     &       R2DRDIREL,RMREL,IRMIND,LMPOTD,NSPOTD,NPOTD,JWSREL,ZREL,
     &       ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
C
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
      
C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U treatment
C
      IF ( IDOLDAU.EQ.1 ) THEN 
C
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         READ (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
C
C -> Calculate Coulomb matrix ULDAU
C    it calculates U matrix only once. Remove the next IF statement 
C    to have U calculated for each iteration anew.
C

c        IF ( ITRUNLDAU.LE.0 ) THEN
            CALL INITLDAU(NSRA,NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU,
     &                    VISP,NSPIN,RMESH,DRDI,ZAT,IPAN,IRCUT,
     &                    PHILDAU,ULDAU)
c        END IF
      END IF
C
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
!
#ifdef CPP_MPI
! MPI:     
      ntot1 = t_inc%NATYP
#endif
C     
C -> no need to recalculate the reference system in SCF decimation case
C
c ITSCF is initialised to 0 in main0
      LREFSYS = .TRUE.
      IF (  OPT('DECIMATE').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF (  OPT('rigid-ef').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF ( TEST('no-neutr').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF (  OPT('no-neutr').AND.(ITSCF.GT.0) ) LREFSYS = .FALSE.
      IF ( TEST('lrefsysf').OR.OPT('lrefsysf') ) LREFSYS = .FALSE.
C
        
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
C
      
#ifdef CPP_MPI
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
     &                 t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                           master,ntot1,ntot_pT,ioff_pT,.true.,
     + .true.)

      i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
      i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+
     &                ntot_pT(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)

      if(.not. (allocated(t_mpi_c_grid%ntot_pT1) .and. 
     &            allocated(t_mpi_c_grid%ioff_pT1))) 
     &     allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_ie-1),
     &              t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_ie-1))
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
C
            CALL CALCTMAT(ICST,INS,IELAST,
     +                 NSRA,ISPIN,NSPIN,
     +                 I1,EZ,
     +                 DRDI(1,I1),RMESH(1,I1),
     +                 VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),
     +                 VISP(1,IPOT),ZAT(I1),IRMIN(I1),IPAN(I1),   ! Added IRMIN 1.7.2014
     +                 IRCUT(0,I1),CLEB,LOFLM,ICLEB,IEND,SOLVER,
     +                 SOCSCL(1,KREL*I1+(1-KREL)),
     +                 CSCL(1,KREL*I1+(1-KREL)),
     +                 VTREL(1,I1),BTREL(1,I1),
     +                 RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
     +                 ZREL(I1),JWSREL(I1),
     +                 IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1),
     &                 LLY,DELTAE) ! LLY
C
          END DO
        END DO
       
     
      ELSE
       
       !for calculation of Jij-tensor: create array for additional t-matrices and
       !set atom-dependent flags which indicate if t-matrix is needed
       call init_t_dtmatJij(t_inc,t_dtmatJij) 
       IF(OPT('XCPL    '))THEN
        call set_Jijcalc_flags(t_dtmatJij,NATYPD,NATOMIMPD,NATOMIMP,
     +    ATOMIMP,IQAT)
       END IF!OPT('XCPL')

       ! nonco angles: defined in mod_wunfiles
       call read_angles(t_params,NATYP,THETA,PHI)

c interpolate potential
       CALL INTERPOLATE_POTEN(LPOTD,IRMD,IRNSD,NATYPD,IPAND,
     +                        NSPOTD,NTOTD,NCHEBD,NTOTD*(NCHEBD+1),
     +                        NSPIN,RMESH,IRMIN,IRWS,IRCUT,VINS,
     +                        VISP,NPAN_LOG,NPAN_EQ,NPAN_TOT,
     +                        RNEW,IPAN_INTERVALL,
     +                        VINSNEW)

       DO I1=i1_start,i1_end

        IPOT=NSPIN*(I1-1)+1

        CALL TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),
     +                      SOCSCALE(I1),EZ,NSRA,CLEB(1,1),ICLEB,IEND,
     &                      NCHEB,NPAN_TOT(I1),
     +                      RPAN_INTERVALL(0,I1),IPAN_INTERVALL(0,I1),
     +                      RNEW(1,I1),VINSNEW,THETA(I1),PHI(I1),I1,IPOT
     &                      ,LLY,DELTAE, ! LLY
     &                      IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1), ! LDAU
     &                      t_dtmatJij(I1))

       ENDDO !I1, atom loop


c      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
c      TTTTTTTTTTTTT TESTOUTPUT   Dij-implementation TTTTTTTTTTTTTTTTTT
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
c      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      ENDIF !NEWSOSOL
C
      IF ( IDOLDAU.EQ.1 ) THEN 
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      END IF
C
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
         call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
     &        mytot,t_mpi_c_grid%mympi_comm_ie, t_mpi_c_grid%nranks_ie)
      end if
      
      if(lly/=0 .and. .not.t_lloyd%dtmat_to_file) then
         call gather_lly_dtmat(t_mpi_c_grid, t_lloyd, lmmaxd, 
     &                         t_mpi_c_grid%mympi_comm_ie)
      endif
               
      !for calculation of Jij-tensor
      if(OPT('XCPL    ').and.OPT('NEWSOSOL'))then
        DO I1=1,t_inc%NATYP
          !initialize t_dtmatJij on other tasks
          !   t_dtmatJij was already allocated for certain atoms within the atom loop
          !   (in tmat_newsolver). This initialization cannot be made before tmat_newsolver,
          !   because the division of the enegry loop (done in there) influences t_dtmatJij.
          call init_t_dtmatJij_at(t_inc,t_mpi_c_grid,t_dtmatJij(I1))

          !communicate
          if(t_dtmatJij(I1)%calculate)then
            iwork = product(shape(t_dtmatJij(I1)%dtmat_xyz))
            allocate(work_jij(iwork,1,1,1),STAT=ierr)
            if(ierr/=0) stop 
     & 'Error allocating work in main1a for communication of dtmat_xyz'
            call MPI_ALLREDUCE(t_dtmatJij(I1)%dtmat_xyz,work_jij,
     &                         iwork,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                         t_mpi_c_grid%mympi_comm_ie,IERR)
            if(ierr/=MPI_SUCCESS) stop 'error communicating t_dtmatJij'
            CALL ZCOPY(iwork,work_jij,1,t_dtmatJij(I1)%dtmat_xyz,1)
            deallocate(work_jij, stat=ierr)
         if(ierr/=0) stop '[main1a] Error deallocating work_jij'
          end if!t_dtmatJij(I1)%calculate

        END DO!I1=1,t_inc%NATYP
c      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
c      TTTTTTTTTTTTT TESTOUTPUT   Dij-implementation TTTTTTTTTTTTTTTTTT
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
c      TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

      end if!OPT('XCPL    ').and.OPT('NEWSOSOL')
#endif


! ------------------------------------------------------------------------
! determine the spin operator, torque operator and spin flux operator
! ------------------------------------------------------------------------
      if (OPT('OPERATOR')) then
#ifdef CPP_TIMING
        call timing_start('main1a - operator')
#endif

        ALLOCATE(PNS_SO_ALL(LMMAXD,LMMAXD,IRMD,2,NATYP))
        PNS_SO_ALL = CZERO

#ifdef CPP_MPI
        i1_start = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)+1
        i1_end   = t_mpi_c_grid%ioff_pT1(t_mpi_c_grid%myrank_ie)+
     &             t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
#else
        i1_start = 1
        i1_end   = NATYP
#endif

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

          !if(myrank==master.and.t_inc%i_write>0) then
          !  WRITE(1337,*) 'atom: ',I1
          !end if
#ifdef CPP_MPI
          ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
          ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
          ie_start = 0 ! offset
          ie_end   = IELAST
#endif

          DO ie_num=1,ie_end

            IE = ie_start+ie_num

            ! make sure only calculated at the Fermi level
            !if(ie_end==1 .or. ie==1) then
            if(ie==1) then
    
               !write(*,*) 'ie,iat', ie,i1
             
               if(t_wavefunctions%Nwfsavemax>0) then ! read wavefunctions?
                 ! read in wavefunction from memory
                 call read_wavefunc(t_wavefunctions,rll,rllleft,sll,
     &                          sllleft,i1,ie,NSRA,LMMAXD,t_inc%IRMDNEW,
     &                          0,1,rll_was_read_in,sll_was_read_in, 
     &                          rllleft_was_read_in,sllleft_was_read_in)
               end if ! t_wavefunctions%Nwfsavemax
               !write(*,*) 'read rll? ', ie, i1, rll_was_read_in
               !write(*,*) 'read sll? ', ie, i1, sll_was_read_in
               !write(*,*) 'read rlll?', ie, i1, rllleft_was_read_in
               !write(*,*) 'read slll?', ie, i1, sllleft_was_read_in

            end if ! ie==1

          END DO !ie_num=1,ie_end

       !if (I1==1) then
       !  write(*,*) 'write RLL before transform'
       !  write(7770,*) RLL
       !endif

          ! transform radial wavefunction back to old mesh
          ALLOCATE(RLLTEMP(t_inc%IRMDNEW,LMMAXD))
          ALLOCATE(PNSTEMP(IRWS(I1),LMMAXD))
          DO LM1=1,LMMAXD
           RLLTEMP=CZERO
           PNSTEMP=CZERO
           IR = 0
           !write(*,*) t_inc%IRMDNEW, IR
           DO IR=1,t_inc%IRMDNEW
            DO LM2=1,LMMAXD
             RLLTEMP(IR,LM2)=RLL(LM1,LM2,IR,0)
            ENDDO
           ENDDO
           CALL CHEB2OLDGRID(IRWS(I1),t_inc%IRMDNEW,LMMAXD,RMESH(:,I1),
     &                       NCHEB, NPAN_TOT(I1),RPAN_INTERVALL(:,I1),
     &                       IPAN_INTERVALL(:,I1),RLLTEMP,PNSTEMP,IRMD)
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
            CALL CHEB2OLDGRID(IRWS(I1),t_inc%IRMDNEW,LMMAXD,RMESH(:,I1),
     &                        NCHEB, NPAN_TOT(I1),RPAN_INTERVALL(:,I1),
     &                        IPAN_INTERVALL(:,I1),RLLTEMP,PNSTEMP,IRMD)
            DO IR=1,IRWS(I1)
             DO LM2=1,LMMAXD
              PNS_SO(LM1,LM2,IR,2)=PNSTEMP(IR,LM2)
             ENDDO
            ENDDO
           ENDDO ! LM1
          ENDIF ! NSRA.EQ.2

       !if (I1==1) then
       !  write(*,*) 'write PNSSO before transform'
       !  write(7771,*) PNS_SO
       !endif
          
          ! rotate radial wavefunction to global frame
          DO IR=1,IRMD
            CALL ROTATEMATRIX(PNS_SO(1,1,IR,1), THETA(I1), PHI(I1),
     &                        LMMAXD/2,0)
            CALL ROTATEMATRIX(PNS_SO(1,1,IR,2), THETA(I1), PHI(I1),
     &                        LMMAXD/2,0)
          ENDDO

       !if (I1==1) then
       !  write(*,*) 'write PNSSO after transform'
       !  write(7772,*) PNS_SO
       !endif

          ! finally collect wavefuncitons in global frame and old mesh for all atoms to be used in normcoeff-routines below
          PNS_SO_ALL(:,:,:,:,i1) = PNS_SO(:,:,:,:)

          DEALLOCATE(RLL,SLL,RLLLEFT,SLLLEFT,PNS_SO)
          DEALLOCATE(RLLTEMP)
          DEALLOCATE(PNSTEMP)

        END DO ! I1=i1_start, i1_end

#ifdef CPP_MPI
        ! finally gather PNS_SO_ALL on master in case of MPI run
        allocate(work(LMMAXD,LMMAXD,IRMD,2,NATYP), stat=ierr)
        if(ierr.ne.0) stop 
     &    'Error allocating work for MPI comm of PNS_SO_ALL in main1a'
        ihelp = LMMAXD*LMMAXD*IRMD*2*NATYP
        call MPI_ALLREDUCE(PNS_SO_ALL, work, ihelp, MPI_DOUBLE_COMPLEX,
     &                  MPI_SUM, t_mpi_c_grid%myMPI_comm_ie,ierr)

        if(ierr.ne.MPI_SUCCESS) stop 
     &    'Error in MPI comm of PNS_SO_ALL in main1a'
        PNS_SO_ALL(:,:,:,:,:) = work(:,:,:,:,:)
        deallocate(work, stat=ierr)
        if(ierr.ne.0) stop 
     &    'Error deallocating work for MPI comm of PNS_SO_ALL in main1a'
#endif


        ! done wf(myrank==master) then
 
        ! done with preparations, call normcoeff routines that construct operators
        !open(888, file='inp_normcoeff.txt', form='formatted')
        !write(8881 ,*) 'IRCUT',shape(IRCUT),IRCUT
        !write(8882 ,*) 'LMMAXD',LMMAXD/2
        !write(8883 ,*) 'PNS_SO',shape(PNS_SO_ALL),PNS_SO_ALL
        !write(8884 ,*) 'THETAS',shape(t_params%THETAS), t_params%THETAS
        !write(8885 ,*) 'NTCELL',t_params%NTCELL
        !write(8886 ,*) 'IFUNM',t_params%IFUNM
        !write(8887 ,*) 'IPAN',IPAN
        !write(8888 ,*) 'LMSP',t_params%LMSP
        !write(8889 ,*) 'KVREL',KORBIT !t_params%KREL
        !write(88810,*) 'CLEB',CLEB
        !write(88811,*) 'ICLEB',t_params%ICLEB
        !write(88812,*) 'IEND',IEND
        !write(88813,*) 'DRDI',DRDI
        !write(88814,*) 'IRWS',IRWS
        !write(88815,*) 'NSPOH',1+KORBIT
        !close(888)!ith preparations, call normcoeff routines that construct operators

           
        if(myrank==master) WRITE(*,*) 'Computing spin operator'
        call MPI_BARRIER(t_mpi_c_grid%myMPI_comm_ie,ierr)
        CALL NORMCOEFF_SO(IRCUT,LMMAXD/2,PNS_SO_ALL,
     +         t_params%THETAS,t_params%NTCELL,t_params%IFUNM,IPAN,
     +         t_params%LMSP,t_inc%KVREL,CLEB,t_params%ICLEB,IEND,
     +         DRDI,IRWS,1+KORBIT)
           
        if(myrank==master) then

           WRITE(*,*) 'Computing torq operator'
           CALL NORMCOEFF_SO_TORQ(IRCUT,LMMAXD/2,PNS_SO_ALL,
     +            t_params%NTCELL,t_params%IFUNM,IPAN,
     +            t_params%LMSP,t_inc%KVREL,CLEB,t_params%ICLEB,IEND,
     +            DRDI,IRWS,VISP,NSPIN,VINS,IRMIN)
           
           WRITE(*,*) 'Computing spinflux operator'
           CALL NORMCOEFF_SO_SPINFLUX(IRCUT,LMMAXD/2,PNS_SO_ALL,
     +            t_inc%KVREL,DRDI)

        end if ! myrank==master
#ifdef CPP_MPI
        ! make sure to wait fo all ranks to finish
        call MPI_BARRIER(t_mpi_c_grid%myMPI_comm_ie,ierr)
#endif

#ifdef CPP_TIMING
        call timing_stop('main1a - operator')
#endif
      end if ! OPERATOR


#ifdef CPP_TIMING
        call timing_start('main1a - tbref')
#endif
        IF ( LREFSYS ) CALL TBREF(EZ,IELAST,ALAT,VREF,IEND,LMAX,NCLS,
     &                            NINEQ,NREF,CLEB,RCLS,ATOM,CLS,ICLEB,
     &                            LOFLM,NACLS,REFPOT,RMTREF,TOLRDIF,
     &                            TMPDIR,ITMPDIR,ILTMP,
     &                            NAEZ,LLY) ! LLY Lloyd
#ifdef CPP_TIMING
        call timing_stop('main1a - tbref')
#endif


      if(t_inc%i_write>0)
     &  WRITE (1337,'(79(1H=),/,30X,"< KKR1a finished >",/,79(1H=),/)')

      deallocate(vins, vinsnew, stat=ierr)
      if(ierr/=0) stop '[main1a] Error deallocating vins, vinsnew'


      END SUBROUTINE !main1a
      
      END MODULE
