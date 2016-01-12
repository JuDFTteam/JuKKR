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
      
      use mod_types, only: t_tgmat, t_inc, t_lloyd
#ifdef CPP_MPI
     &                     ,gather_tmat, gather_gref,
     &                      t_mpi_c_grid, save_t_mpi_c_grid,
     &                      get_ntot_pT_ioff_pT_2D 
#endif
      use mod_mympi, only: myrank, nranks, master
#ifdef CPP_MPI
     &                     ,find_dims_2d,
     &                     distribute_linear_on_tasks
#endif
#ifdef CPP_TIMING
      use mod_timing
#endif
      use mod_wunfiles
     
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
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,PI
      PARAMETER  (PI=4.d0*datan(1.d0))
      INTEGER I1,ICST,IELAST,IEND,INS,IPOT,ISPIN,ITSCF,LMAX,NATYP,NAEZ, ! LLY added NAEZ
     +        NCLS,NINEQ,NREF,NSPIN,NSRA,LM1,IR
      INTEGER LLY ! LLY <> 0: apply Lloyds formula
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION R_LOG,THETA(NATYPD),PHI(NATYPD),TOLRDIF
      LOGICAL          LCHECKANGLES
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
      CHARACTER*80 TMPDIR
      LOGICAL LREFSYS,OPT,TEST,LREAD
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU,OLD
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD)
      DOUBLE PRECISION, allocatable :: ULDAU(:,:,:,:,:) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
#ifdef CPP_MPI
      integer :: myMPI_comm_grid, myMPI_comm_at, myMPI_comm_ie,
     &           myrank_at, myrank_ie, nranks_at,nranks_ie
      integer :: ntot1, dims(2), mytot, ii
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),
     &           ntot_all(0:nranks-1), ioff_all(0:nranks-1)
#endif
      integer :: i1_start, i1_end

C-----------------------------------------------------------------------
C     ..
C     .. Arrays in Common ..
      CHARACTER*8 OPTC(32),TESTC(32)
C     ..
C     .. Common blocks ..
      COMMON /OPTC/OPTC
      COMMON /TESTC/TESTC
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
     &       SOLVER,TMPDIR,DELTAE,TESTC,OPTC,tolrdif,NPAN_LOG,NPAN_EQ,
     &       NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,LMAXD,
     &       NTOTD,NRMAXD,R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,
     &       R2DRDIREL,RMREL,IRMIND,LMPOTD,NSPOTD,NPOTD,JWSREL,ZREL,
     &       ITSCF)
C
      IF ( TEST('Vspher  ') ) VINS(IRMIND:IRMD,2:LMPOTD,1:NSPOTD) = 0.D0
      
      if(IDOLDAU==1) ALLOCATE( ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) )

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
      IF (.NOT.OPT('NEWSOSOL')) THEN
      
#ifdef CPP_MPI
!        write(*,*) 'distr:',ntot1,t_mpi_c_grid%nranks_ie

        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
     &                   t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                             master,ntot1,ntot_pT,ioff_pT,.true.,
     + .true.)

!            write(*,*) 'red myrank_at',myrank,ntot_pT,ioff_PT
!         if(t_mpi_c_grid%nranks_ie>t_mpi_c_grid%dims(2)) then
!            i1_start = ioff_pT(t_mpi_c_grid%myrank_ie
!      &               -(t_mpi_c_grid%nranks_ie-t_mpi_c_grid%dims(2)))+1
!            i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie
!      &               -(t_mpi_c_grid%nranks_ie-t_mpi_c_grid%dims(2)))+
!      &                ntot_pT(t_mpi_c_grid%myrank_ie)
!            t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie
!      &               -(t_mpi_c_grid%nranks_ie-t_mpi_c_grid%dims(2)))
!         else
           i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
           i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+
     &                ntot_pT(t_mpi_c_grid%myrank_ie)
           t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)
!         endif

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
       
!        write(*,*) 'done calctmat'
     
      ELSE
       
      ! nonco angles: defined in mod_wunfiles
      call read_angles(t_params,NATYP,THETA,PHI)

c interpolate potential

       CALL INTERPOLATE_POTEN(NSPIN,RMESH,IRMIN,IRWS,IPAN,IRCUT,VINS,
     +                        VISP,NPAN_LOG,NPAN_EQ,NCHEB,NPAN_TOT,
     +                        RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                        VINSNEW)

#ifdef CPP_MPI

        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
     &                   t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                             master,ntot1,ntot_pT,ioff_pT,.true.,
     &  .true.)

        i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
        i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+
     &             ntot_pT(t_mpi_c_grid%myrank_ie)
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

!        write(*,*) 'i1 loop', myrank, i1_start, i1_end

       DO I1=i1_start,i1_end

        IPOT=NSPIN*(I1-1)+1
        CALL TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,RMESH(1,I1),ZAT(I1),
     +                      SOCSCALE(I1),EZ,NSRA,CLEB(1,1),ICLEB,IEND,
     &                      NCHEB,NPAN_TOT(I1),
     +                      RPAN_INTERVALL(0,I1),IPAN_INTERVALL(0,I1),
     +                      RNEW(1,I1),VINSNEW,THETA(I1),PHI(I1),I1,IPOT
     &                     ,LLY,DELTAE) ! LLY

       ENDDO !I1, atom loop
       CLOSE(10)

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
      if(.not.t_tgmat%tmat_to_file) then ! .and. t_mpi_c_grid%dims(1)>1) then
!       if(.not.t_tgmat%tmat_to_file) then
      
         do ii=0,t_mpi_c_grid%nranks_ie-1
          ntot_all(ii) = t_mpi_c_grid%ntot_pT1(ii)
          ioff_all(ii) = t_mpi_c_grid%ioff_pT1(ii)
         end do
         
          mytot = t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
!      
         call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
     &        mytot,t_mpi_c_grid%mympi_comm_ie, t_mpi_c_grid%nranks_ie)
      end if
#endif


! #ifdef CPP_MPI
!       if(t_mpi_c_grid%myrank_ie==0) then
! #endif
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

! not needed any more, already done in tbref
! #ifdef CPP_MPI
!       end if !(t_mpi_c_grid%myrank_at==0)
!       
!       if(.not.t_tgmat%gref_to_file .and. t_mpi_c_grid%dims(1)>1) then
!       
!          do ii=0,t_mpi_c_grid%nranks_ie-1
!           ntot_all(ii) = t_mpi_c_grid%ntot_pT1(ii)
!           ioff_all(ii) = t_mpi_c_grid%ioff_pT1(ii)
!          end do
!           mytot = t_mpi_c_grid%ntot2 !t_mpi_c_grid%ntot_pT1(t_mpi_c_grid%myrank_ie)
!          call gather_gref(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
!      &        mytot,t_mpi_c_grid%mympi_comm_ie, t_mpi_c_grid%nranks_ie)
!         end if
! #endif


      if(t_inc%i_write>0)
     &  WRITE (1337,'(79(1H=),/,30X,"< KKR1a finished >",/,79(1H=),/)')


      if(IDOLDAU==1) DEALLOCATE( ULDAU )
      deallocate(vins, vinsnew)


      END SUBROUTINE !main1a
      
      END MODULE
