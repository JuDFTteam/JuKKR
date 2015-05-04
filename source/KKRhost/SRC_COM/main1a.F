      module MOD_MAIN1A
      
      
      !use modulename, only rountine1
      
      implicit none
      
      contains
      
      
      subroutine main1a()
      
#ifdef CPP_MPI
      use mpi
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
     &                     distribute_linear_on_tasks,
     &                     create_subarr_comm
#endif
     
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
      DOUBLE PRECISION ALAT
      INTEGER I1,ICST,IELAST,IEND,INS,IPOT,ISPIN,ITSCF,LMAX,NATYP,NAEZ, ! LLY added NAEZ
     +        NCLS,NINEQ,NREF,NSPIN,NSRA,LM1,IR
      INTEGER LLY ! LLY <> 0: apply Lloyds formula
      DOUBLE COMPLEX DELTAE  ! Energy difference for numerical derivative
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION R_LOG,THETA(NATYPD),PHI(NATYPD),TOLRDIF
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     +                 RNEW(NRMAXD,NATYPD),
     +                 VINSNEW(NRMAXD,LMPOTD,NSPOTD)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
C     ..
C     .. Local Arrays ..
CMPI  INTEGER MYRANK,NROFNODES,IERR
CMPI  COMMON /MPI/MYRANK,NROFNODES
CMPI  DOUBLE COMPLEX WORK(IRMD,LMPOTD,NATYPD,NSPIND)
CMPI  DOUBLE COMPLEX WORK1(NACLSD+IEMXD,LMMAXD,LMMAXD,NCLSD)
CMPI  EXTERNAL CINIT,MPI_ALLREDUCE,MPI_COMM_RANK,MPI_COMM_SIZE,
CMPI +         MPI_FINALIZE,MPI_INIT
CT3E  INTEGER ICLKTCK,ICLOCKS,IE,MEND,MSTART
CT3E  DOUBLE PRECISION TIME1,TIME2,TIME3,TIME4
      DOUBLE COMPLEX EZ(IEMXD)
      DOUBLE PRECISION DRDI(IRMD,NATYPD),RMESH(IRMD,NATYPD),ZAT(NATYPD),
     +                 RCLS(3,NACLSD,NCLSD),RMTREF(NREFD),VREF(NREFD),
     +                 VINS(IRMIND:IRMD,LMPOTD,NSPOTD),VISP(IRMD,NPOTD),
     +                 CLEB(NCLEB,2)
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
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
#ifdef CPP_MPI
      integer :: myMPI_comm_grid, myMPI_comm_row, myMPI_comm_col,
     &           myrank_grid, myrank_row, myrank_col, nranks_row,
     &           nranks_col
      integer :: ntot1, dims(2), i1_start, i1_end, mytot
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),
     &           ntot_all(0:nranks-1), ioff_all(0:nranks-1)
#endif

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
      DATA DELTAE /(1.D-3,0.D0)/
C
Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C
CT3E  CALL SYSTEM_CLOCK(MSTART)
CMPI  CALL MPI_INIT(IERR)
CMPI  CALL MPI_COMM_RANK(MPI_COMM_WORLD,MYRANK,IERR)
CMPI  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,NROFNODES,IERR)
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C -------------------------------------------------------------- input1a
C
      OPEN (67,FILE='input1a.unformatted',FORM='unformatted')
      READ (67) NSRA,INS,NAEZ,NATYP,NSPIN,ICST,IPAN,IRCUT,    ! LLY added NAEZ
     &          LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY
C ......................................................................
Consistency check 
C
      IF ( (KREL.EQ.1) .AND. (INS.NE.0) ) THEN
         WRITE(6,*)
     &        ' FULL-POTENTIAL RELATIVISTIC mode not implemented '
         STOP ' set INS = 0 in the input'
      END IF
C     
      IF ( NSRA.LE.2 ) THEN
         IF ( KREL.EQ.1 ) STOP
     &        ' KVREL <= 1 in input, but relativistic program used'
      ELSE
         IF ( KREL.EQ.0 ) STOP
     &        ' KVREL > 1 in input, but non-relativistic program used'
      END IF
C ......................................................................
      READ (67) ALAT,ZAT,DRDI,RMESH,RMTREF,VREF,IEND,CLEB,RCLS,
     &          ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,TESTC,OPTC,
     &          IRWS,IRMIN,TOLRDIF,DELTAE,SOCSCALE
      READ (67) TMPDIR,ITMPDIR,ILTMP
      READ (67) NPAN_LOG,NPAN_EQ,NCHEB,R_LOG,NPAN_TOT
      READ (67) RNEW,RPAN_INTERVALL,IPAN_INTERVALL
      IF (KREL.EQ.1) READ(67) SOLVER,SOCSCL,CSCL
      IF ( IDOLDAU.EQ.1 ) READ(67) NTLDAU,ITLDAU,LOPT,UEFF,JEFF,EREFLDAU
      CLOSE (67)
C ---------------------------------------------------------- energy_mesh
C
      IF (t_inc%i_iteration.eq.0) then
        OPEN (67,FILE='energy_mesh',FORM='unformatted')
      else
        OPEN (67,FILE='new_energy_mesh',FORM='unformatted')
      end if
      READ (67) IELAST,EZ
      CLOSE (67)
C ------------------------------------------------------ input_potential
C
!       OPEN (67,FILE='input_potential',FORM='unformatted')
      IF (t_inc%i_iteration.eq.0) then
        OPEN (67,FILE='input_scf.unformatted',FORM='unformatted')
      else
        OPEN (67,FILE='output_scf.unformatted',FORM='unformatted')
      end if
      READ (67) VINS,VISP
      IF (KREL.EQ.1) THEN
         READ (67) RMREL,DRDIREL,R2DRDIREL
         READ (67) ZREL,JWSREL
         READ (67) VTREL,BTREL
      END IF
      READ (67) ITSCF
      CLOSE (67)
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
! MPI:     
      ntot1 = t_inc%NATYP
!       if(.not.OPT('NEWSOSOL')) ntot1 = t_inc%NSPIN*t_inc%NATYP
!
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
        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_row, 
     &  t_mpi_c_grid%myrank_row,master,ntot1,ntot_pT,ioff_pT,.false.)
!         call distribute_linear_on_tasks(t_mpi_c_grid%nranks_row, 
!      &  t_mpi_c_grid%myrank_row,master,ntot1,ntot_pT,ioff_pT,.true.)
        i1_start = ioff_pT(t_mpi_c_grid%myrank_row)+1
        i1_end   = ioff_pT(t_mpi_c_grid%myrank_row)+
     &             ntot_pT(t_mpi_c_grid%myrank_row)
        t_mpi_c_grid%ntot1  = ntot1 
        if(.not. (allocated(t_mpi_c_grid%ntot_pT1) .and. 
     &            allocated(t_mpi_c_grid%ioff_pT1))) 
     &     allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_row-1),
     &              t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_row-1))
        t_mpi_c_grid%ntot_pT1 = ntot_pT
        t_mpi_c_grid%ioff_pT1 = ioff_pT
#else
      i1_start = 1
      i1_end   = NATYP
#endif
!        DO IPOT = i1_start,i1_end
        DO I1 = i1_start,i1_end
       DO ISPIN = 1,NSPIN
         IPOT=NSPIN*(I1-1)+ispin
!          ISPIN = MOD(IPOT+1,NSPIN) + 1
!          I1 = (IPOT-ISPIN)/NSPIN + 1
C
!          write(456456,*) 'line1',ICST,INS,IELAST,
!      +         'line2',NSRA,ISPIN,NSPIN,I1,EZ,
!      +                 'line3',DRDI(1,I1),RMESH(1,I1),
!      +                 'line4',VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),
!      +                 'line5',VISP(1,IPOT),ZAT(I1),IRMIN(I1),IPAN(I1),   ! Added IRMIN 1.7.2014
!      +              'line6',IRCUT(0,I1),CLEB,LOFLM,ICLEB,IEND,SOLVER,
!      +                 'line7',SOCSCL(1,KREL*I1+(1-KREL)),
!      +                 'line8',CSCL(1,KREL*I1+(1-KREL)),
!      +                 'line9',VTREL(1,I1),BTREL(1,I1),
!      +              'line10',RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
!      +                 'line11',ZREL(I1),JWSREL(I1),
!      +                 'line12',IDOLDAU,LOPT(I1),WLDAU(1,1,1,I1),
!      &                 'line13',LLY,DELTAE
!              write(456456,*) 
!      +            'line6.1',IRCUT(0,I1),
!      +            'line6.2',CLEB,
!      +            'line6.3',LOFLM,
!      +            'line6.4',ICLEB,
!      +            'line6.5',IEND,
!      +            'line6.6',SOLVER

         CALL CALCTMAT(ICST,INS,IELAST,
     +                 NSRA,ISPIN,NSPIN,I1,EZ,
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
       
       THETA(:) = 0.D0
       PHI(:) = 0.D0
       LREAD = .FALSE.
       INQUIRE(file='nonco_angle.dat',EXIST=LREAD)
       IF (LREAD) THEN
          OPEN(UNIT=10,FILE='nonco_angle.dat',FORM='FORMATTED')
          DO I1 = 1,NATYP
             READ(10,*) THETA(I1),PHI(I1)
          ENDDO
       ENDIF

c interpolate potential

       CALL INTERPOLATE_POTEN(NSPIN,RMESH,IRMIN,IRWS,IPAN,IRCUT,VINS,
     +                        VISP,NPAN_LOG,NPAN_EQ,NCHEB,NPAN_TOT,
     +                        RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                        VINSNEW)

#ifdef CPP_MPI

        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_row, 
     &  t_mpi_c_grid%myrank_row,master,ntot1,ntot_pT,ioff_pT,.false.)
!         call distribute_linear_on_tasks(t_mpi_c_grid%nranks_row, 
!      &  t_mpi_c_grid%myrank_row,master,ntot1,ntot_pT,ioff_pT,.true.)
        i1_start = ioff_pT(myrank_row)+1
        i1_end   = ioff_pT(myrank_row)+ntot_pT(myrank_row)
        t_mpi_c_grid%ntot1  = ntot1 
        if(.not. (allocated(t_mpi_c_grid%ntot_pT1) .and. 
     &            allocated(t_mpi_c_grid%ioff_pT1))) 
     &     allocate(t_mpi_c_grid%ntot_pT1(0:t_mpi_c_grid%nranks_row-1),
     &              t_mpi_c_grid%ioff_pT1(0:t_mpi_c_grid%nranks_row-1))
        t_mpi_c_grid%ntot_pT1 = ntot_pT
        t_mpi_c_grid%ioff_pT1 = ioff_pT
#else
      i1_start = 1
      i1_end   = NATYP
#endif
       DO I1=i1_start,i1_end

c read from file theta and phi

        IPOT=NSPIN*(I1-1)+1
        CALL TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,RMESH(1,I1),ZAT(I1),
     +                      SOCSCALE(I1),EZ,NSRA,CLEB(1,1),ICLEB,IEND,
     &                      NCHEB,NPAN_TOT(I1),
     +                      RPAN_INTERVALL(0,I1),IPAN_INTERVALL(0,I1),
     +                      RNEW(1,I1),VINSNEW,THETA(I1),PHI(I1),I1,IPOT
     &                     ,LLY,DELTAE) ! LLY

       ENDDO !I1, atom loop
       CLOSE(10)

      ENDIF
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
      if(.not.t_tgmat%tmat_to_file .and. t_mpi_c_grid%dims(1)>0) then
         mytot = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_col)     
         call get_ntot_pT_ioff_pT_2D(t_mpi_c_grid,ntot_all,ioff_all)
!          call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
!      &                    mytot,MPI_COMM_WORLD)
!      
         call gather_tmat(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
     &                    mytot,t_mpi_c_grid%myMPI_comm_row)
      end if
#endif


CT3E  CALL SYSTEM_CLOCK(MEND)
CT3E  TIME1 = (MEND-MSTART)
CT3E  CALL SYSTEM_CLOCK(MSTART)
#ifdef CPP_MPI
      if(t_mpi_c_grid%myrank_row==0) then
#endif
        IF ( LREFSYS ) CALL TBREF(EZ,IELAST,ALAT,VREF,IEND,LMAX,NCLS,
     &                            NINEQ,NREF,CLEB,RCLS,ATOM,CLS,ICLEB,
     &                            LOFLM,NACLS,REFPOT,RMTREF,TOLRDIF,
     &                            TMPDIR,ITMPDIR,ILTMP,
     &                            NAEZ,LLY) ! LLY Lloyd

#ifdef CPP_MPI
        if(.not.t_tgmat%gref_to_file .and. t_mpi_c_grid%dims(1)>0) then
          call get_ntot_pT_ioff_pT_2D(t_mpi_c_grid,ntot_all,ioff_all)
     
         call gather_gref(t_inc,t_tgmat,t_mpi_c_grid,
     &                    t_mpi_c_grid%ntot_pT2,t_mpi_c_grid%ioff_pT2,
     &                    mytot,t_mpi_c_grid%myMPI_comm_col)
!          call gather_gref(t_inc,t_tgmat,t_mpi_c_grid,ntot_all,ioff_all,
!      &                    mytot,t_mpi_c_grid%myMPI_comm_row)
        end if
      
      end if !(t_mpi_c_grid%myrank_row==0)
#endif


      if(myrank.eq.master) 
     &     WRITE (6,'(79(1H=),/,30X,"< KKR1a finished >",/,79(1H=),/)')


      END SUBROUTINE !main1a
      
      END MODULE
