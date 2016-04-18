      MODULE MOD_MAIN1C
      
      implicit none
      
      contains
      
      
      subroutine main1c()

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
     &                     distribute_linear_on_tasks,
     &                     mympi_main1c_comm,mympi_main1c_comm_newsosol2
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
      INTEGER LMMAXD,LMPOTD
      PARAMETER (LMPOTD= (LPOTD+1)**2)
      PARAMETER (LMMAXD= (KREL+KORBIT+1) * (LMAXD+1)**2)
      INTEGER LMAXD1
      PARAMETER (LMAXD1=LMAXD+1)
      INTEGER LMMAXD1
      PARAMETER (LMMAXD1=LMMAXD+1)
      INTEGER MMAXD
      PARAMETER (MMAXD=2*LMAXD+1)
      INTEGER LM2D
      PARAMETER (LM2D= (2*LMAXD+1)**2)
      INTEGER LMXSPD
      PARAMETER (LMXSPD= (2*LPOTD+1)**2)
      INTEGER NPOTD
      PARAMETER (NPOTD= (2*(KREL+KORBIT) + 
     +           (1-(KREL+KORBIT))*NSPIND)*NATYPD)
      INTEGER IRMIND
      PARAMETER (IRMIND=IRMD-IRNSD)
      INTEGER NRMAXD
      PARAMETER (NRMAXD=NTOTD*(NCHEBD+1))
C     .. itermdir parameter
      INTEGER NMVECMAX
      PARAMETER (NMVECMAX = 4)
C     .. 
      INTEGER LRECTMT,LRECTMT2
      PARAMETER (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
!     ..
      DOUBLE COMPLEX CZERO
      PARAMETER (CZERO=(0.D0,0.D0))
!     ..
      INTEGER LLY ! LLY <> 0 : apply Lloyd's formula

C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION ALAT,DENEF,E1,E2,PI,TK,EFERMI,CHRGSEMICORE
      INTEGER I1,ICELL,ICST,IE,IELAST,IEND,INS,INTERVX,INTERVY,INTERVZ,
     +        IPOT,IPOT1,IR,IS,ISPIN,L,LM,NACLS1,
     +        NATYP,NAEZ,NPOL,NSPIN,NSRA,LMAX,LMAXP1
      INTEGER KMROT,IQ,NSPINPOT,IHOST,ITMPDIR,ILTMP,IESEMICORE
      CHARACTER*5 TEXTNS
      CHARACTER*10 SOLVER
      CHARACTER*80 TMPDIR
      LOGICAL ITERMVDIR,LDORHOEF,LMOMVEC,LREAD
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C LDA+U
      INTEGER IDOLDAU,ITRUNLDAU,NTLDAU
      INTEGER LOPT(NATYPD),ITLDAU(NATYPD)
      DOUBLE PRECISION EREFLDAU(NATYPD),UEFF(NATYPD),JEFF(NATYPD),
     &                 EDC(NATYPD),EU(NATYPD)
      DOUBLE PRECISION ULDAU(MMAXD,MMAXD,MMAXD,MMAXD,NATYPD) 
      DOUBLE PRECISION WLDAU(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE PRECISION WLDAUOLD(MMAXD,MMAXD,NSPIND,NATYPD)
      DOUBLE COMPLEX PHILDAU(IRMD,NATYPD) 
      DOUBLE COMPLEX DENMATC(MMAXD,MMAXD,NPOTD)
      DOUBLE COMPLEX DENMATN(MMAXD,MMAXD,2,2,NATYPD)
C for new spin-orbit solver
      INTEGER NPAN_LOG(NATYPD),NPAN_EQ(NATYPD),NCHEB,NPAN_TOT(NATYPD)
      DOUBLE PRECISION R_LOG,THETA(NATYPD),PHI(NATYPD)
      DOUBLE PRECISION RPAN_INTERVALL(0:NTOTD,NATYPD),
     +                 RNEW(NRMAXD,NATYPD)
      double precision, allocatable :: THETASNEW(:,:,:), vinsnew(:,:,:)
      INTEGER          IPAN_INTERVALL(0:NTOTD,NATYPD)
C LDA+U
C LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
C     ..
C     .. Local Arrays ..
#ifdef CPP_MPI
      DOUBLE COMPLEX,   allocatable :: WORK(:,:)
      integer ierr,idim, nranks_local, myrank_ie_tmp
      double complex, allocatable :: workc(:,:,:,:)
#endif
      DOUBLE COMPLEX EZ(IEMXD),WEZ(IEMXD)
      DOUBLE COMPLEX, allocatable :: den(:,:,:,:), denlm(:,:,:,:) ! DEN(0:LMAXD1,IEMXD,NPOTD,NQDOS),DENLM(LMMAXD1,IEMXD,NPOTD,NQDOS)
      DOUBLE COMPLEX DEN1(0:LMAXD1,IEMXD,2),WEZNORM(IEMXD)
      DOUBLE COMPLEX CDOSAT0(IEMXD),CDOSAT1(IEMXD),CSUM                       ! LLY Lloyd
      DOUBLE COMPLEX CDOS0(IEMXD),CDOS_LLY(IEMXD,NSPIND),CDOS1(IEMXD),
     +               CDOS2(NATYPD),cdos_dum(IEMXD,NSPIND)                     ! LLY Lloyd
      DOUBLE PRECISION CHARGE_LLY(NSPIND)  ! LLY
      DOUBLE COMPLEX EREAD                  !LLY
      DOUBLE PRECISION A(NATYPD),B(NATYPD),CHARGE(0:LMAXD1,NATYPD,2),
     +               CLEB(NCLEB,2),DRDI(IRMD,NATYPD),RMESH(IRMD,NATYPD),
     +                 DOSTOT(0:LMAXD1,2),ECORE(20,NPOTD),
     +                 ESPV(0:LMAXD1,NPOTD),DENEFAT(NATYPD)
      DOUBLE PRECISION ESPV1(0:LMAXD1,2)
C ----------------------------------------------------------------------
C  R2NEF (IRMD,LMPOTD,NATYPD,2)  ! rho at FERMI energy
C  RHO2NS(IRMD,LMPOTD,NATYPD,2)  ! radial density
C   nspin=1            : (*,*,*,1) radial charge density
C   nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
C                               (*,*,*,2) rho(2) - rho(1) -> mag. moment
C  RHOC(IRMD,NPOTD)              ! core charge density
C ----------------------------------------------------------------------
      DOUBLE PRECISION RHOC(IRMD,NPOTD),
     +                 RHO2M1(IRMD,LMPOTD,4),
     +                 RHO2M2(IRMD,LMPOTD,4)
      double precision, allocatable :: rho2n1(:,:,:), rho2n2(:,:,:),
     +                                 rho2ns(:,:,:,:), r2nef(:,:,:,:)
      DOUBLE PRECISION VISP(IRMD,NPOTD),
     +                 ZAT(NATYPD),CONC(NATYPD),
     &                 SOCSCALE(NATYPD)
      double precision, allocatable :: THETAS(:,:,:), VINS(:,:,:)
C----------------------------------- orbital magnetic moment
C     attention: muorb second index means both spins and total
C----------------------------------- orbital density
      DOUBLE PRECISION  MUORB(0:LMAXD1+1,3,NATYPD)
      DOUBLE PRECISION  RHOORB(IRMD*KREL + (1-KREL),NATYPD)
C---------------------------------------------------------------      
      DOUBLE PRECISION SOCSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION CSCL(KREL*LMAXD+1,KREL*NATYPD+(1-KREL))
      DOUBLE PRECISION VTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION BTREL(IRMD*KREL+(1-KREL),NATYPD)
      DOUBLE PRECISION DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 R2DRDIREL(IRMD*KREL+(1-KREL),NATYPD),
     &                 RMREL(IRMD*KREL+(1-KREL),NATYPD)
      INTEGER JWSREL(NATYPD),ZREL(NATYPD),IRSHIFT(NATYPD)
C ===================================================================
C   RELATIVISTIC TREATMENT OF CORE ELECTRONS   July/2002
C   SEE ROUTINE <DRVCORE> FOR A SHORT DESCRIPTION OF THE VARIABLES
C
      DOUBLE PRECISION ECOREREL(KREL*20+(1-KREL),NPOTD) 
      INTEGER NKCORE(20,NATYPD),KAPCORE(20,NPOTD)
C ===================================================================
      INTEGER ICLEB(NCLEB,4),IFUNM1(LMXSPD,NATYPD),IPAN(NATYPD),
     +        IRCUT(0:IPAND,NATYPD),ITITLE(20,NPOTD),IRMIN(NATYPD),
     +        LMSP1(LMXSPD,NATYPD),JEND(LMPOTD,0:LMAXD,0:LMAXD),
     +        LCORE(20,NPOTD),NCORE(NPOTD),LOFLM(LM2D),NTCELL(NATYPD),
     +        IRWS(NATYPD),NFU(NATYPD),LLMSP(NATYPD,NFUND)
      CHARACTER*4 TEXTL(0:6)
      CHARACTER*7 TEXTS(3)
C     .. itermdir variables 
      DOUBLE PRECISION QMTET(NAEZD),QMPHI(NAEZD)
      INTEGER IQAT(NATYPD),NQDOS,LMMAXSO
      PARAMETER (LMMAXSO = 2*LMMAXD)
      DOUBLE COMPLEX MVEVI(NATYPD,3,NMVECMAX),         ! OUTPUT
     &               MVEVIL(0:LMAXD,NATYPD,3,NMVECMAX) 
      DOUBLE COMPLEX MVEVIEF(NATYPD,3,NMVECMAX)        ! OUTPUT
      DOUBLE COMPLEX MVEVIL1(0:LMAXD,3,NMVECMAX)
      DOUBLE COMPLEX MVEVIL2(0:LMAXD,3,NMVECMAX) ! WORK ARRAYS
C     .. qdos and lmlm-dos ..
      CHARACTER*8 QDOSOPT
      DOUBLE COMPLEX DF(IEMXD)!,GFLLE(LMMAXSO,LMMAXSO,IEMXD,100)
      DOUBLE COMPLEX DENTOT ! qdos
      DOUBLE PRECISION, allocatable :: QVEC(:,:)
      INTEGER IREC,LM1,LM2
#ifdef CPP_MPI
      integer :: myMPI_comm_grid, myMPI_comm_at, myMPI_comm_ie,
     &           myrank_at, myrank_ie, nranks_at,nranks_ie
      integer :: ntot1, dims(2), mytot, ii
      integer :: ntot_pT(0:nranks-1), ioff_pT(0:nranks-1),
     &           ntot_all(0:nranks-1), ioff_all(0:nranks-1)
      integer :: ihelp
#endif
      integer :: i1_start, i1_end
      integer :: ie_start, ie_end, ie_num
      
      double precision :: angles_new(2,NATYPD)
C     ..
C     .. External Subroutines ..
      EXTERNAL CINIT,RHOCORE,RHOVAL,WMATLDAU,WRLDAUPOT,WRLDOS,WRMOMS
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,DBLE,DIMAG,DREAL
C     ..
C     .. External Functions ..
      LOGICAL OPT,TEST
      EXTERNAL OPT,TEST
C     ..
C     .. Data statements ..
      DATA TEXTL/' s =',' p =',' d =',' f =',' g =',' h =',' i ='/
      DATA TEXTS/'spin dn','spin up','       '/
      DATA TEXTNS/' ns ='/
      DATA LDORHOEF/.TRUE./
      DATA IHOST / 1 /          ! this is the host program
      
      
      !allocations
      allocate(THETASNEW(NRMAXD,NFUND,NCELLD),THETAS(IRID,NFUND,NCELLD))
      allocate(VINSNEW(NRMAXD,LMPOTD,NSPOTD), 
     +         VINS(IRMIND:IRMD,LMPOTD,NSPOTD))
      allocate(RHO2N1(IRMD,LMPOTD,NPOTD),RHO2NS(IRMD,LMPOTD,NATYPD,2),
     +         RHO2N2(IRMD,LMPOTD,NPOTD),R2NEF(IRMD,LMPOTD,NATYPD,2))


      !initialze to zero
      THETASNEW = 0.0d0
      THETAS = 0.0d0
      VINSNEW = 0.0d0
      VINS = 0.0d0
      RHO2NS = 0.0d0
      R2NEF = 0.0d0
      
C     ..
C
Consistency check
      IF ( (KREL.LT.0) .OR. (KREL.GT.1) )
     &     STOP ' set KREL=0/1 (non/fully) relativistic mode in inc.p'
      IF ( (KREL.EQ.1) .AND. (NSPIND.EQ.2) ) 
     &   STOP ' set NSPIND = 1 for KREL = 1 in inc.p'
C
C ======================================================================
C =             read in variables from unformatted files               =
C ======================================================================
C
C -------------------------------------------------------------- input1c
      call get_params_1c(t_params,KREL,NAEZD,NATYPD,NCLEB,LM2D,
     &          NCHEBD,IPAND,LMPOTD,LMAXD,LMXSPD,NFUND,NPOTD,NCLSD,
     &          NTOTD,MMAXD,IEMXD,IRMD,NSRA,INS,NATYP,NAEZ,NSPIN,NACLS1,
     &          ICST,KMROT,IQAT,IDOLDAU,LMAX,IRWS,IPAN,IRCUT,IEND,ICLEB,
     &          LOFLM,JEND,IFUNM1,LMSP1,NFU,LLMSP,LCORE,NCORE,NTCELL,
     &          IRMIN,ITITLE,INTERVX,INTERVY,INTERVZ,LLY,ITMPDIR,
     &          ILTMP,NPAN_EQ,IPAN_INTERVALL,NPAN_LOG,NCHEB,NPAN_TOT,
     &          NTLDAU,LOPT,ITLDAU,IELAST,IESEMICORE,NPOL,IRSHIFT,JWSREL
     &         ,ZREL,ITRUNLDAU,QMTET,QMPHI,CONC,ALAT,ZAT,DRDI,RMESH,A,B,
     &          CLEB,THETAS,SOCSCALE,RPAN_INTERVALL,CSCL,RNEW,SOCSCL,
     &          THETASNEW,EFERMI,EREFLDAU,UEFF,JEFF,E1,E2,TK,VINS,VISP,
     &          ECORE,DRDIREL,R2DRDIREL,RMREL,VTREL,BTREL,WLDAU,ULDAU,EZ
     &         ,WEZ,PHILDAU,TMPDIR,SOLVER,NSPIND,NSPOTD,
     &          IRMIND,LMAXD1,NCELLD,IRID,R_LOG)
      ! initialization needed due to merging to one executable
      ESPV(:,:) = 0.d0
      RHO2N1(:,:,:) = 0.d0
      RHO2N2(:,:,:) = 0.d0
C ======================================================================
C =                     End read in variables                          =
C ======================================================================
C
      IF ( IDOLDAU.EQ.1 ) THEN 
C
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         READ (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
      ENDIF

      NQDOS=1
      IF (OPT('qdos    ')) THEN                         ! qdos
C        Read BZ path for qdos calculation:
         OPEN(67,FILE='qvec.dat')                       ! qdos
         READ(67,*) NQDOS                               ! qdos
         ALLOCATE(QVEC(3,NQDOS))                        ! qdos
         DO IQ = 1,NQDOS                                ! qdos
            READ(67,*) (QVEC(I1,IQ),I1=1,3)             ! qdos
         ENDDO                                          ! qdos
         CLOSE(67)                                      ! qdos
      END IF
      
      allocate( DEN(0:LMAXD1,IEMXD,NQDOS,NPOTD),
     &          DENLM(LMMAXD,IEMXD,NQDOS,NPOTD))

      CALL CINIT(IEMXD*(LMAXD+2)*NPOTD*NQDOS,DEN)
      CALL CINIT(IEMXD*(LMMAXD)*NPOTD*NQDOS,DENLM)
      DENEF = 0.0D0
      CALL RINIT(NATYPD,DENEFAT)
C
      ITERMVDIR = OPT('ITERMDIR') 
      LMOMVEC = ( ITERMVDIR .OR. ( KMROT.NE.0 ) )
      PI = 4.0D0*ATAN(1.0D0)
      NSPINPOT = KREL*2 + (1-KREL)*NSPIN
C
C --> no need to calculate charge correction if no host program, if 
C     decimation or if no energy contour
C
      LDORHOEF = (IHOST.EQ.1) .AND. (.NOT.OPT('DECIMATE')) 
     &           .AND. (NPOL.NE.0)
C
C LDA+U 
      IF ( IDOLDAU.EQ.1 ) CALL CINIT(MMAXD*MMAXD*NPOTD,DENMATC(1,1,1))
C LDA+U 
C
      if (t_tgmat%gmat_to_file) CALL OPENDAFILE(69,'gmat',
     +                               4,LRECTMT,TMPDIR,ITMPDIR,ILTMP)

!    write parameters file that contains passed parameters for further treatment of gflle
      IF (OPT('lmlm-dos')) THEN                                          ! lmlm-dos
        QDOSOPT = 'n'                                                    ! lmlm-dos
        IF (OPT('qdos    ')) THEN                                        ! lmlm-dos qdos
           QDOSOPT = 'y'                                                 ! lmlm-dos qdos
        ENDIF                                                            ! lmlm-dos qdos
        OPEN(67,FORM='formatted',FILE='parameters.gflle')                ! lmlm-dos
        DF(:)=WEZ(:)/DBLE(NSPIN)                                         ! lmlm-dos
        WRITE(67,*) IELAST,IEMXD,NATYP,NSPIN,LMAX,QDOSOPT,DF(1:IELAST),  ! lmlm-dos
     &              EZ(1:IELAST),KORBIT                                  ! lmlm-dos
        CLOSE(67)                                                        ! lmlm-dos
      ENDIF  ! OPT('lmlm-dos')                                           ! lmlm-dos

! -------------------------------------------------------------------------! LLY Lloyd
      IF (LLY.NE.0) THEN                                                   ! LLY Lloyd
                                                                           ! LLY Lloyd
! Calculate free-space contribution to dos                                 ! LLY Lloyd
         CDOS0(1:IEMXD) = CZERO                                            ! LLY Lloyd
         CDOS1(1:IEMXD) = CZERO                                             ! LLY Lloyd
         CDOS2(1:NAEZ) = CZERO                                             ! LLY Lloyd
         DO I1 = 1,NAEZ                                                    ! LLY Lloyd
            CDOSAT0(1:IEMXD) = CZERO                                       ! LLY Lloyd
            CDOSAT1(1:IEMXD) = CZERO                                       ! LLY Lloyd
            ICELL = NTCELL(I1)                                             ! LLY Lloyd
            DO IE = 1,IELAST                                               ! LLY Lloyd
               CALL RHOVAL0(                                               ! LLY Lloyd
     &           EZ(IE),WEZ(IE),                                           ! LLY Lloyd
     &           DRDI(1,I1),RMESH(1,I1),                                   ! LLY Lloyd
     &           IPAN(I1),IRCUT(0,I1),IRWS(I1),                            ! LLY Lloyd
     &           THETAS(1,1,ICELL),LMAX,                                   ! LLY Lloyd
     &           CDOSAT0(IE),CDOSAT1(IE))  

c calculate contribution from free space

                 CDOS2(I1) = CDOS2(I1) + CDOSAT1(IE)*WEZ(IE)                ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
            CDOS0(1:IEMXD) = CDOS0(1:IEMXD) + CDOSAT0(1:IEMXD)             ! LLY Lloyd
            CDOS1(1:IEMXD) = CDOS1(1:IEMXD) + CDOSAT1(1:IEMXD)             ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CDOS0(:) = -CDOS0(:) / PI                                         ! LLY Lloyd
         CDOS1(:) = -CDOS1(:) / PI                                         ! LLY Lloyd
                                                                           ! LLY Lloyd
         if(myrank==master) then
         OPEN(701,FILE='freedos.dat',FORM='FORMATTED')                     ! LLY Lloyd
         DO IE = 1,IELAST                                                  ! LLY Lloyd
            WRITE(701,FMT='(10E16.8)') EZ(IE),CDOS0(IE),CDOS1(IE)          ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
         OPEN(701,FILE='singledos.dat',FORM='FORMATTED')                   ! LLY Lloyd
         DO I1 = 1,NATYP                                                   ! LLY Lloyd
            WRITE(701,FMT='(I5,10E16.8)') I1,CDOS2(I1)                     ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
         end if ! myrank==master
         CDOS_LLY(1:IEMXD,1:NSPIND) = CZERO                                ! LLY Lloyd
         IF (.NOT.OPT('NEWSOSOL')) THEN
         if(t_lloyd%cdos_diff_lly_to_file) then
            OPEN (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')               ! LLY Lloyd
            DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
               DO IE = 1,IELAST                                               ! LLY Lloyd
                  READ(701,FMT='(10E25.16)') EREAD,CDOS_LLY(IE,ISPIN)         ! LLY Lloyd
               ENDDO                                                          ! LLY Lloyd
            ENDDO                                                             ! LLY Lloyd
            CLOSE(701)                                                        ! LLY Lloyd
         else   !(t_lloyd%cdos_diff_lly_to_file)
            DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
#ifdef CPP_MPI
               ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
               ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
               ie_start = 0
               ie_end = IELAST
#endif
               DO ie_num=1,ie_end
                  IE = ie_start+ie_num
                  CDOS_LLY(IE,ISPIN) = t_lloyd%cdos(ie_num,ispin)             ! LLY Lloyd
               ENDDO  !ie_num                                                 ! LLY Lloyd
            ENDDO    !ispin     
#ifdef CPP_MPI            
            ! MPI gather cdos_lly on all processors
            ihelp      = IELAST*NSPIN  !IELAST*NSPIN
            allocate(work(ielast,nspin))
            work = (0.d0, 0.d0)
            CALL MPI_ALLREDUCE(cdos_lly,work,ihelp,
     &      MPI_DOUBLE_COMPLEX,MPI_SUM,t_mpi_c_grid%myMPI_comm_at,ierr)
            call zcopy(ihelp,work,1,cdos_lly,1)
            deallocate(work)
#endif
         end if   !(t_lloyd%cdos_diff_lly_to_file)
                                                                           ! LLY Lloyd
         ! Add free-space contribution cdos0                               ! LLY Lloyd
         DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
            CDOS_LLY(1:IELAST,ISPIN) = CDOS_LLY(1:IELAST,ISPIN) +          ! LLY Lloyd
     &                                CDOS0(1:IELAST)                      ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
                                                                           ! LLY Lloyd
         CHARGE_LLY(1:NSPIND) = 0.D0                                       ! LLY Lloyd
         DO ISPIN = 1,NSPIN                                                ! LLY Lloyd
            CSUM = CZERO                                                   ! LLY Lloyd
            DO IE = 1,IELAST                                               ! LLY Lloyd
               CSUM = CSUM + CDOS_LLY(IE,ISPIN) * WEZ(IE)                  ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
            CHARGE_LLY(ISPIN) = -DIMAG(CSUM) * PI / NSPINPOT               ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
                                                                           ! LLY Lloyd
                                                                           ! LLY Lloyd
         if(myrank==master) then
         OPEN (701,FILE='cdos_lloyd.dat',FORM='FORMATTED')                 ! LLY Lloyd
         DO ISPIN=1,NSPIN                                                  ! LLY Lloyd
            DO IE=1,IELAST                                                 ! LLY Lloyd
               WRITE(701,FMT='(10E16.8)') EZ(IE),CDOS_LLY(IE,ISPIN)        ! LLY Lloyd
            ENDDO                                                          ! LLY Lloyd
         ENDDO                                                             ! LLY Lloyd
         CLOSE(701)                                                        ! LLY Lloyd
         WRITE(*,*) 'Valence charge from Lloyds formula:',                 ! LLY Lloyd
     &        (CHARGE_LLY(ISPIN),ISPIN=1,NSPIN)                            ! LLY Lloyd
         if(t_inc%i_write>0)                                               ! LLY Lloyd
     &      WRITE(1337,*) 'Valence charge from Lloyds formula:',           ! LLY Lloyd
     &        (CHARGE_LLY(ISPIN),ISPIN=1,NSPIN)                            ! LLY Lloyd
         end if ! myrank==master
         
         ELSE !NEWSOSOL
         
          if(t_lloyd%cdos_diff_lly_to_file) then
            OPEN (701,FILE='cdosdiff_lly.dat',FORM='FORMATTED')            ! LLY
            DO IE = 1,IELAST                                               ! LLY
               READ(701,FMT='(10E25.16)') EREAD,CDOS_LLY(IE,1)             ! LLY
            ENDDO                                                          ! LLY
            CLOSE(701)                                                     ! LLY
          else  !(t_lloyd%cdos_diff_lly_to_file)
#ifdef CPP_MPI
               ie_start = t_mpi_c_grid%ioff_pT2(t_mpi_c_grid%myrank_at)
               ie_end   = t_mpi_c_grid%ntot_pT2(t_mpi_c_grid%myrank_at)
#else
               ie_start = 0
               ie_end = IELAST
#endif
            DO ie_num=1,ie_end
               IE = ie_start+ie_num
               CDOS_LLY(IE,1) = t_lloyd%cdos(ie_num,1)
            ENDDO
#ifdef CPP_MPI            
            ! MPI gather cdos_lly on all processors
            ihelp      = IEMXD*NSPIND  !IELAST*NSPIN
            allocate(work(iemxd,nspind))
            work = (0.d0, 0.d0)
            CALL MPI_ALLREDUCE(cdos_lly,work,ihelp,
     &      MPI_DOUBLE_COMPLEX,MPI_SUM,t_mpi_c_grid%myMPI_comm_at,ierr)
            call zcopy(ihelp,work,1,cdos_lly,1)
            deallocate(work)
#endif
          end if  !(t_lloyd%cdos_diff_lly_to_file)
          
c Add free-space contribution cdos0                                     ! LLY
          CDOS_LLY(1:IEMXD,1) = CDOS_LLY(1:IEMXD,1)+2D0*CDOS0(1:IEMXD)  ! LLY 
          CSUM=CZERO                                                    ! LLY
          DO IE = 1,IELAST                                              ! LLY
           CSUM = CSUM + CDOS_LLY(IE,1) * WEZ(IE)                       ! LLY
          ENDDO                                                         ! LLY
          CHARGE_LLY(1) = -DIMAG(CSUM) * PI                             ! LLY
          if(myrank==master) then
          OPEN (701,FILE='cdos_lloyd.dat',FORM='FORMATTED')             ! LLY 
          DO IE=1,IELAST                                                ! LLY
             WRITE(701,FMT='(10E16.8)') EZ(IE),CDOS_LLY(IE,1)           ! LLY 
          ENDDO                                                         ! LLY 
          CLOSE(701)                                                    ! LLY 
          WRITE(*,*) 'Valence charge from Lloyds formula:',
     &        (CHARGE_LLY(1))                                          ! LLY 
          if(t_inc%i_write>0) 
     &       WRITE(1337,*) 'Valence charge from Lloyds formula:',
     &        (CHARGE_LLY(1))                                          ! LLY 
          end if ! myrank==master
          
         ENDIF ! NEWSOSOL

         ENDIF !LLY<>0                                                     ! LLY
! -------------------------------------------------------------------------! LLY Lloyd


      IF (.NOT.OPT('NEWSOSOL')) THEN
C
C ================================================================ NATYP

#ifdef CPP_MPI
        ntot1 = t_inc%NATYP
        call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
     &                   t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                             master,ntot1,ntot_pT,ioff_pT,.true.)
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
      DO I1=i1_start, i1_end
C ----------------------------------------------------------------- SPIN
         IQ = IQAT(I1)
         DO ISPIN = 1,NSPIN
            ICELL = NTCELL(I1)
            IPOT = (I1-1) * NSPINPOT + ISPIN
            IPOT1 = (I1-1) * NSPINPOT + 1
#ifdef CPP_TIMING
            call timing_start('main1c - rhoval')
#endif
            CALL RHOVAL(IHOST,LDORHOEF,ICST,INS,IELAST,
     &           NSRA,ISPIN,NSPIN,NSPINPOT,I1,EZ,WEZ,
     &           DRDI(1,I1),RMESH(1,I1),
     &           VINS(IRMIND,1,KNOSPH*IPOT+(1-KNOSPH)),VISP(1,IPOT),
     &           ZAT(I1),IPAN(I1),IRCUT(0,I1),IRMIN(I1),
     &           THETAS(1,1,ICELL),IFUNM1(1,ICELL),LMSP1(1,ICELL),
     &           RHO2N1(1,1,IPOT1),RHO2N2(1,1,IPOT1),
     &           RHOORB(1,I1),DEN(0,1,1,IPOT),DENLM(1,1,1,IPOT),
     &           MUORB(0,1,I1),ESPV(0,IPOT1),
     &           CLEB,LOFLM,ICLEB,IEND,JEND,SOLVER,
     &           SOCSCL(1,KREL*I1+(1-KREL)),
     &           CSCL(1,KREL*I1+(1-KREL)),
     &           VTREL(1,I1),BTREL(1,I1),
     &           RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
     &           ZREL(I1),JWSREL(I1),IRSHIFT(I1),
     &           LMOMVEC,QMTET(IQ),QMPHI(IQ),MVEVIL1,MVEVIL2,NMVECMAX,
     &           IDOLDAU,LOPT(I1),PHILDAU(1,I1),WLDAU(1,1,1,I1),
     &           DENMATC(1,1,IPOT),
     &           LLY,NATYP,NQDOS)              ! LLY Lloyd
#ifdef CPP_TIMING
            call timing_pause('main1c - rhoval')
#endif
         END DO
C ----------------------------------------------------------------- SPIN
C
        IPOT1 = (I1-1)*NSPINPOT + 1
C
        DO LM = 1,LMPOTD
           DO IR = 1,IRMD
              RHO2NS(IR,LM,I1,1) = RHO2N1(IR,LM,IPOT1)
              R2NEF(IR,LM,I1,1)  = RHO2N2(IR,LM,IPOT1)
           END DO
        END DO

C     
        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +           DIMAG(DEN(L,IELAST,1,IPOT1))/PI/DBLE(NSPINPOT)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &           DIMAG(DEN(L,IELAST,1,IPOT1))/PI/DBLE(NSPINPOT)
        END DO
C     
        IF (NSPINPOT.EQ.2) THEN
           DO LM = 1,LMPOTD
              DO IR = 1,IRMD
                 RHO2NS(IR,LM,I1,2) = RHO2N1(IR,LM,IPOT1+1)
                 R2NEF(IR,LM,I1,2)  = RHO2N2(IR,LM,IPOT1+1)
              END DO
           END DO
C     
           DO L = 0,LMAXD1
              DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +              DIMAG(DEN(L,IELAST,1,IPOT1+1))/PI/DBLE(NSPINPOT)
              DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &              DIMAG(DEN(L,IELAST,1,IPOT1+1))/PI/DBLE(NSPINPOT) 
           END DO
        END IF

        IF ( TEST('RHOVALW ') ) THEN !Bauer
        open(unit=324234,file='out_rhoval')
        WRITE(324234,*) '#IATOM',I1
        write(324234,'(50000F)') RHO2NS(:,:,I1,1)
        IF (NSPIN==2) write(324234,'(50000F)') RHO2NS(:,:,I1,2)
        END IF


C ----------------------------------------------------------------------
C        
C --> itermdir/kmrot <> 0 
C
        IF (LMOMVEC) THEN
           DO IS = 1, NMVECMAX
              DO LM=1, 3
                 MVEVI(I1,LM,IS) = (0.0D0,0.0D0)
                 MVEVIEF(I1,LM,IS) = (0.0D0,0.0D0) 
                 DO L = 0, LMAXD
                    MVEVIL(L,I1,LM,IS) = MVEVIL1(L,LM,IS)
                    MVEVI(I1,LM,IS) = MVEVI(I1,LM,IS) +
     &                   MVEVIL1(L,LM,IS)
C     
                    MVEVIEF(I1,LM,IS) = MVEVIEF(I1,LM,IS) + 
     &                   MVEVIL2(L,LM,IS)
                 END DO
              END DO
           END DO
        END IF
C ----------------------------------------------------------------------
      END DO
C ================================================================ NATYP
#ifdef CPP_TIMING
      if(ie_end>=ie_start) call timing_stop('main1c - rhoval')
#endif
      CLOSE (69) !gmat file
#ifndef CPP_MPI   
      CLOSE (30) !close lmdos file if no mpi is used
      CLOSE (31) !close qdos file if no mpi is used, otherwise writeout in the following
#endif
      CLOSE (96) !close gflle file
      
      
#ifdef CPP_MPI
      !move writeout of qdos file here
      IF (OPT('qdos    ')) THEN
       ! first communicate den array to write out qdos files
       IDIM = (LMAXD1+1)*IEMXD*NQDOS*NPOTD
       allocate(workc(0:LMAXD1,IEMXD,NQDOS,NPOTD))
       workc = (0.d0, 0.d0)
       
!        if ( (t_mpi_c_grid%nranks_ie>t_mpi_c_grid%dims(1))
!      &   .and. (t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1)) ) then
!          den = 0.d0
!        endif
       
       CALL MPI_ALLREDUCE(DEN,workc,IDIM,MPI_DOUBLE_COMPLEX,MPI_SUM,
     &                    MPI_COMM_WORLD,IERR)
       CALL ZCOPY(IDIM,WORKC,1,DEN,1)
       deallocate(workc)
         
       if(myrank==master) then
       
       DO I1=1,NATYP
      
        DO ISPIN = 1,NSPIN
        
         if(NATYP.ge.100) then
           OPEN(31,
     +     FILE="qdos."//char(48+I1/100)//char(48+mod(I1/10,10))//
     +     char(48+mod(I1,10))//"."//char(48+ISPIN)//".dat")
         else
           OPEN(31,
     +     FILE="qdos."//char(48+I1/10)//char(48+mod(I1,10))//"."//
     +     char(48+ISPIN)//".dat")
         end if
         
         IPOT = (I1-1) * NSPINPOT + ISPIN

         DO IE=1,IELAST
           DO IQ=1,NQDOS
             DENTOT = DCMPLX(0.D0,0.D0)
             DO L = 0,LMAXD1
               DENTOT = DENTOT + DEN(L,IE,IQ,IPOT)
             ENDDO
            WRITE(31,9000) EZ(IE),QVEC(1,IQ),QVEC(2,IQ),QVEC(3,IQ),      ! qdos
     &                    -DIMAG(DENTOT)/PI,
     &                   (-DIMAG(DEN(L,IE,IQ,IPOT))/PI,L=0,LMAXD1)      ! qdos
           END DO ! IQ=1,NQDOS
         END DO ! IE=1,IELAST
         close(31)
        END DO !ISPIN=1,NSPIN
       ENDDO !I1
       end if ! myrank_at==master
      ENDIF ! OPT('qdos    ')
 9000 FORMAT(5F10.6,40E16.8)
 
      !reset NQDOS number to avoid loo large communication which is not needed anyways for qdos run
      NQDOS=1
      
#ifdef CPP_TIMING
      call timing_start('main1c - communication')
#endif
      !MPI: reduce these arrays, so that all processors have the results (MPI_ALLREDUCE)
#ifdef CPP_TIMING
      call timing_start('main1c - communication')
#endif
      call mympi_main1c_comm(IRMD,LMPOTD,NATYPD,LMAXD,LMAXD1,LMMAXD,
     &              NPOTD,IEMXD,MMAXD,IDOLDAU,NATYP,KREL,LMOMVEC,
     &              NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,
     &              denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,
     &              t_mpi_c_grid%mympi_comm_ie)
      call mympi_main1c_comm(IRMD,LMPOTD,NATYPD,LMAXD,LMAXD1,LMMAXD,
     &              NPOTD,IEMXD,MMAXD,IDOLDAU,NATYP,KREL,LMOMVEC,
     &              NMVECMAX,NQDOS,rho2ns,r2nef,espv,den,denlm,denmatc,
     &              denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,
     &              t_mpi_c_grid%mympi_comm_at)
#ifdef CPP_TIMING
      call timing_stop('main1c - communication')
#endif

!     lmdos writeout
      if(myrank==master) then                                            ! lm-dos
c        IF (.not.OPT('qdos    ')) THEN                                   ! lm-dos
        IF (OPT('lmdos    ')) THEN                                   ! lm-dos
          DO I1=1,NATYP                                                  ! lm-dos
            DO ISPIN = 1,NSPIN                                           ! lm-dos
              IPOT = (I1-1) * NSPINPOT + ISPIN                           ! lm-dos
              if(NATYP.ge.100) then                                      ! lm-dos
                OPEN(30,                                                 ! lm-dos
     +               FILE="lmdos."//char(48+I1/100)//                    ! lm-dos
     +                     char(48+mod(I1/10,10))//                      ! lm-dos
     +                     char(48+mod(I1,10))//"."//                    ! lm-dos
     +                     char(48+ISPIN)//".dat")                       ! lm-dos
              else                                                       ! lm-dos
                OPEN(30,                                                 ! lm-dos
     +               FILE="lmdos."//char(48+I1/10)//                     ! lm-dos
     +                     char(48+mod(I1,10))//"."//                    ! lm-dos
     +                     char(48+ISPIN)//".dat")                       ! lm-dos
              end if                                                     ! lm-dos
              WRITE (30,*) ' '                                           ! lm-dos
              WRITE (30,8600) '# ISPIN=',ISPIN,' I1=',I1                 ! lm-dos
 8600         FORMAT (a8,I3,a4,I5)                                       ! lm-dos
              DO IE=1,IELAST                                             ! lm-dos
                WRITE(30,9001) DREAL(EZ(IE)),                            ! lm-dos
     &               (-DIMAG(DENLM(LM,IE,1,IPOT))/PI,LM=1,LMMAXD)        ! lm-dos
              END DO ! IE                                                ! lm-dos
            END DO ! ISPIN                                               ! lm-dos
 9001       FORMAT(30E12.4)                                              ! lm-dos
            close(30)                                                    ! lm-dos
          END DO !I1                                                     ! lm-dos
        ENDIF  ! not qdos option                                         ! lm-dos
      end if ! myrank==master                                            ! lm-dos
#endif

      ELSE ! new spin-orbit solver

      ! nonco angles
      call read_angles(t_params,NATYP,THETA,PHI)


c interpolate potential
      IF ( IDOLDAU.EQ.1 ) 
     +           CALL CINIT(MMAXD*MMAXD*4*NATYPD,DENMATN(1,1,1,1,1))

       CALL INTERPOLATE_POTEN(NSPIN,RMESH,IRMIN,IRWS,IPAN,IRCUT,VINS,
     +                        VISP,NPAN_LOG,NPAN_EQ,NCHEB,NPAN_TOT,
     +                        RNEW,RPAN_INTERVALL,IPAN_INTERVALL,
     +                        VINSNEW)
#ifdef CPP_MPI
        ntot1 = t_inc%NATYP
        
        
!         call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie, 
        if(t_mpi_c_grid%dims(1)>1) then
           nranks_local = t_mpi_c_grid%nranks_ie
           if(t_mpi_c_grid%nranks_ie>t_mpi_c_grid%dims(1)) 
     &              nranks_local = t_mpi_c_grid%dims(1)
        else
           nranks_local = 1
        endif
        call distribute_linear_on_tasks(nranks_local, 
     &                   t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,
     &                             master,ntot1,ntot_pT,ioff_pT,.true.,
     &                             .true.)
        if(t_mpi_c_grid%nranks_ie<=t_mpi_c_grid%dims(1)) then
           i1_start = ioff_pT(t_mpi_c_grid%myrank_ie)+1
           i1_end   = ioff_pT(t_mpi_c_grid%myrank_ie)+
     &                ntot_pT(t_mpi_c_grid%myrank_ie)
           t_mpi_c_grid%ntot1  = ntot_pT(t_mpi_c_grid%myrank_ie)
        else
           myrank_ie_tmp = t_mpi_c_grid%myrank_ie
           if(t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1))
     &       myrank_ie_tmp = myrank_ie_tmp - t_mpi_c_grid%myrank_ie
     +                          + (t_mpi_c_grid%dims(1)-1)
           i1_start = ioff_pT(myrank_ie_tmp)+1
           i1_end   = ioff_pT(myrank_ie_tmp)+ntot_pT(myrank_ie_tmp)
           t_mpi_c_grid%ntot1  = ntot_pT(myrank_ie_tmp)
        endif

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

      !initialize output nonco angles
      angles_new = 0.0d0

      DO I1 = i1_start,i1_end
      
       ICELL = NTCELL(I1)
       IPOT = (I1-1) * NSPIN + 1

#ifdef CPP_TIMING
        call timing_start('main1c - rhovalnew')
#endif

        CALL RHOVALNEW(LDORHOEF,IELAST,NSRA,NSPIN,LMAX,EZ,WEZ,
     &           ZAT(I1),SOCSCALE(I1),CLEB(1,1),ICLEB,IEND,
     &           IFUNM1(1,ICELL),LMSP1(1,ICELL),
     &           NCHEB,NPAN_TOT(I1),NPAN_LOG(I1),
     &           NPAN_EQ(I1),RMESH(1,I1),IRWS(I1),RPAN_INTERVALL(0,I1),
     &           IPAN_INTERVALL(0,I1),RNEW(1,I1),VINSNEW,
     &           THETASNEW(1,1,ICELL),THETA(I1),PHI(I1),I1,IPOT,
     &           DEN1(0,1,1),ESPV1(0,1),RHO2M1,RHO2M2,MUORB(0,1,I1),
     &           angles_new(:,i1),
     &           IDOLDAU,LOPT(I1),PHILDAU(1,I1),WLDAU(1,1,1,I1),   ! LDAU
     &           DENMATN(1,1,1,1,I1))                              ! LDAU
#ifdef CPP_TIMING
        call timing_pause('main1c - rhovalnew')
#endif

        DO L = 0,LMAXD1
         ESPV(L,IPOT)=ESPV1(L,1)
         ESPV(L,IPOT+NSPIN-1)=ESPV1(L,NSPIN)
         DO IE=1,IELAST
          DEN(L,IE,1,IPOT)=DEN1(L,IE,1)
          DEN(L,IE,1,IPOT+NSPIN-1)=DEN1(L,IE,NSPIN)
         ENDDO
        ENDDO
        DO ISPIN=1,NSPIN
         DO LM = 1,LMPOTD
          DO IR = 1,IRMD
           RHO2NS(IR,LM,I1,ISPIN) = RHO2M1(IR,LM,ISPIN)
           R2NEF(IR,LM,I1,ISPIN)  = RHO2M2(IR,LM,ISPIN)
          END DO
         END DO
        END DO

        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +           DIMAG(DEN(L,IELAST,1,IPOT))/PI/DBLE(NSPIN)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &           DIMAG(DEN(L,IELAST,1,IPOT))/PI/DBLE(NSPIN)
        END DO
        DO L = 0,LMAXD1
           DENEF = DENEF - 2.0D0 * CONC(I1) * 
     +            DIMAG(DEN(L,IELAST,1,IPOT+1))/PI/DBLE(NSPIN)
           DENEFAT(I1) = DENEFAT(I1) - 2.0D0 * 
     &            DIMAG(DEN(L,IELAST,1,IPOT+1))/PI/DBLE(NSPIN) 
        END DO
       DO ISPIN=1,NSPIN
        DO L=0,LMAXD1
         MUORB(L,3,I1)=MUORB(L,3,I1)+MUORB(L,ISPIN,I1)
        ENDDO
       ENDDO
       DO ISPIN=1,3
        DO L=0,LMAXD1
         MUORB(LMAXD1+1,ISPIN,I1)=MUORB(LMAXD1+1,ISPIN,I1)+
     &                            MUORB(L,ISPIN,I1)
        ENDDO
       ENDDO
      END DO !I1
      
#ifdef CPP_MPI
      !reset NQDOS to 1 to avoid endless communication
      NQDOS = 1
      call mympi_main1c_comm_newsosol2(LMAXD1,LMMAXD,IEMXD,NQDOS,         &
     &                                 NPOTD,NATYPD,LMPOTD,IRMD,MMAXD,    &
     &                                 den, denlm, muorb, espv, r2nef,    &
     &                                 rho2ns, denefat, denef,denmatn,    &
     &                            angles_new,t_mpi_c_grid%mympi_comm_ie)
#endif
      
#ifdef CPP_TIMING
      call timing_stop('main1c - rhovalnew')
#endif

      CLOSE (69)
      
      
      if(myrank==master) then
c rewrite new theta and phi to nonco_angle_out.dat, noncao_angle.dat is the input
        IF (.NOT.TEST('FIXMOM  ')) THEN
          OPEN(UNIT=13,file='nonco_angle_out.dat',form='formatted')
          DO I1=1,NATYP
            ! save to file in converted units (degrees)
            WRITE(13,*) angles_new(1,I1)/(2.0D0*PI)*360.0D0,
     &                  angles_new(2,I1)/(2.0D0*PI)*360.0D0
            ! use internal units here
            t_params%THETA(I1) = angles_new(1,I1)
            t_params%PHI(I1)   = angles_new(2,I1)
          ENDDO
         CLOSE(13)

        ENDIF ! .not.test('FIXMOM  ')
      end if !(myrank==master)

      CLOSE (29) !close lmdos file
      CLOSE (30) !close lmdos file
      CLOSE (31) !close qdos file
      CLOSE (32) !close qdos file
      CLOSE (91) !close gflle file
      
      ENDIF !NEWSOSOL

#ifdef CPP_MPI
      if(myrank==master) then
#endif
#ifdef CPP_TIMING
      call timing_start('main1c - serial part')
#endif

! In case of Lloyds formula renormalize valence charge                   ! LLY Lloyd
      IF (LLY.GT.0) THEN                                                 ! LLY Lloyd
         LMAXP1 = LMAX                                                   ! LLY Lloyd
         IF (INS.NE.0) LMAXP1 = LMAX + 1                                 ! LLY Lloyd
         CALL RENORM_LLY(
     &     CDOS_LLY,IELAST,NSPIN,NATYP,DEN,LMAXP1,CONC,
     &     1,IELAST,WEZ,IRCUT,IPAN,EZ,ZAT,
     &     RHO2NS,CDOS1,CDOS2,R2NEF,DENEF,DENEFAT,ESPV)


      ENDIF                                                              ! LLY Lloyd

C
C================================================================ NATYP
      CHRGSEMICORE = 0D0
      DO I1 = 1,NATYP
C
C ---> l/m_s/atom-resolved charges 
C
          DO ISPIN = 1,NSPINPOT
             IPOT = (I1-1)*NSPINPOT + ISPIN
             DO L = 0,LMAXD1
                CHARGE(L,I1,ISPIN) = 0.0D0
C
                DO IE = 1,IELAST
                   CHARGE(L,I1,ISPIN) = CHARGE(L,I1,ISPIN) +
     +                  DIMAG(WEZ(IE)*DEN(L,IE,1,IPOT))/
     +                  DBLE(NSPINPOT)
                   IF ( IE.EQ.IESEMICORE ) CHRGSEMICORE = 
     +                  CHRGSEMICORE + CONC(I1)*CHARGE(L,I1,ISPIN)
                END DO
C
             END DO
          END DO
          EU(I1) = 0D0
          EDC(I1) = 0D0
C
C ---> orbital magnetic moments (array initialised to 0.0D0 in rhoval)
C
          IF (KREL.EQ.1) THEN
             DO ISPIN = 1,3
                DO L = 0,LMAXD + 1
                   MUORB(LMAXD1+1,ISPIN,I1) = MUORB(LMAXD1+1,ISPIN,I1)
     &                                       + MUORB(L,ISPIN,I1)
                END DO
             END DO
          END IF
      END DO
C================================================================ NATYP
C
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C LDA+U
      IF ( IDOLDAU.EQ.1 ) THEN
C -> Save old LDA+U interaction matrix for mixing
         CALL DCOPY(MMAXD*MMAXD*NSPIND*NATYPD,WLDAU,1,WLDAUOLD,1)
C
C -> Construct LDA+U interaction matrix for next iteration
C
        IF (.NOT.OPT('NEWSOSOL')) THEN
         CALL WMATLDAU(NTLDAU,ITLDAU,NSPINPOT,DENMATC,LOPT,
     &                 UEFF,JEFF,ULDAU,WLDAU,EU,EDC,MMAXD,NPOTD)
        ELSE
         CALL WMATLDAUSOC(NTLDAU,ITLDAU,NSPINPOT,DENMATN,LOPT,
     &                 UEFF,JEFF,ULDAU,WLDAU,EU,EDC,MMAXD,NPOTD)
        ENDIF
C -> Mix old and new LDA+U interaction matrices
         CALL MIXLDAU(
     >        MMAXD,NSPIND,NATYPD,NATYP,NSPIN,LOPT,WLDAUOLD,
     X        WLDAU)
C
C -> update variables-file
C
         ITRUNLDAU = ITRUNLDAU + 1
         OPEN (67,FILE='ldau.unformatted',FORM='unformatted')
         WRITE (67) ITRUNLDAU,WLDAU,ULDAU,PHILDAU
         CLOSE(67)
C
C -> write full lda+u information in ascii file ldaupot_new
C
         CALL WRLDAUPOT(ITRUNLDAU,LOPT,UEFF,JEFF,EREFLDAU,NATYP,WLDAU,
     &                  ULDAU,PHILDAU,IRMD,NATYPD,NSPIND,MMAXD,IRWS)
      END IF
C
C LDA+U
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

       CALL WRMOMS(KREL+KORBIT,NATYP,NSPINPOT,TEXTS,TEXTL,TEXTNS,CHARGE,
     &             MUORB,LMAXD,LMAXD1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      ITERMDIR
C
      IF ( ( KREL.EQ.1 ) .AND. LMOMVEC ) THEN
         DO I1=1,NATYP
            IQ = IQAT(I1)
            CALL MVECGLOBAL(I1,IQ,NATYP,QMPHI(IQ),QMTET(IQ),
     &                      MVEVI,MVEVIL,MVEVIEF,
     &                      NATYPD,LMAXD,NMVECMAX)
         END DO
      END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
c    TEST BRAHIM
      IF (NPOL.EQ.0 .OR.TEST('DOS     ')) CALL WRLDOS(DEN,EZ,WEZ,
     +    LMAXD1,IEMXD,NPOTD,ITITLE,EFERMI,E1,E2,ALAT,TK,NACLS1,
     +    NSPINPOT,NATYP,CONC,IELAST,INTERVX,INTERVY,INTERVZ,DOSTOT)


C
C ---------------------------------------------------------- CORE STATES
C   RHO_core is calculated only if also RHO_valence was
C
      IF (NPOL.NE.0) THEN 
         if(t_inc%i_write>0) then
         WRITE (1337,*)
         WRITE (1337,'(78(1H#))')
         WRITE (1337,'(33X,A)') 'CORE  STATES'
         WRITE (1337,'(78(1H#))')
         endif
         DO I1 = 1,NATYP
            DO ISPIN = 1,NSPIN
               IPOT = (I1-1) * NSPINPOT + ISPIN
               IPOT1 = (I1-1) * NSPINPOT + 1
C
               CALL RHOCORE(NSRA,ISPIN,NSPIN,I1,
     &              DRDI(1,I1),RMESH(1,I1),VISP(1,IPOT),
     &              A(I1),B(I1),ZAT(I1),
     &              IRCUT(0,I1),RHOC(1,IPOT1),
     &              ECORE(1,IPOT),NCORE(IPOT),LCORE(1,IPOT),
     &              CSCL(1,KREL*I1+(1-KREL)),
     &              VTREL(1,I1),BTREL(1,I1),
     &              RMREL(1,I1),DRDIREL(1,I1),R2DRDIREL(1,I1),
     &              ZREL(I1),JWSREL(I1),IRSHIFT(I1),
     &              ECOREREL(1,IPOT1),NKCORE(1,I1),
     &              KAPCORE(1,IPOT1))
C
            END DO
         END DO
         if(t_inc%i_write>0) then
         WRITE (1337,*)
         WRITE (1337,'(78(1H#))')
         WRITE (1337,*)
         endif
      END IF
C
C ---------------------------------------------------------- CORE STATES
C
      call save_density(t_params,RHO2NS,R2NEF,RHOC,DENEF,DENEFAT,
     &                  ESPV,ECORE,IDOLDAU,LOPT,EU,EDC,
     &                  CHRGSEMICORE,RHOORB,ECOREREL,NKCORE,KAPCORE,
     &                  KREL,NATYPD,NPOTD,IRMD,LMPOTD,LMAXD1)
      
      IF (TEST('den-asci')) THEN
         OPEN (67,FILE='densitydn.ascii',FORM='formatted')
         DO I1 = 1,NATYP
            DO LM = 1,LMPOTD
               DO IR = 1,IRMD
                  WRITE(67,FMT='(I6,2I5,2E25.16)') 
     &               I1,LM,IR,RHO2NS(IR,LM,I1,1),RHO2NS(IR,LM,I1,2)
               ENDDO
            ENDDO
         ENDDO
         CLOSE(67)
      ENDIF

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C      ITERMDIR
C
      IF (ITERMVDIR) THEN
         t_params%MVEVI   = MVEVI
         t_params%MVEVIEF = MVEVIEF
      END IF
C
C      ITERMDIR
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

C
      if(t_inc%i_write>0) 
     &   WRITE (1337,'(79(1H=),/,30X,"< KKR1c finished >",/,79(1H=),/)')
     
#ifdef CPP_TIMING
      call timing_stop('main1c - serial part')
#endif
      if(test('rhoqtest')) then

        open(9999, file='cleb_shapefun.txt')
        write(9999,*) 'params'
        write(9999,*) t_params%ncleb, t_params%lm2d, t_params%lmaxd,
     &                t_params%iend, t_params%nfund, t_params%lmxspd,
     &                nrmaxd
        write(9999,*) 'cleb'
        do lm1=1,t_params%ncleb
          write(9999,'(E16.7,4I9)') t_params%cleb(lm1,1), 
     &                              t_params%icleb(lm1,1:4)
        end do
        write(9999,*) 'loflm'
        do lm1=1,t_params%lm2d
          write(9999,'(I9)') t_params%loflm(lm1)
        end do
        open(8888, file='mu0')
        read(8888,*) i1
        close(8888)
        ICELL = NTCELL(I1)
        write(9999,*) 'ifunm'
        do lm1=1,lmxspd
          write(9999,'(I9)') ifunm1(lm1,icell)
        end do
        write(9999,*) 'lmsp'
        do lm1=1,lmxspd
          write(9999,'(I9)') LMSP1(LM1,icell)
        end do
        write(9999,*) 'thetasnew',nrmaxd,nfund
        do lm1=1,nrmaxd
          do lm2=1,nfund
            write(9999,'(E16.7)') thetasnew(lm1,lm2,icell)
          end do
        end do
        close(9999)
        open(987987, file='shapefun_test')
        write(987987,*) thetasnew
        close(987987)

      endif
#ifdef CPP_MPI
      end if !myrank==master
#endif
      
      deallocate(den,denlm)
      
      deallocate(thetasnew, thetas)
      deallocate(vins, vinsnew)
      deallocate(RHO2N1,RHO2NS,RHO2N2,R2NEF)
      
      
      END SUBROUTINE !MAIN1c
      
      END MODULE
