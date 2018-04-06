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

   implicit none

contains

   !----------------------------------------------------------------------------
   ! SUBROUTINE: main1a
   !> @brief Main subroutine regarding the calculation of the t-matrix
   !> @author Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
   !> and many others ...
   !----------------------------------------------------------------------------
   subroutine main1a(INS,LLY,IRM,LM2D,ICST,IEND,NCLS,LMAX,NREF,NSRA,KREL,NEMB,   &
      NAEZ,NATYP,NCLSD,NPOTD,ITSCF,NTOTD,MMAXD,LMPOT,IPAND,NINEQ,NSPIN,NCHEB,    &
      LMMAXD,IELAST,NRMAXD,IRMIND,NATOMIMP,ALAT,R_LOG,TOLRDIF,DELTAE,CLS,IQAT,   &
      IRWS,NACLS,REFPOT,ATOM,ZAT,VREF,RMTREF,RCLS,SOLVER,SOCSCL,SOCSCALE,CSCL,   &
      NTLDAU,IDOLDAU,ITLDAU,UEFF,JEFF,IPAN,LOFLM,IRMIN,ATOMIMP,ICLEB,IRCUT,      &
      IPAN_INTERVALL,PHI,THETA,CLEB,VISP,DRDI,RNEW,RMESH,RPAN_INTERVALL,VINS,    &
      ZREL,JWSREL,VTREL,BTREL,RMREL,DRDIREL,R2DRDIREL,ITRUNLDAU,LOPT,EREFLDAU,   &
      WLDAU,ULDAU,PHILDAU)

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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For KREL = 1 (relativistic mode)
      !
      !  NPOTD = 2 * NATYPD
      !  LMMAXD = 2 * (LMAX+1)^2
      !  NSPIND = 1
      !  LMGF0D = (LMAX+1)^2 dimension of the reference system Green
      !          function, set up in the spin-independent non-relativstic
      !          (l,m_l)-representation
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !     ..
      !     .. Input variables
      integer, intent(in) :: INS       !< 0 (MT), 1(ASA), 2(Full Potential)
      integer, intent(in) :: LLY       !< LLY <> 0: apply Lloyds formula
      integer, intent(in) :: IRM       !< Maximum number of radial points
      integer, intent(in) :: LM2D      !< (2*LMAX+1)**2
      integer, intent(in) :: ICST      !< Number of Born approximation
      integer, intent(in) :: IEND      !< Number of nonzero gaunt coefficients
      integer, intent(in) :: NCLS      !< Number of reference clusters
      integer, intent(in) :: LMAX      !< Maximum l component in wave function expansion
      integer, intent(in) :: NREF      !< Number of diff. ref. potentials
      integer, intent(in) :: NSRA
      integer, intent(in) :: KREL      !< Switch for non-relativistic/relativistic (0/1) program. Attention: several other parameters depend explicitly on KREL, they are set automatically Used for Dirac solver in ASA
      integer, intent(in) :: NEMB      !< Number of 'embedding' positions
      integer, intent(in) :: NAEZ      !< Number of atoms in unit cell
      integer, intent(in) :: NATYP     !< Number of kinds of atoms in unit cell
      integer, intent(in) :: NCLSD     !< Maximum number of different TB-clusters
      integer, intent(in) :: NPOTD     !< (2*(KREL+KORBIT)+(1-(KREL+KORBIT))*NSPIND)*NATYP)
      integer, intent(in) :: ITSCF
      integer, intent(in) :: NTOTD
      integer, intent(in) :: MMAXD     !< 2*LMAX+1
      integer, intent(in) :: LMPOT     !< (LPOT+1)**2
      integer, intent(in) :: IPAND     !< Number of panels in non-spherical part
      integer, intent(in) :: NINEQ     !< Number of ineq. positions in unit cell
      integer, intent(in) :: NSPIN     !< Counter for spin directions
      integer, intent(in) :: NCHEB     !< Number of Chebychev pannels for the new solver
      integer, intent(in) :: LMMAXD    !< (KREL+KORBIT+1)(LMAX+1)^2
      integer, intent(in) :: IELAST
      integer, intent(in) :: NRMAXD    !< NTOTD*(NCHEB+1)
      integer, intent(in) :: IRMIND    !< IRM-IRNSD
      integer, intent(in) :: NATOMIMP  !< Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
      double precision, intent(in) :: ALAT      !< Lattice constant in a.u.
      double precision, intent(in) :: R_LOG     !< Radius up to which log-rule is used for interval width. Used in conjunction with runopt NEWSOSOL
      double precision, intent(in) :: TOLRDIF   !< For distance between scattering-centers smaller than [<TOLRDIF>], free GF is set to zero. Units are Bohr radii.
      double complex, intent(in) :: DELTAE  !< Energy difference for numerical derivative
      integer, dimension(NAEZ+NEMB), intent(in)          :: CLS   !< Cluster around atomic sites
      integer, dimension(NATYP), intent(in)              :: IQAT  !< The site on which an atom is located on a given site
      integer, dimension(NATYP), intent(in)              :: IRWS  !< R point at WS radius
      integer, dimension(NCLSD), intent(in)              :: NACLS !< Number of atoms in cluster
      integer, dimension(NAEZ+NEMB), intent(in)          :: REFPOT   !< Ref. pot. card  at position
      integer, dimension(NACLSD,NAEZ+NEMB), intent(in)   :: ATOM  !< Atom at site in cluster
      double precision, dimension(NATYP), intent(in)  :: ZAT      !< Nuclear charge
      double precision, dimension(NREF), intent(in)   :: VREF
      double precision, dimension(NREF), intent(in)   :: RMTREF   !< Muffin-tin radius of reference system
      double precision, dimension(3,NACLSD,NCLSD), intent(in) :: RCLS   !< Real space position of atom in cluster

      !-------------------------------------------------------------------------
      !     RELATIVISTIC MODE
      !-------------------------------------------------------------------------
      character(len=10), intent(in) :: SOLVER   !< Type of solver
      double precision, dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(in)   :: SOCSCL
      double precision, dimension(NATYP), intent(in)                             :: SOCSCALE    !< Spin-orbit scaling
      double precision, dimension(KREL*LMAX+1,KREL*NATYP+(1-KREL)), intent(in)   :: CSCL        !< Speed of light scaling
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      integer, intent(in) :: NTLDAU    !< number of atoms on which LDA+U is applied
      integer, intent(in) :: IDOLDAU   !< flag to perform LDA+U
      integer, dimension(NATYP), intent(in) :: ITLDAU !< integer pointer connecting the NTLDAU atoms to heir corresponding index in the unit cell
      double precision, dimension(NATYP), intent(in) :: UEFF   !< input U parameter for each atom
      double precision, dimension(NATYP), intent(in) :: JEFF   !< input J parameter for each atom

      ! .. In/out Variables
      integer, dimension(NATYP), intent(inout)                             :: IPAN !< Number of panels in non-MT-region
      integer, dimension(LM2D), intent(inout)                              :: LOFLM !< l of lm=(l,m) (GAUNT)
      integer, dimension(NATYP), intent(inout)                             :: IRMIN !< Max R for spherical treatment
      integer, dimension(NATOMIMP), intent(inout)                          :: ATOMIMP
      integer, dimension(NCLEB,4), intent(inout)                           :: ICLEB !< Pointer array
      integer, dimension(0:IPAND,NATYP), intent(inout)                     :: IRCUT !< R points of panel borders
      integer, dimension(0:NTOTD,NATYP), intent(inout)                     :: IPAN_INTERVALL
      double precision, dimension(NATYP), intent(inout)                    :: PHI
      double precision, dimension(NATYP), intent(inout)                    :: THETA
      double precision, dimension(NCLEB,2), intent(inout)                  :: CLEB  !< GAUNT coefficients (GAUNT)
      double precision, dimension(IRM,NPOTD), intent(inout)                :: VISP  !< Spherical part of the potential
      double precision, dimension(IRM,NATYP), intent(inout)                :: DRDI  !< Derivative dr/di
      double precision, dimension(NRMAXD,NATYP), intent(inout)             :: RNEW
      double precision, dimension(IRM,NATYP), intent(inout)                :: RMESH !< Radial mesh ( in units a Bohr)
      double precision, dimension(0:NTOTD,NATYP), intent(inout)            :: RPAN_INTERVALL
      double precision, dimension(IRMIND:IRM,LMPOT,NSPOTD), intent(inout)  :: VINS !< Non-spherical part of the potential
      double complex, dimension(IEMXD), intent(inout) :: EZ
      !-------------------------------------------------------------------------
      !     RELATIVISTIC MODE
      !-------------------------------------------------------------------------
      integer, dimension(NATYP), intent(inout) :: ZREL      !< atomic number (cast integer)
      integer, dimension(NATYP), intent(inout) :: JWSREL    !< index of the WS radius
      double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(inout) :: VTREL       !< potential (spherical part)
      double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(inout) :: BTREL       !< magnetic field
      double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(inout) :: RMREL       !< radial mesh
      double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(inout) :: DRDIREL     !< derivative of radial mesh
      double precision, dimension(IRM*KREL+(1-KREL),NATYP), intent(inout) :: R2DRDIREL   !< \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
      !-------------------------------------------------------------------------
      ! LDA+U
      !-------------------------------------------------------------------------
      integer, intent(inout) :: ITRUNLDAU !< Iteration index for LDA+U
      integer, dimension(NATYP), intent(inout) :: LOPT   !< angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
      double precision, dimension(NATYP), intent(inout) :: EREFLDAU  !< the energies of the projector's wave functions (REAL)
      double precision, dimension(MMAXD,MMAXD,NSPIND,NATYP), intent(inout) :: WLDAU !< potential matrix
      double precision, dimension(MMAXD,MMAXD,MMAXD,MMAXD,NATYP), intent(inout) :: ULDAU  !< calculated Coulomb matrix elements (EREFLDAU)
      double complex, dimension(IRM,NATYP), intent(inout) :: PHILDAU

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
      integer, dimension(NATYP) :: NPAN_EQ
      integer, dimension(NATYP) :: NPAN_LOG
      integer, dimension(NATYP) :: NPAN_TOT
      double precision, dimension(:,:,:), allocatable :: VINSNEW

      ! Assignment of values to parameters
      parameter (LRECTMT=WLENGTH*4*LMMAXD*LMMAXD)
      parameter (LRECTRA=WLENGTH*4)

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
      external TBREF,CALCTMAT,OPT
      !-------------------------------------------------------------------------
      !     .. Intrinsic Functions ..
      !-------------------------------------------------------------------------
      intrinsic ATAN,MOD
      !     ..
      data TOLRDIF /1.5D0/ ! Set free GF to zero if R<TOLRDIF in case of virtual atoms
      data LLY /0/
      !

      allocate(VINSNEW(NRMAXD,LMPOT,NSPOTD),stat=i_stat)
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
      call get_params_1a(t_params,IPAND,NATYP,IRM,NACLSD,IELAST,   &
         NCLSD,NREF,NCLEB,NEMB,NAEZ,LM2D,NSRA,INS,NAEZ,NATYP,     &
         NSPIN,ICST,IPAN,IRCUT,LMAX,NCLS,NINEQ,NREF,IDOLDAU,LLY,     &
         KREL,ATOM,CLS,ICLEB,LOFLM,NACLS,REFPOT,IRWS,IEND,EZ,VINS,   &
         IRMIN,ITMPDIR,ILTMP,ALAT,DRDI,RMESH,ZAT,RCLS,IEMXD,VISP,    &
         RMTREF,VREF,CLEB,CSCL,SOCSCALE,SOCSCL,EREFLDAU,UEFF,JEFF,   &
         SOLVER,TMPDIR,DELTAE,tolrdif,NPAN_LOG,NPAN_EQ,              &
         NCHEB,NPAN_TOT,IPAN_INTERVALL,RPAN_INTERVALL,RNEW,LMAX,    &
         NTOTD,NRMAXD,R_LOG,NTLDAU,ITLDAU,LOPT,VTREL,BTREL,DRDIREL,  &
         R2DRDIREL,RMREL,IRMIND,LMPOT,NSPOTD,NPOTD,JWSREL,ZREL,     &
         ITSCF,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
      !
      if ( TEST('Vspher  ') ) VINS(IRMIND:IRM,2:LMPOT,1:NSPOTD) = 0.D0

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

      !skip this part with GREENIMP option
      if(opt('GREENIMP')) then
         if(myrank==master) write(*,*) 'Skipping atom loop in main1a'
         i1_start = 1
         i1_end = 1
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
            call set_Jijcalc_flags(t_dtmatJij,NATYPD,NATOMIMPD,NATOMIMP,ATOMIMP,IQAT)
         end if !OPT('XCPL')

         ! nonco angles: defined in mod_wunfiles
         call read_angles(t_params,NATYP,THETA,PHI)

         ! Interpolate potential
         call INTERPOLATE_POTEN(LPOTD,IRM,IRNSD,NATYPD,IPAND,NSPOTD,NTOTD,   &
            NCHEBD,NTOTD*(NCHEBD+1),NSPIN,RMESH,IRMIN,IRWS,IRCUT,VINS,VISP,   &
            NPAN_LOG,NPAN_EQ,NPAN_TOT,RNEW,IPAN_INTERVALL,VINSNEW)

         do I1=i1_start,i1_end

            IPOT=NSPIN*(I1-1)+1

            call TMAT_NEWSOLVER(IELAST,NSPIN,LMAX,ZAT(I1),SOCSCALE(I1),EZ, &
               NSRA,CLEB(1,1),ICLEB,IEND,NCHEB,NPAN_TOT(I1),               &
               RPAN_INTERVALL(0,I1),IPAN_INTERVALL(0,I1),RNEW(1,I1),       &
               VINSNEW,THETA(I1),PHI(I1),I1,IPOT,LLY,DELTAE,IDOLDAU,       &
               LOPT(I1),WLDAU(1,1,1,I1),t_dtmatJij(I1))

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
      !-------------------------------------------------------------------------
      ! End of calculation of Jij-tensor
      !-------------------------------------------------------------------------
#endif


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

      if(t_inc%i_write>0) write (1337,'(79(1H=),/,30X,"< KKR1a finished >",/,79(1H=),/)')

      ! Deallocate leftover arrays
      if (allocated(VINSNEW)) then
         i_all=-product(shape(VINSNEW))*kind(VINSNEW)
         deallocate(VINSNEW,stat=i_stat)
         call memocc(i_stat,i_all,'VINSNEW','main1a')
      endif

   end subroutine main1a

end module MOD_MAIN1A
