!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the calculation of the density for the JM-KKR package
!> Author: Philipp Ruessmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,        
!> and many others ...
!> Wrapper module for the calculation of the density for the JM-KKR package.
!> The code uses the information obtained in the main0 module, this is
!> mostly done via the `get_params_1c()` call, that obtains parameters of the type
!> `t_params` and passes them to local variables
!------------------------------------------------------------------------------------
module mod_main1c

  use :: mod_profiling
  use :: constants
  use :: global_variables
  use :: mod_datatypes

  use :: mod_mvecglobal
  use :: mod_mixldau
  use :: mod_interpolate_poten
  use :: mod_rhoval
  use :: mod_rhoval0
  use :: mod_rhocore
  use :: mod_renorm_lly
  use :: mod_getscratch, only: opendafile
  use :: mod_rhovalnew
  use :: mod_wmatldau
  use :: mod_wmatldausoc
  use :: mod_wrldaupot
  use :: mod_wrldos
  use :: mod_wrmoms

  implicit none

contains

  !------------------------------------------------------------------------------- 
  !> Summary: Main subroutine regarding the calculation of the electronic density
  !> Author: Philipp Rüssmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,
  !> and many others ...
  !> Category: communication, physical-observables, KKRhost 
  !> Deprecated: False 
  !> Main subroutine regarding the calculation of the electronic density
  !------------------------------------------------------------------------------- 
  subroutine main1c()

#ifdef CPP_MPI
    use :: mpi
#endif
#ifdef CPP_TIMING
    use :: mod_timing
#endif
    use :: mod_types, only: t_tgmat, t_inc, t_lloyd
#ifdef CPP_MPI
    use :: mod_types, only: gather_tmat, t_mpi_c_grid, save_t_mpi_c_grid, get_ntot_pt_ioff_pt_2d
#endif
    use :: mod_mympi, only: myrank, master
#ifdef CPP_MPI
    use :: mod_mympi, only: nranks, find_dims_2d, distribute_linear_on_tasks, mympi_main1c_comm, mympi_main1c_comm_newsosol2
#endif
    use :: mod_wunfiles
    use :: mod_version_info

    use :: mod_main0

    ! .. Parameters
    integer :: nmvecmax
    parameter (nmvecmax=4)
    ! ..
    ! .. Local variables
    integer :: nqdos
    integer :: lmaxp1
    integer :: nacls1
    integer :: lmaxd1
    integer :: lrectmt
    integer :: lmmaxd1
    integer :: itmpdir
    integer :: nspinpot

    integer :: i1_start, i1_end, i_stat
    integer :: ie_start, ie_end, ie_num
    integer :: l, i1, ie, ir, is, lm, iq, ipot
    integer :: icell, ipot1, ispin, ihost, iltmp
    real (kind=dp) :: denef
    real (kind=dp) :: chrgsemicore
    character (len=5) :: textns
    character (len=80) :: tmpdir
    character (len=8) :: qdosopt
    logical :: lmomvec
    logical :: ldorhoef
    logical :: itermvdir
    complex (kind=dp) :: csum      ! LLY
    complex (kind=dp) :: eread     ! LLY
    integer, dimension (20, natyp) :: nkcore !! Number of KAPPA values for a given (n,l) core state
    integer, dimension (20, npotd) :: kapcore !! The (maximum 2) values of KAPPA
    real (kind=dp), dimension (natypd) :: eu
    real (kind=dp), dimension (natypd) :: edc
    real (kind=dp), dimension (natypd) :: phi
    real (kind=dp), dimension (natypd) :: theta
    real (kind=dp), dimension (natypd) :: denefat
    real (kind=dp), dimension (nspind) :: charge_lly ! LLY
    real (kind=dp), dimension (0:lmaxd+1, npotd) :: espv
    real (kind=dp), dimension (0:lmaxd+1, 2) :: espv1
    real (kind=dp), dimension (0:lmaxd+1, 2) :: dostot
    real (kind=dp), dimension (krel*20+(1-krel), npotd) :: ecorerel !! for a given (n,l) state the core energies corresponding first/second KAPPA value, AVERAGED over \mu's  These values are written out to the  potential file (routine <RITES>), but the read in (routine <STARTB1>) updates the ECORE array
    real (kind=dp), dimension (2, natypd) :: angles_new
    real (kind=dp), dimension (0:lmaxd+1, natypd, 2) :: charge
    real (kind=dp), dimension (mmaxd, mmaxd, nspind, natypd) :: wldauold
    complex (kind=dp), dimension (iemxd) :: df
    complex (kind=dp), dimension (0:lmaxd+1, ielast, 2) :: den1
    complex (kind=dp), dimension (natypd, 3, nmvecmax) :: mvevi ! OUTPUT
    complex (kind=dp), dimension (natypd, 3, nmvecmax) :: mvevief ! OUTPUT
    complex (kind=dp), dimension (mmaxd, mmaxd, npotd) :: denmatc
    complex (kind=dp), dimension (mmaxd, mmaxd, 2, 2, natypd) :: denmatn
    complex (kind=dp), dimension (0:lmaxd, 3, nmvecmax) :: mvevil1
    complex (kind=dp), dimension (0:lmaxd, 3, nmvecmax) :: mvevil2 ! WORK ARRAYS
    complex (kind=dp), dimension (0:lmaxd, natypd, 3, nmvecmax) :: mvevil
    character (len=7), dimension (3) :: texts
    character (len=4), dimension (0:6) :: textl
    ! -------------------------------------------------------------------------
    !> @note attention: muorb second index means both spins and total
    ! -------------------------------------------------------------------------
    real (kind=dp), dimension (irmd*krel+(1-krel), natypd) :: rhoorb !! orbital density
    real (kind=dp), dimension (0:lmaxd+1+1, 3, natypd) :: muorb !! orbital magnetic moment
    ! ----------------------------------------------------------------------
    ! R2NEF (IRMD,LMPOTD,NATYP,2)  ! rho at FERMI energy
    ! RHO2NS(IRMD,LMPOTD,NATYP,2)  ! radial density
    ! nspin=1            : (*,*,*,1) radial charge density
    ! nspin=2 or krel=1  : (*,*,*,1) rho(2) + rho(1) -> charge
    ! (*,*,*,2) rho(2) - rho(1) -> mag. moment
    ! RHOC(IRMD,NPOTD)              ! core charge density
    ! ----------------------------------------------------------------------
    real (kind=dp), dimension (irmd, npotd) :: rhoc !! core charge density
    real (kind=dp), dimension (irmd, lmpotd, 4) :: rho2m1
    real (kind=dp), dimension (irmd, lmpotd, 4) :: rho2m2
    real (kind=dp), dimension (nrmaxd, lmpotd, nspotd) :: vinsnew
    ! ---------------------------------------------------------------

    ! .. Alocatable arrays
    real (kind=dp), dimension (:, :), allocatable :: qvec
    real (kind=dp), dimension (:, :, :), allocatable :: rho2n1
    real (kind=dp), dimension (:, :, :), allocatable :: rho2n2
    real (kind=dp), dimension (:, :, :, :), allocatable :: r2nef !! rho at FERMI energy
    real (kind=dp), dimension (:, :, :, :), allocatable :: rho2ns !! radial density
    complex (kind=dp), dimension (:, :, :, :), allocatable :: den ! DEN(0:LMAXD1,IELAST,NPOTD,NQDOS)
    complex (kind=dp), dimension (:, :, :, :), allocatable :: denlm ! DENLM(LMMAXD1,IELAST,NPOTD,NQDOS)
    complex (kind=dp), dimension (:), allocatable :: cdos2 ! LLY
    complex (kind=dp), dimension (:), allocatable :: cdos0
    complex (kind=dp), dimension (:), allocatable :: cdos1
    complex (kind=dp), dimension (:), allocatable :: cdosat0
    complex (kind=dp), dimension (:), allocatable :: cdosat1
    complex (kind=dp), dimension (:, :), allocatable :: cdos_lly

    ! -------------------------------------------------------------------------
    ! MPI parameters
    ! -------------------------------------------------------------------------
#ifdef CPP_MPI
    integer :: ierr, i_all
    integer :: ntot1
    integer :: ihelp
    integer :: idim, nranks_local, myrank_ie_tmp
    complex (kind=dp) :: dentot    ! qdos
    integer, dimension (0:nranks-1) :: ntot_pt
    integer, dimension (0:nranks-1) :: ioff_pt
    complex (kind=dp), dimension (:, :), allocatable :: work
    complex (kind=dp), dimension (:, :, :, :), allocatable :: workc
#endif
    ! .. External Functions ..
    logical :: opt, test
    external :: opt, test
    ! ..
    ! .. Data statements ..
    data textl/' s =', ' p =', ' d =', ' f =', ' g =', ' h =', ' i ='/
    data texts/'spin dn', 'spin up', '       '/
    data textns/' ns ='/
    data ldorhoef/.true./
    data ihost/1/                  ! this is the host program

    ! .. Calculate parameters
    lmaxd1 = lmax + 1
    lmmaxd1 = lmmaxd + 1
    lrectmt = wlength*4*lmmaxd*lmmaxd


    ! Allocate arrays
    allocate (rho2n1(irmd,lmpotd,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(rho2n1))*kind(rho2n1), 'RHO2N1', 'main1c')
    allocate (rho2ns(irmd,lmpotd,natypd,2), stat=i_stat)
    call memocc(i_stat, product(shape(rho2ns))*kind(rho2ns), 'RHO2NS', 'main1c')
    allocate (rho2n2(irmd,lmpotd,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(rho2n2))*kind(rho2n2), 'RHO2N2', 'main1c')
    allocate (r2nef(irmd,lmpotd,natypd,2), stat=i_stat)
    call memocc(i_stat, product(shape(r2nef))*kind(r2nef), 'R2NEF', 'main1c')

    ! Initialze to zero
    vins = 0.0d0
    r2nef = 0.0d0
    muorb = 0.0d0
    rho2ns = 0.0d0
    thetas = 0.0d0
    vinsnew = 0.0d0
    thetasnew = 0.0d0
    angles_new = 0.0d0

    ! Consistency check
    if ((krel<0) .or. (krel>1)) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
    if ((krel==1) .and. (nspind==2)) stop ' set NSPIND = 1 for KREL = 1 in the inputcard'
    ! -------------------------------------------------------------------------
    ! This routine previously used to read from unformatted files created by
    ! the main0 module, now  instead of unformatted files take parameters from
    ! types defined in wunfiles.F90
    ! -------------------------------------------------------------------------
    call get_params_1c(t_params,krel,naezd,natypd,ncleb,lm2d,ncheb,ipand,lmpotd,    &
      lmaxd,lmxspd,nfund,npotd,ntotd,mmaxd,iemxd,irmd,nsra,ins,nspin, nacls1, icst, &
      kmrot,iqat,idoldau,irws,ipan,ircut,iend,icleb,loflm,jend,ifunm1,lmsp1,nfu,    &
      llmsp,lcore,ncore,ntcell,irmin,ititle,intervx,intervy,intervz,lly,itmpdir,    &
      iltmp,npan_eq_at,ipan_intervall,npan_log_at,npan_tot,ntldau,lopt,itldau,      &
      ielast,iesemicore,npol,irshift,jwsrel,zrel,itrunldau,qmtet,qmphi,conc,alat,   &
      zat,drdi,rmesh,a,b,cleb,thetas,socscale,rpan_intervall,cscl,rnew,socscl,      &
      thetasnew,efermi,erefldau,ueff,jeff,emin,emax,tk,vins,visp,ecore,drdirel,     &
      r2drdirel,rmrel,vtrel,btrel,wldau,uldau,ez,wez,phildau,tmpdir,solver,nspind,  &
      nspotd, irmind, lmaxd1, ncelld, irid, r_log, naez, natyp, lmax)

    ! Initialization needed due to merging to one executable
    espv(:, :) = 0.d0
    rho2n1(:, :, :) = 0.d0
    rho2n2(:, :, :) = 0.
    ! -------------------------------------------------------------------------
    ! End read in variables
    ! -------------------------------------------------------------------------

    if (idoldau==1) then

      open (67, file='ldau.unformatted', form='unformatted')
      read (67) itrunldau, wldau, uldau, phildau
      close (67)
    end if

    nqdos = 1
    if (opt('qdos    ')) then      ! qdos
      ! Read BZ path for qdos calculation:                             ! qdos
      open (67, file='qvec.dat')   ! qdos
      read (67, *) nqdos           ! qdos
      allocate (qvec(3,nqdos), stat=i_stat) ! qdos
      call memocc(i_stat, product(shape(qvec))*kind(qvec), 'QVEC', 'main1c') ! qdos
      do iq = 1, nqdos             ! qdos
        read (67, *)(qvec(i1,iq), i1=1, 3) ! qdos
      end do                       ! qdos
      close (67)                   ! qdos
    end if

    allocate (den(0:lmaxd1,ielast,nqdos,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(den))*kind(den), 'DEN', 'main1c')
    allocate (denlm(lmmaxd,ielast,nqdos,npotd), stat=i_stat)
    call memocc(i_stat, product(shape(denlm))*kind(denlm), 'DENLM', 'main1c')

    call cinit(ielast*(lmax+2)*npotd*nqdos, den)
    call cinit(ielast*(lmmaxd)*npotd*nqdos, denlm)
    denef = 0.0d0
    call rinit(natyp, denefat)

    itermvdir = opt('ITERMDIR')
    lmomvec = (itermvdir .or. (kmrot/=0))
    nspinpot = krel*2 + (1-krel)*nspin
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! no need to calculate charge correction if no host program, if decimation
    ! or if no energy contour
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ldorhoef = (ihost==1) .and. (.not. opt('DECIMATE')) .and. (npol/=0)
    ! -------------------------------------------------------------------------
    ! LDA+U
    ! -------------------------------------------------------------------------
    if (idoldau==1) call cinit(mmaxd*mmaxd*npotd, denmatc(1,1,1))
    ! -------------------------------------------------------------------------
    ! LDA+U
    ! -------------------------------------------------------------------------
    if (t_tgmat%gmat_to_file) then
      call opendafile(69, 'gmat', 4, lrectmt, tmpdir, itmpdir, iltmp)
    end if

    ! write parameters file that contains passed parameters for further treatment of gflle
    if (opt('lmlm-dos')) then      ! lmlm-dos
      qdosopt = 'n'                ! lmlm-dos
      if (opt('qdos    ')) then    ! lmlm-dos qdos
        qdosopt = 'y'              ! lmlm-dos qdos
      end if                       ! lmlm-dos qdos
      open (67, form='formatted', file='parameters.gflle') ! lmlm-dos
      df(:) = wez(:)/dble(nspin)   ! lmlm-dos
      write (67, *) ielast, iemxd, natyp, nspin, lmax, qdosopt, df(1:ielast), & ! lmlm-dos
        ez(1:ielast), korbit       ! lmlm-dos
      close (67)                   ! lmlm-dos
    end if                         ! OPT('lmlm-dos')                                                   ! lmlm-dos

    ! ------------------------------------------------------------------------
    ! LLY
    ! ------------------------------------------------------------------------
    if (lly/=0) then               ! LLY

      allocate (cdos0(ielast), stat=i_stat)
      call memocc(i_stat, product(shape(cdos0))*kind(cdos0), 'CDOS0', 'main1c')
      allocate (cdos1(ielast), stat=i_stat)
      call memocc(i_stat, product(shape(cdos1))*kind(cdos1), 'CDOS1', 'main1c')
      allocate (cdos2(natypd), stat=i_stat)
      call memocc(i_stat, product(shape(cdos2))*kind(cdos2), 'CDOS2', 'main1c')
      allocate (cdosat0(ielast), stat=i_stat)
      call memocc(i_stat, product(shape(cdosat0))*kind(cdosat0), 'CDOSAT0', 'main1c')
      allocate (cdosat1(ielast), stat=i_stat)
      call memocc(i_stat, product(shape(cdosat1))*kind(cdosat1), 'CDOSAT1', 'main1c')
      allocate (cdos_lly(ielast,nspind), stat=i_stat)
      call memocc(i_stat, product(shape(cdos_lly))*kind(cdos_lly), 'CDOS_LLY', 'main1c')

      ! LLY
      ! Calculate free-space contribution to dos                              ! LLY
      cdos0(1:ielast) = czero      ! LLY
      cdos1(1:ielast) = czero      ! LLY
      cdos2(1:naez) = czero        ! LLY
      do i1 = 1, naez              ! LLY
        cdosat0(1:ielast) = czero  ! LLY
        cdosat1(1:ielast) = czero  ! LLY
        icell = ntcell(i1)         ! LLY
        do ie = 1, ielast          ! LLY
          call rhoval0(ez(ie), drdi(1,i1), rmesh(1,i1), ipan(i1), & ! LLY
            ircut(0,i1), irws(i1), thetas(1,1,icell), cdosat0(ie) & ! LLY
            , cdosat1(ie), irmd, lmax) ! LLY

          ! calculate contribution from free space

          cdos2(i1) = cdos2(i1) + cdosat1(ie)*wez(ie) ! LLY
        end do                     ! LLY
        cdos0(1:ielast) = cdos0(1:ielast) + cdosat0(1:ielast) ! LLY
        cdos1(1:ielast) = cdos1(1:ielast) + cdosat1(1:ielast) ! LLY
      end do                       ! LLY
      cdos0(:) = -cdos0(:)/pi      ! LLY
      cdos1(:) = -cdos1(:)/pi      ! LLY
      ! LLY
      if (myrank==master) then
        open (701, file='freedos.dat', form='FORMATTED') ! LLY
        do ie = 1, ielast          ! LLY
          write (701, fmt='(10E16.8)') ez(ie), cdos0(ie), cdos1(ie) ! LLY
        end do                     ! LLY
        close (701)                ! LLY
        open (701, file='singledos.dat', form='FORMATTED') ! LLY
        do i1 = 1, natyp           ! LLY
          write (701, fmt='(I5,10E16.8)') i1, cdos2(i1) ! LLY
        end do                     ! LLY
        close (701)                ! LLY
      end if                       ! myrank==master

      ! energy integration to get cdos_lly
      cdos_lly(1:ielast, 1:nspin) = czero ! LLY

      if (t_lloyd%cdos_diff_lly_to_file) then
        if (myrank==master) then
          open (701, file='cdosdiff_lly.dat', form='FORMATTED') ! LLY
          do ispin = 1, nspin/(1+korbit) ! LLY
            do ie = 1, ielast      ! LLY
              read (701, fmt='(10e25.16)') eread, cdos_lly(ie, ispin) ! lly
            end do                 ! LLY
          end do                   ! LLY
          close (701)              ! LLY
        end if
#ifdef CPP_MPI
        call mpi_bcast(cdos_lly, ielast*nspin, mpi_double_complex, master, t_mpi_c_grid%mympi_comm_at, ierr)
#endif
      else                         ! (t_lloyd%cdos_diff_lly_to_file)
#ifdef CPP_MPI
        ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
        ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
        ie_start = 0
        ie_end = ielast
#endif
        do ispin = 1, nspin/(1+korbit) ! lly
          do ie_num = 1, ie_end
            ie = ie_start + ie_num
            cdos_lly(ie, ispin) = t_lloyd%cdos(ie_num, ispin)
          end do                   ! ie_num                                                 ! lly
        end do                     ! ispin
#ifdef CPP_MPI
        ! MPI gather cdos_lly on all processors
        ihelp = ielast*nspin       ! IELAST*NSPIN
        allocate (work(ielast,nspin), stat=i_stat)
        call memocc(i_stat, product(shape(work))*kind(work), 'work', 'main1c')
        work = czero
        call mpi_allreduce(cdos_lly(1:ielast,1:nspin), work, ihelp, mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_at, ierr)
        call zcopy(ihelp, work, 1, cdos_lly(1:ielast,1:nspin), 1)

        i_all = -product(shape(work))*kind(work)
        deallocate (work, stat=i_stat)
        call memocc(i_stat, i_all, 'work', 'main1c')
#endif
      end if                       ! (t_lloyd%cdos_diff_lly_to_file)            ! LLY

      ! Add free-space contribution cdos0                                  ! LLY
      do ispin = 1, nspin/(1+korbit) ! LLY
        cdos_lly(1:ielast, ispin) = cdos_lly(1:ielast, ispin) + real(korbit+1, kind=dp)*cdos0(1:ielast) ! LLY
      end do                       ! LLY
      do ispin = 1, nspin/(1+korbit) ! LLY
        csum = czero               ! LLY
        do ie = 1, ielast          ! LLY
          csum = csum + cdos_lly(ie, ispin)*wez(ie) ! LLY
        end do                     ! LLY
        charge_lly(ispin) = -aimag(csum)*pi/nspinpot*(1+korbit) ! LLY
      end do                       ! LLY
      if (myrank==master) then
        open (701, file='cdos_lloyd.dat', form='formatted') ! LLY
        do ispin = 1, nspin/(1+korbit) ! LLY
          do ie = 1, ielast        ! LLY
            write (701, fmt='(10e16.8)') ez(ie), cdos_lly(ie, ispin) ! LLY
          end do                   ! LLY
        end do                     ! LLY
        close (701)                ! LLY
        write (*, *) 'valence charge from lloyds formula:', (charge_lly(ispin), ispin=1, nspin) ! LLY
        if (t_inc%i_write>0) then  ! LLY
          write (1337, *) 'valence charge from lloyds formula:', (charge_lly(ispin), ispin=1, nspin) ! LLY
        end if
      end if                       ! myrank==master

    end if                         ! LLY<>0                                                     ! LLY
    ! -------------------------------------------------------------------------
    ! LLY
    ! -------------------------------------------------------------------------

    if (.not. opt('NEWSOSOL')) then
      ! ----------------------------------------------------------------------
      ! NATYP
      ! ----------------------------------------------------------------------
#ifdef CPP_MPI
      ntot1 = t_inc%natyp
      call distribute_linear_on_tasks(t_mpi_c_grid%nranks_ie,                       &
        t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at,master,ntot1,ntot_pt,ioff_pt, &
        .true.)

      i1_start = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
      i1_end = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)
      t_mpi_c_grid%ntot1 = ntot_pt(t_mpi_c_grid%myrank_ie)

      if (.not. (allocated(t_mpi_c_grid%ntot_pt1) .and. allocated(t_mpi_c_grid%ioff_pt1))) then
        allocate (t_mpi_c_grid%ntot_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
        call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt1))*kind(t_mpi_c_grid%ntot_pt1), 't_mpi_c_grid%ntot_pT1', 'main1c')
        allocate (t_mpi_c_grid%ioff_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
        call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt1))*kind(t_mpi_c_grid%ioff_pt1), 't_mpi_c_grid%ioff_pT1', 'main1c')
      end if
      t_mpi_c_grid%ntot_pt1 = ntot_pt
      t_mpi_c_grid%ioff_pt1 = ioff_pt
#else
      i1_start = 1
      i1_end = natyp
#endif
      do i1 = i1_start, i1_end
        ! -------------------------------------------------------------------
        ! SPIN
        ! -------------------------------------------------------------------
        iq = iqat(i1)
        do ispin = 1, nspin
          icell = ntcell(i1)
          ipot = (i1-1)*nspinpot + ispin
          ipot1 = (i1-1)*nspinpot + 1
#ifdef CPP_TIMING
          call timing_start('main1c - rhoval')
#endif
          call rhoval(ihost,ldorhoef,icst,ins,ielast,nsra,ispin,nspin,nspinpot,i1,  &
            ez,wez,drdi(1,i1),rmesh(1,i1),vins(irmind,1,knosph*ipot+(1-knosph)),    &
            visp(1,ipot),zat(i1),ipan(i1),ircut(0,i1),irmin(i1),thetas(1,1,icell),  &
            ifunm1(1,icell),lmsp1(1,icell),rho2n1(1,1,ipot1),rho2n2(1,1,ipot1),     &
            rhoorb(1,i1),den(0,1,1,ipot),denlm(1,1,1,ipot),muorb(0,1,i1),           &
            espv(0,ipot1),cleb,loflm,icleb,iend,jend,solver,                        &
            socscl(1,krel*i1+(1-krel)),cscl(1,krel*i1+(1-krel)),vtrel(1,i1),        &
            btrel(1,i1),rmrel(1,i1),drdirel(1,i1),r2drdirel(1,i1),zrel(i1),         &
            jwsrel(i1),irshift(i1),lmomvec,qmtet(iq),qmphi(iq),mvevil1,mvevil2,     &
            nmvecmax,idoldau,lopt(i1),phildau(1,i1),wldau(1,1,1,i1),                &
            denmatc(1,1,ipot),natyp,nqdos,lmax)
#ifdef CPP_TIMING
          call timing_pause('main1c - rhoval')
#endif
        end do
        ! -------------------------------------------------------------------
        ! SPIN
        ! -------------------------------------------------------------------
        ipot1 = (i1-1)*nspinpot + 1

        do lm = 1, lmpotd
          do ir = 1, irmd
            rho2ns(ir, lm, i1, 1) = rho2n1(ir, lm, ipot1)
            r2nef(ir, lm, i1, 1) = rho2n2(ir, lm, ipot1)
          end do
        end do

        do l = 0, lmaxd1
          denef = denef - 2.0d0*conc(i1)*aimag(den(l,ielast,1,ipot1))/pi/dble(nspinpot)
          denefat(i1) = denefat(i1) - 2.0d0*aimag(den(l,ielast,1,ipot1))/pi/dble(nspinpot)
        end do

        if (nspinpot==2) then
          do lm = 1, lmpotd
            do ir = 1, irmd
              rho2ns(ir, lm, i1, 2) = rho2n1(ir, lm, ipot1+1)
              r2nef(ir, lm, i1, 2) = rho2n2(ir, lm, ipot1+1)
            end do
          end do

          do l = 0, lmaxd1
            denef = denef - 2.0d0*conc(i1)*aimag(den(l,ielast,1,ipot1+1))/pi/dble(nspinpot)
            denefat(i1) = denefat(i1) - 2.0d0*aimag(den(l,ielast,1,ipot1+1))/pi/dble(nspinpot)
          end do
        end if

        if (test('RHOVALW ')) then ! Bauer
          open (unit=324234, file='out_rhoval')
          write (324234, *) '#IATOM', i1
          write (324234, '(50000F14.7)') rho2ns(:, :, i1, 1)
          if (nspin==2) write (324234, '(50000F14.7)') rho2ns(:, :, i1, 2)
        end if

        ! -----------------------------------------------------------------------
        ! itermdir/kmrot <> 0
        ! -----------------------------------------------------------------------
        if (lmomvec) then
          do is = 1, nmvecmax
            do lm = 1, 3
              mvevi(i1, lm, is) = (0.0d0, 0.0d0)
              mvevief(i1, lm, is) = (0.0d0, 0.0d0)
              do l = 0, lmax
                mvevil(l, i1, lm, is) = mvevil1(l, lm, is)
                mvevi(i1, lm, is) = mvevi(i1, lm, is) + mvevil1(l, lm, is)

                mvevief(i1, lm, is) = mvevief(i1, lm, is) + mvevil2(l, lm, is)
              end do
            end do
          end do
        end if
      end do
      ! ----------------------------------------------------------------------
      ! NATYP
      ! ----------------------------------------------------------------------
#ifdef CPP_TIMING
      ! call timing_stop('main1c - rhoval')
      if (i1_end>=i1_start .and. t_inc%i_time>0) call timing_stop('main1c - rhoval')
#endif
      close (69)                   ! gmat file
#ifndef CPP_MPI
      close (30)                   ! close lmdos file if no mpi is used
      close (31)                   ! close qdos file if no mpi is used, otherwise writeout in the following
#endif
      close (96)                   ! close gflle file

#ifdef CPP_MPI
      ! move writeout of qdos file here                                           ! qdos
      if (opt('qdos    ')) then    ! qdos
        ! first communicate den array to write out qdos files                   ! qdos
        idim = (lmaxd1+1)*ielast*nqdos*npotd ! qdos
        allocate (workc(0:lmaxd1,ielast,nqdos,npotd)) ! qdos
        call memocc(i_stat, product(shape(workc))*kind(workc), 'workc', 'main1c') ! qdos
        workc = (0.d0, 0.d0)       ! qdos
        ! qdos
        call mpi_allreduce(den, workc, idim, mpi_double_complex, mpi_sum, & ! qdos
          mpi_comm_world, ierr)    ! qdos
        call zcopy(idim, workc, 1, den, 1) ! qdos
        i_all = -product(shape(workc))*kind(workc) ! qdos
        deallocate (workc, stat=i_stat) ! qdos
        call memocc(i_stat, i_all, 'workc', 'main1c') ! qdos
        ! qdos
        if (myrank==master) then   ! qdos
          ! qdos
          do i1 = 1, natyp         ! qdos
            ! qdos
            do ispin = 1, nspin    ! qdos
              ! qdos
              if (natyp>=100) then ! qdos
                open (31, &        ! qdos
                  file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10))// & ! qdos
                  char(48+mod(i1,10))//'.'//char(48+ispin)//'.dat') ! qdos
              else                 ! qdos
                open (31, &        ! qdos
                  file='qdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'// & ! qdos
                  char(48+ispin)//'.dat') ! qdos
              end if               ! qdos
              call version_print_header(31) ! qdos
              write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', & ! qdos
                'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! qdos
              ! qdos
              ipot = (i1-1)*nspinpot + ispin ! qdos
              ! qdos
              do ie = 1, ielast    ! qdos
                do iq = 1, nqdos   ! qdos
                  dentot = cmplx(0.d0, 0.d0, kind=dp) ! qdos
                  do l = 0, lmaxd1 ! qdos
                    dentot = dentot + den(l, ie, iq, ipot) ! qdos
                  end do           ! qdos
                  write (31, 100) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), -aimag(dentot)/pi, & ! qdos
                    (-aimag(den(l,ie,iq,ipot))/pi, l=0, lmaxd1) ! qdos
                end do             ! IQ=1,NQDOS                                         ! qdos
              end do               ! IE=1,IELAST                                           ! qdos
              close (31)           ! qdos
              ! qdos
              if (test('compqdos')) then ! complex qdos
                if (natyp>=100) then ! complex qdos
                  open (31, &      ! complex qdos
                    file='cqdos.'//char(48+i1/100)//char(48+mod(i1/10,10))// & ! complex qdos
                    char(48+mod(i1,10))//'.'//char(48+ispin)//'.dat') ! complex qdos
                else               ! complex qdos
                  open (31, &      ! complex qdos
                    file='cqdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'// & ! complex qdos
                    char(48+ispin)//'.dat') ! complex qdos
                end if             ! complex qdos
                call version_print_header(31) ! complex qdos
                write (31, '(A)') '#   lmax, natyp, nspin, nqdos, ielast:' ! complex qdos
                write (31, '(5I9)') lmax, natyp, nspin, nqdos, ielast ! complex qdos
                write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', & ! complex qdos
                  'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! complex qdos
                ! complex qdos
                ipot = (i1-1)*nspinpot + ispin ! complex qdos
                ! complex qdos
                do ie = 1, ielast  ! complex qdos
                  do iq = 1, nqdos ! complex qdos
                    dentot = cmplx(0.d0, 0.d0, kind=dp) ! complex qdos
                    do l = 0, lmaxd1 ! complex qdos
                      den(l, ie, iq, ipot) = -2.0d0/pi*den(l, ie, iq, ipot) ! complex qdos
                      dentot = dentot + den(l, ie, iq, ipot) ! complex qdos
                    end do         ! complex qdos
                    write (31, 110) ez(ie), qvec(1, iq), qvec(2, iq), & ! complex qdos
                      qvec(3, iq), dentot, (den(l,ie,iq,ipot), l=0, lmaxd1) ! complex qdos
                  end do           ! IQ=1,NQDOS                                            ! complex qdos
                end do             ! IE=1,IELAST                                              ! complex qdos
                close (31)         ! complex qdos
              end if               ! qdos
              ! qdos
            end do                 ! ISPIN=1,NSPIN                                             ! qdos
          end do                   ! I1                                                            ! qdos
        end if                     ! myrank_at==master                                              ! qdos
      end if                       ! OPT('qdos    ')                                                    ! qdos
100   format (5f10.6, 40e16.8)     ! qdos
110   format (5f10.6, 80e16.8)     ! complex qdos

      ! reset NQDOS number to avoid loo large communication which is not needed anyways for qdos run
      nqdos = 1

#ifdef CPP_TIMING
      call timing_start('main1c - communication')
#endif
      ! MPI: reduce these arrays, so that master processor has the results (MPI_REDUCE)
#ifdef CPP_TIMING
      call timing_start('main1c - communication')
#endif
      call mympi_main1c_comm(irmd,lmpotd,natyp,lmax,lmaxd1,lmmaxd,npotd,ielast,     &
        mmaxd,idoldau,natyp,krel,lmomvec,nmvecmax,nqdos,rho2ns,r2nef,espv,den,denlm,&
        denmatc,denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,                    &
        t_mpi_c_grid%mympi_comm_ie)
      call mympi_main1c_comm(irmd,lmpotd,natyp,lmax,lmaxd1,lmmaxd,npotd,ielast,     &
        mmaxd,idoldau,natyp,krel,lmomvec,nmvecmax,nqdos,rho2ns,r2nef,espv,den,denlm,&
        denmatc,denef,denefat,rhoorb,muorb,mvevi,mvevil,mvevief,                    &
        t_mpi_c_grid%mympi_comm_at)
#ifdef CPP_TIMING
      call timing_stop('main1c - communication')
#endif

      ! lmdos writeout
      if (myrank==master) then     ! lm-dos
        ! IF (.not.OPT('qdos    ')) THEN                              ! lm-dos
        if (opt('lmdos    ')) then ! lm-dos
          do i1 = 1, natyp         ! lm-dos
            do ispin = 1, nspin    ! lm-dos
              ipot = (i1-1)*nspinpot + ispin ! lm-dos
              if (natyp>=100) then ! lm-dos
                open (30, file='lmdos.'//char(48+i1/100)// & ! lm-dos
                  char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'// & ! lm-dos
                  char(48+ispin)//'.dat') ! lm-dos
              else                 ! lm-dos
                open (30, file='lmdos.'//char(48+i1/10)// & ! lm-dos
                  char(48+mod(i1,10))//'.'//char(48+ispin)//'.dat') ! lm-dos
              end if               ! lm-dos
              call version_print_header(30)
              write (30, *) ' '    ! lm-dos
              write (30, 120) '# ISPIN=', ispin, ' I1=', i1 ! lm-dos
120           format (a8, i3, a4, i5) ! lm-dos
              do ie = 1, ielast    ! lm-dos
                write (30, 130) real(ez(ie), kind=dp), & ! lm-dos
                  (-aimag(denlm(lm,ie,1,ipot))/pi, lm=1, lmmaxd) ! lm-dos
              end do               ! IE                                                 ! lm-dos
            end do                 ! ISPIN                                                 ! lm-dos
130         format (30e12.4)       ! lm-dos
            close (30)             ! lm-dos
          end do                   ! I1                                                        ! lm-dos
        end if                     ! not qdos option                                             ! lm-dos
      end if                       ! myrank==master                                                 ! lm-dos
#endif

    else                           ! new spin-orbit solver

      ! nonco angles
      call read_angles(t_params, natyp, theta, phi)

      ! interpolate potential
      if (idoldau==1) then
        call cinit(mmaxd*mmaxd*4*natyp, denmatn(1,1,1,1,1))
      end if

      call interpolate_poten(lpotd,irmd,irnsd,natyp,ipand,lmpotd,nspotd,ntotd,      &
        ntotd*(ncheb+1),nspin,rmesh,irmin,irws,ircut,vins,visp,npan_log_at,         &
        npan_eq_at,npan_tot,rnew,ipan_intervall,vinsnew)


#ifdef CPP_MPI
      ntot1 = t_inc%natyp

      if (t_mpi_c_grid%dims(1)>1) then
        nranks_local = t_mpi_c_grid%nranks_ie
        if (t_mpi_c_grid%nranks_ie>t_mpi_c_grid%dims(1)) then
          nranks_local = t_mpi_c_grid%dims(1)
        end if
      else
        nranks_local = 1
      end if
      call distribute_linear_on_tasks(nranks_local,t_mpi_c_grid%myrank_ie+t_mpi_c_grid%myrank_at, master, ntot1, ntot_pt, ioff_pt, .true., .true.)
      if (t_mpi_c_grid%nranks_ie<=t_mpi_c_grid%dims(1)) then
        i1_start = ioff_pt(t_mpi_c_grid%myrank_ie) + 1
        i1_end = ioff_pt(t_mpi_c_grid%myrank_ie) + ntot_pt(t_mpi_c_grid%myrank_ie)
        t_mpi_c_grid%ntot1 = ntot_pt(t_mpi_c_grid%myrank_ie)
      else
        myrank_ie_tmp = t_mpi_c_grid%myrank_ie
        if (t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1)) then
          myrank_ie_tmp = myrank_ie_tmp - t_mpi_c_grid%myrank_ie + (t_mpi_c_grid%dims(1)-1)
        end if
        i1_start = ioff_pt(myrank_ie_tmp) + 1
        i1_end = ioff_pt(myrank_ie_tmp) + ntot_pt(myrank_ie_tmp)
        t_mpi_c_grid%ntot1 = ntot_pt(myrank_ie_tmp)
      end if

      if (.not. (allocated(t_mpi_c_grid%ntot_pt1) .and. allocated(t_mpi_c_grid%ioff_pt1))) then
        allocate (t_mpi_c_grid%ntot_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
        call memocc(i_stat, product(shape(t_mpi_c_grid%ntot_pt1))*kind(t_mpi_c_grid%ntot_pt1), 't_mpi_c_grid%ntot_pT1', 'main1c')
        allocate (t_mpi_c_grid%ioff_pt1(0:t_mpi_c_grid%nranks_ie-1), stat=i_stat)
        call memocc(i_stat, product(shape(t_mpi_c_grid%ioff_pt1))*kind(t_mpi_c_grid%ioff_pt1), 't_mpi_c_grid%ioff_pT1', 'main1c')
      end if
      t_mpi_c_grid%ntot_pt1 = ntot_pt
      t_mpi_c_grid%ioff_pt1 = ioff_pt
#else
      i1_start = 1
      i1_end = natyp
#endif

      do i1 = i1_start, i1_end

        icell = ntcell(i1)
        ipot = (i1-1)*nspin + 1

#ifdef CPP_TIMING
        call timing_start('main1c - rhovalnew')
#endif

        call rhovalnew(ldorhoef,ielast,nsra,nspin,lmax,ez,wez,zat(i1),socscale(i1), &
          cleb(1,1),icleb,iend,ifunm1(1,icell),lmsp1(1,icell),ncheb,npan_tot(i1),   &
          npan_log_at(i1),npan_eq_at(i1),rmesh(1,i1),irws(i1),rpan_intervall(0,i1), &
          ipan_intervall(0,i1),rnew(1,i1),vinsnew,thetasnew(1,1,icell),theta(i1),   &
          phi(i1),i1,ipot,den1(0,1,1),espv1(0,1),rho2m1,rho2m2,muorb(0,1,i1),       &
          angles_new(:,i1), idoldau, lopt(i1), wldau(1,1,1,i1),                     & ! LDAU
          denmatn(1,1,1,1,i1), natyp)                                                 ! LDAU

#ifdef CPP_TIMING
        call timing_pause('main1c - rhovalnew')
#endif

        do l = 0, lmaxd1
          espv(l, ipot) = espv1(l, 1)
          espv(l, ipot+nspin-1) = espv1(l, nspin)
          do ie = 1, ielast
            den(l, ie, 1, ipot) = den1(l, ie, 1)
            den(l, ie, 1, ipot+nspin-1) = den1(l, ie, nspin)
          end do
        end do
        do ispin = 1, nspin
          do lm = 1, lmpotd
            do ir = 1, irmd
              rho2ns(ir, lm, i1, ispin) = rho2m1(ir, lm, ispin)
              r2nef(ir, lm, i1, ispin) = rho2m2(ir, lm, ispin)
            end do
          end do
        end do

        do l = 0, lmaxd1
          denef = denef - 2.0d0*conc(i1)*aimag(den(l,ielast,1,ipot))/pi/dble(nspin)
          denefat(i1) = denefat(i1) - 2.0d0*aimag(den(l,ielast,1,ipot))/pi/dble(nspin)
        end do
        do l = 0, lmaxd1
          denef = denef - 2.0d0*conc(i1)*aimag(den(l,ielast,1,ipot+1))/pi/dble(nspin)
          denefat(i1) = denefat(i1) - 2.0d0*aimag(den(l,ielast,1,ipot+1))/pi/dble(nspin)
        end do
        do ispin = 1, nspin
          do l = 0, lmaxd1
            muorb(l, 3, i1) = muorb(l, 3, i1) + muorb(l, ispin, i1)
          end do
        end do
        do ispin = 1, 3
          do l = 0, lmaxd1
            muorb(lmaxd1+1, ispin, i1) = muorb(lmaxd1+1, ispin, i1) + muorb(l, ispin, i1)
          end do
        end do
      end do                       ! I1

#ifdef CPP_MPI
      ! reset NQDOS to 1 to avoid endless communication
      nqdos = 1
      call mympi_main1c_comm_newsosol2(lmaxd1,lmmaxd,ielast,nqdos,npotd,natyp,      &
        lmpotd,irmd,mmaxd,den,denlm,muorb,espv,r2nef,rho2ns,denefat,denef,denmatn,  &
        angles_new,t_mpi_c_grid%mympi_comm_ie)
#endif

#ifdef CPP_TIMING
      ! call timing_stop('main1c - rhovalnew')
      if (i1_end>=i1_start) call timing_stop('main1c - rhovalnew')
#endif

      close (69)

      if (myrank==master) then
        ! rewrite new theta and phi to nonco_angle_out.dat, nonco_angle.dat is the input
        if (.not. test('FIXMOM  ')) then
          open (unit=13, file='nonco_angle_out.dat', form='formatted')
          call version_print_header(13)
          do i1 = 1, natyp
            ! save to file in converted units (degrees)
            write (13, *) angles_new(1, i1)/(2.0d0*pi)*360.0d0, angles_new(2, i1)/(2.0d0*pi)*360.0d0
            ! use internal units here
            t_params%theta(i1) = angles_new(1, i1)
            t_params%phi(i1) = angles_new(2, i1)
          end do
          close (13)

        end if                     ! .not.test('FIXMOM  ')
      end if                       ! (myrank==master)

      close (29)                   ! close lmdos file
      close (30)                   ! close lmdos file
      close (31)                   ! close qdos file
      close (32)                   ! close qdos file
      close (91)                   ! close gflle file

    end if                         ! NEWSOSOL

#ifdef CPP_MPI
    if (myrank==master) then
#endif
#ifdef CPP_TIMING
      call timing_start('main1c - serial part')
#endif

      ! In case of Lloyds formula renormalize valence charge                  ! LLY
      if (lly>0) then              ! LLY
        lmaxp1 = lmax              ! LLY
        if (ins/=0) lmaxp1 = lmax + 1 ! LLY

        call renorm_lly(cdos_lly,ielast,nspin,natyp,den(:,:,1,:),lmaxp1,conc,1,     &
          ielast,wez,ircut,ipan,ez,zat,rho2ns,r2nef,denef,denefat,espv)

      end if                       ! LLY

      ! ----------------------------------------------------------------------
      ! NATYP
      ! ----------------------------------------------------------------------
      chrgsemicore = 0d0
      do i1 = 1, natyp
        !----------------------------------------------------------------------------
        ! l/m_s/atom-resolved charges
        !----------------------------------------------------------------------------

        do ispin = 1, nspinpot
          ipot = (i1-1)*nspinpot + ispin
          do l = 0, lmaxd1
            charge(l, i1, ispin) = 0.0d0

            do ie = 1, ielast
              charge(l, i1, ispin) = charge(l, i1, ispin) + aimag(wez(ie)*den(l,ie,1,ipot))/dble(nspinpot)
              if (ie==iesemicore) then
                chrgsemicore = chrgsemicore + conc(i1)*charge(l, i1, ispin)
              end if
            end do

          end do
        end do
        eu(i1) = 0d0
        edc(i1) = 0d0
        !----------------------------------------------------------------------------
        ! Orbital magnetic moments (array initialised to 0.0D0 in rhoval)
        !----------------------------------------------------------------------------
        if (krel==1) then
          do ispin = 1, 3
            do l = 0, lmax + 1
              muorb(lmaxd1+1, ispin, i1) = muorb(lmaxd1+1, ispin, i1) + muorb(l, ispin, i1)
            end do
          end do
        end if
      end do
      ! ----------------------------------------------------------------------
      ! NATYP
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! LDA+U
      ! ----------------------------------------------------------------------
      if (idoldau==1) then
        ! -------------------------------------------------------------------
        ! Save old LDA+U interaction matrix for mixing
        ! -------------------------------------------------------------------
        call dcopy(mmaxd*mmaxd*nspind*natyp, wldau, 1, wldauold, 1)
        ! -------------------------------------------------------------------
        ! Construct LDA+U interaction matrix for next iteration
        ! -------------------------------------------------------------------
        if (.not. opt('NEWSOSOL')) then
          call wmatldau(ntldau,itldau,nspinpot,denmatc,lopt,ueff,jeff,uldau,wldau,  &
            eu,edc,mmaxd,npotd,natyp,nspin,lmax)
        else
          call wmatldausoc(ntldau,itldau,nspinpot,denmatn,lopt,ueff,jeff,uldau,     &
            wldau,eu,edc,mmaxd,natyp,nspin,lmax)
        end if
        ! -> Mix old and new LDA+U interaction matrices
        call mixldau(mmaxd, nspind, natyp, natyp, nspin, lopt, wldauold, wldau)
        ! -------------------------------------------------------------------
        ! Update variables-file
        ! -------------------------------------------------------------------
        itrunldau = itrunldau + 1
        open (67, file='ldau.unformatted', form='unformatted')
        write (67) itrunldau, wldau, uldau, phildau
        close (67)
        ! -------------------------------------------------------------------
        ! Write full lda+u information in ascii file ldaupot_new
        ! -------------------------------------------------------------------
        call wrldaupot(itrunldau,lopt,ueff,jeff,erefldau,natyp,wldau,uldau,phildau, &
          irmd,natyp,nspind,mmaxd,irws)
      end if
      ! ----------------------------------------------------------------------
      ! LDA+U
      ! ----------------------------------------------------------------------

      call wrmoms(krel+korbit,natyp,nspinpot,texts,textl,textns,charge,muorb,lmax,  &
        lmaxd1)
      ! ----------------------------------------------------------------------
      ! ITERMDIR
      ! ----------------------------------------------------------------------
      if ((krel==1) .and. lmomvec) then
        do i1 = 1, natyp
          iq = iqat(i1)
          call mvecglobal(i1,iq,natyp,qmphi(iq),qmtet(iq),mvevi,mvevil,mvevief,     &
            natyp,lmax,nmvecmax)
        end do
      end if
      ! ----------------------------------------------------------------------
      ! TERMDIR
      ! ----------------------------------------------------------------------
      ! ----------------------------------------------------------------------
      ! TEST BRAHIM
      ! ----------------------------------------------------------------------
      if (npol==0 .or. test('DOS     ')) then
        call wrldos(den,ez,wez,lmaxd1,iemxd,npotd,ititle,efermi,emin,emax,alat,tk,  &
          nacls1,nspinpot,natyp,conc,ielast,intervx,intervy,intervz,dostot)
      end if

      ! ----------------------------------------------------------------------
      ! CORE STATES
      ! ----------------------------------------------------------------------
      !------------------------------------------------------------------------------
      ! RHO_core is calculated only if also RHO_valence was
      !------------------------------------------------------------------------------
      if (npol/=0) then
        if (t_inc%i_write>0) then
          write (1337, *)
          write (1337, '(78("#"))')
          write (1337, '(33X,A)') 'CORE  STATES'
          write (1337, '(78("#"))')
        end if
        do i1 = 1, natyp
          do ispin = 1, nspin
            ipot = (i1-1)*nspinpot + ispin
            ipot1 = (i1-1)*nspinpot + 1

            call rhocore(nsra,ispin,nspin,i1,drdi(1,i1),rmesh(1,i1),visp(1,ipot),   &
              a(i1),b(i1),zat(i1),ircut(0,i1),rhoc(1,ipot1),ecore(1,ipot),          &
              ncore(ipot),lcore(1,ipot),cscl(1,krel*i1+(1-krel)),vtrel(1,i1),       &
              btrel(1,i1),rmrel(1,i1),drdirel(1,i1),r2drdirel(1,i1),zrel(i1),       &
              jwsrel(i1),irshift(i1),ecorerel(1,ipot1),nkcore(1,i1),kapcore(1,ipot1))

          end do
        end do
        if (t_inc%i_write>0) then
          write (1337, *)
          write (1337, '(78("#"))')
          write (1337, *)
        end if
      end if
      ! ----------------------------------------------------------------------
      ! CORE STATES
      ! ----------------------------------------------------------------------
      call save_density(t_params,rho2ns,r2nef,rhoc,denef,denefat,espv,ecore,idoldau,&
        lopt,eu,edc,chrgsemicore,rhoorb,ecorerel,nkcore,kapcore,krel,natyp,npotd,   &
        irmd,lmpotd,lmaxd1)

      if (test('den-asci')) then
        open (67, file='densitydn.ascii', form='formatted')
        do i1 = 1, natyp
          do lm = 1, lmpotd
            do ir = 1, irmd
              write (67, fmt='(I6,2I5,2E25.16)') i1, lm, ir, rho2ns(ir, lm, i1, 1), rho2ns(ir, lm, i1, 2)
            end do
          end do
        end do
        close (67)
      end if

      ! ----------------------------------------------------------------------
      ! TERMDIR
      ! ----------------------------------------------------------------------
      if (itermvdir) then
        t_params%mvevi = mvevi
        t_params%mvevief = mvevief
      end if
      ! ----------------------------------------------------------------------
      ! TERMDIR
      ! ----------------------------------------------------------------------


      if (t_inc%i_write>0) then
        write (1337, '(79("="),/,30X,"< KKR1c finished >",/,79("="),/)')
      end if

#ifdef CPP_TIMING
      call timing_stop('main1c - serial part')
#endif
#ifdef CPP_MPI
    end if                         ! myrank==master
#endif


    if (lly/=0) then
      ! cleanup allocations
      deallocate (cdos0, stat=i_stat)
      call memocc(i_stat, -product(shape(cdos0))*kind(cdos0), 'CDOS0', 'main1c')
      deallocate (cdos1, stat=i_stat)
      call memocc(i_stat, -product(shape(cdos1))*kind(cdos1), 'CDOS1', 'main1c')
      deallocate (cdos2, stat=i_stat)
      call memocc(i_stat, -product(shape(cdos2))*kind(cdos2), 'CDOS2', 'main1c')
      deallocate (cdosat0, stat=i_stat)
      call memocc(i_stat, -product(shape(cdosat0))*kind(cdosat0), 'CDOSAT0', 'main1c')
      deallocate (cdosat1, stat=i_stat)
      call memocc(i_stat, -product(shape(cdosat1))*kind(cdosat1), 'CDOSAT1', 'main1c')
      deallocate (cdos_lly, stat=i_stat)
      call memocc(i_stat, -product(shape(cdos_lly))*kind(cdos_lly), 'CDOS_LLY', 'main1c')
    end if

  end subroutine main1c

end module mod_main1c
