!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Wrapper module for the calculation of the structural Greens function `gmat`
!> Author: Philipp Ruessmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,       
!> and many others ... 
!> Wrapper module for the calculation of the structural Greens function `gmat` 
!------------------------------------------------------------------------------------
module mod_main1b

  private
  public :: main1b

contains

  !-------------------------------------------------------------------------------  
  !> Summary: Main subroutine regarding the claculation of the structural Greens function `gmat`
  !> Author: Philipp Ruessmann, Bernd Zimmermann, Phivos Mavropoulos, R. Zeller,     
  !> and many others ... 
  !> Category: structural-greensfunction, k-points, reference-system, KKRhost 
  !> Deprecated: False 
  !> Main subroutine regarding the claculation of the structural Green's function `gmat`
  !-------------------------------------------------------------------------------  
  subroutine main1b() 

#ifdef CPP_MPI
    use :: mpi
    use :: mod_mympi, only: find_dims_2d, distribute_linear_on_tasks, mpiadapt
    use :: mod_types, only: t_mpi_c_grid, save_t_mpi_c_grid, get_ntot_pt_ioff_pt_2d, init_params_t_imp, init_t_imp, bcast_t_imp_scalars, &
      bcast_t_imp_arrays
#endif
    use :: mod_mympi, only: myrank, master
    use :: mod_datatypes, only: dp
    use :: mod_runoptions, only: calc_exchange_couplings, formatted_files, set_gmat_to_zero, use_Chebychev_solver, &
      use_qdos, use_readcpa, write_deci_tmat, write_gmat_plain, write_green_host, write_green_imp, write_kkrimp_input, &
      write_pkkr_input, write_pkkr_operators, write_rhoq_input
    use :: mod_constants, only: czero, cone, pi, nsymaxd
    use :: mod_profiling, only: memocc
    use :: mod_operators_for_fscode, only: operators_for_fscode
    use :: mod_getscratch, only: opendafile
    use :: mod_kloopz1, only: kloopz1_qdos
    use :: mod_greenimp, only: greenimp
    use :: mod_changerep, only: changerep
    use :: mod_tmatimp_newsolver, only: tmatimp_newsolver
    use :: mod_setfactl, only: setfactl
    use :: mod_calctref13, only: calctref13
    use :: mod_rotatespinframe, only: rotatematrix
    use :: mod_types, only: t_tgmat, t_inc, t_lloyd, t_cpa, init_t_cpa, t_imp
    use :: mod_timing, only: timing_start, timing_pause, timing_stop, timings_1b, print_time_and_date
    use :: mod_wunfiles, only: get_params_1b, t_params, read_angles
    use :: mod_tbxccpljij, only: tbxccpljij
    use :: mod_tbxccpljijdij, only: tbxccpljijdij
    use :: mod_rhoqtools, only: rhoq_save_refpot
    use :: mod_cinit, only: cinit
    ! array dimensions
    use :: global_variables, only: maxmshd, iemxd, natypd, naezd, kpoibz, lmmaxd, lmgf0d, lmaxd, nrefd, nsheld, wlength, nofgij, &
      naclsd, nspind, nclsd, nembd, krel, korbit, natomimpd, nrd, nembd1, nspindd, nprincd, lmmaxso, irmind, nspotd, irmd, lpotd, &
      ncleb, ipand, irnsd, lmpotd, irid, nfund, ntotd
    ! stuff defined in main0 already
    use :: mod_main0, only: natyp, ielast, npol, nref, naez, nsra, ins, nspin, ncls, lly, atom, cls, nacls, refpot, ez, alat, rmtref, &
      vref, atomimp, icc, igf, nlbasis, nrbasis, ncpa, icpa, itcpamax, cpatol, rbasis, rr, ezoa, nshell, kmrot, kaoez, ish, jsh, nsh1, &
      nsh2, noq, iqat, natomimp, conc, kmesh, maxmesh, nsymat, nqcalc, ratom, rrot, drotq, ijtabcalc, ijtabcalc_i, ijtabsym, ijtabsh, &
      iqcalc, dsymll, invmod, icheck, symunitary, rc, crel, rrel, srrel, nrrel, irrel, lefttinvll, righttinvll, wez, rclsimp, vacflag, &
      iend, lmax, r_log, vins, visp, ipan, irmin, icleb, zat, rmesh, cleb, ncheb, ircut, rcls

    implicit none

    ! .. Local variables
    integer :: nspin1
    integer :: l
    integer :: i
    integer :: i1
    integer :: ie
    integer :: iq
    integer :: ix
    integer :: ic
    integer :: l1
    integer :: ilm
    integer :: lm1
    integer :: lm2
    integer :: irec
    integer :: iltmp
    integer :: nmesh
    integer :: nqdos               !! number of qdos points
    integer :: isite               ! qdos ruess
    integer :: ideci
    integer :: ispin
    integer :: iprint
    integer :: itmpdir
    integer :: lrecgrf
    integer :: lrectmt
    integer :: lrectra             ! LLY Lloyd
    integer :: iqdosrun            !! counter to organise qdos run
    integer :: naclsmax
    integer :: lrecgrf1
    integer :: ncpafail
    integer :: icpaflag
    integer :: reclength
    integer :: lrecgreen
    real (kind=dp) :: phi
    real (kind=dp) :: theta
    real (kind=dp) :: rfctor       !! rfctor=a/(2*pi) conversion factor to p.u.
    complex (kind=dp) :: eryd
    complex (kind=dp) :: tread     ! qdos ruess
    complex (kind=dp) :: cfctor
    complex (kind=dp) :: tralpha1  ! LLY Lloyd
    complex (kind=dp) :: cfctorinv
    logical :: opt
    logical :: test
    logical :: lcpaij

    character (len=80) :: tmpdir
    character (len=80) :: text                             ! qdos ruess

    ! .. Local arrays
    integer, dimension (maxmshd) :: nofks
    integer, dimension (iemxd) :: iecpafail
    integer, dimension (natypd, naezd) :: itoq
    real (kind=dp), dimension (natypd) :: phi_at
    real (kind=dp), dimension (maxmshd) :: volbz
    real (kind=dp), dimension (natypd) :: theta_at
    real (kind=dp), dimension (kpoibz, maxmshd) :: volcub
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: w1
    complex (kind=dp), dimension (lmgf0d, lmgf0d) :: wn1
    complex (kind=dp), dimension (lmgf0d, lmgf0d) :: wn2 ! LLY
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: tmat
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gmat0
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: factl
    complex (kind=dp), dimension (0:lmaxd, nrefd) :: alpharef
    complex (kind=dp), dimension (0:lmaxd, nrefd) :: dalpharef !! LLY Lloyd Alpha matrix and deriv.
    complex (kind=dp), dimension (lmmaxd, lmmaxd, natypd) :: msst
    complex (kind=dp), dimension (lmmaxd, lmmaxd, natypd) :: tsst
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naezd) :: tqdos ! qdos ruess
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nrefd) :: trefll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsheld) :: gmatll !! GMATLL = diagonal elements of the G matrix (system)
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nrefd) :: dtrefll !! LLY Lloyd dtref/dE
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naezd) :: dtmatll !! LLY Lloyd  dt/dE
    complex (kind=dp), dimension (lmmaxd*lmmaxd) :: gimp !!  Cluster GF (ref. syst.)
    character (len=35), dimension (0:3), parameter :: invalg = [ &
      'FULL MATRIX                        ', &
      'BANDED MATRIX (slab)               ', &
      'BANDED + CORNERS MATRIX (supercell)', &
      'GODFRIN: nonuniform block partition' ]

    ! .. Allocatable local arrays
    real (kind=dp), dimension (:, :), allocatable :: qvec !! qdos ruess, q-vectors for qdos
    real (kind=dp), dimension (:, :, :), allocatable :: bzkp
    complex (kind=dp), dimension (:, :), allocatable :: dtmtrx !! For GREENIMP
    complex (kind=dp), dimension (:, :, :), allocatable :: ginp !! Cluster GF (ref syst.) GINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
    complex (kind=dp), dimension (:, :, :), allocatable :: dginp !! LLY Lloyd Energy derivative of GINP DGINP(NACLSD*LMGF0D,LMGF0D,NCLSD)
    complex (kind=dp), dimension (:), allocatable :: lly_g0tr !! LLY Lloyd  Trace[ X ], Eq.5.27 PhD Thiess
    complex (kind=dp), dimension (:), allocatable :: tralpharef ! LLY Lloyd
    complex (kind=dp), dimension (:), allocatable :: cdosref_lly ! LLY Lloyd
    complex (kind=dp), dimension (:, :), allocatable :: tracet !! Tr[ (t-tref)^-1 d(t-tref)/dE ]  ! LLY Lloyd
    complex (kind=dp), dimension (:, :), allocatable :: tralpha
    complex (kind=dp), dimension (:, :), allocatable :: lly_grtr !! LLY Lloyd  Trace[ M^-1 dM/dE ], Eq.5.38 PhD Thiess
    complex (kind=dp), dimension (:, :), allocatable :: cdos_lly

#ifdef CPP_MPI
    integer :: ihelp
    complex (kind=dp), allocatable :: work(:, :)
#endif
    integer :: ie_start, ntot2
    integer :: ie_num, ie_end, ierr, i_stat, i_all

    ! for OPERATOR option
    logical :: lexist, operator_imp
    ! -------------------------------------------------------------------------
    ! for conductivity calculation
    ! INTEGER NCPAIRD
    ! PARAMETER(NCPAIRD=10)
    ! INTEGER IATCONDL(NCPAIRD),IATCONDR(NCPAIRD),NCONDPAIR
    ! -------------------------------------------------------------------------
    ! .. Intrinsic Functions ..
    intrinsic :: atan


    ! .. Set the parameters
    lrectra = wlength*4            ! LLY Lloyd
    lrecgrf = wlength*4*naclsd*lmgf0d*lmgf0d*nclsd ! 4 words = 16 bytes / complex number (in ifort 4; in gfort 16) word/byte distiction moved to subroutine opendafile to be the same for all unformatted files
    lrectmt = wlength*4*lmmaxd*lmmaxd
    lrecgreen = wlength*2*natomimp*lmmaxd*natomimp*lmmaxd
    ! ..
    ! .. Data statements
    ! ..
    iprint = 0
    if (t_inc%i_write>0) iprint = 1

    ! allocatable arrays
    allocate (bzkp(3,kpoibz,maxmshd), stat=i_stat)
    call memocc(i_stat, product(shape(bzkp))*kind(bzkp), 'BZKP', 'main1b')
    allocate (lly_g0tr(ielast), stat=i_stat)
    call memocc(i_stat, product(shape(lly_g0tr))*kind(lly_g0tr), 'LLY_G0TR', 'main1b')
    allocate (tralpharef(ielast), stat=i_stat)
    call memocc(i_stat, product(shape(tralpharef))*kind(tralpharef), 'TRALPHAREF', 'main1b')
    allocate (cdosref_lly(ielast), stat=i_stat)
    call memocc(i_stat, product(shape(cdosref_lly))*kind(cdosref_lly), 'CDOSREF_LLY', 'main1b')
    allocate (tracet(ielast,nspind), stat=i_stat)
    call memocc(i_stat, product(shape(tracet))*kind(tracet), 'TRACET', 'main1b')
    allocate (tralpha(ielast,nspind), stat=i_stat)
    call memocc(i_stat, product(shape(tralpha))*kind(tralpha), 'TRALPHA', 'main1b')
    allocate (lly_grtr(ielast,nspind), stat=i_stat)
    call memocc(i_stat, product(shape(lly_grtr))*kind(lly_grtr), 'LLY_GRTR', 'main1b')
    allocate (cdos_lly(ielast,nspind), stat=i_stat)
    call memocc(i_stat, product(shape(cdos_lly))*kind(cdos_lly), 'CDOS_LLY', 'main1b')

    ! Consistency check
    if ((krel<0) .or. (krel>1)) stop ' set KREL=0/1 (non/fully) relativistic mode in the inputcard'
    if ((krel==1) .and. (nspind==2)) stop ' set NSPIND = 1 for KREL = 1 in the inputcard'

    ! -------------------------------------------------------------------------
    ! This routine previously used to read from unformatted files created by
    ! the main0 module, now  instead of unformatted files take parameters from
    ! types defined in wunfiles.F90
    ! -------------------------------------------------------------------------
    call get_params_1b(t_params,natypd,naezd,natyp,naclsd,ielast,npol,nclsd,nrefd,  &
      nref,nembd,naez,nsra,ins,nspin,lmaxd,ncls,lly,krel,atom,cls,nacls,refpot,     &
      ez, itmpdir, iltmp, alat, rcls, iemxd, rmtref, vref, tmpdir, nsheld, nprincd, &
      kpoibz,atomimp,natomimpd,icc,igf,nlbasis,nrbasis,ncpa,icpa,itcpamax,cpatol,   &
      nrd,ideci,rbasis,rr,ezoa,nshell,kmrot,kaoez,ish,jsh,nsh1,nsh2,noq,iqat,       &
      nofgij,natomimp,conc,kmesh,maxmesh,nsymat,nqcalc,ratom,rrot,drotq,ijtabcalc,  &
      ijtabcalc_i,ijtabsym,ijtabsh,iqcalc,dsymll,invmod,icheck,symunitary,rc,crel,  &
      rrel,srrel,nrrel,irrel,lefttinvll,righttinvll,vacflag,nofks,volbz,bzkp,volcub,&
      wez, nembd1, lmmaxd, nsymaxd, nspindd, maxmshd, rclsimp)

    if (write_rhoq_input) then
      open (9889, access='direct', file='tau0_k', form='unformatted', recl=(lmmaxd*lmmaxd+1)*4) ! lm blocks
    end if


    if (test('gmatasci')) open (298347, file='gmat.ascii', form='formatted')
    ! -------------------------------------------------------------------------
    ! End of reading the variables
    ! -------------------------------------------------------------------------

    ! -------------------------------------------------------------------------       !fswrt
    ! open file to store the output for the (external) Fermi-surface program          !fswrt
    ! this file is already partly filled with data by main0. More data                !fswrt
    ! will be stored in                                                               !fswrt
    ! -------------------------------------------------------------------------       !fswrt
    if (write_pkkr_input .and. myrank==master) then                                    !fswrt
      open (6801, file='TBkkr_container.txt', form='formatted', position='append')    !fswrt
    end if                                                                            !fswrt
    ! -------------------------------------------------------------------------       !fswrt
    ! open file for WRTGREEN option (writes green_host file for                       !fswrt
    ! GMATLL_GES creation in zulapi part) file is filled in ROTGLL called in kloopz   !fswrt
    ! -------------------------------------------------------------------------       !fswrt
    if (write_green_host .and. myrank==master) then
      open (58, file='green_host', form='formatted')
    end if
    !--------------------------------------------------------------------------------
    ! If qdos option is used set IQDOSRUN so that in a first run the
    ! (t(E)-t_ref(E))^-1 matrix (tmat.qdos) and the gref matrix can be
    ! written out for one k point, in a second run these matrices are
    ! read in to continue the calculation with the k points specified by
    ! the user in the qvec.dat file
    !--------------------------------------------------------------------------------
    if (use_qdos) then                                                         ! qdos ruess
      iqdosrun = 0                                                                    ! qdos ruess
    else                                                                              ! qdos ruess
      iqdosrun = -1                                                                   ! qdos ruess
    end if                                                                            ! qdos ruess
    ! Jump back here to continue with second run if qdos option is selected           ! qdos ruess
100 continue                                                                          ! qdos ruess
    ! Reset GMATLL for calculation in second run                                      ! qdos ruess
    if (iqdosrun==1) then                                                             ! qdos ruess
      do i1 = 1, nshell(0)                                                            ! qdos ruess
        gmatll(1:lmmaxd, 1:lmmaxd, i1) = czero                                        ! qdos ruess
      end do                                                                          ! qdos ruess
    end if                                                                            ! qdos ruess

    if ((use_qdos) .and. (write_deci_tmat)) then
      stop 'ERROR: qdos and deci-out cannot be used simultaniously'
    else if (use_qdos) then
      if (.not. formatted_files) then
        ! wlength needs to take double complex values
        open (37, access='direct', recl=wlength*16, file='tmat.qdos', form='unformatted')
      else
        open (37, file='tmat.qdos', form='formatted') ! qdos ruess
      end if
    else if (write_deci_tmat) then
      open (37, file='decifile', form='formatted', position='append') ! ruess: needed in case of deci-out option to prepare decifile
    end if

    do i = 1, naez
      do l = 1, noq(i)
        itoq(l, i) = kaoez(l, i)
      end do
    end do
    rfctor = alat/(2*pi)           ! = ALAT/(2*PI)
    cfctor = cone*rfctor
    cfctorinv = cone/rfctor

    call setfactl(factl, lmax, krel, lmmaxd)

    if (t_inc%i_write>0) then
      write (1337, '(79("="))')
      write (1337, '(2A)') '      Inversion algorithm used : ', invalg(invmod)
      write (1337, '(79("="))')
    end if

    naclsmax = 1
    do ic = 1, ncls
      if (nacls(ic)>naclsmax) naclsmax = nacls(ic)
    end do
    lrecgrf1 = wlength*4*naclsmax*lmgf0d*lmgf0d*ncls

    if (.not. allocated(ginp)) then
      allocate (ginp(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
      call memocc(i_stat, product(shape(ginp))*kind(ginp), 'GINP', 'main1b')
    end if
    if (.not. allocated(dginp)) then
      allocate (dginp(naclsmax*lmgf0d,lmgf0d,ncls), stat=i_stat)
      call memocc(i_stat, product(shape(dginp))*kind(dginp), 'DGINP', 'main1b')
    end if

    if (t_tgmat%gref_to_file) then
      call opendafile(68, 'gref', 4, lrecgrf1, tmpdir, itmpdir, iltmp)
    end if
    if (t_tgmat%tmat_to_file) then
      call opendafile(69, 'tmat', 4, lrectmt, tmpdir, itmpdir, iltmp)
    end if
    if (t_tgmat%gmat_to_file) then
      call opendafile(70, 'gmat', 4, lrectmt, tmpdir, itmpdir, iltmp)
    end if
    if (lly/=0) then               ! LLY Lloyd
      if (t_lloyd%dgref_to_file) then
        call opendafile(681, 'dgrefde', 7, lrecgrf1, tmpdir, itmpdir, iltmp) ! LLY Lloyd: derivative of Gref
      end if
      if (t_lloyd%g0tr_to_file) then
        open (682, file='lly_g0tr_ie.ascii', form='FORMATTED') ! LLY Lloyd: trace eq.5.27 PhD Thiess
      end if
      if (t_lloyd%dtmat_to_file) then
        call opendafile(691, 'dtmatde', 7, lrectmt, tmpdir, itmpdir, iltmp) ! LLY Lloyd: derivative of t-matrix
      end if
      if (t_lloyd%tralpha_to_file) then
        call opendafile(692, 'tralpha', 7, lrectra, tmpdir, itmpdir, iltmp) ! LLY Lloyd: Tr[alpha^{-1} dalpha/dE]
      end if
    end if                         ! LLY Lloyd

    lcpaij = .false.
    if ((ncpa/=0) .and. (nshell(0)>natyp)) lcpaij = .true.

    if (lcpaij) then
      if (t_cpa%dmatproj_to_file) then
        open (71, access='direct', recl=2*lrectmt, file='dmatproj.unformatted', form='unformatted')
      else
#ifdef CPP_MPI
        ntot2 = t_mpi_c_grid%ntot2
#else
        ntot2 = ielast
#endif
        call init_t_cpa(t_inc, t_cpa, ntot2)
      end if                       ! t_cpa%dmatproj_to_file
    end if                         ! LCPAIJ

    if (igf/=0) then

      if ((write_gmat_plain)) then
        open (8888, file='kkrflex_green.dat')
      end if

      if ((write_kkrimp_input)) then
        ! ! Green functions has (lmmaxd*natomimp)**2 double complex (i.e. factor '4') values
        ! RECLENGTH = WLENGTH*4*NATOMIMP*LMMAXD*NATOMIMP*LMMAXD
        ! at the moment kkrflex_green file is only written with single precision (factor'2')
        reclength = wlength*2*natomimp*lmmaxd*natomimp*lmmaxd
        ! sometimes (lmax=2) the record length might be too small to store the parameters, then reclength needs to be bigger
        if (reclength<8*ielast+6) then
          stop '[main1b] record length for kkrflex_green is too small to store parameters, use either more atoms in the cluster, a higher lmax or less energy points'
        end if
        open (888, access='direct', recl=reclength, file='kkrflex_green', form='unformatted')
      end if

      !------------------------------------------------------------------------------
      ! the following write-out has been disabled, because it was assumed to be       !no-green
      ! obsolete with the implementation of the MPI-communicated arrays. If I am      !no-green
      ! wrong and the write-out is needed in subsequent parts, construct a            !no-green
      ! test-option around it so that it is only written out in this case.            !no-green
      ! OPEN (88,ACCESS='direct',RECL=LRECGREEN,                                      !no-green
      ! &             FILE='green',FORM='unformatted')                                !no-green
      !------------------------------------------------------------------------------
      irec = 1

      if ((write_kkrimp_input)) then
        if (myrank==master) then
          write (888, rec=irec) ielast, nspin, natomimp, natomimp, lmmaxd, korbit,  &
            (ez(ie), ie=1, ielast), (wez(ie), ie=1, ielast)
          if ((write_gmat_plain)) then
            ! WRITE(8888,'(I5,50000F)') IELAST,NSPIN,NATOMIMP,NATOMIMP,&
            write (8888, *) ielast, nspin, natomimp, natomimp, (lmax+1)**2, (ez(ie), &
              ie=1, ielast), (wez(ie), ie=1, ielast)
          end if
        end if
#ifdef CPP_MPI
        call mpi_barrier(mpi_comm_world, ierr)
#endif
      end if
      ! IF ( (.not. write_kkrimp_input ) ) THEN                     !no-green
      ! WRITE(88,REC=IREC) IELAST,NSPIN,                         !no-green
      ! &         (EZ(IE),IE=1,IELAST),(WEZ(IE),IE=1,IELAST),    !no-green
      ! &         NATOMIMPD*LMMAXD                               !no-green
      ! END IF                                                   !no-green
    end if

    ! Value of NQDOS changes to a read-in value if option qdos is applied, otherwise:
    nqdos = 1                                                                       ! qdos ruess
    if (use_qdos .and. (iqdosrun==1)) then                                   ! qdos ruess
      ! Read BZ path for qdos calculation:                                          ! qdos ruess
      open (67, file='qvec.dat')                                                    ! qdos ruess
      read (67, *) nqdos                                                            ! qdos ruess
      i_all = -product(shape(qvec))*kind(qvec)                                      ! qdos ruess
      ! deallocate in first run allocated array to change it                        ! qdos ruess
      deallocate (qvec, stat=i_stat)                                                ! qdos ruess
      call memocc(i_stat, i_all, 'QVEC', 'main1b')                                  ! qdos ruess
                                                                                    ! qdos ruess  
      allocate (qvec(3,nqdos), stat=i_stat)                                         ! qdos ruess
      call memocc(i_stat, product(shape(qvec))*kind(qvec), 'QVEC', 'main1b')        ! qdos ruess
      do iq = 1, nqdos                                                              ! qdos ruess
        read (67, *)(qvec(ix,iq), ix=1, 3)                                          ! qdos ruess
      end do                                                                        ! qdos ruess
      close (67)                                                                    ! qdos ruess
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ! qdos ruess
      ! Prepare k-mesh information to be appropriate for qdos calculation.          ! qdos ruess
      ! The idea is that subr. KLOOPZ1 is called for only one point at a time,      ! qdos ruess
      ! with weight equal to the full BZ; in this way we avoid changing the         ! qdos ruess
      ! calling list or the contents of kloopz1.                                    ! qdos ruess
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      ! qdos ruess
      kmesh(1:ielast) = 1                                                           ! qdos ruess
      nofks(1) = 1                                                                  ! qdos ruess
      volcub(1, 1) = volbz(1)                                                       ! qdos ruess
      nsymat = 1                                                                    ! qdos ruess
    else if (use_qdos .and. (iqdosrun==0)) then                              ! qdos ruess
      ! Call the k loop just once with one k point to write out the tmat.qdos file  ! qdos ruess
      allocate (qvec(3,nqdos), stat=i_stat)                                         ! qdos ruess
      call memocc(i_stat, product(shape(qvec))*kind(qvec), 'QVEC', 'main1b')        ! qdos ruess
      if (i_stat/=0) stop '[main1b] Error allocating qvec'                          ! qdos ruess
      qvec(1:3, 1) = 0.d0                                                           ! qdos ruess
      kmesh(1:ielast) = 1                                                           ! qdos ruess
      nofks(1) = 1                                                                  ! qdos ruess
      volcub(1, 1) = volbz(1)                                                       ! qdos ruess
    end if                                                                          ! qdos ruess

    ncpafail = 0

    ! Initialize trace for Lloyd formula
    lly_grtr(:, :) = czero ! 1:IELAST,1:NSPIND

    ! determine extend of spin loop
    nspin1 = nspin/(1+korbit) ! factor (1+korbit) takes care of NOSOC option
    if (use_Chebychev_solver) then
      ! nonco angles
      call read_angles(t_params, natyp, theta_at, phi_at)
    end if

#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    ie_start = 0
    ie_end = ielast
#endif
    if (write_rhoq_input) then
      ie_start = 1
      ie_end = 1
    end if

    ! ----------------------------------------------------------------------
    ! BEGIN do loop over spins and energies
    ! ----------------------------------------------------------------------
    do ispin = 1, nspin1

      do ie_num = 1, ie_end
        ie = ie_start + ie_num

#ifdef CPP_MPI
        ! start timing measurement for this ie, needed for MPIadapt
        if (mpiadapt>0 .and. t_mpi_c_grid%myrank_ie==0) then
          call timing_start('time_1b_ie')
        end if
#endif

        ! write energy into green_host file
        if (write_green_host .and. myrank==master) then
          write (58, '(2e17.9)') ez(ie)
        end if

        ! read in Green function of reference system
        if (t_tgmat%gref_to_file) then
          read (68, rec=ie) ginp
        else
          ginp(:, :, :) = t_tgmat%gref(:, :, :, ie_num)
        end if
        if (t_lloyd%dgref_to_file) then
          if (lly/=0) read (681, rec=ie) dginp ! LLY Lloyd
        else
          if (lly/=0) dginp(:, :, :) = t_lloyd%dgref(:, :, :, ie_num)
        end if

        eryd = ez(ie)
        nmesh = kmesh(ie)
        if (t_inc%i_write>0) then
          write (1337, '(A,I3,A,2(1X,F10.6),A,I3)') ' ************ IE = ', ie, ' ENERGY =', ez(ie), ' KMESH = ', nmesh
        end if

        ! ----------------------------------------------------------------
        ! I1 = 1,NREF
        ! calculate t(ll') of the reference system (on clusters)
        ! ----------------------------------------------------------------
#ifdef CPP_TIMING
        call timing_start('main1b - calctref13')
#endif
        trefll(:, :, :) = czero
        dtrefll(:, :, :) = czero
        if (krel==0) then
          do i1 = 1, nref
            if (.not. use_Chebychev_solver) then
              call calctref13(eryd, vref(i1), rmtref(i1), lmax, lm1, wn1, wn2, & ! LLY Lloyd
                alpharef(0,i1), dalpharef(0,i1), lmax+1, lmmaxd)                 ! LLY Lloyd
            else
              call calctref13(eryd, vref(i1), rmtref(i1), lmax, lm1, wn1, wn2, & ! LLY
                alpharef(0,i1), dalpharef(0,i1), lmax+1, lmgf0d)                 ! LLY
            end if
            do i = 1, lm1
              trefll(i, i, i1) = wn1(i, i)
              if (use_Chebychev_solver .and. .not.test('NOSOC   ')) trefll(lm1+i, lm1+i, i1) = wn1(i, i)
              dtrefll(i, i, i1) = wn2(i, i)                              ! LLY
              if (use_Chebychev_solver .and. .not.test('NOSOC   ')) dtrefll(lm1+i, lm1+i, i1) = wn2(i, i) ! LLY
            end do

            if (write_rhoq_input) then
              call rhoq_save_refpot(ielast,i1,nref,natyp,refpot(1:natyp),wlength,   &
                lmmaxd,ie,trefll)
            end if                 ! rhoqtest

          end do                   ! I1
        else
          do i1 = 1, nref
            call calctref13(eryd, vref(i1), rmtref(i1), lmax, lm1, &     ! LLY Lloyd
              wn1, wn2, alpharef(0,i1), dalpharef(0,i1), lmax+1, lmgf0d) ! LLY Lloyd
            ! -------------------------------------------------------
            ! add second spin-block for relativistic calculation and transform
            ! from NREL to REL representation
            ! -------------------------------------------------------
            call cinit(lmmaxd*lmmaxd, w1)
            if (lmmaxd/=lm1*2) stop 'LMMAXD <> LM1*2 '
            do i = 1, lm1
              w1(i, i) = wn1(i, i)
              w1(lm1+i, lm1+i) = wn1(i, i)
            end do
            call changerep(w1, 'RLM>REL', trefll(1,1,i1), lmmaxd, lmmaxd, rc, crel, rrel, 'TREFLL', 0)
          end do
        end if
#ifdef CPP_TIMING
        call timing_pause('main1b - calctref13')
#endif
        ! ----------------------------------------------------------------
        tralpha(ie, ispin) = czero ! LLY
        tralpharef(ie) = czero     ! LLY
        ! read in t matrix
        do i1 = 1, natyp
          ! read in t-matrix from file
          if (t_tgmat%tmat_to_file) then
            irec = ie + ielast*(ispin-1) + ielast*nspin1*(i1-1)
            read (69, rec=irec) tmat
          else
            irec = ie_num + ie_end*(ispin-1) + ie_end*nspin1*(i1-1)
            tmat(:, :) = t_tgmat%tmat(:, :, irec)
          end if

          if (use_Chebychev_solver) then
            ! read in theta and phi for noncolinear
            theta = theta_at(i1)
            phi = phi_at(i1)

            ! rotate t-matrix from local to global frame
            call rotatematrix(tmat, theta, phi, lmgf0d, 0)
          end if

          tsst(1:lmmaxd, 1:lmmaxd, i1) = tmat(1:lmmaxd, 1:lmmaxd)

          if (lly/=0) then         ! LLY
            if (t_lloyd%dtmat_to_file) then
              irec = ie + ielast*(ispin-1) + ielast*nspin1*(i1-1)
              read (691, rec=irec) tmat ! LLY dt/dE
            else
              irec = ie_num + ie_end*(ispin-1) + ie_end*nspin1*(i1-1)
              tmat(:, :) = t_lloyd%dtmat(:, :, irec)
            end if

            if (use_Chebychev_solver) call rotatematrix(tmat, theta, phi, lmgf0d, 0) ! LLY

            dtmatll(1:lmmaxd, 1:lmmaxd, i1) = tmat(1:lmmaxd, 1:lmmaxd) ! LLY
            if (t_lloyd%dtmat_to_file) then
              irec = ie + ielast*(ispin-1) + ielast*nspin1*(i1-1)
              read (692, rec=irec) tralpha1 ! LLY
            else
              irec = ie_num + ie_end*(ispin-1) + ie_end*nspin1*(i1-1)
              tralpha1 = t_lloyd%tralpha(irec)
            end if

            tralpha(ie, ispin) = tralpha(ie, ispin) + tralpha1 ! LLY Tr[ alpha^{-1} dalpha/dE]

            if (ispin==1) then     ! Ref. system is spin-independent     ! LLY
              tralpha1 = czero     ! LLY
              do l1 = 0, lmax      ! LLY
                tralpha1 = tralpha1 + (2*l1+1)*dalpharef(l1, refpot(i1))/alpharef(l1, refpot(i1)) ! LLY
              end do
              tralpharef(ie) = tralpharef(ie) + tralpha1 ! LLY Tr[ alpharef^{-1} dalpharef/dE
            end if

          end if                   ! LLY

        end do                     ! i1 = 1,natyp

        if (t_lloyd%g0tr_to_file) then
          if (lly/=0 .and. ispin==1) read (682, fmt='(2E24.16)') lly_g0tr(ie) ! LLY
        else
          if (lly/=0 .and. ispin==1) lly_g0tr(ie) = t_lloyd%g0tr(ie_num) ! LLY
        end if

        !----------------------------------------------------------------------------
        ! qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos qdos
        !----------------------------------------------------------------------------
        if (use_readcpa .or. (use_qdos .and. (iqdosrun==1))) then ! qdos ruess: read in cpa t-matrix
          do isite = 1, naez                                                 ! qdos ruess
            tqdos(:, :, isite) = czero                                       ! qdos ruess

            if ( .not. (formatted_files .or. write_deci_tmat) ) then
              do lm1 = 1, lmmaxd
                do lm2 = 1, lmmaxd
                  irec = lm2 + (lm1-1)*lmmaxd + lmmaxd**2*(isite-1) + lmmaxd**2*naez*(ie-1) + lmmaxd**2*ielast*naez*(ispin-1)
                  read (37, rec=irec) tread
                  if ((lm1+lm2)/=0) then                      ! qdos ruess
                    tqdos(lm1, lm2, isite) = tread/cfctorinv  ! qdos ruess
                  end if                                      ! qdos ruess
                end do
              end do
            else ! 'filverb'
              read (37, *) text                               ! qdos ruess
              read (37, *) text                               ! qdos ruess
110           continue                                        ! qdos ruess
              read (37, fmt='(2i5,2d22.14)') lm1, lm2, tread  ! qdos ruess
              if ((lm1+lm2)/=0) then                          ! qdos ruess
                tqdos(lm1, lm2, isite) = tread/cfctorinv      ! qdos ruess
                if ((lm1+lm2)<2*lmmaxd) go to 110             ! qdos ruess
              end if                                          ! qdos ruess
            endif ! 'filverb'

          end do ! isite                                      ! qdos ruess
        end if                                                ! qdos ruess
        ! -------------------------------------------------------------------
        ! Loop over all QDOS points and change volume for KLOOPZ run accordingly
        ! -------------------------------------------------------------------
        do iq = 1, nqdos           ! qdos ruess
          if (use_qdos) bzkp(:, 1, 1) = qvec(:, iq) ! qdos ruess: Set q-point x,y,z

#ifdef CPP_TIMING
          call timing_start('main1b - kloopz')
#endif
          call kloopz1_qdos(eryd, gmatll, ins, alat, ie, igf, nshell, naez, nofks(nmesh), volbz(nmesh), bzkp(1,1,nmesh), volcub(1,nmesh), cls, nacls, naclsmax, ncls, rr, rbasis, &
            ezoa, atom, rcls, icc, ginp, ideci, lefttinvll(1,1,1,1,ie), righttinvll(1,1,1,1,ie), vacflag, nlbasis, nrbasis, factl, natomimp, nsymat, dsymll, ratom, rrot, nsh1, &
            nsh2, ijtabsym, ijtabsh, icheck, invmod, refpot, trefll, tsst, msst, cfctor, cfctorinv, crel, rc, rrel, srrel, irrel, nrrel, drotq, symunitary, kmrot, natyp, ncpa, &
            icpa, itcpamax, cpatol, noq, iqat, itoq, conc, iprint, icpaflag, ispin, nspindd, tqdos, iqdosrun, & ! qdos
            dtrefll, dtmatll, dginp, lly_grtr(ie,ispin), & ! LLY Lloyd
            tracet(ie,ispin), lly) ! LLY Lloyd

#ifdef CPP_TIMING
          call timing_pause('main1b - kloopz')
#endif
          ! -------------------------------------------------------------
          ! Skip this part if first part of the qdos is running
          ! -------------------------------------------------------------
          if (.not. (use_qdos .and. (iqdosrun==0))) then
            if (ncpa/=0) then
              if (icpaflag/=0) then
                ncpafail = ncpafail + 1
                iecpafail(ncpafail) = ie
              end if
            end if                 ! (NCPA.NE.0)

            do i1 = 1, nshell(0)
              gmat0(1:lmmaxd, 1:lmmaxd) = gmatll(1:lmmaxd, 1:lmmaxd, i1)
              irec = iq + nqdos*(ie-1) + nqdos*ielast*(ispin-1) + nqdos*ielast*nspin1*(i1-1) ! qdos ruess
              if (t_tgmat%gmat_to_file) then
                write (70, rec=irec) gmat0
                ! human readable writeout if test option is hit
                if (formatted_files) then
                  write (707070, '(i9,200000F15.7)') irec, gmat0
                end if
              else
                irec = iq + nqdos*(ie_num-1) + nqdos*ie_end*(ispin-1) + nqdos*ie_end*nspin1*(i1-1)
                t_tgmat%gmat(:, :, irec) = gmat0
              end if
            end do
            if (test('gmatasci')) then
              write (*, *) 'Writing out gmat.ascii'
              do i1 = 1, nshell(0)
                do lm1 = 1, lmmaxd
                  do lm2 = 1, lmmaxd
                    write (298347, fmt='(3I5,2E25.16)') i1, lm1, lm2, gmatll(lm1, lm2, i1)
                  end do
                end do
              end do
            end if

            ! writeout of host green function for impurity code for single-atom cluster (not captured in rotgll)
            if (natomimp==1) then
              i1 = atomimp(1)
              if (write_kkrimp_input) then
                ilm = 0
                gimp = czero ! complex*8
                do lm2 = 1, lmmaxd
                  do lm1 = 1, lmmaxd
                    ilm = ilm + 1
                    gimp(ilm) = gmatll(lm1, lm2, i1)
                  end do
                end do
                irec = ielast*(ispin-1) + ie + 1
                write (888, rec=irec) gimp
                if (write_gmat_plain) then
                  ! write(8888,'(50000E)') GIMP
                  write (8888, *) gimp
                end if
              end if               ! KKRFLEX
              if (write_green_host .and. myrank==master) then
                do lm2 = 1, lmmaxd
                  do lm1 = 1, lmmaxd
                    ! writeout of green_host for WRTGREEN option
                    write (58, '((2I5),(2e17.9))') lm2, lm1, gmatll(lm1, lm2, i1)
                  end do
                end do
              end if               ! WRTGREEN
            end if                 ! ( NATOMIMP==1 )

            if (lcpaij) then
              if (t_cpa%dmatproj_to_file) then
                do i1 = 1, natyp
                  gmat0(:, :) = tsst(:, :, i1)
                  w1(:, :) = msst(:, :, i1)
                  irec = ie + ielast*(ispin-1) + ielast*nspin1*(i1-1)
                  write (71, rec=irec) gmat0, w1
                end do             ! I1
              else                 ! t_cpa%dmatproj_to_file
                irec = ie_num + ie_end*(ispin-1)
                t_cpa%dmatts(:, :, :, irec) = tsst(:, :, :)
                t_cpa%dtilts(:, :, :, irec) = msst(:, :, :)
              end if               ! t_cpa%dmatproj_to_file
            end if                 ! ( LCPAIJ )

          end if                   ! ( .NOT.(use_qdos.AND.(IQDOSRUN.EQ.0)) )
        end do                     ! IQ = 1,NQDOS                                       ! qdos ruess


        if (lly/=0) then           ! LLY

          if (use_Chebychev_solver .and. .not.test('NOSOC   ')) then
            cdos_lly(ie, ispin) = tralpha(ie, ispin) - lly_grtr(ie, ispin)/volbz(1) + 2.0_dp*lly_g0tr(ie) ! LLY
          else
            if (lly/=2) then       ! LLY Lloyd
              cdos_lly(ie, ispin) = tralpha(ie, ispin) - lly_grtr(ie, ispin)/volbz(1) + lly_g0tr(ie) ! LLY Lloyd
            else                   ! LLY Lloyd
              cdos_lly(ie, ispin) = tracet(ie, ispin) + tralpharef(ie) - lly_grtr(ie, ispin)/volbz(1) + lly_g0tr(ie) ! LLY Lloyd
            end if                 ! LLY Lloyd
          end if

          if (set_gmat_to_zero) then ! LLY Lloyd
            cdos_lly(ie, ispin) = tralpha(ie, ispin) ! LLY Lloyd
            if (lly==2) then       ! LLY Lloyd
              cdos_lly(ie, ispin) = tracet(ie, ispin) + tralpharef(ie) ! LLY Lloyd
            end if                 ! LLY Lloyd
          end if                   ! LLY Lloyd

          cdos_lly(ie, ispin) = cdos_lly(ie, ispin)/pi ! LLY Lloyd

          if (ispin==1) cdosref_lly(ie) = tralpharef(ie) - lly_g0tr(ie) ! LLY Lloyd

        end if                     ! LLY

        ! ---------------------------------------------------------------

#ifdef CPP_MPI
        ! stop timing measurement for this ie, needed for MPIadapt
        if (mpiadapt>0 .and. t_mpi_c_grid%myrank_ie==0) then
          call timing_stop('time_1b_ie', save_out=timings_1b(ie))
        end if
#endif

      end do                       ! IE = 1,IELAST
#ifdef CPP_TIMING
      if (.not. write_green_imp) then
        if (t_inc%i_time>0) call timing_stop('main1b - calctref13')
        if (t_inc%i_time>0) call timing_stop('main1b - fourier')
        if (t_inc%i_time>0) call timing_stop('main1b - inversion')
        if (t_inc%i_time>0) call timing_stop('main1b - kloopz')
      end if
      if (t_inc%i_time>0 .and. write_rhoq_input) then
        call timing_stop('main1b - kkrmat01 - writeout_rhoq')
      end if
#endif

      if (ncpafail/=0) then
        if (t_inc%i_write>0) then
          write (1337, *)
          write (1337, '(1X,79(''*''),/)')
          write (1337, '(1X,79(''*''))')
          write (1337, '(" tolerance for CPA-cycle:",F15.7)') cpatol
          write (1337, '(" CPA not converged for",I3," energies:")') ncpafail
          write (1337, '(3(" E:",I3,F7.4,:,2X))')(iecpafail(ie), dble(ez(iecpafail(ie))), ie=1, ncpafail)
          write (1337, '(1X,79(''*''),/)')
          write (1337, *)
        end if
      else
        if (ncpa/=0) then
          if (t_inc%i_write>0) then
            write (1337, *)
            write (1337, 120)
            write (1337, *)
          end if
        end if
      end if

    end do                         ! ISPIN = 1,NSPIN1
    ! -------------------------------------------------------------------------
    ! END of do loop over spins and energies
    ! -------------------------------------------------------------------------

    ! close green_host file after write out
    if (write_green_host .and. myrank==master) then
      close (58)
    end if

    close (68)
    ! IF ( IGF.NE.0 ) CLOSE (88)  !no-green


    if (lly/=0) then               ! LLY Lloyd

      if (t_lloyd%cdos_diff_lly_to_file) then
        if (myrank==master) then
          open (701, file='cdosdiff_lly.dat', form='FORMATTED') ! LLY Lloyd
        end if
#ifdef CPP_MPI
        ihelp = ielast*nspin       ! IELAST*NSPIN
        allocate (work(ielast,nspin))
        work = czero
        call mpi_allreduce(cdos_lly, work, ihelp, mpi_double_complex, mpi_sum, t_mpi_c_grid%mympi_comm_at, ierr)
        call zcopy(ihelp, work, 1, cdos_lly, 1)
        deallocate (work)
#endif
      end if                       ! t_lloyd%cdos_diff_lly_to_file


      do ispin = 1, nspin/(1+korbit) ! LLY
#ifdef CPP_MPI
        ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
        ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
        ie_start = 0
        ie_end = ielast
#endif
        do ie_num = 1, ie_end
          ie = ie_start + ie_num
          if (t_lloyd%cdos_diff_lly_to_file .and. myrank==master) then
            write (701, fmt='(10E25.16)') ez(ie), cdos_lly(ie, ispin), tralpha(ie, ispin), lly_grtr(ie, ispin) ! LLY
          else
            t_lloyd%cdos(ie_num, ispin) = cdos_lly(ie, ispin)
          end if
        end do                     ! IE                                                     ! LLY
      end do                       ! ISPIN                                                     ! LLY
      if (t_lloyd%cdos_diff_lly_to_file .and. myrank==master) close (701) ! LLY
    end if                         ! LLY/=0                                                       ! LLY



    if ((calc_exchange_couplings) .and. (icc<=0)) then
#ifdef CPP_TIMING
      call timing_start('main1b - tbxccpl')
#endif
      if (nqdos/=1) stop 'QDOS option not compatible with XCPL'
      if (.not. use_Chebychev_solver) then
        call tbxccpljij(69,ielast,ez,wez,nspindd,ncpa,naez,natyp,noq,itoq,iqat,     &
          nshell,natomimp,atomimp,ratom,nofgij,nqcalc,iqcalc,ijtabcalc,ijtabsym,    &
          ijtabsh,ish, jsh, dsymll, iprint, natyp, nsheld, lmmaxd, npol)
      else                         ! .NOT.use_Chebychev_solver)
        call tbxccpljijdij(naez,natyp,lmmaxd,lmgf0d,natomimpd,iemxd,theta_at,phi_at,&
          natomimp,atomimp,nofgij,iqat,rclsimp,ijtabcalc,ijtabcalc_i,ijtabsh,       &
          ijtabsym,ielast, ez, wez, npol, dsymll, noq, itoq, ncpa)
      end if
#ifdef CPP_TIMING
      call timing_stop('main1b - tbxccpl')
#endif
    end if
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    close (69)
    close (70)
    if (write_pkkr_input .and. myrank==master) close (6801) ! fswrt
    if (lcpaij .and. t_cpa%dmatproj_to_file) close (71)

    close (37)                     ! qdos ruess
    !--------------------------------------------------------------------------------
    ! Finished first qdos run. Now re-run the whole kkr1b program to
    ! calculate the GF for every energy (defined in inputcard) and
    ! kpoint (defined in qvec.dat)
    !--------------------------------------------------------------------------------
    iqdosrun = iqdosrun + 1        ! qdos ruess
    if (t_lloyd%dgref_to_file) close (681)
    if (t_lloyd%g0tr_to_file) close (682)
    if (t_lloyd%dtmat_to_file) close (691)
    if (t_lloyd%tralpha_to_file) close (692)
    if (iqdosrun==1) go to 100     ! qdos ruess


    if (write_rhoq_input) then
#ifdef CPP_MPI
      call mpi_barrier(mpi_comm_world, ierr)
#endif
      close (9889)                 ! tau0_k file
      close (99992)
      call timing_stop('Time in Iteration')
      if (myrank==master) call print_time_and_date('Done w. rhoq!!!')
#ifdef CPP_MPI
      call mpi_barrier(mpi_comm_world, ierr)
      call mpi_finalize(ierr)
#endif
      stop 'finished with rhoq output'

    end if

    if (write_pkkr_operators .and. myrank==master) then
      ! check if impurity files are present (otherwise no imp.
      ! wavefunctions can be calculated)
      operator_imp = .true.
      inquire (file='potential_imp', exist=lexist)
      if (.not. lexist) operator_imp = .false.
      inquire (file='shapefun_imp', exist=lexist)
      if (.not. lexist) operator_imp = .false.
      inquire (file='scoef', exist=lexist)
      if (.not. lexist) operator_imp = .false.
    else
      operator_imp = .false.
    end if
#ifdef CPP_MPI
    call mpi_bcast(operator_imp, 1, mpi_logical, master, mpi_comm_world, ierr)
    if (ierr/=mpi_success) stop 'error broadcasting operator_imp'
#endif

    ! Do stuff for GREENIMP option (previously done in zulapi code)
    ! run for OPERATOR option to precalculate impurity wavefunctions
    ! that are then stored in t_imp (used only if potential_imp,
    ! scoef, shapefun_imp)
    if (write_green_imp .or. operator_imp) then

      ! consistency checks
      if (.not. (ielast==1 .or. ielast==3)) stop 'Error: GREENIMP option only possible with 1 () or 3 () energy points in contour'
      if (ielast==1 .and. abs(aimag(ez(1)))>1e-10) stop 'Error: T>0 for GREENIMP (DTMTRX writeout, IELAST==3)'
      if (ielast==3 .and. abs(aimag(ez(1)))<1e-10) stop 'Error: T==0 for GREENIMP (GMATLL_GES writeout, IELAST==3)'
      ! end consistency checks

#ifdef CPP_MPI
      ! init arrays and communicate parameters of t_imp for all ranks
      ! that are not the master
      call bcast_t_imp_scalars(t_imp, master)
      if (myrank/=master) call init_t_imp(t_inc, t_imp)
      call bcast_t_imp_arrays(t_imp, t_inc, master)
#endif

      do ie = 1, ielast            ! big ie loop (use only for GMATLL output)
        if (ielast==3 .and. ie==1 .and. myrank==master) then
          open (unit=59, file='GMATLL_GES', form='FORMATTED')
          open (unit=60, file='green_host', form='FORMATTED')
        end if
        allocate (dtmtrx(lmmaxd*t_imp%natomimp,lmmaxd*t_imp%natomimp), stat=i_stat)
        call memocc(i_stat, product(shape(dtmtrx))*kind(dtmtrx), 'DTMTRX', 'main1b')
        dtmtrx = czero

        ! find DTMTRX (written out for IELAST==1), parallelized with
        ! mpi over atoms
        call tmatimp_newsolver(irmd,nsra-1,lmax,iend,irid,lpotd,natyp,ncleb,ipand,  &
          irnsd,nfund,t_imp%ihost,ntotd,nspin,lmpotd,ncheb,lmmaxd/(1+korbit),korbit,&
          nspotd,ielast,irmind,t_params%npan_eq,t_params%npan_log,t_imp%natomimp,   &
          r_log, vins, visp, ipan, irmin, t_imp%hostimp(1:t_imp%ihost),          &
          t_imp%ipanimp(1:t_imp%natomimp), t_imp%irwsimp(1:t_imp%natomimp),         &
          atomimp(1:t_imp%natomimp), t_imp%irminimp(1:t_imp%natomimp), icleb, ircut,&
          t_imp%ircutimp(0:ipand,1:t_imp%natomimp),zat,t_imp%zimp(1:t_imp%natomimp),&
          rmesh,cleb(1,1),t_imp%rimp(1:irmd,1:t_imp%natomimp),rclsimp,ez(ie),       &
          t_imp%vispimp,t_imp%vinsimp, dtmtrx, lmmaxso)

        ! compute GMATLL_GES, on master rank only
        if (ielast==3 .and. myrank==master) then
          call greenimp(t_imp%natomimp, dtmtrx, ez(ie))
        end if

        i_all = -product(shape(dtmtrx))*kind(dtmtrx)
        deallocate (dtmtrx, stat=i_stat)
        call memocc(i_stat, i_all, 'DTMTRX', 'main1b')

      end do                       ! ie-loop

      ! done with GREENIMP option, stopping now
      if (.not. write_pkkr_operators) then
        if (myrank==master) write (*, *) 'done with GREENIMP, stop here!'
#ifdef CPP_MPI
        call mpi_finalize(ierr)
#endif
        stop
      end if

    end if                         ! GREENIMP .or. OPERATOR

    ! ------------------------------------------------------------------------
    ! determine the spin operator, torque operator and spin flux operator
    ! used in FScode do compute spin expectation values etc. within Boltzmann
    ! formalism
    ! ------------------------------------------------------------------------
    if (write_pkkr_operators) then
#ifdef CPP_TIMING
      call timing_start('main1b - operator')
#endif

      call operators_for_fscode(korbit, operator_imp)

#ifdef CPP_TIMING
      call timing_stop('main1b - operator')
#endif
    end if                         ! OPERATOR

    ! ----------------------------------------------------------------------

    if (write_rhoq_input .and. (myrank==master)) then
      open (9999, file='params.txt')
      write (9999, *) 2*lmmaxd, t_params%natyp
      write (9999, *) t_params%naez, t_params%nclsd, t_params%nr, t_params%nembd1, t_params%lmax
      write (9999, *) t_params%alat, naclsmax
      close (9999)

      open (9999, file='host.txt')
      write (9999, *) t_params%rbasis(1:3, 1:t_params%natyp)
      write (9999, *) t_params%rcls(1:3, 1:t_params%nclsd, 1:t_params%nclsd), t_params%rr(1:3, 0:t_params%nr), t_params%atom(1:t_params%nclsd, 1:t_params%naez+t_params%nembd1)
      write (9999, *) t_params%cls(1:t_params%naez+t_params%nembd1), t_params%ezoa(1:t_params%nclsd, 1:t_params%naez+t_params%nembd1), t_params%nacls(1:t_params%nclsd)
      close (9999)
    end if

    ! deallocate arrays
    deallocate (bzkp, stat=i_stat)
    call memocc(i_stat, -product(shape(bzkp))*kind(bzkp), 'BZKP', 'main1b')
    deallocate (lly_g0tr, stat=i_stat)
    call memocc(i_stat, -product(shape(lly_g0tr))*kind(lly_g0tr), 'LLY_G0TR', 'main1b')
    deallocate (tralpharef, stat=i_stat)
    call memocc(i_stat, -product(shape(tralpharef))*kind(tralpharef), 'TRALPHAREF', 'main1b')
    deallocate (cdosref_lly, stat=i_stat)
    call memocc(i_stat, -product(shape(cdosref_lly))*kind(cdosref_lly), 'CDOSREF_LLY', 'main1b')
    deallocate (tracet, stat=i_stat)
    call memocc(i_stat, -product(shape(tracet))*kind(tracet), 'TRACET', 'main1b')
    deallocate (tralpha, stat=i_stat)
    call memocc(i_stat, -product(shape(tralpha))*kind(tralpha), 'TRALPHA', 'main1b')
    deallocate (lly_grtr, stat=i_stat)
    call memocc(i_stat, -product(shape(lly_grtr))*kind(lly_grtr), 'LLY_GRTR', 'main1b')
    deallocate (cdos_lly, stat=i_stat)
    call memocc(i_stat, -product(shape(cdos_lly))*kind(cdos_lly), 'CDOS_LLY', 'main1b')

    if (t_inc%i_write>0) write (1337, '(79("="),/,30X,"< KKR1b finished >",/,79("="),/)')

    ! 99019 FORMAT('(/,1X,79(*),/," tolerance for CPA-cycle:",F15.7,/," CPA not converged for",I3," energies:",/,3(" E:",I3,F7.4,:,2X))')
120 format ('(/,1X,79(*),/,25X,"no problems with","  CPA-cycle ",/,1X,79(*),/)')

  end subroutine main1b

end module mod_main1b
