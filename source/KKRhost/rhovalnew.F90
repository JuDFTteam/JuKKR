
!-----------------------------------------------------------------------------------------!
! Copyright (c) 2018 Peter Grünberg Institut, Forschungszentrum Jülich, Germany           !
! This file is part of Jülich KKR code and available as free software under the conditions!
! of the MIT license as expressed in the LICENSE.md file in more detail.                  !
!-----------------------------------------------------------------------------------------!

!------------------------------------------------------------------------------------
!> Summary: Calculation of the density for the new solver 
!> Author:
!> Calculation of the density for the new solver
!------------------------------------------------------------------------------------
module mod_rhovalnew

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the density for the new solver
  !> Author: 
  !> Category: physical-observables, KKRhost 
  !> Deprecated: False 
  !> Calculation of the density for the new solver
  !-------------------------------------------------------------------------------
  subroutine rhovalnew(ldorhoef, ielast, nsra, nspin, lmax, ez, wez, zat, socscale, cleb, icleb, iend, ifunm, lmsp, ncheb, &
    npan_tot, npan_log, npan_eq, rmesh, irws, rpan_intervall, ipan_intervall, rnew, vinsnew, thetasnew, theta, phi, fixdir, i1, ipot, &
    den_out, espv, rho2ns, r2nef, muorb, angles_new, totmoment, idoldau, lopt, wldau, denmatn, natyp, ispin)

#ifdef CPP_OMP
    use :: omp_lib
#endif
#ifdef CPP_MPI
    use :: mpi
    use :: mod_types, only: gather_tmat, t_mpi_c_grid, save_t_mpi_c_grid, get_ntot_pt_ioff_pt_2d
    use :: mod_mympi, only: find_dims_2d, distribute_linear_on_tasks, mympi_main1c_comm_newsosol
#endif
#ifdef CPP_TIMING
    use :: mod_timing
#endif
    use :: mod_types, only: t_tgmat, t_inc
    use :: mod_mympi, only: myrank, master
    use :: mod_save_wavefun, only: t_wavefunctions, read_wavefunc
    use :: mod_runoptions, only: calc_exchange_couplings, calc_gmat_lm_full, disable_tmat_sratrick, fix_nonco_angles, &
                                 use_qdos, write_complex_qdos, write_pkkr_operators, write_DOS_lm, set_cheby_nospeedup, &
                                 decouple_spins_cheby, disable_print_serialnumber
    use :: mod_version_info, only: version_print_header
    use :: global_variables, only: lmmaxd, iemxd, ncleb, lmxspd, irmd, ntotd, nrmaxd, lmpotd, nspotd, nfund, korbit, mmaxd, nspind
    use :: mod_constants, only: czero, cvlight, cone, pi, ci
    use :: mod_profiling, only: memocc 
    use :: mod_datatypes, only: dp
    use :: mod_rhooutnew, only:  rhooutnew
    use :: mod_calcsph, only: calcsph
    use :: mod_cheb2oldgrid, only: cheb2oldgrid
    use :: mod_rll_global_solutions, only: rll_global_solutions
    use :: mod_intcheb_cell, only: intcheb_cell
    use :: mod_rllsllsourceterms, only: rllsllsourceterms
    use :: mod_rllsll, only: rllsll
    use :: mod_spinorbit_ham, only: spinorbit_ham
    use :: mod_sll_global_solutions, only: sll_global_solutions
    use :: mod_rotatespinframe, only: rotatematrix, rotatevector
    use :: mod_vllmatsra, only: vllmatsra
    use :: mod_vllmat, only: vllmat
    use :: mod_wunfiles, only: t_params
    use :: mod_bfield, only: add_bfield
    use :: mod_torque, only: calc_torque

    implicit none

    integer, intent (in) :: ispin
    integer, intent (in) :: i1
    integer, intent (in) :: nsra
    integer, intent (in) :: lmax   !! Maximum l component in wave function expansion
    integer, intent (in) :: iend   !! Number of nonzero gaunt coefficients
    integer, intent (in) :: ipot
    integer, intent (in) :: irws   !! R point at WS radius for a given atom
    integer, intent (in) :: lopt   !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: ncheb  !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: ielast
    integer, intent (in) :: idoldau !! flag to perform LDA+U
    integer, intent (in) :: npan_eq !! Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
    integer, intent (in) :: npan_tot
    integer, intent (in) :: npan_log !! Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
    real (kind=dp), intent (in) :: zat !! Nuclear charge for a given atom
    real (kind=dp), intent (in) :: socscale !! Spin-orbit scaling for a given atom
    logical, intent (in) :: ldorhoef
    integer, dimension (lmxspd), intent (in) :: lmsp !! 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (lmxspd), intent (in) :: ifunm
    integer, dimension (0:ntotd), intent (in) :: ipan_intervall
    integer, dimension (ncleb, 4), intent (in) :: icleb
    real (kind=dp), dimension (*), intent (in) :: cleb !! GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irmd), intent (in) :: rmesh
    real (kind=dp), dimension (mmaxd, mmaxd, nspind), intent (in) :: wldau !! potential matrix

    ! .. In/Out variables
    real (kind=dp), intent (inout) :: phi
    real (kind=dp), intent (inout) :: theta
    logical, intent(in) :: fixdir
    real (kind=dp), dimension (nrmaxd), intent (inout) :: rnew
    real (kind=dp), dimension (0:ntotd), intent (inout) :: rpan_intervall
    real (kind=dp), dimension (0:lmax+2, 3), intent (inout) :: muorb
    real (kind=dp), dimension (nrmaxd, nfund), intent (inout) :: thetasnew
    real (kind=dp), dimension (nrmaxd, lmpotd, nspotd), intent (inout) :: vinsnew !! Non-spherical part of the potential
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    complex (kind=dp), dimension (iemxd), intent (inout) :: wez
    ! .. Output variables
    real (kind=dp), dimension (2), intent (out) :: angles_new
    real (kind=dp), dimension (0:lmax+1, 2/(nspin-korbit)), intent (out) :: espv
    real (kind=dp), dimension (irmd, lmpotd, nspin/(nspin-korbit)*(1+korbit)), intent (out) :: r2nef
    real (kind=dp), dimension (irmd, lmpotd, nspin/(nspin-korbit)*(1+korbit)), intent (out) :: rho2ns
    real (kind=dp), intent(out) :: totmoment
    complex (kind=dp), dimension (0:lmax+1, ielast, nspin/(nspin-korbit)), intent (out) :: den_out

    ! .. Local variables
    integer :: lmmax0d !! (lmax+1)**2
    integer :: lmaxd1
    integer :: ir, irec, use_sratrick, nvec, lm1, lm2, ie, irmdnew, imt1, jspin, idim, iorb, l1
    integer :: i_stat, i_all
    integer :: use_fullgmat !! use (l,m,s) coupled matrices or not for 'NOSOC' test option (1/0)
    integer :: iq, nqdos           ! qdos ruess: number of qdos points
    integer :: lrecgflle, ierr     ! lmlm-dos
    integer :: lmlo, lmhi, is, js, mmax ! LDAU
    integer :: ix, m1              ! qdos ruess

    real (kind=dp) :: thetanew, phinew
    real (kind=dp) :: totxymoment
    complex (kind=dp) :: ek
    complex (kind=dp) :: df
    complex (kind=dp) :: eryd
    complex (kind=dp) :: temp1
    complex (kind=dp) :: dentemp
    complex (kind=dp) :: gmatprefactor
    integer, dimension (nsra*lmmaxd) :: jlk_index
    real (kind=dp), dimension (3) :: moment
    real (kind=dp), dimension (3) :: denorbmom
    real (kind=dp), dimension (3) :: denorbmomns
    real (kind=dp), dimension (2, 3) :: denorbmomsp
    real (kind=dp), dimension (0:lmax, 3) :: denorbmomlm
    complex (kind=dp), dimension (nspin/(nspin-korbit)*(1+korbit)) :: rho2
    complex (kind=dp), dimension (nspin/(nspin-korbit)*(1+korbit)) :: rho2int
    complex (kind=dp), dimension (nspin/(nspin-korbit)*(lmax+1)) :: alphasph
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gmat0
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gldau ! LDAU
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: tmatll
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: alphall ! LLY
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: tmattemp
    complex (kind=dp), dimension (2, 2) :: rho2ns_temp
    complex (kind=dp), dimension (lmmaxd, lmmaxd, iemxd) :: gmatll
    complex (kind=dp), dimension (mmaxd, mmaxd, 2, 2) :: denmatn ! LDAU

    ! .. Local allocatable arrays
    real (kind=dp), dimension (:, :), allocatable :: qvec ! qdos ruess: q-vectors for qdos
    real (kind=dp), dimension (:, :, :), allocatable :: vins  ! vins in the new mesh!! corresponds to the vinsnew in main1c but only for the two spin channels
    complex (kind=dp), dimension (:, :), allocatable :: tmatsph
    complex (kind=dp), dimension (:, :), allocatable :: cdentemp
    complex (kind=dp), dimension (:, :), allocatable :: rhotemp
    complex (kind=dp), dimension (:, :), allocatable :: rhonewtemp
    complex (kind=dp), dimension (:, :, :), allocatable :: hlk
    complex (kind=dp), dimension (:, :, :), allocatable :: jlk
    complex (kind=dp), dimension (:, :, :), allocatable :: hlk2
    complex (kind=dp), dimension (:, :, :), allocatable :: jlk2
    complex (kind=dp), dimension (:, :, :), allocatable :: cdenns
    complex (kind=dp), dimension (:, :, :), allocatable :: r2nefc
    complex (kind=dp), dimension (:, :, :), allocatable :: rho2nsc
    complex (kind=dp), dimension (:, :, :), allocatable :: vnspll0
    complex (kind=dp), dimension (:, :, :), allocatable :: r2nefnew
    complex (kind=dp), dimension (:, :, :), allocatable :: rho2nsnew
    complex (kind=dp), dimension (:, :, :), allocatable :: gflle_part
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rll
    complex (kind=dp), dimension (:, :, :, :), allocatable :: sll
    complex (kind=dp), dimension (:, :, :, :), allocatable :: den
    complex (kind=dp), dimension (:, :, :, :), allocatable :: cden
    complex (kind=dp), dimension (:, :, :, :), allocatable :: gflle
    complex (kind=dp), dimension (:, :, :, :), allocatable :: denlm
    complex (kind=dp), dimension (:, :, :, :), allocatable :: cdenlm
    complex (kind=dp), dimension (:, :, :, :), allocatable :: vnspll
    complex (kind=dp), dimension (:, :, :, :), allocatable :: r2orbc
    complex (kind=dp), dimension (:, :, :, :), allocatable :: vnspll1
    complex (kind=dp), dimension (:, :, :), allocatable :: vnspll2
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rllleft
    complex (kind=dp), dimension (:, :, :, :), allocatable :: sllleft
    complex (kind=dp), dimension (:, :, :, :), allocatable :: r2nefc_loop
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rho2nsc_loop

#ifdef CPP_MPI
    complex (kind=dp), dimension (2) :: dentot ! qdos ruess
    ! communication
    complex (kind=dp), dimension (:, :, :, :), allocatable :: workc
#endif
    ! OMP - number of threads, thread id
    integer :: nth, ith
    integer :: ie_start, ie_end, ie_num
    integer :: i1_myrank           ! lmlm-dos, needed for MPI with more than one rank per energy point (nranks_ie>1)
    ! read in wavefunctions
    logical :: rll_was_read_in, sll_was_read_in, rllleft_was_read_in, sllleft_was_read_in

    ! lmmax0d is original lm-size (without enhancement through soc etc.)
    lmmax0d = lmmaxd/(1+korbit)

    ! determine if omp is used
    ith = 0
    nth = 1
#ifdef CPP_OMP
    ! $omp parallel shared(nth,ith)
    ! $omp single
    nth = omp_get_num_threads()
    if (t_inc%i_write>0) write (1337, *) 'nth =', nth
    ! $omp end single
    ! $omp end parallel
#endif

    ! .. Parameters
    lmaxd1 = lmax + 1

    if (nsra==2) then
      use_sratrick = 1
      if (disable_tmat_sratrick) use_sratrick = 0
    else
      use_sratrick = 0
    end if

    irmdnew = npan_tot*(ncheb+1)
    imt1 = ipan_intervall(npan_log+npan_eq) + 1
    allocate (vins(irmdnew,lmpotd,nspin/(nspin-korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(vins))*kind(vins), 'VINS', 'RHOVALNEW')
    vins = 0.0_dp
    vins(1:irmdnew, 1:lmpotd, 1) = vinsnew(1:irmdnew, 1:lmpotd, ipot)
    if (.not.decouple_spins_cheby)  vins(1:irmdnew, 1:lmpotd, nspin) = vinsnew(1:irmdnew, 1:lmpotd, ipot+nspin-1)

    ! set up the non-spherical ll' matrix for potential VLL'
    allocate (vnspll0(lmmaxd,lmmaxd,irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'RHOVALNEW')
    vnspll0 = czero
    allocate (vnspll1(lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'RHOVALNEW')
    vnspll1 = czero
    allocate (vnspll2(lmmaxd,lmmaxd,irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(vnspll2))*kind(vnspll2), 'VNSPLL2', 'RHOVALNEW')
    vnspll2 = czero

    call vllmat(1, nrmaxd, irmdnew, lmmax0d, lmmaxd, vnspll0, vins, lmpotd, cleb, icleb, iend, nspin/(nspin-korbit), zat, rnew, use_sratrick, ncleb)
    !--------------------------------------------------------------------------------
    ! LDAU
    !--------------------------------------------------------------------------------
    if (idoldau==1) then
      lmlo = lopt**2 + 1
      lmhi = (lopt+1)**2
      do ir = 1, irmdnew
        vnspll0(lmlo:lmhi, lmlo:lmhi, ir) = vnspll0(lmlo:lmhi, lmlo:lmhi, ir) + wldau(1:mmaxd, 1:mmaxd, 1)
      end do
      lmlo = lmlo + lmmax0d
      lmhi = lmhi + lmmax0d
      do ir = 1, irmdnew
        vnspll0(lmlo:lmhi, lmlo:lmhi, ir) = vnspll0(lmlo:lmhi, lmlo:lmhi, ir) + wldau(1:mmaxd, 1:mmaxd, 2)
      end do
    end if
    !--------------------------------------------------------------------------------
    ! LDAU
    !--------------------------------------------------------------------------------

    ! initial allocate
    if (nsra==2) then
      allocate (vnspll(nsra*lmmaxd,nsra*lmmaxd,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'RHOVALNEW')
    else
      allocate (vnspll(lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'RHOVALNEW')
    end if
    vnspll(:,:,:,:) = czero

    allocate (hlk(nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(hlk))*kind(hlk), 'HLK', 'RHOVALNEW')
    hlk(:,:,:) = czero
    allocate (jlk(nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(jlk))*kind(jlk), 'JLK', 'RHOVALNEW')
    jlk(:,:,:) = czero
    allocate (hlk2(nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(hlk2))*kind(hlk2), 'HLK2', 'RHOVALNEW')
    hlk2(:,:,:) = czero
    allocate (jlk2(nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(jlk2))*kind(jlk2), 'JLK2', 'RHOVALNEW')
    jlk2(:,:,:) = czero
    allocate (tmatsph(nspin/(nspin-korbit)*(lmax+1),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(tmatsph))*kind(tmatsph), 'TMATSPH', 'RHOVALNEW')
    tmatsph(:,:) = czero
    allocate (rll(nsra*lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(rll))*kind(rll), 'RLL', 'RHOVALNEW')
    rll(:,:,:,:) = czero
    allocate (sll(nsra*lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(sll))*kind(sll), 'SLL', 'RHOVALNEW')
    sll(:,:,:,:) = czero
    allocate (rllleft(nsra*lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'RHOVALNEW')
    rllleft(:,:,:,:) = czero
    allocate (sllleft(nsra*lmmaxd,lmmaxd,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'RHOVALNEW')
    sllleft(:,:,:,:) = czero
    allocate (cden(irmdnew,0:lmax,nspin/(nspin-korbit)*(1+korbit),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cden))*kind(cden), 'CDEN', 'RHOVALNEW')
    cden(:,:,:,:) = czero
    allocate (cdenlm(irmdnew,lmmax0d,nspin/(nspin-korbit)*(1+korbit),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdenlm))*kind(cdenlm), 'CDENLM', 'RHOVALNEW')
    cdenlm(:,:,:,:) = czero
    allocate (cdenns(irmdnew,nspin/(nspin-korbit)*(1+korbit),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdenns))*kind(cdenns), 'CDENNS', 'RHOVALNEW')
    cdenns(:,:,:) = czero
    allocate (rho2nsc(irmdnew,lmpotd,nspin/(nspin-korbit)*(1+korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsc))*kind(rho2nsc), 'RHO2NSC', 'RHOVALNEW')
    rho2nsc(:,:,:) = czero
    allocate (rho2nsc_loop(irmdnew,lmpotd,nspin/(nspin-korbit)*(1+korbit),ielast), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsc_loop))*kind(rho2nsc_loop), 'RHO2NSC_loop', 'RHOVALNEW')
    rho2nsc_loop(:,:,:,:) = czero
    allocate (rho2nsnew(irmd,lmpotd,nspin/(nspin-korbit)*(1+korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsnew))*kind(rho2nsnew), 'RHO2NSNEW', 'RHOVALNEW')
    rho2nsnew(:,:,:) = czero
    allocate (r2nefc(irmdnew,lmpotd,nspin/(nspin-korbit)*(1+korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefc))*kind(r2nefc), 'R2NEFC', 'RHOVALNEW')
    r2nefc(:,:,:) = czero
    allocate (r2nefc_loop(irmdnew,lmpotd,nspin/(nspin-korbit)*(1+korbit),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefc_loop))*kind(r2nefc_loop), 'R2NEFC_loop', 'RHOVALNEW')
    r2nefc_loop(:,:,:,:) = czero
    allocate (r2nefnew(irmd,lmpotd,nspin/(nspin-korbit)*(1+korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefnew))*kind(r2nefnew), 'R2NEFNEW', 'RHOVALNEW')
    r2nefnew(:,:,:) = czero
    allocate (r2orbc(irmdnew,lmpotd,nspin/(nspin-korbit)*(1+korbit),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(r2orbc))*kind(r2orbc), 'R2ORBC', 'RHOVALNEW')
    r2orbc(:,:,:,:) = czero
    allocate (cdentemp(irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdentemp))*kind(cdentemp), 'CDENTEMP', 'RHOVALNEW')
    cdentemp(:,:) = czero
    allocate (gflle_part(lmmaxd,lmmaxd,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(gflle_part))*kind(gflle_part), 'GFLLE_PART', 'RHOVALNEW')
    gflle_part(:,:,:) = czero
    allocate (gflle(lmmaxd,lmmaxd,ielast,1), stat=i_stat)
    call memocc(i_stat, product(shape(gflle))*kind(gflle), 'GFLLE', 'RHOVALNEW')
    gflle(:,:,:,:) = czero
    allocate (den(0:lmaxd1,ielast,1,nspin/(nspin-korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(den))*kind(den), 'DEN', 'RHOVALNEW')
    den(:,:,:,:) = czero
    allocate (denlm(lmmax0d,ielast,1,nspin/(nspin-korbit)), stat=i_stat)
    call memocc(i_stat, product(shape(den))*kind(den), 'DENLM', 'RHOVALNEW')
    denlm(:,:,:,:) = czero

    rho2ns = 0.0_dp                  ! fivos 19.7.2014, this was CZERO
    r2nef = 0.0_dp                   ! fivos 19.7.2014, this was CZERO
    espv = 0.0_dp
    rho2int = czero
    denorbmom = 0.0_dp
    denorbmomsp = 0.0_dp
    denorbmomlm = 0.0_dp
    denorbmomns = 0.0_dp
    thetanew = 0.0_dp
    phinew = 0.0_dp
    gldau = czero

    ! DO IR=1,3
    ! DO LM1=0,LMAXD1+1
    ! MUORB(LM1,IR)=0_dp  !zimmer: initialization shifted to main1c
    ! ENDDO
    ! ENDDO

    nqdos = 1                      ! qdos ruess
    if (use_qdos) then      ! qdos ruess
      ! Read BZ path for qdos calculation:                                ! qdos ruess
      open (67, file='qvec.dat', status='old', iostat=ierr, err=100) ! qdos ruess
      read (67, *) nqdos           ! qdos ruess
      allocate (qvec(3,nqdos), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(qvec))*kind(qvec), 'QVEC', 'RHOVALNEW') ! qdos ruess
      do iq = 1, nqdos             ! qdos ruess
        read (67, *)(qvec(ix,iq), ix=1, 3) ! qdos ruess
      end do                       ! qdos ruess
      close (67)                   ! qdos ruess
      ! Change allocation for GFLLE to be suitabel for qdos run           ! qdos ruess
      i_all = -product(shape(gflle))*kind(gflle) ! qdos ruess
      deallocate (gflle, stat=i_stat) ! qdos ruess
      call memocc(i_stat, i_all, 'GFLLE', 'RHOVALNEW') ! qdos ruess
      i_all = -product(shape(den))*kind(den) ! qdos ruess
      deallocate (den, stat=i_stat) ! qdos ruess
      call memocc(i_stat, i_all, 'DEN', 'RHOVALNEW') ! qdos ruess
      i_all = -product(shape(denlm))*kind(denlm) ! qdos ruess
      deallocate (denlm, stat=i_stat) ! qdos ruess
      call memocc(i_stat, i_all, 'DENLM', 'RHOVALNEW') ! qdos ruess
      ! ! qdos ruess
      allocate (gflle(lmmaxd,lmmaxd,ielast,nqdos), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(gflle))*kind(gflle), 'GFLLE', 'RHOVALNEW') ! qdos ruess
      allocate (den(0:lmaxd1,ielast,nqdos,nspin/(nspin-korbit)), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(den))*kind(den), 'DEN', 'RHOVALNEW') ! qdos ruess
      allocate (denlm(lmmax0d,ielast,nqdos,nspin/(nspin-korbit)), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(denlm))*kind(qvec), 'DENLM', 'RHOVALNEW') ! qdos ruess
100   if (ierr/=0) stop 'ERROR READING ''qvec.dat''' ! qdos ruess
    end if                         ! use_qdos                                                     ! qdos ruess

#ifdef CPP_MPI
    i1_myrank = i1 - t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) ! lmlm-dos ruess
#else
    i1_myrank = i1                 ! lmlm-dos ruess
#endif
    if ((calc_gmat_lm_full) .and. (i1_myrank==1)) then ! lmlm-dos ruess
      lrecgflle = nspin*(1+korbit)*lmmaxd*lmmaxd*ielast*nqdos ! lmlm-dos ruess
      open (91, access='direct', recl=lrecgflle, file='gflle' & ! lmlm-dos ruess
        , form='unformatted', status='replace', err=110, iostat=ierr) ! lmlm-dos ruess
110   if (ierr/=0) stop 'ERROR CREATING ''gflle''' ! lmlm-dos ruess
    end if                         ! lmlm-dos ruess

    ! initialize to zero
    den = czero
    denlm = czero
    ! energy loop
    if (myrank==master .and. t_inc%i_write>0) write (1337, *) 'atom: ', i1
#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    ie_start = 0                   ! offset
    ie_end = ielast
#endif

#ifdef CPP_OMP
    ! omp: start parallel region here
    ! $omp parallel do default(none)                                         &
    ! $omp private(eryd,ie,ir,irec,lm1,lm2,gmatprefactor,nvec)               &
    ! $omp private(jlk_index,tmatll,ith)                                     &
    ! $omp private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)              &
    ! $omp private(rho2ns_temp,rho2,temp1,jspin)                             &
    ! $omp private(alphasph,alphall,ie_num)                                  &
    ! $omp private(rll_was_read_in, sll_was_read_in)                         &
    ! $omp private(rllleft_was_read_in, sllleft_was_read_in)                 &
    ! !$omp firstprivate(t_inc)                                              &
    ! $omp shared(t_inc)                                                     &
    ! $omp shared(ldorhoef,nqdos,wez,lmsp,imt1,ifunm)      &
    ! $omp shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)             &
    ! $omp shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)           &
    ! $omp shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)    &
    ! $omp shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)         &
    ! $omp shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)               &
    ! $omp shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)         &
    ! $omp shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)        &
    ! $omp shared(t_tgmat,ie_end, ie_start, t_wavefunctions)                 &
    ! $omp shared(lmmaxd,lmmax0d,lmpotd,NRMAXD,NTOTD,LMAXD1)                  &
    ! $omp reduction(+:rho2int,espv) reduction(-:muorb)                      &
    ! $omp reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
#endif
    do ie_num = 1, ie_end
      ie = ie_start + ie_num

#ifdef CPP_OMP
      ith = omp_get_thread_num()
#else
      ith = 0
#endif

      eryd = ez(ie)
      ek = sqrt(eryd)

      ! set energy integration weight
      df = wez(ie)/real(nspin, kind=dp)

      if (nsra==2) then
        ek = sqrt(eryd+eryd*eryd/(cvlight*cvlight))*(1.0_dp+eryd/(cvlight*cvlight))
      end if
#ifdef CPP_OMP
      ! $omp critical
#endif
      if (t_inc%i_write>0) write (1337, *) 'energy:', ie, '', eryd
#ifdef CPP_OMP
      ! $omp end critical
#endif

      if (t_wavefunctions%nwfsavemax>0) then ! read wavefunctions?
        ! read in wavefunction from memory
        call read_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,i1,ie,nsra,      &
          lmmaxd,irmdnew,ith,nth,rll_was_read_in,sll_was_read_in,                  &
          rllleft_was_read_in,sllleft_was_read_in)
      end if

      ! recalculate wavefuntions, also include left solution
      ! contruct the spin-orbit coupling hamiltonian and add to potential
      if ( .not. decouple_spins_cheby) then
        call spinorbit_ham(lmax, lmmax0d, vins, rnew, eryd, zat, cvlight, socscale, nspin, lmpotd, theta, phi, ipan_intervall, &
          rpan_intervall, npan_tot, ncheb, irmdnew, nrmaxd, vnspll0, vnspll2(:,:,:), '1')
      else
        vnspll2(:,:,:) = vnspll0(:,:,:)
      end if
  
      ! Add magnetic field
      if ( t_params%bfield%lbfield .or. t_params%bfield%lbfield_constr ) then
        imt1 = ipan_intervall(t_params%npan_log+t_params%npan_eq) + 1
        call add_bfield(t_params%bfield,i1,lmax,nspin,irmdnew,imt1,iend,ncheb,theta,phi,t_params%ifunm1(:,t_params%ntcell(i1)),&
                        t_params%icleb,t_params%cleb(:,1),t_params%thetasnew(1:irmdnew,:,t_params%ntcell(i1)),'1',vnspll2(:,:,:), &
                        vnspll1(:,:,:,ith),t_params%bfield%thetallmat(:,:,1:irmdnew,t_params%ntcell(i1)))
      else
        vnspll1(:,:,:,ith) = vnspll2(:,:,:)
      end if

      ! extend matrix for the SRA treatment
      vnspll(:, :, :, ith) = czero
      if (nsra==2) then
        if (use_sratrick==0) then
          call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxd,irmdnew, &
            nrmaxd,eryd,lmax,0,'Ref=0')
        else if (use_sratrick==1) then
          call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxd,irmdnew, &
            nrmaxd,eryd,lmax,0,'Ref=Vsph')
        end if
      else
        vnspll(:, :, :, ith) = vnspll1(:, :, :, ith)
      end if

      if ((t_wavefunctions%nwfsavemax>0 .and. .not. rll_was_read_in) .or. (t_wavefunctions%nwfsavemax==0)) then ! read/recalc wavefunctions

        ! calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        hlk(:, :, ith) = czero
        jlk(:, :, ith) = czero
        hlk2(:, :, ith) = czero
        jlk2(:, :, ith) = czero
        gmatprefactor = czero
        jlk_index = 0
        if (decouple_spins_cheby) then
          use_fullgmat = 0
        else
          use_fullgmat = 1
        end if
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxd, use_fullgmat, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)
        ! using spherical potential as reference
        if (use_sratrick==1) then
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin/(nspin-korbit), zat, eryd, lmpotd, lmmaxd, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk(:,:,ith), jlk(:,:,ith), &
            hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, tmatsph(:,ith), alphasph, use_sratrick)
        end if

        ! calculate the tmat and wavefunctions
        rll(:, :, :, ith) = czero
        sll(:, :, :, ith) = czero

        !----------------------------------------------------------------------------
        ! Right solutions
        !----------------------------------------------------------------------------
        tmatll = czero
        alphall = czero
        ! faster calculation of RLL.
        ! no irregular solutions SLL are needed in self-consistent iterations
        ! because the density depends only on RLL, RLLLEFT and SLLLEFT
        if (.not.set_cheby_nospeedup .and. .not. (calc_exchange_couplings .or. write_pkkr_operators)) then
          call rll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), tmatll, ncheb, npan_tot, lmmaxd, nvec*lmmaxd, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', use_sratrick, alphall)
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), sll(:,:,:,ith), tmatll, ncheb, npan_tot, lmmaxd, nvec*lmmaxd, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alphall)
        end if
        if (nsra==2) then
          rll(lmmaxd+1:nvec*lmmaxd, :, :, ith) = rll(lmmaxd+1:nvec*lmmaxd, :, :, ith)/cvlight
          sll(lmmaxd+1:nvec*lmmaxd, :, :, ith) = sll(lmmaxd+1:nvec*lmmaxd, :, :, ith)/cvlight
        end if

      end if                       ! read/recalc wavefunctions

      !------------------------------------------------------------------------------
      ! Left solutions
      !------------------------------------------------------------------------------
      if ((t_wavefunctions%nwfsavemax>0 .and. (.not. (rllleft_was_read_in .and. sllleft_was_read_in))) .or. (t_wavefunctions%nwfsavemax==0)) then
        ! read/recalc wavefunctions left contruct the TRANSPOSE spin-orbit coupling hamiltonian and add to potential
        if ( .not. decouple_spins_cheby) then
          call spinorbit_ham(lmax, lmmax0d, vins, rnew, eryd, zat, cvlight, socscale, nspin, lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, irmdnew, nrmaxd, &
            vnspll0, vnspll2(:,:,:), 'transpose')
        else
          vnspll2(:,:,:) = vnspll0(:,:,:)
        end if
        
        ! Add magnetic field 
        if ( t_params%bfield%lbfield .or. t_params%bfield%lbfield_constr ) then
          call add_bfield(t_params%bfield,i1,lmax,nspin,irmdnew,imt1,iend,ncheb,theta,phi,t_params%ifunm1(:,t_params%ntcell(i1)),&
                          t_params%icleb,t_params%cleb(:,1),t_params%thetasnew(1:irmdnew,:,t_params%ntcell(i1)),'transpose',vnspll2(:,:,:), &
                          vnspll1(:,:,:,ith),t_params%bfield%thetallmat(:,:,1:irmdnew,t_params%ntcell(i1)))
        else
          vnspll1(:,:,:,ith) = vnspll2(:,:,:)
        end if
        ! extend matrix for the SRA treatment
        vnspll(:, :, :, ith) = czero
        if (nsra==2) then
          if (use_sratrick==0) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxd,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=0')
          else if (use_sratrick==1) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxd,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=Vsph')
          end if
        else
          vnspll(:, :, :, ith) = vnspll1(:, :, :, ith)
        end if

        ! calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        hlk(:, :, ith) = czero
        jlk(:, :, ith) = czero
        hlk2(:, :, ith) = czero
        jlk2(:, :, ith) = czero
        gmatprefactor = czero
        jlk_index = 0
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxd, use_fullgmat, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)

        ! using spherical potential as reference
        ! notice that exchange the order of left and right hankel/bessel functions
        if (use_sratrick==1) then
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin/(nspin-korbit), zat, eryd, lmpotd, lmmaxd, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), &
            hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, alphasph, tmatsph(:,ith), use_sratrick)
        end if

        ! calculate the tmat and wavefunctions
        rllleft(:, :, :, ith) = czero
        sllleft(:, :, :, ith) = czero

        ! left solutions
        ! notice that exchange the order of left and right hankel/bessel functions
        tmattemp = czero
        alphall = czero
        ! faster calculation of RLLLEFT and SLLLEFT.
        if (.not.set_cheby_nospeedup .and. .not. (calc_exchange_couplings .or. write_pkkr_operators)) then
          call rll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), rllleft(:,:,:,ith), tmattemp, ncheb, npan_tot, lmmaxd, nvec*lmmaxd, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', use_sratrick, alphall)
          call sll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), sllleft(:,:,:,ith), ncheb, npan_tot, lmmaxd, nvec*lmmaxd, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', use_sratrick)
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rllleft(:,:,:,ith), sllleft(:,:,:,ith), tmattemp, ncheb, npan_tot, lmmaxd, nvec*lmmaxd, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alphall)
        end if
        if (nsra==2) then
          rllleft(lmmaxd+1:nvec*lmmaxd, :, :, ith) = rllleft(lmmaxd+1:nvec*lmmaxd, :, :, ith)/cvlight
          sllleft(lmmaxd+1:nvec*lmmaxd, :, :, ith) = sllleft(lmmaxd+1:nvec*lmmaxd, :, :, ith)/cvlight
        end if
      end if                       ! read/recalc wavefunctions left

      do iq = 1, nqdos             ! qdos
        ! read in GF
#ifdef CPP_OMP
        ! $omp critical
#endif
        if (t_tgmat%gmat_to_file) then
          irec = iq + nqdos*(ie-1) + nqdos*ielast*(ispin-1) + nqdos*ielast*nspin/(1+korbit)*(i1-1)
          read (70, rec=irec) gmat0
        else
          irec = iq + nqdos*(ie_num-1) + nqdos*ie_end*(ispin-1) + nqdos*ie_end*nspin/(1+korbit)*(i1-1)
          gmat0(:, :) = t_tgmat%gmat(:, :, irec)
        end if
#ifdef CPP_OMP
        ! $omp end critical
#endif

        if ( .not. decouple_spins_cheby) then
          ! rotate gmat from global frame to local frame
          call rotatematrix(gmat0, theta, phi, lmmax0d, 1)
        end if

        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            gmatll(lm1, lm2, ie) = gmat0(lm1, lm2)
          end do
        end do
        ! calculate density
        call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, df, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), commented out since sll is not used in rhooutnew
          rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), rho2nsc_loop(:,:,:,ie), 0, gflle(:,:,ie,iq), rpan_intervall, ipan_intervall, nspin/(nspin-korbit))

        do jspin = 1, nspin/(nspin-korbit)*(1+korbit)
          do lm1 = 0, lmax
            cdentemp(:, ith) = czero
            dentemp = czero
            do ir = 1, irmdnew
              cdentemp(ir, ith) = cden(ir, lm1, jspin, ith)
            end do
            call intcheb_cell(cdentemp(:,ith),dentemp,rpan_intervall,ipan_intervall,&
              npan_tot,ncheb,irmdnew)
            rho2(jspin) = dentemp
            rho2int(jspin) = rho2int(jspin) + rho2(jspin)*df
            if (jspin<=2) then
              den(lm1, ie, iq, jspin) = rho2(jspin)
            end if
          end do

          if (jspin<=2) then
            do lm1 = 1, lmmax0d
              cdentemp(:, ith) = czero
              dentemp = czero
              do ir = 1, irmdnew
                cdentemp(ir, ith) = cdenlm(ir, lm1, jspin, ith)
              end do
              call intcheb_cell(cdentemp(:,ith),dentemp,rpan_intervall,             &
                ipan_intervall,npan_tot,ncheb,irmdnew)
              denlm(lm1, ie, iq, jspin) = dentemp
            end do
          end if
          cdentemp(:, ith) = czero
          dentemp = czero
          cdentemp(1:irmdnew, ith) = cdenns(1:irmdnew, jspin, ith)
          call intcheb_cell(cdentemp(:,ith), dentemp, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
          rho2int(jspin) = rho2int(jspin) + dentemp*df
          if (jspin<=2) then
            den(lmaxd1, ie, iq, jspin) = dentemp
            espv(0:lmaxd1, jspin) = espv(0:lmaxd1, jspin) + aimag(eryd*den(0:lmaxd1,ie,iq,jspin)*df)
          end if
        end do                     ! JSPIN

      end do                       ! IQ = 1,NQDOS

      !------------------------------------------------------------------------------
      ! Get charge at the Fermi energy (IELAST)
      !------------------------------------------------------------------------------
      if (ie==ielast .and. ldorhoef) then
        call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, cone, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
          rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), r2nefc_loop(:,:,:,ith), 0, gflle_part(:,:,ith), rpan_intervall, &
          ipan_intervall, nspin/(nspin-korbit))
      end if

      !------------------------------------------------------------------------------
      ! Get orbital moment
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not. decouple_spins_cheby) then
        do iorb = 1, 3
          call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, cone, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
            rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), r2orbc(:,:,:,ith), iorb, gflle_part(:,:,ith), rpan_intervall, ipan_intervall, nspin)
          do jspin = 1, nspin*(1+korbit)
            if (jspin<=2) then
              do lm1 = 0, lmax
                cdentemp(:, ith) = czero
                dentemp = czero
                do ir = 1, irmdnew
                  cdentemp(ir, ith) = cden(ir, lm1, jspin, ith)
                end do
                call intcheb_cell(cdentemp(:,ith), dentemp, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
        
                rho2(jspin) = dentemp
                muorb(lm1, jspin) = muorb(lm1, jspin) - aimag(rho2(jspin)*df)
                denorbmom(iorb) = denorbmom(iorb) - aimag(rho2(jspin)*df)
                denorbmomsp(jspin, iorb) = denorbmomsp(jspin, iorb) - aimag(rho2(jspin)*df)
                denorbmomlm(lm1, iorb) = denorbmomlm(lm1, iorb) - aimag(rho2(jspin)*df)
                cdentemp(:, ith) = czero
        
                do ir = 1, irmdnew
                  cdentemp(ir, ith) = cdenns(ir, jspin, ith)
                end do
                call intcheb_cell(cdentemp(:,ith), temp1, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
                denorbmomns(iorb) = denorbmomns(iorb) - aimag(temp1*df)
              end do ! lm1
            end if
          end do ! jspin
          ! fill summed values of orbital moment
          do jspin=1, nspin
            ! fill muorb(:,3,:) with sum of both spin channels
            muorb(0:lmaxd1, 3) = muorb(0:lmaxd1, 3) + muorb(0:lmaxd1, jspin)
          end do
          ! sum over l-channels
          do lm1 = 0, lmaxd1
            muorb(lmaxd1+1, 1:3) = muorb(lmaxd1+1, 1:3) + muorb(lm1, 1:3)
          end do
        end do ! IORB
      end if ! .not. decouple_spins_cheby

    end do                         ! IE loop
#ifdef CPP_OMP
    ! $omp end parallel do
#endif

    ! omp: move sum from rhooutnew here after parallel calculation
    do ie = 1, ielast
      rho2nsc(:, :, :) = rho2nsc(:, :, :) + rho2nsc_loop(:, :, :, ie)
    end do
    ! omp: don't forget to do the same with density at fermi energy:
    do ith = 0, nth - 1
      r2nefc(:, :, :) = r2nefc(:, :, :) + r2nefc_loop(:, :, :, ith)
    end do

#ifdef CPP_MPI
    if (use_qdos) then                                                                                    ! qdos
      ! first communicate den array to write out qdos files                                                      ! qdos
      idim = (lmaxd1+1)*ielast*nspin/(nspin-korbit)*nqdos                                                                       ! qdos
      allocate (workc(0:lmaxd1,ielast,nspin/(nspin-korbit),nqdos), stat=i_stat)                                                 ! qdos
      call memocc(i_stat, product(shape(workc))*kind(workc), 'workc', 'RHOVALNEW')                               ! qdos
      workc = czero                                                                                              ! qdos
      call mpi_reduce(den, workc, idim, mpi_double_complex, mpi_sum, master, t_mpi_c_grid%mympi_comm_at, ierr)   ! qdos
      call zcopy(idim, workc, 1, den, 1)                                                                         ! qdos
      i_all = -product(shape(workc))*kind(workc)                                                                 ! qdos
      deallocate (workc, stat=i_stat)                                                                            ! qdos
      call memocc(i_stat, i_all, 'workc', 'RHOVALNEW')                                                           ! qdos
                                                                                                                 
      if (t_mpi_c_grid%myrank_at==master) then                                                                   ! qdos
        ie_start = 0                                                                                             ! qdos
        ie_end = ielast                                                                                          ! qdos
        do ie_num = 1, ie_end                                                                                    ! qdos
          ie = ie_start + ie_num                                                                                 ! qdos
          do iq = 1, nqdos                                                                                       ! qdos
            if ((iq==1) .and. (ie_num==1)) then                                                                  ! qdos
              if (natyp>=100) then                                                                               ! qdos
                open (31, file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! qdos
                open (32, file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! qdos
              else                                                                                               ! qdos
                open (31, file='qdos.'//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+1)//'.dat')    ! qdos
                open (32, file='qdos.'//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+2)//'.dat')    ! qdos
              end if                                                                                             ! qdos
              call version_print_header(31, disable_print=disable_print_serialnumber)                            ! qdos
              write (31, *) ' '                                                                                  ! qdos
              write (31, 150) '# ISPIN=', 1, ' I1=', i1                                                          ! qdos
              if (write_DOS_lm) then 
                write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s, px, pz, py, dx2-y2, dxz, dz2, dyz, dxy, ..., ns'  ! lm-qdos
              else
                write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...,ns'   ! qdos
              end if
              if (nspin>1) then                                                                                  ! qdos
                call version_print_header(32, disable_print=disable_print_serialnumber)                          ! qdos
                write (32, *) ' '                                                                                ! qdos
                write (32, 150) '# ISPIN=', 2, ' I1=', i1                                                        ! qdos
                if (write_DOS_lm) then 
                  write (32, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s, px, pz, py, dx2-y2, dxz, dz2, dyz, dxy, ..., ns'  ! lm-qdos
                else
                  write (32, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...,ns' ! qdos
                end if
              end if                                                                                             ! qdos
            end if ! IQ.EQ.1                                                                                     ! qdos
            do jspin = 1, nspin/(nspin-korbit)                                                                   ! qdos
              dentot(jspin) = czero                                                                              ! qdos
              do l1 = 0, lmaxd1                                                                                  ! qdos
                dentot(jspin) = dentot(jspin) + den(l1, ie, iq, jspin)                                           ! qdos
              end do                                                                                             ! qdos
            end do                                                                                               ! qdos
            ! write qdos.nn.s.dat                                                                                ! qdos
            if (write_DOS_lm) then
              write (31, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), -aimag(dentot(1))/pi, (-aimag(denlm(l1,ie,iq,1))/pi, l1=1, lmmax0d), -aimag(den(lmaxd1,ie,iq,1))/pi ! lm-dos
              write (32, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), -aimag(dentot(2))/pi, (-aimag(denlm(l1,ie,iq,2))/pi, l1=1, lmmax0d), -aimag(den(lmaxd1,ie,iq,2))/pi ! lm-dos
            else
              write (31, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), -aimag(dentot(1))/pi, (-aimag(den(l1,ie,iq,1))/pi, l1=0, lmaxd1) ! qdos
              write (32, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), -aimag(dentot(2))/pi, (-aimag(den(l1,ie,iq,2))/pi, l1=0, lmaxd1) ! qdos
            end if
120         format (5f10.6, 40e16.8)                                                                             ! qdos

            if (write_complex_qdos) then                                                                                            ! complex qdos
              if ((iq==1) .and. (ie_num==1)) then                                                                                 ! complex qdos
                if (natyp>=100) then                                                                                              ! complex qdos
                  open (31, file='cqdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! complex qdos
                  open (32, file='cqdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! complex qdos
                else                                                                                                              ! complex qdos
                  open (31, file='cqdos.'//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+1)//'.dat')                  ! complex qdos
                  open (32, file='cqdos.'//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+2)//'.dat')                  ! complex qdos
                end if                                                                                                            ! complex qdos
                call version_print_header(31, disable_print=disable_print_serialnumber)                                           ! complex qdos
                write (31, *) ' '                                                                                                 ! complex qdos
                write (31, '(A)') '#   lmax, natyp, nspin, nqdos, ielast:'                                                        ! complex qdos
                write (31, '(5I9)') lmax, natyp, nspin, nqdos, ielast                                                             ! complex qdos
                write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...'                       ! complex qdos
                if (nspin>1) then                                                                                                 ! complex qdos
                  call version_print_header(32, disable_print=disable_print_serialnumber)                                         ! complex qdos
                  write (32, *) ' '                                                                                               ! complex qdos
                  write (32, '(A)') '# lmax, natyp, nspin, nqdos, ielast:'                                                        ! complex qdos
                  write (32, '(5I9)') lmax, natyp, nspin, nqdos, ielast                                                           ! complex qdos
                  write (32, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...'                     ! complex qdos
                end if                                                                                                            ! complex qdos
              end if ! IQ.EQ.1                                                                                                    ! complex qdos
              do jspin = 1, nspin/(nspin-korbit)                                                                                  ! complex qdos
                dentot(jspin) = czero                                                                                             ! complex qdos
                do l1 = 0, lmaxd1                                                                                                 ! complex qdos
                  dentot(jspin) = dentot(jspin) + den(l1, ie, iq, jspin)                                                          ! complex qdos
                end do                                                                                                            ! complex qdos
              end do                                                                                                              ! complex qdos
              ! write qdos.nn.s.dat                                                                                               ! complex qdos
              write (31, 130) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), dentot(1), (den(l1,ie,iq,1), l1=0, lmaxd1)           ! complex qdos
              write (32, 130) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), dentot(2), (den(l1,ie,iq,2), l1=0, lmaxd1)           ! complex qdos
130           format (6f10.6, 80e16.8)                                                                                            ! complex qdos
            end if                                                                                                                ! complex qdos
                     
          end do ! IQ             ! qdos
        end do ! IE               ! qdos
      end if ! myrank_at==master  ! qdos
    end if ! use_qdos      ! qdos
#endif

#ifdef CPP_MPI
    ! do communication only when compiled with MPI
#ifdef CPP_TIMING
    call timing_start('main1c - communication')
#endif
    ! reset NQDOS to avoid endless communication
    if (.not. calc_gmat_lm_full) then
      nqdos = 1
    else
      if (myrank==master) write (*, *) '< calc_gmat_lm_full > option, communcation might take a while!', ielast, nqdos
    end if
    ! set these arrays to zero to avoid double counting in cases where extra ranks are used
    if (t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1)) then
      den = czero
      denlm = czero
      gflle = czero
      r2nefc = czero
      rho2nsc = czero
      rho2int = czero
      muorb = 0.0_dp
      espv = 0.0_dp
      denorbmom = 0.0_dp
      denorbmomsp = 0.0_dp
      denorbmomlm = 0.0_dp
      denorbmomns = 0.0_dp
    end if
    call mympi_main1c_comm_newsosol(nspin/(nspin-korbit), korbit, irmdnew, lmpotd, lmax, lmaxd1, lmmax0d, lmmaxd, ielast, nqdos, den, denlm, &
      gflle, rho2nsc, r2nefc, rho2int, espv, muorb, denorbmom, denorbmomsp, denorbmomlm, denorbmomns, t_mpi_c_grid%mympi_comm_at)
#ifdef CPP_TIMING
    call timing_pause('main1c - communication')
#endif

    ! MPI: do these writeout/data collection steps only on master and broadcast important results afterwards
    if (t_mpi_c_grid%myrank_at==master) then
!    if (myrank==master) then
#endif
      ! CPP_MPI
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Magnetic torques
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(t_params%bfield%ltorque) then
        call calc_torque(i1,lmax,irmdnew,nspin,rpan_intervall,ipan_intervall,npan_tot,ncheb,theta,phi,rho2nsc,vins, &
                         t_params%ifunm1(:,t_params%ntcell(i1)), iend, t_params%icleb,t_params%cleb(:,1),&
                         t_params%thetasnew(1:irmdnew,:,t_params%ntcell(i1)))
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (idoldau==1) then
        ! calculate WLDAU
        do ie = 1, ielast
          do lm1 = 1, lmmaxd
            do lm2 = 1, lmmaxd
              gldau(lm1, lm2) = gldau(lm1, lm2) + gflle(lm1, lm2, ie, 1)*wez(ie)/real(nspin, kind=dp)
            end do
          end do
        end do
        ! calculate occupation matrix
        mmax = 2*lopt + 1
        do is = 1, 2
          do js = 1, 2
            lmlo = lopt**2 + 1 + (is-1)*lmmax0d
            lmhi = (lopt+1)**2 + (js-1)*lmmax0d
            lm2 = lopt**2 + 1 + (js-1)*lmmax0d
            do m1 = 1, mmax
              lm1 = lmlo - 1 + m1
              denmatn(1:mmax, m1, js, is) = (1.0/(2.0*ci))*(gldau(lm2:lmhi,lm1)-conjg(gldau(lm1,lm2:lmhi)))
            end do
          end do
        end do
      end if                       ! LDAU
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (.not. use_qdos) then
        ! omp: moved write-out of dos files out of parallel energy loop
        ! Write out lm-dos:                                                     ! lm-dos
        if (write_DOS_lm) then ! qdos ruess
          do ie = 1, ielast        ! lm-dos
            iq = 1                 ! lm-dos
            if (ie==1) then        ! lm-dos
              if (natyp>=100) then ! lm-dos
                open (29, file='lmdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! lm-dos
                if (nspin==2) open (30, file='lmdos.'//char(48+i1/100)//char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! lm-dos
              else                 ! lm-dos
                open (29, file='lmdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! lm-dos
                if (nspin==2) open (30, file='lmdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! lm-dos
              end if               ! lm-dos
              call version_print_header(29, disable_print=disable_print_serialnumber) ! lm-dos
              write (29, *) ' '    ! lm-dos
              write (29, 150) '# ISPIN=', 1, ' I1=', i1 ! lm-dos
              if (nspin==2) call version_print_header(30, disable_print=disable_print_serialnumber) ! lm-dos
              if (nspin==2) write (30, *) ' '    ! lm-dos
              if (nspin==2) write (30, 150) '# ISPIN=', 2, ' I1=', i1 ! lm-dos
            end if                 ! IE==1                                                      ! lm-dos
            write (29, 140) ez(ie), (-aimag(denlm(l1,ie,iq,1))/pi, l1=1, lmmax0d) ! lm-dos
            if (nspin==2) write (30, 140) ez(ie), (-aimag(denlm(l1,ie,iq,nspin))/pi, l1=1, lmmax0d) ! lm-dos
140         format (30e12.4)       ! lm-dos
150         format (a8, i3, a4, i5) ! lm-dos/qdos ruess
          end do                   ! IE
        end if
      end if                       ! .not. use_qdos

      ! write gflle to file                                                 ! lmlm-dos
      if (calc_gmat_lm_full) then    ! lmlm-dos
        if (t_inc%i_write>0) then  ! lmlm-dos
          write (1337, *) 'gflle:', shape(gflle), shape(gflle_part), lrecgflle ! lmlm-dos
        end if                     ! lmlm-dos
        write (91, rec=i1) gflle   ! lmlm-dos
      end if                       ! lmlm-dos

      allocate (rhotemp(irmdnew,lmpotd), stat=i_stat)
      call memocc(i_stat, product(shape(rhotemp))*kind(rhotemp), 'RHOTEMP', 'RHOVALNEW')
      allocate (rhonewtemp(irws,lmpotd), stat=i_stat)
      call memocc(i_stat, product(shape(rhonewtemp))*kind(rhonewtemp), 'RHONEWTEMP', 'RHOVALNEW')

      do jspin = 1, nspin/(nspin-korbit)*(1+korbit)
        rhotemp = czero
        rhonewtemp = czero
        do lm1 = 1, lmpotd
          do ir = 1, irmdnew
            rhotemp(ir, lm1) = rho2nsc(ir, lm1, jspin)
          end do
        end do
        call cheb2oldgrid(irws,irmdnew,lmpotd,rmesh,ncheb,npan_tot,rpan_intervall,  &
          ipan_intervall,rhotemp,rhonewtemp,irmd)
        do lm1 = 1, lmpotd
          do ir = 1, irws
            rho2nsnew(ir, lm1, jspin) = rhonewtemp(ir, lm1)
          end do
        end do

        rhotemp = czero
        rhonewtemp = czero
        rhotemp(1:irmdnew, 1:lmpotd) = r2nefc(1:irmdnew, 1:lmpotd, jspin)
        call cheb2oldgrid(irws, irmdnew, lmpotd, rmesh, ncheb, npan_tot, rpan_intervall, ipan_intervall, rhotemp, rhonewtemp, irmd)
        r2nefnew(1:irws, 1:lmpotd, jspin) = rhonewtemp(1:irws, 1:lmpotd)
      end do

      i_all = -product(shape(rhotemp))*kind(rhotemp)
      deallocate (rhotemp, stat=i_stat)
      call memocc(i_stat, i_all, 'RHOTEMP', 'RHOVALNEW')
      i_all = -product(shape(rhonewtemp))*kind(rhonewtemp)
      deallocate (rhonewtemp, stat=i_stat)
      call memocc(i_stat, i_all, 'RHONEWTEMP', 'RHOVALNEW')
      rho2ns(1:irmd, 1:lmpotd, 1:nspin/(nspin-korbit)) = aimag(rho2nsnew(1:irmd, 1:lmpotd,1:nspin/(nspin-korbit)))
      r2nef(1:irmd, 1:lmpotd, 1:nspin/(nspin-korbit)) = aimag(r2nefnew(1:irmd, 1:lmpotd,1:nspin/(nspin-korbit)))
      ! MdSD: Should this also be corrected if the angles change?
      den_out(0:lmaxd1, 1:ielast, 1:nspin/(nspin-korbit)) = den(0:lmaxd1, 1:ielast, 1, 1:nspin/(nspin-korbit))
      ! calculate new THETA and PHI for non-colinear
      ! if (.not. fix_nonco_angles .and. .not.decouple_spins_cheby) then
      ! MdSD: now the new directions are always calculated, which can be useful for information purposes
      if (.not.decouple_spins_cheby) then
        rho2ns_temp(1, 1) = rho2int(1)
        rho2ns_temp(2, 2) = rho2int(2)
        rho2ns_temp(1, 2) = rho2int(3)
        rho2ns_temp(2, 1) = rho2int(4)

        call rotatematrix(rho2ns_temp, theta, phi, 1, 0)

        rho2int(1) = rho2ns_temp(1, 1)
        rho2int(2) = rho2ns_temp(2, 2)
        rho2int(3) = rho2ns_temp(1, 2)
        rho2int(4) = rho2ns_temp(2, 1)

        moment(1) = aimag(rho2int(3)+rho2int(4))
        moment(2) = -real(rho2int(3)-rho2int(4))
        moment(3) = aimag(-rho2int(1)+rho2int(2))
        
        totmoment = sqrt(moment(1)**2+moment(2)**2+moment(3)**2)
        totxymoment = sqrt(moment(1)**2+moment(2)**2)

        ! MdSD: theta not 0 or pi
        if (abs(totxymoment)>1e-05_dp) then
          thetanew = acos(moment(3)/totmoment)
          phinew = atan2(moment(2), moment(1))
        ! MdSD: theta is 0 or pi
        else
          if (moment(3) < 0.0_dp .and. abs(moment(3)) > 1e-14_dp) then
            thetanew = pi
          else
            thetanew = 0.0_dp
          end if
          phinew = 0.0_dp
        end if

        if (t_inc%i_write>0) then
          write (1337, '(A,i5,3es16.7)') 'moment', myrank, moment(1), moment(2), moment(3)
          write (1337, '(2es16.7)') thetanew/(2.0_dp*pi)*360.0_dp, phinew/(2.0_dp*pi)*360.0_dp
        end if
        ! only on master different from zero:
        angles_new(1) = thetanew
        angles_new(2) = phinew
        ! MdSD: use new angles to correct local frame, which defines the z-component of the spin density
        if (.not.fixdir) then
          call rotatevector(rho2nsnew,rho2ns,irws,lmpotd,thetanew,phinew,theta,phi,irmd)
          call rotatevector(r2nefnew,r2nef,irws,lmpotd,thetanew,phinew,theta,phi,irmd)
        end if
      end if

#ifdef CPP_MPI
    end if                         ! (myrank==master)
    
    ! communicate den_out to all processors with the same atom number
    idim = (lmax+2)*ielast*nspin/(nspin-korbit)
    call mpi_bcast(den_out, idim, mpi_double_complex, master, t_mpi_c_grid%mympi_comm_at, ierr)
    if (ierr/=mpi_success) stop 'error bcast den_out in rhovalnew'
    idim = 2
    call mpi_bcast(angles_new,idim,mpi_double_precision,master,t_mpi_c_grid%mympi_comm_at,ierr)
    if (ierr/=mpi_success) stop 'error bcast angles_new in rhovalnew'
    call mpi_bcast(totmoment,1,mpi_double_precision,master,t_mpi_c_grid%mympi_comm_at,ierr)
    if (ierr/=mpi_success) stop 'error bcast totmoment in rhovalnew'
#endif

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Deallocate arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    i_all = -product(shape(vins))*kind(vins)
    deallocate (vins, stat=i_stat)
    call memocc(i_stat, i_all, 'VINS', 'RHOVALNEW')
    i_all = -product(shape(vnspll0))*kind(vnspll0)
    deallocate (vnspll0, stat=i_stat)
    call memocc(i_stat, i_all, 'VNSPLL0', 'RHOVALNEW')
    i_all = -product(shape(vnspll1))*kind(vnspll1)
    deallocate (vnspll1, stat=i_stat)
    i_all = -product(shape(vnspll2))*kind(vnspll2)
    deallocate (vnspll2, stat=i_stat)
    call memocc(i_stat, i_all, 'VNSPLL1', 'RHOVALNEW')
    i_all = -product(shape(vnspll))*kind(vnspll)
    deallocate (vnspll, stat=i_stat)
    call memocc(i_stat, i_all, 'VNSPLL', 'RHOVALNEW')
    i_all = -product(shape(hlk))*kind(hlk)
    deallocate (hlk, stat=i_stat)
    call memocc(i_stat, i_all, 'HLK', 'RHOVALNEW')
    i_all = -product(shape(jlk))*kind(jlk)
    deallocate (jlk, stat=i_stat)
    call memocc(i_stat, i_all, 'JLK', 'RHOVALNEW')
    i_all = -product(shape(hlk2))*kind(hlk2)
    deallocate (hlk2, stat=i_stat)
    call memocc(i_stat, i_all, 'HLK2', 'RHOVALNEW')
    i_all = -product(shape(jlk2))*kind(jlk2)
    deallocate (jlk2, stat=i_stat)
    call memocc(i_stat, i_all, 'JLK2', 'RHOVALNEW')
    i_all = -product(shape(tmatsph))*kind(tmatsph)
    deallocate (tmatsph, stat=i_stat)
    call memocc(i_stat, i_all, 'TMATSPH', 'RHOVALNEW')
    i_all = -product(shape(rll))*kind(rll)
    deallocate (rll, stat=i_stat)
    call memocc(i_stat, i_all, 'RLL', 'RHOVALNEW')
    i_all = -product(shape(sll))*kind(sll)
    deallocate (sll, stat=i_stat)
    call memocc(i_stat, i_all, 'SLL', 'RHOVALNEW')
    i_all = -product(shape(rllleft))*kind(rllleft)
    deallocate (rllleft, stat=i_stat)
    call memocc(i_stat, i_all, 'RLLLEFT', 'RHOVALNEW')
    i_all = -product(shape(sllleft))*kind(sllleft)
    deallocate (sllleft, stat=i_stat)
    call memocc(i_stat, i_all, 'SLLLEFT', 'RHOVALNEW')
    i_all = -product(shape(cden))*kind(cden)
    deallocate (cden, stat=i_stat)
    call memocc(i_stat, i_all, 'CDEN', 'RHOVALNEW')
    i_all = -product(shape(cdenlm))*kind(cdenlm)
    deallocate (cdenlm, stat=i_stat)
    call memocc(i_stat, i_all, 'CDENLM', 'RHOVALNEW')
    i_all = -product(shape(cdenns))*kind(cdenns)
    deallocate (cdenns, stat=i_stat)
    call memocc(i_stat, i_all, 'CDENNS', 'RHOVALNEW')
    i_all = -product(shape(rho2nsc))*kind(rho2nsc)
    deallocate (rho2nsc, stat=i_stat)
    call memocc(i_stat, i_all, 'RHO2NSC', 'RHOVALNEW')
    i_all = -product(shape(rho2nsc_loop))*kind(rho2nsc_loop)
    deallocate (rho2nsc_loop, stat=i_stat)
    call memocc(i_stat, i_all, 'RHO2NSC_loop', 'RHOVALNEW')
    i_all = -product(shape(rho2nsnew))*kind(rho2nsnew)
    deallocate (rho2nsnew, stat=i_stat)
    call memocc(i_stat, i_all, 'RHO2NSNEW', 'RHOVALNEW')
    i_all = -product(shape(r2nefc))*kind(r2nefc)
    deallocate (r2nefc, stat=i_stat)
    call memocc(i_stat, i_all, 'R2NEFC', 'RHOVALNEW')
    i_all = -product(shape(r2nefc_loop))*kind(r2nefc_loop)
    deallocate (r2nefc_loop, stat=i_stat)
    call memocc(i_stat, i_all, 'R2NEFC_loop', 'RHOVALNEW')
    i_all = -product(shape(r2nefnew))*kind(r2nefnew)
    deallocate (r2nefnew, stat=i_stat)
    call memocc(i_stat, i_all, 'R2NEFNEW', 'RHOVALNEW')
    i_all = -product(shape(r2orbc))*kind(r2orbc)
    deallocate (r2orbc, stat=i_stat)
    call memocc(i_stat, i_all, 'R2ORBC', 'RHOVALNEW')
    i_all = -product(shape(cdentemp))*kind(cdentemp)
    deallocate (cdentemp, stat=i_stat)
    call memocc(i_stat, i_all, 'CDENTEMP', 'RHOVALNEW')
    i_all = -product(shape(gflle_part))*kind(gflle_part)
    deallocate (gflle_part, stat=i_stat)
    call memocc(i_stat, i_all, 'GFLLE_PART', 'RHOVALNEW')
    i_all = -product(shape(gflle))*kind(gflle)
    deallocate (gflle, stat=i_stat)
    call memocc(i_stat, i_all, 'GFLLE', 'RHOVALNEW')
    i_all = -product(shape(den))*kind(den)
    deallocate (den, stat=i_stat)
    call memocc(i_stat, i_all, 'DEN', 'RHOVALNEW')
    i_all = -product(shape(denlm))*kind(denlm)
    deallocate (denlm, stat=i_stat)
    call memocc(i_stat, i_all, 'DENLM', 'RHOVALNEW')

  end subroutine rhovalnew

end module mod_rhovalnew
