module mod_rhovalnew

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: RHOVALNEW
  ! > @note Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
  ! -------------------------------------------------------------------------------
  subroutine rhovalnew(ldorhoef, ielast, nsra, nspin, lmax, ez, wez, zat, socscale, cleb, icleb, iend, ifunm, lmsp, ncheb, npan_tot, npan_log, npan_eq, rmesh, irws, rpan_intervall, &
    ipan_intervall, rnew, vinsnew, thetasnew, theta, phi, i1, ipot, den_out, espv, rho2ns, r2nef, muorb, angles_new, idoldau, lopt, wldau, denmatn, natyp)

#ifdef CPP_OMP
    use :: omp_lib
#endif
#ifdef CPP_MPI
    use :: mpi
#endif
#ifdef CPP_TIMING
    use :: mod_timing
#endif

    use :: mod_types, only: t_tgmat, t_inc
#ifdef CPP_MPI
    use :: mod_types, only: gather_tmat, t_mpi_c_grid, save_t_mpi_c_grid, get_ntot_pt_ioff_pt_2d
#endif
    use :: mod_mympi, only: myrank, master
#ifdef CPP_MPI
    use :: mod_mympi, only: find_dims_2d, distribute_linear_on_tasks, mympi_main1c_comm_newsosol
#endif
    use :: mod_save_wavefun, only: t_wavefunctions, read_wavefunc
    use :: mod_version_info
    use :: global_variables
    use :: constants
    use :: mod_profiling
    use :: mod_datatypes, only: dp
    use :: mod_rhooutnew
    use :: mod_calcsph
    use :: mod_cheb2oldgrid
    use :: mod_rll_global_solutions
    use :: mod_intcheb_cell
    use :: mod_rllsllsourceterms
    use :: mod_rllsll
    use :: mod_spinorbit_ham
    use :: mod_sll_global_solutions
    use :: mod_rotatespinframe, only: rotatematrix, rotatevector
    use :: mod_vllmatsra
    use :: mod_vllmat

    implicit none

    integer, intent (in) :: i1
    integer, intent (in) :: nsra
    integer, intent (in) :: lmax   ! < Maximum l component in wave function expansion
    integer, intent (in) :: iend   ! < Number of nonzero gaunt coefficients
    integer, intent (in) :: ipot
    integer, intent (in) :: irws   ! < R point at WS radius for a given atom
    integer, intent (in) :: lopt   ! < angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, intent (in) :: natyp  ! < Number of kinds of atoms in unit cell
    integer, intent (in) :: nspin  ! < Counter for spin directions
    integer, intent (in) :: ncheb  ! < Number of Chebychev pannels for the new solver
    integer, intent (in) :: ielast
    integer, intent (in) :: idoldau ! < flag to perform LDA+U
    integer, intent (in) :: npan_eq ! < Number of intervals from [R_LOG] to muffin-tin radius Used in conjunction with runopt NEWSOSOL
    integer, intent (in) :: npan_tot
    integer, intent (in) :: npan_log ! < Number of intervals from nucleus to [R_LOG] Used in conjunction with runopt NEWSOSOL
    real (kind=dp), intent (in) :: zat ! < Nuclear charge for a given atom
    real (kind=dp), intent (in) :: socscale ! < Spin-orbit scaling for a given atom
    logical, intent (in) :: ldorhoef
    integer, dimension (lmxspd), intent (in) :: lmsp ! < 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
    integer, dimension (lmxspd), intent (in) :: ifunm
    integer, dimension (0:ntotd), intent (in) :: ipan_intervall
    integer, dimension (ncleb, 4), intent (in) :: icleb
    real (kind=dp), dimension (*), intent (in) :: cleb ! < GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irmd), intent (in) :: rmesh
    real (kind=dp), dimension (mmaxd, mmaxd, nspind), intent (in) :: wldau ! < potential matrix

    ! .. In/Out variables
    real (kind=dp), intent (inout) :: phi
    real (kind=dp), intent (inout) :: theta
    real (kind=dp), dimension (nrmaxd), intent (inout) :: rnew
    real (kind=dp), dimension (0:ntotd), intent (inout) :: rpan_intervall
    real (kind=dp), dimension (0:lmax+1, 3), intent (inout) :: muorb
    real (kind=dp), dimension (nrmaxd, nfund), intent (inout) :: thetasnew
    real (kind=dp), dimension (nrmaxd, lmpotd, nspotd), intent (inout) :: vinsnew ! < Non-spherical part of the potential
    complex (kind=dp), dimension (iemxd), intent (inout) :: ez
    complex (kind=dp), dimension (iemxd), intent (inout) :: wez
    ! .. Output variables
    real (kind=dp), dimension (2), intent (out) :: angles_new
    real (kind=dp), dimension (0:lmax+1, 2), intent (out) :: espv
    real (kind=dp), dimension (irmd, lmpotd, 4), intent (out) :: r2nef
    real (kind=dp), dimension (irmd, lmpotd, 4), intent (out) :: rho2ns
    complex (kind=dp), dimension (0:lmax+1, ielast, 2), intent (out) :: den_out

    ! .. Local variables
    integer :: lmsize
    integer :: lmaxd1
    integer :: ir, irec, use_sratrick, nvec, lm1, lm2, ie, irmdnew, imt1, jspin, idim, iorb, l1
    integer :: i_stat, i_all
    integer :: iq, nqdos           ! qdos ruess: number of qdos points
    integer :: lrecgflle, ierr     ! lmlm-dos
    integer :: lmlo, lmhi, is, js, mmax ! LDAU
    integer :: ix, m1              ! qdos ruess

    real (kind=dp) :: thetanew, phinew
    real (kind=dp) :: totmoment
    real (kind=dp) :: totxymoment
    complex (kind=dp) :: ek
    complex (kind=dp) :: df
    complex (kind=dp) :: eryd
    complex (kind=dp) :: temp1
    complex (kind=dp) :: dentemp
    complex (kind=dp) :: gmatprefactor
    integer, dimension (4) :: lmshift1
    integer, dimension (4) :: lmshift2
    integer, dimension (2*lmmaxso) :: jlk_index
    real (kind=dp), dimension (3) :: moment
    real (kind=dp), dimension (3) :: denorbmom
    real (kind=dp), dimension (3) :: denorbmomns
    real (kind=dp), dimension (2, 4) :: denorbmomsp
    real (kind=dp), dimension (0:lmax, 3) :: denorbmomlm
    complex (kind=dp), dimension (4) :: rho2
    complex (kind=dp), dimension (4) :: rho2int
    complex (kind=dp), dimension (2*(lmax+1)) :: alphasph
    complex (kind=dp), dimension (lmmaxso, lmmaxso) :: gmat0
    complex (kind=dp), dimension (lmmaxso, lmmaxso) :: gldau ! LDAU
    complex (kind=dp), dimension (lmmaxso, lmmaxso) :: tmatll
    complex (kind=dp), dimension (lmmaxso, lmmaxso) :: alphall ! LLY
    complex (kind=dp), dimension (lmmaxso, lmmaxso) :: tmattemp
    complex (kind=dp), dimension (2, 2) :: rho2ns_temp
    complex (kind=dp), dimension (lmmaxso, lmmaxso, iemxd) :: gmatll
    complex (kind=dp), dimension (mmaxd, mmaxd, 2, 2) :: denmatn ! LDAU

    ! .. Local allocatable arrays
    real (kind=dp), dimension (:, :), allocatable :: qvec ! qdos ruess: q-vectors for qdos
    real (kind=dp), dimension (:, :, :), allocatable :: vins
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
    ! ..
    logical :: test, opt
    external :: test, opt

    ! lmsize is original lm-size (without enhancement through soc etc.)
    lmsize = lmmaxd/(1+korbit)

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

    irmdnew = npan_tot*(ncheb+1)
    imt1 = ipan_intervall(npan_log+npan_eq) + 1
    allocate (vins(irmdnew,lmpotd,nspin), stat=i_stat)
    call memocc(i_stat, product(shape(vins))*kind(vins), 'VINS', 'RHOVALNEW')
    vins = 0d0
    do lm1 = 1, lmpotd
      do ir = 1, irmdnew
        vins(ir, lm1, 1) = vinsnew(ir, lm1, ipot)
        vins(ir, lm1, nspin) = vinsnew(ir, lm1, ipot+nspin-1)
      end do
    end do

    ! ! set up the non-spherical ll' matrix for potential VLL'
    if (nsra==2) then
      use_sratrick = 1
      if (test('nosph   ')) use_sratrick = 0
    else
      use_sratrick = 0
    end if
    allocate (vnspll0(lmmaxso,lmmaxso,irmdnew), stat=i_stat)
    call memocc(i_stat, product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'RHOVALNEW')
    vnspll0 = czero
    allocate (vnspll1(lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'RHOVALNEW')
    vnspll0 = czero

    call vllmat(1, nrmaxd, irmdnew, lmsize, lmmaxso, vnspll0, vins, lmpotd, cleb, icleb, iend, nspin, zat, rnew, use_sratrick, ncleb)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDAU
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (idoldau==1) then
      lmlo = lopt**2 + 1
      lmhi = (lopt+1)**2
      do ir = 1, irmdnew
        vnspll0(lmlo:lmhi, lmlo:lmhi, ir) = vnspll0(lmlo:lmhi, lmlo:lmhi, ir) + wldau(1:mmaxd, 1:mmaxd, 1)
      end do
      lmlo = lmlo + lmsize
      lmhi = lmhi + lmsize
      do ir = 1, irmdnew
        vnspll0(lmlo:lmhi, lmlo:lmhi, ir) = vnspll0(lmlo:lmhi, lmlo:lmhi, ir) + wldau(1:mmaxd, 1:mmaxd, 2)
      end do
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDAU
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! initial allocate
    if (nsra==2) then
      allocate (vnspll(2*lmmaxso,2*lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'RHOVALNEW')
    else
      allocate (vnspll(lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'RHOVALNEW')
    end if

    allocate (hlk(4*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(hlk))*kind(hlk), 'HLK', 'RHOVALNEW')
    allocate (jlk(4*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(jlk))*kind(jlk), 'JLK', 'RHOVALNEW')
    allocate (hlk2(4*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(hlk2))*kind(hlk2), 'HLK2', 'RHOVALNEW')
    allocate (jlk2(4*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(jlk2))*kind(jlk2), 'JLK2', 'RHOVALNEW')
    allocate (tmatsph(2*(lmax+1),0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(tmatsph))*kind(tmatsph), 'TMATSPH', 'RHOVALNEW')
    allocate (rll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(rll))*kind(rll), 'RLL', 'RHOVALNEW')
    allocate (sll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(sll))*kind(sll), 'SLL', 'RHOVALNEW')
    allocate (rllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'RHOVALNEW')
    allocate (sllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'RHOVALNEW')
    allocate (cden(irmdnew,0:lmax,4,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cden))*kind(cden), 'CDEN', 'RHOVALNEW')
    allocate (cdenlm(irmdnew,lmsize,4,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdenlm))*kind(cdenlm), 'CDENLM', 'RHOVALNEW')
    allocate (cdenns(irmdnew,4,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdenns))*kind(cdenns), 'CDENNS', 'RHOVALNEW')
    allocate (rho2nsc(irmdnew,lmpotd,4), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsc))*kind(rho2nsc), 'RHO2NSC', 'RHOVALNEW')
    allocate (rho2nsc_loop(irmdnew,lmpotd,4,ielast), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsc_loop))*kind(rho2nsc_loop), 'RHO2NSC_loop', 'RHOVALNEW')
    allocate (rho2nsnew(irmd,lmpotd,4), stat=i_stat)
    call memocc(i_stat, product(shape(rho2nsnew))*kind(rho2nsnew), 'RHO2NSNEW', 'RHOVALNEW')
    allocate (r2nefc(irmdnew,lmpotd,4), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefc))*kind(r2nefc), 'R2NEFC', 'RHOVALNEW')
    allocate (r2nefc_loop(irmdnew,lmpotd,4,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefc_loop))*kind(r2nefc_loop), 'R2NEFC_loop', 'RHOVALNEW')
    allocate (r2nefnew(irmd,lmpotd,4), stat=i_stat)
    call memocc(i_stat, product(shape(r2nefnew))*kind(r2nefnew), 'R2NEFNEW', 'RHOVALNEW')
    allocate (r2orbc(irmdnew,lmpotd,4,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(r2orbc))*kind(r2orbc), 'R2ORBC', 'RHOVALNEW')
    allocate (cdentemp(irmdnew,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(cdentemp))*kind(cdentemp), 'CDENTEMP', 'RHOVALNEW')
    allocate (gflle_part(lmmaxso,lmmaxso,0:nth-1), stat=i_stat)
    call memocc(i_stat, product(shape(gflle_part))*kind(gflle_part), 'GFLLE_PART', 'RHOVALNEW')
    allocate (gflle(lmmaxso,lmmaxso,ielast,1), stat=i_stat)
    call memocc(i_stat, product(shape(gflle))*kind(gflle), 'GFLLE', 'RHOVALNEW')
    allocate (den(0:lmaxd1,ielast,1,2), denlm(lmsize,ielast,1,2), stat=i_stat)
    call memocc(i_stat, product(shape(den))*kind(den), 'DEN', 'RHOVALNEW')
    rho2nsc = czero
    rho2nsc_loop = czero
    r2nefc = czero
    r2nefc_loop = czero
    r2orbc = czero
    rho2ns = 0.d0                  ! fivos 19.7.2014, this was CZERO
    r2nef = 0.d0                   ! fivos 19.7.2014, this was CZERO
    rho2nsnew = czero
    r2nefnew = czero
    den = czero
    denlm = czero
    espv = 0d0
    rho2int = czero
    denorbmom = 0d0
    denorbmomsp = 0d0
    denorbmomlm = 0d0
    denorbmomns = 0d0
    thetanew = 0d0
    phinew = 0d0
    gflle_part = czero
    gflle = czero
    gldau = czero
    ! LM shifts for correct density summation
    lmshift1(1) = 0                ! qdos ruess
    lmshift1(2) = lmsize           ! qdos ruess
    lmshift1(3) = 0                ! qdos ruess
    lmshift1(4) = lmsize           ! qdos ruess
    lmshift2(1) = 0                ! qdos ruess
    lmshift2(2) = lmsize           ! qdos ruess
    lmshift2(3) = lmsize           ! qdos ruess
    lmshift2(4) = 0                ! qdos ruess

    ! DO IR=1,3
    ! DO LM1=0,LMAXD1+1
    ! MUORB(LM1,IR)=0d0  !zimmer: initialization shifted to main1c
    ! ENDDO
    ! ENDDO

    nqdos = 1                      ! qdos ruess
    if (opt('qdos    ')) then      ! qdos ruess
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
      allocate (gflle(lmmaxso,lmmaxso,ielast,nqdos), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(gflle))*kind(gflle), 'GFLLE', 'RHOVALNEW') ! qdos ruess
      allocate (den(0:lmaxd1,ielast,nqdos,2), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(den))*kind(den), 'DEN', 'RHOVALNEW') ! qdos ruess
      allocate (denlm(lmsize,ielast,nqdos,2), stat=i_stat) ! qdos ruess
      call memocc(i_stat, product(shape(denlm))*kind(qvec), 'DENLM', 'RHOVALNEW') ! qdos ruess
100   if (ierr/=0) stop 'ERROR READING ''qvec.dat''' ! qdos ruess
    end if                         ! OPT('qdos    ')                                                     ! qdos ruess

#ifdef CPP_MPI
    i1_myrank = i1 - t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie) ! lmlm-dos ruess
#else
    i1_myrank = i1                 ! lmlm-dos ruess
#endif
    if ((opt('lmlm-dos')) .and. (i1_myrank==1)) then ! lmlm-dos ruess
      lrecgflle = 4*lmmaxso*lmmaxso*ielast*nqdos ! lmlm-dos ruess
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
    ! $omp shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)      &
    ! $omp shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)             &
    ! $omp shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)           &
    ! $omp shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)    &
    ! $omp shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)         &
    ! $omp shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)               &
    ! $omp shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)         &
    ! $omp shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)        &
    ! $omp shared(t_tgmat,ie_end, ie_start, t_wavefunctions)                 &
    ! $omp shared(LMMAXSO,lmsize,lmpotd,NRMAXD,NTOTD,LMAXD1)                  &
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
      df = wez(ie)/dble(nspin)
      if (nsra==2) then
        ek = sqrt(eryd+eryd*eryd/(cvlight*cvlight))*(1d0+eryd/(cvlight*cvlight))
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
        call read_wavefunc(t_wavefunctions, rll, rllleft, sll, sllleft, i1, ie, nsra, lmmaxso, irmdnew, ith, nth, rll_was_read_in, sll_was_read_in, rllleft_was_read_in, &
          sllleft_was_read_in)
      end if

      ! recalculate wavefuntions, also include left solution
      ! contruct the spin-orbit coupling hamiltonian and add to potential
      call spinorbit_ham(lmax, lmsize, vins, rnew, eryd, zat, cvlight, socscale, nspin, lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, irmdnew, nrmaxd, &
        vnspll0, vnspll1(:,:,:,ith), '1')

      ! extend matrix for the SRA treatment
      vnspll(:, :, :, ith) = czero
      if (nsra==2) then
        if (use_sratrick==0) then
          call vllmatsra(vnspll1(:,:,:,ith), vnspll(:,:,:,ith), rnew, lmmaxso, irmdnew, nrmaxd, eryd, lmax, 0, 'Ref=0')
        else if (use_sratrick==1) then
          call vllmatsra(vnspll1(:,:,:,ith), vnspll(:,:,:,ith), rnew, lmmaxso, irmdnew, nrmaxd, eryd, lmax, 0, 'Ref=Vsph')
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
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxso, 1, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)

        ! using spherical potential as reference
        if (use_sratrick==1) then
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin, zat, eryd, lmpotd, lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk(:,:,ith), jlk(:,:,ith), &
            hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, tmatsph(:,ith), alphasph, use_sratrick)
        end if

        ! calculate the tmat and wavefunctions
        rll(:, :, :, ith) = czero
        sll(:, :, :, ith) = czero

        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Right solutions
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        tmatll = czero
        ! faster calculation of RLL.
        ! no irregular solutions SLL are needed in self-consistent iterations
        ! because the density depends only on RLL, RLLLEFT and SLLLEFT
        if (opt('RLL-SLL ') .and. .not. (opt('XCPL    ') .or. opt('OPERATOR'))) then
          call rll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), tmatll, ncheb, npan_tot, lmmaxso, nvec*lmmaxso, 4*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', use_sratrick, alphall)
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), sll(:,:,:,ith), tmatll, ncheb, npan_tot, lmmaxso, nvec*lmmaxso, 4*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alphall)
        end if
        if (nsra==2) then
          rll(lmmaxso+1:nvec*lmmaxso, :, :, ith) = rll(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
          sll(lmmaxso+1:nvec*lmmaxso, :, :, ith) = sll(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
        end if

      end if                       ! read/recalc wavefunctions

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Left solutions
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((t_wavefunctions%nwfsavemax>0 .and. (.not. (rllleft_was_read_in .and. sllleft_was_read_in))) .or. (t_wavefunctions%nwfsavemax==0)) then
        ! read/recalc wavefunctions left contruct the TRANSPOSE spin-orbit coupling hamiltonian and add to potential
        call spinorbit_ham(lmax, lmsize, vins, rnew, eryd, zat, cvlight, socscale, nspin, lmpotd, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, irmdnew, nrmaxd, &
          vnspll0, vnspll1(:,:,:,ith), 'transpose')
        ! extend matrix for the SRA treatment
        vnspll(:, :, :, ith) = czero
        if (nsra==2) then
          if (use_sratrick==0) then
            call vllmatsra(vnspll1(:,:,:,ith), vnspll(:,:,:,ith), rnew, lmmaxso, irmdnew, nrmaxd, eryd, lmax, 0, 'Ref=0')
          else if (use_sratrick==1) then
            call vllmatsra(vnspll1(:,:,:,ith), vnspll(:,:,:,ith), rnew, lmmaxso, irmdnew, nrmaxd, eryd, lmax, 0, 'Ref=Vsph')
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
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxso, 1, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)

        ! using spherical potential as reference
        ! notice that exchange the order of left and right hankel/bessel functions
        if (use_sratrick==1) then
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin, zat, eryd, lmpotd, lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), &
            hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, alphasph, tmatsph(:,ith), use_sratrick)
        end if

        ! calculate the tmat and wavefunctions
        rllleft(:, :, :, ith) = czero
        sllleft(:, :, :, ith) = czero

        ! left solutions
        ! notice that exchange the order of left and right hankel/bessel functions
        tmattemp = czero
        ! faster calculation of RLLLEFT and SLLLEFT.
        if (opt('RLL-SLL ') .and. .not. (opt('XCPL    ') .or. opt('OPERATOR'))) then
          call rll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), rllleft(:,:,:,ith), tmattemp, ncheb, npan_tot, lmmaxso, nvec*lmmaxso, 4*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', use_sratrick, alphall)
          call sll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), sllleft(:,:,:,ith), ncheb, npan_tot, lmmaxso, nvec*lmmaxso, 4*(lmax+1), irmdnew, nsra, jlk_index, &
            hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', use_sratrick)
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rllleft(:,:,:,ith), sllleft(:,:,:,ith), tmattemp, ncheb, npan_tot, lmmaxso, nvec*lmmaxso, 4*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alphall)
        end if
        if (nsra==2) then
          rllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith) = rllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
          sllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith) = sllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
        end if
      end if                       ! read/recalc wavefunctions left

      do iq = 1, nqdos             ! qdos
        ! read in GF
        irec = iq + nqdos*(ie-1) + nqdos*ielast*(i1-1) ! qdos
#ifdef CPP_OMP
        ! $omp critical
#endif
        if (t_tgmat%gmat_to_file) then
          read (69, rec=irec) gmat0
        else
          irec = iq + nqdos*(ie_num-1) + nqdos*ie_end*(i1-1)
          gmat0(:, :) = t_tgmat%gmat(:, :, irec)
        end if
#ifdef CPP_OMP
        ! $omp end critical
#endif

        ! rotate gmat from global frame to local frame
        call rotatematrix(gmat0, theta, phi, lmsize, 1)

        do lm1 = 1, lmmaxso
          do lm2 = 1, lmmaxso
            gmatll(lm1, lm2, ie) = gmat0(lm1, lm2)
          end do
        end do
        ! calculate density
        call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, df, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), commented out since sll is not used in rhooutnew
          rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), rho2nsc_loop(:,:,:,ie), 0, gflle(:,:,ie,iq), rpan_intervall, ipan_intervall)

        do jspin = 1, 4
          do lm1 = 0, lmax
            cdentemp(:, ith) = czero
            dentemp = czero
            do ir = 1, irmdnew
              cdentemp(ir, ith) = cden(ir, lm1, jspin, ith)
            end do
            call intcheb_cell(cdentemp(:,ith), dentemp, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
            rho2(jspin) = dentemp
            rho2int(jspin) = rho2int(jspin) + rho2(jspin)*df
            if (jspin<=2) then
              den(lm1, ie, iq, jspin) = rho2(jspin)
            end if
          end do

          if (jspin<=2) then
            do lm1 = 1, lmsize
              cdentemp(:, ith) = czero
              dentemp = czero
              do ir = 1, irmdnew
                cdentemp(ir, ith) = cdenlm(ir, lm1, jspin, ith)
              end do
              call intcheb_cell(cdentemp(:,ith), dentemp, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
              denlm(lm1, ie, iq, jspin) = dentemp
            end do
            cdentemp(:, ith) = czero
            dentemp = czero
            do ir = 1, irmdnew
              cdentemp(ir, ith) = cdenns(ir, jspin, ith)
            end do
            call intcheb_cell(cdentemp(:,ith), dentemp, rpan_intervall, ipan_intervall, npan_tot, ncheb, irmdnew)
            den(lmaxd1, ie, iq, jspin) = dentemp
            rho2int(jspin) = rho2int(jspin) + den(lmaxd1, ie, iq, jspin)*df
          end if
        end do                     ! JSPIN

        do jspin = 1, 4
          if (jspin<=2) then
            do lm1 = 0, lmaxd1
              espv(lm1, jspin) = espv(lm1, jspin) + aimag(eryd*den(lm1,ie,iq,jspin)*df)
            end do
          end if
        end do
      end do                       ! IQ = 1,NQDOS

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get charge at the Fermi energy (IELAST)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (ie==ielast .and. ldorhoef) then
        call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, cone, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
          rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), r2nefc_loop(:,:,:,ith), 0, gflle_part(:,:,ith), rpan_intervall, &
          ipan_intervall)
      end if

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Get orbital moment
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do iorb = 1, 3
        call rhooutnew(nsra, lmax, gmatll(1,1,ie), ek, lmpotd, cone, npan_tot, ncheb, cleb, icleb, iend, irmdnew, thetasnew, ifunm, imt1, lmsp, rll(:,:,:,ith), & ! SLL(:,:,:,ith), ! commented out since sll is not used in rhooutnew
          rllleft(:,:,:,ith), sllleft(:,:,:,ith), cden(:,:,:,ith), cdenlm(:,:,:,ith), cdenns(:,:,ith), r2orbc(:,:,:,ith), iorb, gflle_part(:,:,ith), rpan_intervall, ipan_intervall)
        do jspin = 1, 4
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
            end do
          end if
        end do
      end do                       ! IORB
    end do                         ! IE loop
#ifdef CPP_OMP
    ! $omp end parallel do
#endif

    ! omp: move sum from rhooutnew here after parallel calculation
    do ir = 1, irmdnew
      do lm1 = 1, lmpotd
        do jspin = 1, 4
          do ie = 1, ielast
            rho2nsc(ir, lm1, jspin) = rho2nsc(ir, lm1, jspin) + rho2nsc_loop(ir, lm1, jspin, ie)
          end do
        end do
      end do
    end do
    ! omp: don't forget to do the same with density at fermi energy:
    do ith = 0, nth - 1
      r2nefc(:, :, :) = r2nefc(:, :, :) + r2nefc_loop(:, :, :, ith)
    end do

#ifdef CPP_MPI
    if (opt('qdos    ')) then      ! qdos
      ! first communicate den array to write out qdos files                      ! qdos
      idim = (lmaxd1+1)*ielast*2*nqdos ! qdos
      allocate (workc(0:lmaxd1,ielast,2,nqdos), stat=i_stat) ! qdos
      call memocc(i_stat, product(shape(workc))*kind(workc), 'workc', 'RHOVALNEW') ! qdos
      workc = czero                ! qdos
      call mpi_reduce(den, workc, idim, mpi_double_complex, mpi_sum, master, & ! qdos
        t_mpi_c_grid%mympi_comm_at, ierr) ! qdos
      call zcopy(idim, workc, 1, den, 1) ! qdos
      i_all = -product(shape(workc))*kind(workc) ! qdos
      deallocate (workc, stat=i_stat) ! qdos
      call memocc(i_stat, i_all, 'workc', 'RHOVALNEW') ! qdos
      ! ! qdos
      if (t_mpi_c_grid%myrank_at==master) then ! qdos
        ie_start = 0               ! qdos
        ie_end = ielast            ! qdos
        do ie_num = 1, ie_end      ! qdos
          ie = ie_start + ie_num   ! qdos
          do iq = 1, nqdos         ! qdos
            if ((iq==1) .and. (ie_num==1)) then ! qdos
              if (natyp>=100) then ! qdos
                open (31, &        ! qdos
                  file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10)) & ! qdos
                  //char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! qdos
                open (32, &        ! qdos
                  file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10)) & ! qdos
                  //char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! qdos
              else                 ! qdos
                open (31, &        ! qdos
                  file='qdos.'//char(48+mod(i1/10,10))// & ! qdos
                  char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! qdos
                open (32, &        ! qdos
                  file='qdos.'//char(48+mod(i1/10,10))// & ! qdos
                  char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! qdos
              end if               ! qdos
              call version_print_header(31) ! qdos
              write (31, *) ' '    ! qdos
              write (31, 150) '# ISPIN=', 1, ' I1=', i1 ! qdos
              write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', 'k_x', 'k_y', & ! qdos
                'k_z', 'DEN_tot', 'DEN_s,p,...' ! qdos
              if (nspin>1) then    ! qdos
                call version_print_header(32) ! qdos
                write (32, *) ' '  ! qdos
                write (32, 150) '# ISPIN=', 2, ' I1=', i1 ! qdos
                write (32, '(7(A,3X))') '#   Re(E)', 'Im(E)', & ! qdos
                  'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! qdos
              end if               ! qdos
            end if                 ! IQ.EQ.1                                                 ! qdos
            do jspin = 1, 2        ! qdos
              dentot(jspin) = cmplx(0.d0, 0.d0, kind=dp) ! qdos
              do l1 = 0, lmaxd1    ! qdos
                dentot(jspin) = dentot(jspin) + den(l1, ie, iq, jspin) ! qdos
              end do               ! qdos
            end do                 ! qdos
            ! write qdos.nn.s.dat                                          ! qdos
            write (31, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! qdos
              -aimag(dentot(1))/pi, (-aimag(den(l1,ie,iq,1))/pi, l1=0, lmaxd1) ! qdos
            write (32, 120) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! qdos
              -aimag(dentot(2))/pi, (-aimag(den(l1,ie,iq,2))/pi, l1=0, lmaxd1) ! qdos
120         format (5f10.6, 40e16.8) ! qdos
            ! ! qdos
            if (test('compqdos')) then ! complex qdos
              if ((iq==1) .and. (ie_num==1)) then ! complex qdos
                if (natyp>=100) then ! complex qdos
                  open (31, &      ! complex qdos
                    file='cqdos.'//char(48+i1/100)// & ! complex qdos
                    char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.' & ! complex qdos
                    //char(48+1)//'.dat') ! complex qdos
                  open (32, &      ! complex qdos
                    file='cqdos.'//char(48+i1/100)// & ! complex qdos
                    char(48+mod(i1/10,10))//char(48+mod(i1,10))//'.' & ! complex qdos
                    //char(48+2)//'.dat') ! complex qdos
                else               ! complex qdos
                  open (31, file='cqdos.'//char(48+mod(i1/10,10))// & ! complex qdos
                    char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! complex qdos
                  open (32, file='cqdos.'//char(48+mod(i1/10,10))// & ! complex qdos
                    char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! complex qdos
                end if             ! complex qdos
                call version_print_header(31) ! complex qdos
                write (31, *) ' '  ! complex qdos
                write (31, '(A)') '#   lmax, natyp, nspin, nqdos, ielast:' ! complex qdos
                write (31, '(5I9)') lmax, natyp, nspin, nqdos, ielast ! complex qdos
                write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)', & ! complex qdos
                  'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! complex qdos
                if (nspin>1) then  ! complex qdos
                  call version_print_header(32) ! complex qdos
                  write (32, *) ' ' ! complex qdos
                  write (32, '(A)') '# lmax, natyp, nspin, nqdos, ielast:' ! complex qdos
                  write (32, '(5I9)') lmax, natyp, nspin, nqdos, ielast ! complex qdos
                  write (32, '(7(A,3X))') '#   Re(E)', 'Im(E)', & ! complex qdos
                    'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! complex qdos
                end if             ! complex qdos
              end if               ! IQ.EQ.1                                              ! complex qdos
              do jspin = 1, 2      ! complex qdos
                dentot(jspin) = cmplx(0.d0, 0.d0, kind=dp) ! complex qdos
                do l1 = 0, lmaxd1  ! complex qdos
                  dentot(jspin) = dentot(jspin) + den(l1, ie, iq, jspin) ! complex qdos
                end do             ! complex qdos
              end do               ! complex qdos
              ! write qdos.nn.s.dat                                       ! complex qdos
              write (31, 130) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! complex qdos
                dentot(1), (den(l1,ie,iq,1), l1=0, lmaxd1) ! complex qdos
              write (32, 130) ez(ie), qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! complex qdos
                dentot(2), (den(l1,ie,iq,2), l1=0, lmaxd1) ! complex qdos
130           format (6f10.6, 80e16.8) ! complex qdos
            end if                 ! complex qdos
            ! qdos
          end do                   ! IQ                                                            ! qdos
        end do                     ! IE                                                               ! qdos
      end if                       ! myrank_at==master                                                  ! qdos
    end if                         ! OPT('qdos    ')                                                  ! qdos
#endif

#ifdef CPP_MPI
    ! do communication only when compiled with MPI
#ifdef CPP_TIMING
    call timing_start('main1c - communication')
#endif
    ! reset NQDOS to avoid endless communication
    if (.not. opt('lmdos    ')) then
      nqdos = 1
    else
      if (myrank==master) write (*, *) 'lmlm-dos option, communcation might take a while!', ielast, nqdos
    end if
    ! set these arrays to zero to avoid double counting in cases where extra ranks are used
    if (t_mpi_c_grid%myrank_ie>(t_mpi_c_grid%dims(1)-1)) then
      den = czero
      denlm = czero
      gflle = czero
      r2nefc = czero
      rho2nsc = czero
      rho2int = czero
      muorb = 0.0d0
      espv = 0.0d0
      denorbmom = 0.0d0
      denorbmomsp = 0.0d0
      denorbmomlm = 0.0d0
      denorbmomns = 0.0d0
    end if
    call mympi_main1c_comm_newsosol(irmdnew, lmpotd, lmax, lmaxd1, lmsize, lmmaxso, ielast, nqdos, den, denlm, gflle, rho2nsc, r2nefc, rho2int, espv, muorb, denorbmom, denorbmomsp, &
      denorbmomlm, denorbmomns, t_mpi_c_grid%mympi_comm_at)
#ifdef CPP_TIMING
    call timing_pause('main1c - communication')
#endif

    ! MPI: do these writeout/data collection steps only on master and broadcast important results afterwards
    if (t_mpi_c_grid%myrank_at==master) then
#endif
      ! CPP_MPI
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDAU
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (idoldau==1) then
        ! calculate WLDAU
        do ie = 1, ielast
          do lm1 = 1, lmmaxso
            do lm2 = 1, lmmaxso
              gldau(lm1, lm2) = gldau(lm1, lm2) + gflle(lm1, lm2, ie, 1)*wez(ie)/dble(nspin)
            end do
          end do
        end do
        ! calculate occupation matrix
        mmax = 2*lopt + 1
        do is = 1, 2
          do js = 1, 2
            lmlo = lopt**2 + 1 + (is-1)*lmsize
            lmhi = (lopt+1)**2 + (js-1)*lmsize
            lm2 = lopt**2 + 1 + (js-1)*lmsize
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

      if (.not. opt('qdos    ')) then
        ! omp: moved write-out of dos files out of parallel energy loop
        ! Write out lm-dos:                                                     ! lm-dos
        if (opt('lmdos    ')) then ! qdos ruess
          do ie = 1, ielast        ! lm-dos
            iq = 1                 ! lm-dos
            if (ie==1) then        ! lm-dos
              if (natyp>=100) then ! lm-dos
                open (29, &        ! lm-dos
                  file='lmdos.'//char(48+i1/100)// & ! lm-dos
                  char(48+mod(i1/10,10))// & ! lm-dos
                  char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! lm-dos
                open (30, &        ! lm-dos
                  file='lmdos.'//char(48+i1/100)// & ! lm-dos
                  char(48+mod(i1/10,10))// & ! lm-dos
                  char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! lm-dos
              else                 ! lm-dos
                open (29, &        ! lm-dos
                  file='lmdos.'//char(48+i1/10)// & ! lm-dos
                  char(48+mod(i1,10))//'.'//char(48+1)//'.dat') ! lm-dos
                open (30, file='lmdos.'//char(48+i1/10)// & ! lm-dos
                  char(48+mod(i1,10))//'.'//char(48+2)//'.dat') ! lm-dos
              end if               ! lm-dos
              call version_print_header(29) ! lm-dos
              write (29, *) ' '    ! lm-dos
              write (29, 150) '# ISPIN=', 1, ' I1=', i1 ! lm-dos
              call version_print_header(30) ! lm-dos
              write (30, *) ' '    ! lm-dos
              write (30, 150) '# ISPIN=', 2, ' I1=', i1 ! lm-dos
            end if                 ! IE==1                                                      ! lm-dos
            write (29, 140) ez(ie), (-aimag(denlm(l1,ie,iq,1))/pi, l1=1, lmsize) ! lm-dos
            write (30, 140) ez(ie), (-aimag(denlm(l1,ie,iq,2))/pi, l1=1, lmsize) ! lm-dos
140         format (30e12.4)       ! lm-dos
150         format (a8, i3, a4, i5) ! lm-dos/qdos ruess
          end do                   ! IE
        end if
      end if                       ! .not. OPT('qdos    ')

      ! write gflle to file                                                 ! lmlm-dos
      if (opt('lmlm-dos')) then    ! lmlm-dos
        if (t_inc%i_write>0) then  ! lmlm-dos
          write (1337, *) 'gflle:', shape(gflle), shape(gflle_part), lrecgflle ! lmlm-dos
        end if                     ! lmlm-dos
        write (91, rec=i1) gflle   ! lmlm-dos
      end if                       ! lmlm-dos

      allocate (rhotemp(irmdnew,lmpotd), stat=i_stat)
      call memocc(i_stat, product(shape(rhotemp))*kind(rhotemp), 'RHOTEMP', 'RHOVALNEW')
      allocate (rhonewtemp(irws,lmpotd), stat=i_stat)
      call memocc(i_stat, product(shape(rhonewtemp))*kind(rhonewtemp), 'RHONEWTEMP', 'RHOVALNEW')

      do jspin = 1, 4
        rhotemp = czero
        rhonewtemp = czero
        do lm1 = 1, lmpotd
          do ir = 1, irmdnew
            rhotemp(ir, lm1) = rho2nsc(ir, lm1, jspin)
          end do
        end do
        call cheb2oldgrid(irws, irmdnew, lmpotd, rmesh, ncheb, npan_tot, rpan_intervall, ipan_intervall, rhotemp, rhonewtemp, irmd)
        do lm1 = 1, lmpotd
          do ir = 1, irws
            rho2nsnew(ir, lm1, jspin) = rhonewtemp(ir, lm1)
          end do
        end do

        rhotemp = czero
        rhonewtemp = czero
        do lm1 = 1, lmpotd
          do ir = 1, irmdnew
            rhotemp(ir, lm1) = r2nefc(ir, lm1, jspin)
          end do
        end do
        call cheb2oldgrid(irws, irmdnew, lmpotd, rmesh, ncheb, npan_tot, rpan_intervall, ipan_intervall, rhotemp, rhonewtemp, irmd)
        do lm1 = 1, lmpotd
          do ir = 1, irws
            r2nefnew(ir, lm1, jspin) = rhonewtemp(ir, lm1)
          end do
        end do
      end do
      i_all = -product(shape(rhotemp))*kind(rhotemp)
      deallocate (rhotemp, stat=i_stat)
      call memocc(i_stat, i_all, 'RHOTEMP', 'RHOVALNEW')
      i_all = -product(shape(rhonewtemp))*kind(rhonewtemp)
      deallocate (rhonewtemp, stat=i_stat)
      call memocc(i_stat, i_all, 'RHONEWTEMP', 'RHOVALNEW')
      ! calculate new THETA and PHI for non-colinear
      if (.not. test('FIXMOM  ')) then
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

        if (abs(totxymoment)>1d-05) then
          if (abs(moment(3))<1d-05) then
            thetanew = pi/2d0
          else
            thetanew = acos(moment(3)/totmoment)
          end if
          if (totxymoment<1d-05) then
            phinew = 0d0
          else
            phinew = atan2(moment(2), moment(1))
          end if
        end if

        if (t_inc%i_write>0) then
          write (1337, *) 'moment', myrank, moment(1), moment(2), moment(3)
          write (1337, *) thetanew/(2.0d0*pi)*360.0d0, phinew/(2.0d0*pi)*360.0d0
        end if
        ! only on master different from zero:
        angles_new(1) = thetanew
        angles_new(2) = phinew
        call rotatevector(rho2nsnew, rho2ns, irws, lmpotd, thetanew, phinew, theta, phi, irmd)
        call rotatevector(r2nefnew, r2nef, irws, lmpotd, thetanew, phinew, theta, phi, irmd)
      else
        rho2ns(:, :, :) = aimag(rho2nsnew(:,:,:))
        r2nef(:, :, :) = aimag(r2nefnew(:,:,:))
      end if

      idim = irmd*lmpotd
      call dscal(idim, 2.d0, rho2ns(1,1,1), 1)
      call daxpy(idim, -0.5d0, rho2ns(1,1,1), 1, rho2ns(1,1,2), 1)
      call daxpy(idim, 1.0d0, rho2ns(1,1,2), 1, rho2ns(1,1,1), 1)
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Do the same at the Fermi energy
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call dscal(idim, 2.d0, r2nef(1,1,1), 1)
      call daxpy(idim, -0.5d0, r2nef(1,1,1), 1, r2nef(1,1,2), 1)
      call daxpy(idim, 1.0d0, r2nef(1,1,2), 1, r2nef(1,1,1), 1)

      do lm1 = 0, lmaxd1
        do ie = 1, ielast
          do jspin = 1, nspin
            den_out(lm1, ie, jspin) = den(lm1, ie, 1, jspin)
          end do
        end do
      end do

#ifdef CPP_MPI
    end if                         ! (myrank==master)

    ! communicate den_out to all processors with the same atom number
    idim = (lmax+2)*ielast*2
    call mpi_bcast(den_out, idim, mpi_double_complex, master, t_mpi_c_grid%mympi_comm_at, ierr)
    if (ierr/=mpi_success) stop 'error bcast den_out in rhovalnew'
    idim = 2
    call mpi_bcast(angles_new, idim, mpi_double_precision, master, t_mpi_c_grid%mympi_comm_at, ierr)
    if (ierr/=mpi_success) stop 'error bcast angles_new in rhovalnew'
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
