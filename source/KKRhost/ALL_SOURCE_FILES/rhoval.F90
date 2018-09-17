module mod_rhoval

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: RHOVAL

  ! > @note - Average WLDAU for spherical wavefunctions:
  ! > The spherical part of the d or f wavefunction is found by adding
  ! > the average interaction potential WLDAUAV to the spherical
  ! > potential. Then the non-spherical parts are found by using only
  ! > the deviation of WLDAU from the average. This speeds up the
  ! > convergence of the Born series. See also subroutines
  ! > regsol, pnstmat and pnsqns
  ! >
  ! > @note -LDA+U implementation Mar. 2002-Dec.2004 Ph. Mavropoulos, H. Ebert, V. Popescu
  ! > @note -Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
  ! -------------------------------------------------------------------------------
  subroutine rhoval(ihost, ldorhoef, icst, ins, ielast, nsra, ispin, nspin, nspinpot, i1, ez, wez, drdi, r, vins, visp, zat, ipan, ircut, irmin, thetas, ifunm, lmsp, rho2ns, r2nef, &
    rhoorb, den, denlm, muorb, espv, cleb, loflm, icleb, iend, jend, solver, soctl, ctl, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, irshift, itermvdir, qmtet, qmphi, &
    mvevil, mvevilef, nmvecmax, idoldau, lopt, phildau, wldau, denmatc, natyp, nqdos, lmax)

#ifdef CPP_MPI
    use :: mpi
#endif
    use :: mod_mympi, only: myrank, master
#ifdef CPP_MPI
    use :: mod_mympi, only: distribute_linear_on_tasks
#endif
    use :: mod_types, only: t_tgmat, t_inc, t_mpi_c_grid, init_tgmat
    use :: constants
    use :: mod_profiling
    use :: mod_version_info
    use :: global_variables
    use :: mod_datatypes
    use :: mod_pnsqns
    use :: mod_cradwf
    use :: mod_rholm
    use :: mod_rhons
    use :: mod_wfmesh
    use :: mod_drvrho

    implicit none

    ! .. Input variables
    integer, intent (in) :: i1
    integer, intent (in) :: ins
    integer, intent (in) :: lmax   ! < Maximum l component in wave function expansion
    integer, intent (in) :: iend   ! < Number of nonzero gaunt coefficients
    integer, intent (in) :: ipan   ! < Number of panels in non-MT-region
    integer, intent (in) :: icst   ! < Number of Born approximation
    integer, intent (in) :: nsra
    integer, intent (in) :: zrel   ! < atomic number (cast integer)
    integer, intent (in) :: lopt   ! < angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, intent (in) :: ispin
    integer, intent (in) :: nspin  ! < Counter for spin directions
    integer, intent (in) :: natyp  ! < Number of kinds of atoms in unit cell
    integer, intent (in) :: ihost
    integer, intent (in) :: irmin  ! < Max R for spherical treatment
    integer, intent (in) :: ielast
    integer, intent (in) :: jwsrel ! < index of the WS radius
    integer, intent (in) :: irshift ! < shift of the REL radial mesh with respect no NREL
    integer, intent (in) :: idoldau ! < flag to perform LDA+U
    integer, intent (in) :: nspinpot
    integer, intent (in) :: nmvecmax
    integer, intent (inout) :: nqdos
    real (kind=dp), intent (in) :: zat ! < Nuclear charge
    logical, intent (in) :: ldorhoef
    character (len=10), intent (in) :: solver
    real (kind=dp), dimension (irmd), intent (in) :: r
    real (kind=dp), dimension (krel*lmax+1), intent (in) :: ctl
    real (kind=dp), dimension (irmd), intent (in) :: drdi ! < Derivative dr/di
    real (kind=dp), dimension (irmd), intent (in) :: visp ! < Spherical part of the potential
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (in) :: vtrel ! < potential (spherical part)
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (in) :: btrel ! < magnetic field
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (in) :: rmrel ! < radial mesh
    real (kind=dp), dimension (krel*lmax+1), intent (in) :: soctl
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (in) :: drdirel ! < derivative of radial mesh
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (in) :: r2drdirel ! < \f$ r^2 \frac{\partial}{\partial \mathbf{r}}\frac{\partial}{\partial i}\f$ (r**2 * drdi)
    real (kind=dp), dimension (irmind:irmd, lmpotd), intent (in) :: vins ! < Non-spherical part of the potential
    real (kind=dp), dimension (ncleb, 2), intent (in) :: cleb ! < GAUNT coefficients (GAUNT)
    real (kind=dp), dimension (irid, nfund), intent (in) :: thetas ! < shape function THETA=0 outer space THETA =1 inside WS cell in spherical harmonics expansion
    complex (kind=dp), dimension (iemxd), intent (in) :: ez
    complex (kind=dp), dimension (iemxd), intent (in) :: wez
    complex (kind=dp), dimension (irmd), intent (in) :: phildau
    complex (kind=dp), dimension (0:lmax+1, ielast*(1+krel), nqdos), intent (in) :: den
    complex (kind=dp), dimension (lmmaxd, ielast*(1+krel), nqdos), intent (in) :: denlm
    ! .. In/Out variables
    real (kind=dp), dimension (mmaxd, mmaxd, nspind), intent (inout) :: wldau ! < potential matrix
    ! ---------------------------------------------------------------------------
    ! IHOST = 1   < -- this routine is called by the HOST tbkkr-program
    ! IHOST <> 1  < --                 called by the IMPURITY program
    ! ---------------------------------------------------------------------------
    ! .. Output variables
    real (kind=dp), dimension (irmd*krel+(1-krel)), intent (out) :: rhoorb
    real (kind=dp), dimension (0:lmax+1+1, 3), intent (out) :: muorb ! < orbital magnetic moment
    real (kind=dp), dimension (0:lmax+1, 2), intent (out) :: espv ! < changed for REL case
    real (kind=dp), dimension (irmd, lmpotd, 2), intent (out) :: r2nef ! < rho at FERMI energy
    real (kind=dp), dimension (irmd, lmpotd, 2), intent (out) :: rho2ns ! < radial density
    complex (kind=dp), dimension (mmaxd, mmaxd), intent (out) :: denmatc
    ! ----------------------------------------------------------------------------
    ! ITERMDIR variables
    ! ----------------------------------------------------------------------------
    logical, intent (in) :: itermvdir
    real (kind=dp), intent (in) :: qmtet ! < \f$ \theta\f$ angle of the agnetization with respect to the z-axis
    real (kind=dp), intent (in) :: qmphi ! < \f$ \phi\f$ angle of the agnetization with respect to the z-axis
    complex (kind=dp), dimension (0:lmax, 3, nmvecmax), intent (out) :: mvevil ! OUTPUT
    complex (kind=dp), dimension (0:lmax, 3, nmvecmax), intent (out) :: mvevilef ! OUTPUT
    ! ----------------------------------------------------------------------------
    ! ITERMDIR variables
    ! ----------------------------------------------------------------------------
    integer, dimension (lmxspd), intent (in) :: lmsp
    integer, dimension (lmxspd), intent (in) :: ifunm
    integer, dimension (0:ipand), intent (in) :: ircut ! < R points of panel borders
    integer, dimension (lm2d), intent (in) :: loflm ! < l of lm=(l,m) (GAUNT)
    integer, dimension (ncleb, 4), intent (in) :: icleb ! < Pointer array
    integer, dimension (lmpotd, 0:lmax, 0:lmax), intent (in) :: jend ! < Pointer array for icleb()
    ! .. Local Scalars
    ! .. Parameters
    integer :: lmaxd1
    integer :: i_stat, i_all
    real (kind=dp) :: wldauav
    complex (kind=dp) :: df, eryd, ek
#ifndef CPP_MPI
    complex (kind=dp) :: dentot    ! qdos
#endif
    integer :: idim, ie, ir, l, lm1, lm2, lmhi, lmlo, irec, ispinpot, lastez, m1, mmax
    integer :: iq                  ! NQDOS ! qdos number of qdos points
    integer :: ix                  ! qdos
    integer :: lrecgflle           ! lmlm-dos
    integer, dimension (4) :: lmshift1 ! lmlm-dos
    integer, dimension (4) :: lmshift2 ! lmlm-dos

    ! .. Local Arrays
    real (kind=dp), dimension (0:lmax) :: s
    real (kind=dp), dimension (irmd) :: cutoff
    real (kind=dp), dimension (irmd, 0:lmax) :: rs
    complex (kind=dp), dimension (0:lmax) :: ekl
    complex (kind=dp), dimension (0:lmax) :: tmat
    complex (kind=dp), dimension (0:lmax) :: alpha
    complex (kind=dp), dimension (0:lmax+1) :: dendum
    complex (kind=dp), dimension (lmmaxd) :: dum_denlm
    complex (kind=dp), dimension (irmd, 0:lmax) :: qz
    complex (kind=dp), dimension (irmd, 0:lmax) :: sz
    complex (kind=dp), dimension (irmd, 0:lmax) :: pz
    complex (kind=dp), dimension (irmd, 0:lmax) :: fz
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: ar
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: cr
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: dr
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gmat0
    complex (kind=dp), dimension (lmmaxd, lmmaxd, ielast) :: gmatll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd, 2) :: pns
    complex (kind=dp), dimension (lmmaxd, lmmaxd, irmind:irmd, 2) :: qns
    ! .. first 2 indices in dmuorb are the spin-resolved contributions,
    ! .. the 3rd one should be the sum of them
    complex (kind=dp), dimension (0:krel*lmax+(1-krel), 3) :: dmuorb
    ! .. Local allocatable arrays
    complex (kind=dp), dimension (:, :), allocatable :: qvec ! < qdos, q-vectors for qdos
    complex (kind=dp), dimension (:, :), allocatable :: gldau
    complex (kind=dp), dimension (:, :), allocatable :: dum_gflle ! lmlm-dos
    complex (kind=dp), dimension (:, :, :, :), allocatable :: gflle ! qdos

    ! This routine needs irregular wavefunctions
    logical :: lirrsol
    parameter (lirrsol=.true.)

#ifdef CPP_MPI
    integer :: ie_start
#endif
    integer :: ie_end, ie_num
    ! .. External Functions
    logical :: opt, test           ! qdos
    external :: opt, test          ! qdos

    lmaxd1 = lmax + 1

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDAU
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (idoldau==1) then
      wldauav = 0.d0
      lmlo = lopt*lopt + 1
      lmhi = (lopt+1)*(lopt+1)
      mmax = lmhi - lmlo + 1
      do m1 = 1, mmax
        wldauav = wldauav + wldau(m1, m1, ispin)
      end do
      wldauav = wldauav/dble(mmax)
      ! -------------------------------------------------------------------------
      ! Note: Application if WLDAU makes the potential discontinuous.
      ! A cutoff can be used if necessary to make the potential continuous
      ! for example (array bounds should be adjusted):

      ! CUTOFF(IR) = ( 1.D0 + DEXP( 20.D0*(R(IR)-R(349)) ) ) * &
      ! ( 1.D0 + DEXP( 20.D0*(R(276)-R(IR)) ) )
      ! CUTOFF(IR) = 1D0/CUTOFF(IR)
      ! -------------------------------------------------------------------------
      do m1 = 1, irmd
        cutoff(m1) = 1.d0
      end do
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDAU
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise variables
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (krel==0) then
      do lm1 = 1, lmpotd
        do ir = 1, irmd
          rho2ns(ir, lm1, ispin) = 0.0d0
        end do
      end do

      do l = 0, lmaxd1
        espv(l, ispin) = 0.0d0
      end do
      ! -------------------------------------------------------------------------
    else
      ! -------------------------------------------------------------------------
      do ispinpot = 1, 2
        do lm1 = 1, lmpotd
          do ir = 1, irmd
            rho2ns(ir, lm1, ispinpot) = 0.0d0
            r2nef(ir, lm1, ispinpot) = 0.0d0
          end do
        end do

        do l = 0, lmaxd1
          espv(l, ispinpot) = 0.0d0
        end do

      end do

      do ir = 1, irmd
        rhoorb(ir) = 0.0d0
      end do

      do ir = 1, 3
        do l = 0, lmax
          dmuorb(l, ir) = czero
        end do
      end do
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ITERMDIR
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (itermvdir) then
        do lm1 = 1, 3
          do lm2 = 1, nmvecmax
            do l = 0, lmax
              mvevil(l, lm1, lm2) = czero
              mvevilef(l, lm1, lm2) = czero
            end do
          end do
        end do
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! ITERMDIR
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if                         ! KREL = 0/1
    ! ----------------------------------------------------------------------------
    lastez = ielast

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End initialise variables
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef CPP_MPI
    if (opt('qdos    ')) then      ! qdos
      if (natyp>=100) then         ! qdos
        open (31, file='qdos.'//char(48+i1/100)//char(48+mod(i1/10,10))// & ! qdos
          char(48+mod(i1,10))//'.'//char(48+ispin)//'.dat') ! qdos
      else                         ! qdos
        open (31, file='qdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'// & ! qdos
          char(48+ispin)//'.dat')  ! qdos
      end if                       ! qdos
      call version_print_header(31) ! qdos
      write (31, '(7(A,3X))') '#   Re(E)', 'Im(E)' & ! qdos
        , 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! qdos
    end if                         ! qdos

    open (30, file='lmdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'// & ! lmdos
      char(48+ispin)//'.dat')      ! lmdos
    call version_print_header(31)  ! lmdos
    write (30, *) ' '              ! lmdos
    write (30, 100) '# ISPIN=', ispin, ' I1=', i1 ! lmdos
100 format (a8, i3, a4, i5)        ! lmdos

    ! write out complex qdos for interpolation to the real axis                   ! complex qdos
    if (test('compqdos')) then     ! complex qdos
      if (natyp>=100) then         ! complex qdos
        open (3031, file='cqdos.'//char(48+i1/100)//char(48+mod(i1/10,10))// & ! complex qdos
          char(48+mod(i1,10))//'.'//char(48+ispin)//'.dat') ! complex qdos
      else                         ! complex qdos
        open (3031, file='cqdos.'//char(48+i1/10)//char(48+mod(i1,10))//'.'// & ! complex qdos
          char(48+ispin)//'.dat')  ! complex qdos
      end if                       ! complex qdos
      call version_print_header(3031) ! complex qdos
      write (3031, '(A)') '# lmax, natyp, nspin, nqdos, ielast:' ! complex qdos
      write (3031, '(5I9)') lmax, natyp, nspin, nqdos, ielast ! complex qdos
      write (3031, '(7(A,3X))') '#   Re(E)', 'Im(E)' & ! complex qdos
        , 'k_x', 'k_y', 'k_z', 'DEN_tot', 'DEN_s,p,...' ! complex qdos
    end if                         ! qdos
#endif

    ! set LMSHIFT value which is need to construct dentmp
    lmshift1(1) = 0
    lmshift1(2) = lmmaxd
    lmshift1(3) = 0
    lmshift1(4) = lmmaxd
    lmshift2(1) = 0
    lmshift2(2) = lmmaxd
    lmshift2(3) = lmmaxd
    lmshift2(4) = 0
    nqdos = 1                      ! qdos
    if (opt('qdos    ')) then      ! qdos
      ! Read BZ path for qdos calculation:                                       ! qdos
      open (67, file='qvec.dat')   ! qdos
      read (67, *) nqdos           ! qdos
      allocate (qvec(3,nqdos), stat=i_stat) ! qdos
      call memocc(i_stat, product(shape(qvec))*kind(qvec), 'QVEC', 'RHOVAL') ! qdos
      do iq = 1, nqdos             ! qdos
        read (67, *)(qvec(ix,iq), ix=1, 3) ! qdos
      end do                       ! qdos
      close (67)                   ! qdos
    end if

    allocate (gflle(lmmaxd,lmmaxd,ielast,nqdos), stat=i_stat)
    call memocc(i_stat, product(shape(gflle))*kind(gflle), 'GFLLE', 'RHOVAL')
    allocate (dum_gflle(lmmaxd,lmmaxd), stat=i_stat)
    call memocc(i_stat, product(shape(dum_gflle))*kind(dum_gflle), 'DUM_GFLLE', 'RHOVAL')
    if (idoldau==1) then
      allocate (gldau(lmmaxd,lmmaxd), stat=i_stat)
      call memocc(i_stat, product(shape(gldau))*kind(gldau), 'GLDAU', 'RHOVAL')
      gldau = czero
    end if

    if (myrank==master) write (1337, *) 'atom', i1
#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)

    do ie_num = 1, ie_end
      ie = ie_start + ie_num
#else
      ! ! omp: start parallel region here
      ! !$omp parallel do default(none)
      ! !$omp& private(eryd,ie,ir,irec,lm1,lm2)
      ! !$omp& private(jlk_index,tmatll,ith)
      ! !$omp& shared(nspin,nsra,iend,ipot,ielast,npan_tot,ncheb,lmax)
      ! !$omp& shared(zat,socscale,ez,rmesh,cleb,rnew,nth,icleb,thetasnew,i1)
      ! !$omp& shared(rpan_intervall,vinsnew,ipan_intervall,r2nefc_loop)
      ! !$omp& shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)
      ! !$omp& shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,cdentemp)
      ! !$omp& shared(tmatsph,den,denlm,gflle,gflle_part,rllleft,sllleft)
      ! !$omp& private(iq,df,ek,tmattemp,gmatll,gmat0,iorb,dentemp)
      ! !$omp& private(rho2ns_temp,dentot,dentmp,rho2,temp1)
      ! !$omp& shared(ldorhoef,nqdos,lmshift1,lmshift2,wez,lmsp,imt1,ifunm)
      ! !$omp& shared(r2orbc,r2nefc,cden,cdenlm,cdenns,rho2nsc_loop)
      ! !$omp& reduction(+:rho2int,espv) reduction(-:muorb)
      ! !$omp& reduction(-:denorbmom,denorbmomsp,denorbmomlm,denorbmomns)
      ! !$omp& shared(t_tgmat)
      ! !$omp& private(alphasph,alphall)
      do ie = 1, ielast
        ie_num = ie
        ie_end = ielast
#endif
        if (t_inc%i_write>0) write (1337, *) 'energy', ie, ez(ie)

        eryd = ez(ie)
        df = wez(ie)/dble(nspinpot)
        ! -------------------------------------------------------------------------
        ! non/scalar-relativistic OR relativistic
        ! -------------------------------------------------------------------------
        if (krel==0) then
          call wfmesh(eryd, ek, cvlight, nsra, zat, r, s, rs, ircut(ipan), irmd, lmax)
          call cradwf(eryd, ek, nsra, alpha, ipan, ircut, cvlight, rs, s, pz, fz, qz, sz, tmat, visp, drdi, r, zat, lirrsol, idoldau, lopt, wldauav, cutoff)
          ! -----------------------------------------------------------------------
          ! Non-spherical
          ! -----------------------------------------------------------------------
          if (ins>0) then
            call pnsqns(ar, cr, dr, drdi, ek, icst, pz, qz, fz, sz, pns, qns, nsra, vins, ipan, irmin, ircut, & ! Added IRMIN 1.7.2014
              cleb, icleb, iend, loflm, lmax, idoldau, lopt, lmlo, lmhi, wldau(1,1,ispin), wldauav, cutoff, lmax)
          end if

          do l = 0, lmax
            ekl(l) = ek*dble(2*l+1)
          end do

          do iq = 1, nqdos         ! qdos
            ! -------------------------------------------------------------------
            ! Read in Green function
            ! -------------------------------------------------------------------
            irec = iq + nqdos*(ie-1) + nqdos*ielast*(ispin-1) + & ! qdos (without qdos, IQ=NQDOS=1)
              nqdos*ielast*nspin*(i1-1) ! qdos
            if (t_tgmat%gmat_to_file) then
              read (69, rec=irec) gmat0
            else
              irec = iq + nqdos*(ie_num-1) + nqdos*t_mpi_c_grid%ntot2*(ispin-1) + nqdos*t_mpi_c_grid%ntot2*nspin*(i1-1)
              gmat0(:, :) = t_tgmat%gmat(:, :, irec)
            end if
            if (test('GMAT=0  ')) then
              write (*, *) 'TEST GMAT=0, setting GMAT to zero'
              gmat0 = czero
            end if
            ! -------------------------------------------------------------------
            ! Spherical/non-spherical input potential
            ! -------------------------------------------------------------------
            if (ins==0) then
              call rholm(den(0,ie,iq), df, gmat0, nsra, rho2ns(1,1,ispin), drdi, ipan, ircut, pz, fz, qz, sz, cleb(1,1), icleb, iend, jend, ekl)
            else
              call rhons(den(0,ie,iq), df, drdi, gmat0, ek, rho2ns(1,1,ispin), ipan, ircut, irmin, thetas, ifunm, lmsp, & ! Added IRMIN 1.7.2014
                nsra, qns, pns, ar, cr, pz, fz, qz, sz, cleb(1,1), icleb, jend, iend, ekl, denlm(1,ie,iq), gflle(:,:,ie,iq))
            end if
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! LDA+U
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if (idoldau==1) then
              do lm1 = 1, lmmaxd
                do lm2 = 1, lmmaxd
                  gldau(lm1, lm2) = gldau(lm1, lm2) + df*gflle(lm1, lm2, ie, iq)
                end do
              end do
            end if
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! LDA+U
            ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifndef CPP_MPI
            ! Write out qdos:
            if (opt('qdos    ')) then ! qdos
              dentot = cmplx(0.d0, 0.d0, kind=dp) ! qdos
              do l = 0, lmaxd1     ! qdos
                dentot = dentot + den(l, ie, iq) ! qdos
              end do               ! qdos
              write (30, 110) eryd, qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! lmdos
                -aimag(dentot)/pi, (-aimag(denlm(l,ie,iq))/pi, l=1, lmmaxd) ! lmdos
              write (31, 110) eryd, qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! qdos
                -aimag(dentot)/pi, (-aimag(denlm(l,ie,iq))/pi, l=1, lmmaxd) ! qdos
110           format (5f10.6, 40e16.8) ! qdos
              ! writeout complex qdos for interpolation                         ! complex qdos
              if (test('compqdos')) then ! complex qdos
                write (3031, 120) eryd, qvec(1, iq), qvec(2, iq), qvec(3, iq), & ! complex qdos
                  dentot, (denlm(l,ie,iq), l=1, lmmaxd) ! complex qdos
              end if               ! complex qdos
120           format (6f10.6, 80e16.8) ! qdos
            end if                 ! qdos
#endif
          end do                   ! IQ                                                             ! qdos
          ! ----------------------------------------------------------------------
          do l = 0, lmaxd1
            espv(l, ispin) = espv(l, ispin) + aimag(eryd*den(l,ie,1)*df)
          end do
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! get the charge at the Fermi energy (IELAST)
          ! call RHOLM/RHONS with the energy weight CONE --> not overwrite DF
          ! with the dummy DENDUM       --> not overwrite DEN
          ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if ((ie==ielast) .and. ldorhoef) then
            if (ins==0) then
              call rholm(dendum, cone, gmat0, nsra, r2nef(1,1,ispin), drdi, ipan, ircut, pz, fz, qz, sz, cleb(1,1), icleb, iend, jend, ekl)
            else
              call rhons(dendum, cone, drdi, gmat0, ek, r2nef(1,1,ispin), ipan, ircut, irmin, thetas, ifunm, lmsp, & ! Added IRMIN 1.7.2014
                nsra, qns, pns, ar, cr, pz, fz, qz, sz, cleb(1,1), icleb, jend, iend, ekl, dum_denlm, dum_gflle)
            end if
          end if
          ! ----------------------------------------------------------------------
        else                       ! ( KREL.EQ.0 )
          ! ----------------------------------------------------------------------
          iq = 1                   ! reset IQ to zero, problem with qdos!!!

          ! !!!! PROBLEM WITH ARRAY DIMENSIONS FOR VTREL ETC. !!!!
#ifdef CPP_MPI
          ! call MPI_FINALIZE(L)
#endif
          ! stop '[rhoval] ERROR array dimensions need to be checked!'
          call drvrho_qdos(ldorhoef, rho2ns, r2nef, den, dmuorb, rhoorb, ie, eryd, df, lastez, gmatll, vtrel, btrel, rmrel, drdirel, r2drdirel, zrel, jwsrel, irshift, solver, &
            soctl, ctl, qmtet, qmphi, itermvdir, mvevil, mvevilef, lmmaxd, lmax, irmd, lmpotd, ielast, nmvecmax, i1, nqdos) ! qdos

          do l = 0, lmaxd1
            espv(l, 1) = espv(l, 1) + aimag(eryd*den(l,ie,iq)*df)
            espv(l, 2) = espv(l, 2) + aimag(eryd*den(l,ie+ielast,iq)*df)
          end do

          do ir = 1, 3
            do l = 0, lmax
              muorb(l, ir) = muorb(l, ir) + aimag(dmuorb(l,ir)*df)
            end do
          end do
        end if
        ! -------------------------------------------------------------------------
        ! Non/scalar-relativistic OR relativistic
        ! -------------------------------------------------------------------------
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LDAU
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LDA+U calculation
        ! IF ( ( IDOLDAU.EQ.1 ).AND.( LOPT.GE.0 ) )
        ! &        CALL DENSITYMAT(DF,PZ,QZ,PNS,QNS,AR,CR,DR,GMATLL(1,1,IE),
        ! &                        IPAN,IRCUT,DRDI,EK,
        ! &                        IRMIN,LOPT,MMAX,LMLO,LMHI,PHILDAU,DENMATC
        ! &        ,den,ie) ! test fivos 19.9.08
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! LDAU
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      end do                       ! IE = 1,IELAST

      ! LDA+U
      if (idoldau==1) then
        do m1 = 1, mmax
          lm1 = lmlo - 1 + m1
          denmatc(1:mmax, m1) = (1.0/(2.0*ci))*(gldau(lmlo:lmhi,lm1)-conjg(gldau(lm1,lmlo:lmhi)))
        end do
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Write out gflle
      if (opt('lmlm-dos')) then    ! lmlm-dos
        if (ispin==1) then         ! lmlm-dos
          ! 4 words = 16 bytes / complex number                                   ! lmlm-dos
          lrecgflle = wlength*4*lmmaxd*lmmaxd*ielast*nqdos ! lmlm-dos
          open (96, access='direct', recl=lrecgflle, file='gflle', form='unformatted') ! lmlm-dos
        end if                     ! lmlm-dos
        irec = i1 + natyp*(ispin-1) ! lmlm-dos
        write (96, rec=irec) gflle(:, :, :, :) ! lmlm-dos
      end if

      i_all = -product(shape(gflle))*kind(gflle)
      deallocate (gflle, stat=i_stat)
      call memocc(i_stat, i_all, 'GFLLE', 'RHOVAL')
      i_all = -product(shape(dum_gflle))*kind(dum_gflle)
      deallocate (dum_gflle, stat=i_stat)
      call memocc(i_stat, i_all, 'DUM_GFLLE', 'RHOVAL')

      if (idoldau==1) then
        i_all = -product(shape(gldau))*kind(gldau)
        deallocate (gldau, stat=i_stat)
        call memocc(i_stat, i_all, 'GLDAU', 'RHOVAL')
      end if

      if (ihost/=1) return

      ! Transformation of ISPIN=1,2 from (spin-down,spin-up) to (charge-density,spin-density)
      if (ispin==2) then
        idim = irmd*lmpotd
        call dscal(idim, 2.d0, rho2ns(1,1,1), 1)
        call daxpy(idim, -0.5d0, rho2ns(1,1,1), 1, rho2ns(1,1,2), 1)
        call daxpy(idim, 1.0d0, rho2ns(1,1,2), 1, rho2ns(1,1,1), 1)
        ! -------------------------------------------------------------------------
        ! Do the same at the Fermi energy
        ! -------------------------------------------------------------------------
        call dscal(idim, 2.d0, r2nef(1,1,1), 1)
        call daxpy(idim, -0.5d0, r2nef(1,1,1), 1, r2nef(1,1,2), 1)
        call daxpy(idim, 1.0d0, r2nef(1,1,2), 1, r2nef(1,1,1), 1)
      end if

    end subroutine rhoval

  end module mod_rhoval
