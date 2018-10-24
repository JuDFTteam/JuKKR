!------------------------------------------------------------------------------------
!> Summary: Calculation of the t-matrix for the new solver
!> Author: 
!> Calculation of the t-matrix for the new solver
!------------------------------------------------------------------------------------
module mod_tmatnewsolver

  implicit none

contains

  !-------------------------------------------------------------------------------
  !> Summary: Calculation of the t-matrix for the new solver
  !> Author: 
  !> Category: solver, single-site, KKRhost
  !> Deprecated: False 
  !> Constructs potential matrix (2x2 for SOC) adding SOC potential with proper form 
  !> of small-component in the case of a scalar-relativistic calculation.
  !> Then creates source terms needed to solve Lippmann-Schwinger equations as described 
  !> in the PhD thesis of David Bauer.
  !-------------------------------------------------------------------------------
  !> @note 
  !> - in the case of using the `SRATRICK` (default behavior) first the spherical solutions are computed (diagonal)
  !>   which are then used to compute the non-spherical solutions in a second step
  !> - Jonathan Chico Apr. 2019: Removed inc.p dependencies and rewrote to Fortran90
  !> @endnote
  !-------------------------------------------------------------------------------
  subroutine tmat_newsolver(ielast,nspin,lmax,zat,socscale,ez,nsra,cleb,icleb,iend, &
    ncheb,npan_tot,rpan_intervall,ipan_intervall,rnew,vinsnew,theta,phi,i1,ipot,    &
    lmpot,lly,deltae,idoldau,lopt,wldau,t_dtmatjij_at,ispin)

#ifdef CPP_OMP
    use :: omp_lib ! necessary for omp functions
#endif
#ifdef CPP_MPI
    use :: mpi ! necessary for MPI functions
    use :: mod_mympi, only: mpiadapt, nranks
    use :: mod_timing, only: timing_start, timing_stop, timings_1a
#endif
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero, cone, cvlight
    use :: global_variables, only: ntotd, ncleb, nrmaxd, mmaxd, nspind, nspotd, iemxd, lmmaxd, lmmaxso, korbit
    use :: mod_wunfiles, only: t_params
    use :: mod_profiling, only: memocc
    use :: mod_mympi, only: myrank, master, distribute_work_energies
    use :: mod_types, only: t_tgmat, t_inc, t_mpi_c_grid, init_tgmat, t_lloyd, init_tlloyd, type_dtmatjijdij, init_t_dtmatjij_at
    use :: mod_save_wavefun, only: t_wavefunctions, find_isave_wavefun, save_wavefunc
    use :: mod_jijhelp, only: calc_dtmatjij
    use :: mod_calcsph, only: calcsph
    use :: mod_rll_global_solutions, only: rll_global_solutions
    use :: mod_rllsllsourceterms, only: rllsllsourceterms
    use :: mod_rllsll, only: rllsll
    use :: mod_spinorbit_ham, only: spinorbit_ham
    use :: mod_vllmat, only: vllmat
    use :: mod_vllmatsra, only: vllmatsra
    use :: mod_regns, only: zgeinv1
#ifdef CPP_BdG
    use :: mod_ioinput, only: ioinput ! to read in something from inputcard
#endif

    implicit none

    ! inputs
    integer, intent (in) :: i1               !! atom index
    integer, intent (in) :: ispin            !! spin index, only used for 'NOSOC' test option where external spin loop is used in main1a 
    integer, intent (in) :: lly              !! LLY /= 0: apply Lloyds formula
    integer, intent (in) :: lopt             !! angular momentum QNUM for the atoms on which LDA+U should be applied (-1 to switch it OFF)
    integer, intent (in) :: lmax             !! Maximum l component in wave function expansion
    integer, intent (in) :: nsra             !! use scalar-relativistic (nsra=2) or non-relativistic (nsra=1) wavefunctions  
    integer, intent (in) :: iend             !! Number of nonzero gaunt coefficients
    integer, intent (in) :: ipot             !! potential index (ipot=(iatom-1)*nspin+ispin)
    integer, intent (in) :: ncheb            !! Number of Chebychev pannels for the new solver
    integer, intent (in) :: nspin            !! Number of spin directions
    integer, intent (in) :: lmpot            !! maximal LM-value of potential expansion: (LPOT+1)**2
    integer, intent (in) :: ielast           !! number of energy points in contour
    integer, intent (in) :: npan_tot         !! total number of panels for Chebychev radial mesh
    integer, intent (in) :: idoldau          !! flag to perform LDA+U
    real (kind=dp), intent (in) :: zat       !! Nuclear charge for a given atom
    real (kind=dp), intent (in) :: phi       !! phi of local spin frame
    real (kind=dp), intent (in) :: theta     !! theta of local spin frame, relative to z-axis
    real (kind=dp), intent (in) :: socscale  !! Spin-orbit scaling for a given atom
    complex (kind=dp), intent (in) :: deltae !! Energy difference for numerical derivative
    integer, dimension (0:ntotd), intent (in) :: ipan_intervall !! indices where panels start in radial Chebychev mesh
    integer, dimension (ncleb, 4), intent (in) :: icleb         !! index array of nonzero Gaunt coefficients [mapping of (lm1, lm2) to lm3]
    real (kind=dp), dimension (ncleb), intent (in) :: cleb      !! values of GAUNT coefficients 
    real (kind=dp), dimension (nrmaxd), intent (in) :: rnew     !! radial mesh points in Chebychev mesh
    real (kind=dp), dimension (0:ntotd), intent (in) :: rpan_intervall        !! radial meshpoints of panel boundaries
    real (kind=dp), dimension (mmaxd, mmaxd, nspind), intent (in) :: wldau    !! potential matrix for LDA+U
    real (kind=dp), dimension (nrmaxd, lmpot, nspotd), intent (in) :: vinsnew !! potential interpolated to Chebychev radial mesh
    complex (kind=dp), dimension (iemxd), intent (in) :: ez  !! list of complex energy points in contour
    type (type_dtmatjijdij), intent (inout) :: t_dtmatjij_at !! derived Data type to store \[\Delta t\]-matrix for Jij calculation

    ! .. Local variables
    integer :: ir, irec, use_sratrick, nvec, lm1, lm2, ie, irmdnew, i11
    integer :: i_stat, lmsize
    integer :: use_fullgmat !! use (l,m,s) coupled matrices or not for 'NOSOC' test option (1/0)
    complex (kind=dp) :: eryd !! energy in Ry
    complex (kind=dp), dimension (nspin*(lmax+1)) :: alphasph !! spherical part of alpha-matrix
    ! .. Local allocatable arrays
    integer, dimension (:), allocatable :: jlk_index
    real (kind=dp), dimension (:, :, :), allocatable :: vins     !! Non-spherical part of the potential
    complex (kind=dp), dimension (:, :), allocatable :: aux      ! LLY
    complex (kind=dp), dimension (:, :), allocatable :: tmat0    !
    complex (kind=dp), dimension (:, :), allocatable :: alpha0   ! LLY
    complex (kind=dp), dimension (:, :), allocatable :: tmatll   !! t-matrix
    complex (kind=dp), dimension (:, :), allocatable :: dtmatll  !! derivative of t-matrix for Lloyd
    complex (kind=dp), dimension (:, :), allocatable :: tmatsph  !! spherical part of t-matrix
    complex (kind=dp), dimension (:, :), allocatable :: alphall  !! alpha matrix for Lloyd
    complex (kind=dp), dimension (:, :), allocatable :: dalphall !! derivatve of alpha matrix for Lloyd
    complex (kind=dp), dimension (:, :, :), allocatable :: hlk
    complex (kind=dp), dimension (:, :, :), allocatable :: jlk
    complex (kind=dp), dimension (:, :, :), allocatable :: hlk2
    complex (kind=dp), dimension (:, :, :), allocatable :: jlk2
    complex (kind=dp), dimension (:, :, :), allocatable :: vnspll0
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rll !! regular solution of radial equation
    complex (kind=dp), dimension (:, :, :, :), allocatable :: sll !! irregular solution of radial equation
    complex (kind=dp), dimension (:, :, :, :), allocatable :: vnspll
    complex (kind=dp), dimension (:, :, :, :), allocatable :: vnspll1
    complex (kind=dp), dimension (:, :, :, :), allocatable :: rllleft !! regular left solution of radial equation
    complex (kind=dp), dimension (:, :, :, :), allocatable :: sllleft !! irregular left solution of radial equation

    ! .. LDAU local variables
    integer :: lmlo, lmhi

    ! .. LLoyd local variables
    integer :: ideriv, signde
    complex (kind=dp) :: tralpha
    complex (kind=dp) :: gmatprefactor
    integer, dimension (:), allocatable :: ipiv ! LLY

    ! .. OMP local variables
    integer :: nth, ith !! total number of threads and thread id for openmp

#ifdef CPP_MPI
    integer, dimension (0:nranks-1) :: ntot_pt, ioff_pt !! auxiliary arrays for MPI communication
#endif
    integer :: ie_end, ie_num, ie_start

    ! rhoqtest
    logical, external :: test, opt
    integer :: mu0, nscoef

    ! BdG
    character (len=100) :: filename
    complex (kind=dp) :: e_shift
    integer :: ier, KBdG
    character (len=256) :: uio
    !-------------------------------------------------------------------------------

#ifdef CPP_OMP
    ! determine if omp parallelisation is used (compiled with -openmp flag and OMP_NUM_THREADS>1)
    !$omp parallel shared(nth,ith)
    !$omp single
    nth = omp_get_num_threads()
    if (t_inc%i_write>0) write (1337, *) 'nth =', nth
    !$omp end single
    !$omp end parallel
#else
    nth = 1
    ith = 0
#endif

    lmsize = lmmaxd/(1+korbit)
    irmdnew = npan_tot*(ncheb+1)

    if (nsra==2) then
      use_sratrick = 1
      if (test('nosph   ')) then
        if (myrank==master .and. ith==0 .and. i1==1 .and. ispin==1) then
          write (*, *) 'Found test option "nosph   ", deactivate SRATRICK'
        end if
        use_sratrick = 0
      end if
    else if (nsra==1) then
      use_sratrick = 0
    end if


    ! .. allocate and initialize arrays
    call allocate_locals_tmat_newsolver(1, irmdnew, lmpot, nspin/(2-korbit), vins, aux, ipiv, tmat0, tmatll, alpha0, dtmatll, alphall, dalphall, jlk_index, nsra, lmmaxso, nth, lmax, vnspll, &
      vnspll0, vnspll1, hlk, jlk, hlk2, jlk2, tmatsph, rll, sll, rllleft, sllleft)

    vins(1:irmdnew, 1:lmpot, 1) = vinsnew(1:irmdnew, 1:lmpot, ipot)
    if (.not.test('NOSOC   ')) vins(1:irmdnew, 1:lmpot, nspin) = vinsnew(1:irmdnew, 1:lmpot, ipot+nspin-1)

    KBdG = 0
#ifdef CPP_BdG
    ! shift potential by EF to change referece point of energy to Fermi level
    ! should later be done automatically in main0
    if (test('BdG_dev ')) then
      ! e_shift = complex(0.723775735132693_dp, 0.0_dp)
      ! e_shift = complex(0.724775735132693_dp, 0.0_dp)

      call ioinput('eshift          ', uio, 1, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) e_shift
        write (*, *) 'e_shift=', e_shift
      else
        e_shift = czero
      end if
    else
      e_shift = czero
    end if
#endif

    ! set up the non-spherical ll' matrix for potential VLL' (done in VLLMAT)
    call vllmat(1, nrmaxd, irmdnew, lmsize, lmmaxso, vnspll0, vins, lmpot, cleb, icleb, iend, nspin/(2-korbit), zat, rnew, use_sratrick, ncleb)
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

    ! start energy loop
    if (myrank==master .and. (t_inc%i_write>0)) write (1337, *) 'atom: ', i1, ' NSRA:', nsra

    ! handle mpi parallelization
    call distribute_work_energies(ielast)
#ifdef CPP_MPI
    ie_start = t_mpi_c_grid%ioff_pt2(t_mpi_c_grid%myrank_at)
    ie_end = t_mpi_c_grid%ntot_pt2(t_mpi_c_grid%myrank_at)
#else
    ie_start = 0
    ie_end = ielast
#endif
    ! Now initialize arrays for tmat, gmat, and gref
    call init_tgmat(t_inc, t_tgmat, t_mpi_c_grid)
    if (lly/=0) call init_tlloyd(t_inc, t_lloyd, t_mpi_c_grid)

    ! consistency check
    if (test('rhoqtest')) then
      if (ielast/=3) stop 'Error: wrong energy contour for rhoqtest'
      ie_start = 1
      ie_end = 1
    end if

    ! For Jij-tensor calculation: allocate array to hold additional t-matrices
    call init_t_dtmatjij_at(t_inc, t_mpi_c_grid, t_dtmatjij_at)

    ! Initialize wfsave
    if (t_inc%i_iteration==0) then
      call find_isave_wavefun(t_wavefunctions)
      ! Reset Nwfsavemax to 0 if test option 'STOP1B  ' is found
      ! to prevent unnessesary storing of wavefunctions
      if (test('STOP1B  ') .and. .not. opt('OPERATOR')) then
        t_wavefunctions%nwfsavemax = 0
      end if
    end if

#ifdef CPP_OMP
    !$omp parallel do default(none)                                            &
    !$omp private(eryd,ie,ir,nvec,lm1,lm2,gmatprefactor)                       &
    !$omp private(jlk_index,tmatll,ith,irec, ie_num)                           &
    !$omp private(tralpha, aux, ideriv, ipiv)                                  &
    !$omp private(alpha0)                                                      &
    !$omp private(alphall)                                                     &
    !$omp private(tmat0)                                                       &
    !$omp private(alphasph)                                                    &
    !$omp private(dtmatll)                                                     &
    !$omp private(dalphall)                                                    &
    !$omp shared(t_inc)                                                        &
    !$omp shared(nspin,nsra,lmax,lmsize,iend,ipot,ielast,npan_tot,ncheb)       &
    !$omp shared(zat,socscale,ez,cleb,rnew,nth,LMPOT,NRMAXD,LMMAXSO,NTOTD)     &
    !$omp shared(rpan_intervall,vinsnew,ipan_intervall,NCLEB)                  &
    !$omp shared(use_sratrick,irmdnew,theta,phi,vins,vnspll0)                  &
    !$omp shared(vnspll1,vnspll,hlk,jlk,hlk2,jlk2,rll,sll,rllleft,sllleft)     &
    !$omp shared(tmatsph, ie_end,t_tgmat,t_lloyd, ie_start, t_dtmatjij_at)     &
    !$omp shared(lly,deltae,i1,t_mpi_c_grid, t_wavefunctions, icleb)           &
    !$omp shared(mu0, nscoef, e_shift, filename, use_fullgmat)
#endif

    do ie_num = 1, ie_end

      ie = ie_start + ie_num

#ifdef CPP_MPI
      ! start timing measurement for this pair of ie and i1, needed for MPIadapt
      call timing_start('time_1a_ieiatom')
#endif

      ! get current thread
      if (nth>=1) then
#ifdef CPP_OMP
        ith = omp_get_thread_num()
#endif
      else
        ith = 0
      end if

      ! In case of Lloyds formula the derivative of t is needed.
      ! Then calculate t at E+dE, E-dE and average for t, subtract for dt/dE
      tmatll = czero
      alphall = czero              ! LLY
      dtmatll = czero              ! LLY
      dalphall = czero             ! LLY
      ideriv = 0
      if (lly/=0) ideriv = 1
      do signde = -ideriv, ideriv, 2
        eryd = ez(ie) + real(signde, kind=dp)*deltae/2.d0 ! LLY

#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        if (test('BdG_dev ')) then
          write (*, '(A,4ES21.7)') 'shifting energy by e_fermi:', eryd, e_shift
          eryd = eryd + e_shift
        end if
#endif

        if (t_inc%i_write>0) write (1337, *) 'energy:', ie, '', eryd
#ifdef CPP_OMP
        if (ie==1 .and. (t_inc%i_write>0)) write (1337, *) 'nested omp?', omp_get_nested()
        !$omp end critical
#endif

        if ( .not. test('NOSOC   ')) then
          ! Contruct the spin-orbit coupling hamiltonian and add to potential
          call spinorbit_ham(lmax, lmsize, vins, rnew, eryd, zat, cvlight, socscale, nspin, lmpot, theta, phi, ipan_intervall, rpan_intervall, npan_tot, ncheb, irmdnew, nrmaxd, &
            vnspll0(:,:,:), vnspll1(:,:,:,ith), '1')
        else
          vnspll1(:,:,:,ith) = vnspll0(:,:,:)
        end if

#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        ! test writeout of VNSPLL1
        if (test('BdG_dev ')) then
          open (7352834, file='vnspll_SOC.txt', form='formatted')
          write (7352834, '(A,3I9)') '# LMMAXSO,LMMAXSO,IRMDNEW=', lmmaxso, lmmaxso, irmdnew
          write (7352834, '(2F25.14)') vnspll1(:, :, :, ith)
          close (7352834)
        end if
#endif
#ifdef CPP_OMP
        !$omp end critical
#endif

        ! now extend matrix for the SRA treatment
        vnspll(:, :, :, ith) = czero

        if (nsra==2) then
          if (use_sratrick==0) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxso,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=0')
          else if (use_sratrick==1) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxso,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=Vsph')
          end if
        else
          vnspll(:, :, :, ith) = vnspll1(:, :, :, ith)
        end if

#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        ! test writeout of VNPSLL
        if (test('BdG_dev ')) then
          open (7352834, file='vnspll_sra.txt', form='formatted')
          if (nsra==2) then
            write (7352834, '(A,3I9)') '# 2*LMMAXSO,2*LMMAXSO,IRMDNEW=', 2*lmmaxso, 2*lmmaxso, irmdnew
          else
            write (7352834, '(A,3I9)') '# LMMAXSO,LMMAXSO,IRMDNEW=', lmmaxso, lmmaxso, irmdnew
          end if
          write (7352834, '(2F25.14)') vnspll(:, :, :, ith)
          close (7352834)
        end if
#endif
#ifdef CPP_OMP
        !$omp end critical
#endif

        ! Calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        hlk(:, :, ith) = czero
        jlk(:, :, ith) = czero
        hlk2(:, :, ith) = czero
        jlk2(:, :, ith) = czero
        gmatprefactor = czero
        if (test('NOSOC   ')) then
          use_fullgmat = 0
        else
          use_fullgmat = 1
        end if
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxso, use_fullgmat, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)

#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        if (test('BdG_dev ')) then
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') jlk(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') hlk(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_jlk2_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') jlk2(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_source_hlk2_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') hlk2(:, :, ith)
          close (888888)
        end if
#endif
#ifdef CPP_OMP
        !$omp end critical
#endif

        ! Using spherical potential as reference
        if (use_sratrick==1) then
          tmatsph(:, ith) = czero
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin/(2-korbit), zat, eryd, lmpot, lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk(:,:,ith), jlk(:,:,ith), &
            hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, tmatsph(:,ith), alphasph, use_sratrick)
#ifdef CPP_BdG
        if (test('BdG_dev ')) then
          write (filename, '(A,I0.3,A,I0.3,A)') 'tmatsph_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso=', lmmaxso, ' lmmaxso=', lmmaxso
          write (888888, '(2ES21.9)') tmatsph(:, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_jlk_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') jlk(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_hlk_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') hlk(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_jlk2_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') jlk2(:, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_sph_hlk2_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,2ES15.7)') '# dimension: 4*(LMAX+1)=', 4*(lmax+1), ' IRMDNEW=', irmdnew, ' ; ERYD=', eryd
          write (888888, '(2ES21.9)') hlk2(:, :, ith)
          close (888888)
        end if
#endif
        end if
        ! Calculate the tmat and wavefunctions
        rll(:, :, :, ith) = czero
        sll(:, :, :, ith) = czero

        ! Right solutions
        tmat0 = czero
        alpha0 = czero             ! LLY
        ! faster calculation of RLL.
        ! no irregular solutions are needed in self-consistent iterations
        ! because the t-matrix depends only on RLL
        if (opt('RLL-SLL ') .and. .not. (opt('XCPL    ') .or. opt('OPERATOR'))) then
          call rll_global_solutions(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), tmat0(:,:), ncheb, npan_tot, lmmaxso, nvec*lmmaxso, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', use_sratrick, alpha0(:,:))
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rll(:,:,:,ith), sll(:,:,:,ith), tmat0(:,:), ncheb, npan_tot, lmmaxso, nvec*lmmaxso, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alpha0(:,:))
        end if

        if (nsra==2) then
          rll(lmmaxso+1:nvec*lmmaxso, :, :, ith) = rll(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
          sll(lmmaxso+1:nvec*lmmaxso, :, :, ith) = sll(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
        end if
#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        if (test('BdG_dev ')) then
          write (filename, '(A,I0.3,A,I0.3,A)') 'rll_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso*nvec=', nvec*lmmaxso, ' lmmaxso=', lmmaxso, ' irmdnew=', irmdnew
          write (888888, '(2ES21.9)') rll(:, :, :, ith)
          close (888888)
          write (filename, '(A,I0.3,A,I0.3,A)') 'sll_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso*nvec=', nvec*lmmaxso, ' lmmaxso=', lmmaxso, ' irmdnew=', irmdnew
          write (888888, '(2ES21.9)') sll(:, :, :, ith)
          close (888888)
        end if
#endif
#ifdef CPP_OMP
        !$omp end critical
#endif

        ! add spherical contribution of tmatrix
        if (use_sratrick==1) then
          do lm1 = 1, lmmaxso
            tmat0(lm1, lm1) = tmat0(lm1, lm1) + tmatsph(jlk_index(lm1), ith)
          end do
          if (lly/=0) then
            do lm2 = 1, lmmaxso
              do lm1 = 1, lmmaxso
                ! alphasph is multiplied not added
                alpha0(lm1, lm2) = alphasph(jlk_index(lm1))*alpha0(lm1, lm2) ! LLY
              end do
            end do
          end if                   ! LLY
        end if
        tmatll(:, :) = tmatll(:, :) + tmat0(:, :)
        if (lly/=0) then
          alphall(:, :) = alphall(:, :) + alpha0(:, :)
          dtmatll(:, :) = dtmatll(:, :) + real(signde, kind=dp)*tmat0(:, :) ! LLY
          dalphall(:, :) = dalphall(:, :) + real(signde, kind=dp)*alpha0(:, :) ! LLY
        end if

#ifdef CPP_OMP
        !$omp critical
#endif
#ifdef CPP_BdG
        if (test('BdG_dev ')) then
          write (filename, '(A,I0.3,A,I0.3,A)') 'tmat_atom_', i1, '_energ_', ie, '.dat'
          open (888888, file=trim(filename), form='formatted')
          write (888888, '(A,I9,A,I9,A,I9)') '# dimension: lmmaxso=', lmmaxso, ' lmmaxso=', lmmaxso
          write (888888, '(2ES21.9)') tmatll(:, :)
          close (888888)
        end if
#endif
#ifdef CPP_OMP
        !$omp end critical
#endif

      end do                       ! signde=-ideriv,ideriv,2 ! lly

      ! Average values of t-matrix and alpha at e+de and e-de
      tmatll(:, :) = tmatll(:, :)/real(1+ideriv, kind=dp) ! LLY
      if (lly/=0) then
        alphall(:, :) = alphall(:, :)/real(1+ideriv, kind=dp) ! LLY
        ! Contruct derivative of t-matrix and alpha
        dtmatll(:, :) = dtmatll(:, :)/deltae ! LLY
        dalphall(:, :) = dalphall(:, :)/deltae ! LLY
      end if
      if (lly/=0) then
        ! calculate Tr[alpha^-1*dalpha/de] for LLoyd's formula
        alpha0 = czero             ! LLY
        aux = czero                ! LLY
        call zgeinv1(alphall, alpha0, aux, ipiv, lmmaxso)
        call zgemm('N','N',lmmaxso,lmmaxso,lmmaxso,cone,alpha0,lmmaxso,dalphall,    &
          lmmaxso,czero,aux,lmmaxso) ! LLY
        ! Trace of AUX
        tralpha = czero            ! LLY
        do lm1 = 1, lmmaxso
          tralpha = tralpha + aux(lm1, lm1) ! LLY
        end do
      end if                       ! LLY

      if (test('rhoqtest') .and. ie==2) then
        ! read in mu0 atom index
        open (9999, file='mu0')
        read (9999, *) mu0, nscoef
        close (9999)
      end if

      ! Calculate additional t-matrices for Jij-tensor calculation
      if (t_dtmatjij_at%calculate .or. (t_wavefunctions%isave_wavefun(i1,ie)>0 .and.&
         (t_wavefunctions%save_rllleft .or. t_wavefunctions%save_sllleft)) .or.     &
         ((test('rhoqtest') .and. ie==2) .and. (i1==mu0))) then ! rhoqtest
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Calculate the left-hand side solution this needs to be done for the
        ! calculation of t-matrices for Jij tensor or if wavefunctions should be saved
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Contruct the spin-orbit coupling hamiltonian and add to potential
        call spinorbit_ham(lmax,lmsize,vins,rnew,eryd,zat,cvlight,socscale,nspin,   &
          lmpot,theta,phi,ipan_intervall,rpan_intervall,npan_tot,ncheb,irmdnew,     &
          nrmaxd,vnspll0(:,:,:),vnspll1(:,:,:,ith),'transpose')

        ! Extend matrix for the SRA treatment
        vnspll(:, :, :, ith) = czero
        if (nsra==2) then
          if (use_sratrick==0) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxso,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=0')
          else if (use_sratrick==1) then
            call vllmatsra(vnspll1(:,:,:,ith),vnspll(:,:,:,ith),rnew,lmmaxso,       &
              irmdnew,nrmaxd,eryd,lmax,0,'Ref=Vsph')
          end if
        else
          vnspll(:, :, :, ith) = vnspll1(:, :, :, ith)
        end if

        ! Calculate the source terms in the Lippmann-Schwinger equation
        ! these are spherical hankel and bessel functions
        hlk(:, :, ith) = czero
        jlk(:, :, ith) = czero
        hlk2(:, :, ith) = czero
        jlk2(:, :, ith) = czero
        gmatprefactor = czero
        jlk_index = 0
        call rllsllsourceterms(nsra, nvec, eryd, rnew, irmdnew, nrmaxd, lmax, lmmaxso, use_fullgmat, jlk_index, hlk(:,:,ith), jlk(:,:,ith), hlk2(:,:,ith), jlk2(:,:,ith), gmatprefactor)

        ! Using spherical potential as reference
        ! notice that exchange the order of left and right hankel/bessel functions
        if (use_sratrick==1) then
          tmatsph(:, ith) = czero
          call calcsph(nsra, irmdnew, nrmaxd, lmax, nspin/(2-korbit), zat, eryd, lmpot, lmmaxso, rnew, vins, ncheb, npan_tot, rpan_intervall, jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), &
            hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, alphasph, tmatsph(:,ith), use_sratrick)
        end if

        ! Calculate the tmat and wavefunctions
        rllleft(:, :, :, ith) = czero
        sllleft(:, :, :, ith) = czero

        ! Left solutions
        ! notice that exchange the order of left and right hankel/bessel functions
        tmat0 = czero
        alpha0 = czero             ! LLY
        ! faster calculation of RLL.
        ! no left solutions are needed in self-consistent iterations
        ! because the t-matrix depends only on RLL
        if (opt('RLL-SLL ') .and. .not. (opt('XCPL    ') .or. opt('OPERATOR'))) then
          ! do nothing
        else
          call rllsll(rpan_intervall, rnew, vnspll(:,:,:,ith), rllleft(:,:,:,ith), sllleft(:,:,:,ith), tmat0, ncheb, npan_tot, lmmaxso, nvec*lmmaxso, nsra*(1+korbit)*(lmax+1), irmdnew, nsra, &
            jlk_index, hlk2(:,:,ith), jlk2(:,:,ith), hlk(:,:,ith), jlk(:,:,ith), gmatprefactor, '1', '1', '0', use_sratrick, alpha0)
        end if
        if (nsra==2) then
          rllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith) = rllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
          sllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith) = sllleft(lmmaxso+1:nvec*lmmaxso, :, :, ith)/cvlight
        end if

        if (test('rhoqtest')) then
#ifdef CPP_OMP
          write (*, *) 'rhoqtest does not work in OMP version!!'
          write (*, *) 'please use hybrid compilation mode'
          stop
#else
          open (9999, file='params.txt')
          write (9999, *) lmmaxso, t_params%natyp
          write (9999, *) t_params%naez, t_params%nclsd, t_params%nr, t_params%nembd1 - 1, t_params%lmax
          write (9999, *) t_params%alat
          close (9999)

          open (9999, file='host.txt')
          write (9999, *) t_params%rbasis(1:3, 1:t_params%natyp)
          write (9999, *) t_params%rcls(1:3, 1:t_params%nclsd, 1:t_params%nclsd),t_params%rr(1:3, 0:t_params%nr),t_params%atom(1:t_params%nclsd,1:t_params%naez+t_params%nembd1-1)
          write (9999, *) t_params%cls(1:t_params%naez+t_params%nembd1-1), t_params%ezoa(1:t_params%nclsd, 1:t_params%naez+t_params%nembd1-1), t_params%nacls(1:t_params%nclsd)
          close (9999)

          open (9999, file='wavefunctions.txt')
          write (9999, '(100I9)') ntotd, npan_tot, ncheb, nsra, irmdnew
          write (9999, '(1000E26.17)') rnew(1:irmdnew)
          do ir = 1, irmdnew
            do lm1 = 1, nsra*lmmaxso
              do lm2 = 1, lmmaxso
                write (9999, '(20000E16.7)') rll(lm1, lm2, ir, ith), rllleft(lm1, lm2, ir, ith)
              end do
            end do
          end do
          do lm1 = 0, npan_tot
            write (9999, '(E16.7,I9)') rpan_intervall(lm1), ipan_intervall(lm1)
          end do
          close (9999)
#endif
        end if                     ! test('rhoqtest')

        !----------------------------------------------------------------------------
        ! Calculate the left-hand side solution
        !----------------------------------------------------------------------------

      end if                       ! t_dtmatJij_at%calculate .or. t_wavefunctions%Nwfsavemax>0

      ! save_wavefuncions
      if (t_wavefunctions%nwfsavemax>0) then
        ! here all four (left, right, regular and irregular) are stored, the memory demand cound be reduced by a factor 2 if only the right solution would be computed here and saved and the left solution would be calculated later in main1c
        call save_wavefunc(t_wavefunctions,rll,rllleft,sll,sllleft,i1,ie,nsra,      &
          lmmaxso,irmdnew,ith)
      end if

      if (t_dtmatjij_at%calculate) then
        call calc_dtmatjij(lmsize,lmmaxso,lmpot,ntotd,nrmaxd,nsra,irmdnew,nspin,    &
          vins,rllleft(:,:,:,ith),rll(:,:,:,ith),rpan_intervall,ipan_intervall,     &
          npan_tot,ncheb,cleb,icleb,iend,ncleb,rnew,t_dtmatjij_at%dtmat_xyz(:,:,:,ie_num))

      end if                       ! t_dtmatJij_at%calculate

      ! writeout
#ifdef CPP_OMP
      !$omp critical
#endif
      if (t_tgmat%tmat_to_file) then
#ifndef CPP_OMP
        if (test('rhoqtest')) then

          if (ie_num==1 .and. i1==1) then
            write (*, *)           ! status bar
            write (*, *) 'rhoq: write-out t-mat', ie_end, t_params%natyp
            write (*, '("Loop over points:|",5(1X,I2,"%",5X,"|"),1X,I3,"%")') 0, 20, 40, 60, 80, 100
            write (*, fmt=100, advance='no') ! beginning of statusbar
          end if

          if (myrank==master) then
            if (t_params%natyp*ie_end>=50) then
              if (mod(i1+t_params%natyp*(ie_num-1),(t_params%natyp*ie_end/50))==0) write (6, fmt=110, advance='no')
            else
              write (6, fmt=110, advance='no')
            end if
          end if

          ! write(*,*) 'rotating with', theta, phi
          ! lmGF0D= (LMAXD+1)**2
          ! caLL ROTATEMATRIX(TMATLL,THETA,PHI,LMGF0D,0)

        end if                     ! test('rhoqtest')
#endif
        irec = ie + ielast*(ispin-1) + ielast*nspin/(1+korbit)*(i1-1)
        write (69, rec=irec) tmatll(:, :)
        ! human readable writeout if test option is hit
        if (test('fileverb')) then
          write (696969, '(i9,20000F15.7)') irec, tmatll(:, :)
        end if
      else
#ifdef CPP_MPI
        i11 = i1-t_mpi_c_grid%ioff_pt1(t_mpi_c_grid%myrank_ie)
#else
        i11 = i1
#endif
        irec = ie_num + ie_end*(ispin-1) + ie_end*nspin/(1+korbit)*(i11-1)
        t_tgmat%tmat(:, :, irec) = tmatll(:, :)
      end if
      if (lly/=0) then
        if (t_lloyd%dtmat_to_file) then
          irec = ie + ielast*(ispin-1) + ielast*nspin/(1+korbit)*(i1-1)
          write (691, rec=irec) dtmatll(:, :) ! LLY
          if (test('fileverb')) then
            write (691691691, '(i9,20000F15.7)') irec, dtmatll(:, :)
          end if
        else
          irec = ie_num + ie_end*(ispin-1) + ie_end*nspin/(1+korbit)*(i1-1)
          t_lloyd%dtmat(:, :, irec) = dtmatll(:, :)
        end if
        if (t_lloyd%tralpha_to_file) then
          irec = ie + ielast*(ispin-1) + ielast*nspin/(1+korbit)*(i1-1)
          write (692, rec=irec) tralpha ! LLY
          if (test('fileverb')) then
            write (692692692, '(i9,20000F15.7)') irec, tralpha
          end if
        else
          irec = ie_num + ie_end*(ispin-1) + ie_end*nspin/(1+korbit)*(i1-1)
          t_lloyd%tralpha(irec) = tralpha
        end if
      end if
#ifdef CPP_OMP
      !$omp end critical
#endif

#ifdef CPP_MPI
      ! stop timing measurement for this pair of ie and i1, needed for MPIadapt
      if (mpiadapt>0) call timing_stop('time_1a_ieiatom', save_out=timings_1a(ie,i1))
#endif

    end do                         ! IE loop
#ifdef CPP_OMP
    !$omp end parallel do
#endif

100 format ('                 |')  ! status bar
110 format ('|')                   ! status bar
    if (test('rhoqtest') .and. i1==t_params%natyp .and. myrank==master) write (6, *) ! status bar
    ! finished kpts status bar

    ! deallocate arrays
    call allocate_locals_tmat_newsolver(-1,irmdnew,lmpot,nspin,vins,aux,ipiv,tmat0, &
      tmatll,alpha0,dtmatll,alphall,dalphall,jlk_index,nsra,lmmaxso,nth,lmax,vnspll,&
      vnspll0,vnspll1,hlk,jlk,hlk2,jlk2,tmatsph,rll,sll,rllleft,sllleft)

  end subroutine tmat_newsolver


  !-------------------------------------------------------------------------------
  !> Summary: Wrapper routine for the allocation/deallocation of the arrays for the t-matrix
  !> calculation 
  !> Author: Philipp Ruessmann
  !> Category: solver, single-site, memory-management, KKRhost
  !> Deprecated: False 
  !> Wrapper routine for the allocation/deallocation of the arrays for the t-matrix
  !> calculation.
  !-------------------------------------------------------------------------------
  subroutine allocate_locals_tmat_newsolver(allocmode,irmdnew,lmpot,nspin,vins,aux, &
    ipiv,tmat0,tmatll,alpha0,dtmatll,alphall,dalphall,jlk_index,nsra,lmmaxso,nth,   &
    lmax,vnspll,vnspll0,vnspll1,hlk,jlk,hlk2,jlk2,tmatsph,rll,sll,rllleft,sllleft)
    use :: mod_datatypes, only: dp
    use :: mod_constants, only: czero
    use :: global_variables, only: korbit
    use :: mod_profiling, only: memocc
    use :: mod_save_wavefun, only: t_wavefunctions
    implicit none

    ! array dimensions
    integer, intent (in) :: allocmode !! allocation mode (1: allocate and initialize, other: deallocate)
    integer, intent (in) :: irmdnew   !! number of radial points in Chebycheb mesh
    integer, intent (in) :: lmpot     !! lm-cutoff of potential expansion
    integer, intent (in) :: nspin     !! number of spin channels
    integer, intent (in) :: nsra      !! scalar-relativistic (nsra=2) or non-relativistic (nsra=1)
    integer, intent (in) :: lmmaxso   !! cutoff of combined (l,m,s) index
    integer, intent (in) :: nth       !! number of OpenMP threads
    integer, intent (in) :: lmax      !! lmax cutoff

    ! allocatable arrays
    real (kind=dp), allocatable, intent (inout) :: vins(:, :, :)
    complex (kind=dp), allocatable, intent (inout) :: aux(:, :), tmat0(:, :), tmatll(:, :), alpha0(:, :), dtmatll(:, :), alphall(:, :), dalphall(:, :)
    integer, allocatable, intent (inout) :: ipiv(:), jlk_index(:)
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: vnspll
    complex (kind=dp), allocatable, dimension (:, :, :), intent (inout) :: vnspll0
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: vnspll1
    complex (kind=dp), allocatable, dimension (:, :, :), intent (inout) :: jlk
    complex (kind=dp), allocatable, dimension (:, :, :), intent (inout) :: hlk
    complex (kind=dp), allocatable, dimension (:, :, :), intent (inout) :: jlk2
    complex (kind=dp), allocatable, dimension (:, :, :), intent (inout) :: hlk2
    complex (kind=dp), allocatable, dimension (:, :), intent (inout) :: tmatsph
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: rll
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: sll
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: rllleft
    complex (kind=dp), allocatable, dimension (:, :, :, :), intent (inout) :: sllleft

    ! local
    integer :: i_stat

    logical, external :: test, opt


    if (allocmode==1) then ! allocate and initialize

      ! potential arrays
      allocate (vnspll(nsra*lmmaxso,nsra*lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'allocate_locals_tmat_newsolver')
      vnspll = czero
      allocate (vnspll0(lmmaxso,lmmaxso,irmdnew), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'allocate_locals_tmat_newsolver')
      vnspll0 = czero
      allocate (vnspll1(lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'allocate_locals_tmat_newsolver')
      vnspll1 = czero

      ! source terms (bessel and hankel functions)
      allocate (hlk(1:nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(hlk))*kind(hlk), 'HLK', 'allocate_locals_tmat_newsolver')
      hlk = czero
      allocate (jlk(1:nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(jlk))*kind(jlk), 'JLK', 'allocate_locals_tmat_newsolver')
      jlk = czero
      allocate (hlk2(1:nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(hlk2))*kind(hlk2), 'HLK2', 'allocate_locals_tmat_newsolver')
      hlk2 = czero
      allocate (jlk2(1:nsra*(1+korbit)*(lmax+1),irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(jlk2))*kind(jlk2), 'JLK2', 'allocate_locals_tmat_newsolver')
      jlk2 = czero

      ! Spherical part of tmatrix (used with SRATRICK)
      allocate (tmatsph(nspin*(lmax+1),0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(tmatsph))*kind(tmatsph), 'TMATSPH', 'allocate_locals_tmat_newsolver')
      tmatsph = czero

      ! Regular and irregular wavefunctions
      allocate (rll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(rll))*kind(rll), 'RLL', 'allocate_locals_tmat_newsolver')
      rll = czero
      allocate (sll(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
      call memocc(i_stat, product(shape(sll))*kind(sll), 'SLL', 'allocate_locals_tmat_newsolver')
      sll = czero

      ! Left regular and irregular wavefunctions (used here only in case of XCPL or saving of left wavefunctions)
      if (opt('XCPL    ') .or. (t_wavefunctions%save_rllleft .or. t_wavefunctions%save_sllleft .or. test('rhoqtest'))) then
        allocate (rllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
        call memocc(i_stat, product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'allocate_locals_tmat_newsolver')
        rllleft = czero
        allocate (sllleft(nsra*lmmaxso,lmmaxso,irmdnew,0:nth-1), stat=i_stat)
        call memocc(i_stat, product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'allocate_locals_tmat_newsolver')
        sllleft = czero
      else
        allocate (rllleft(1,1,1,0:nth-1), stat=i_stat)
        call memocc(i_stat, product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'allocate_locals_tmat_newsolver')
        rllleft = czero
        allocate (sllleft(1,1,1,0:nth-1), stat=i_stat)
        call memocc(i_stat, product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'allocate_locals_tmat_newsolver')
        sllleft = czero
      end if                       ! ( opt('XCPL    ') .or. ... )

      allocate (vins(irmdnew,lmpot,nspin), stat=i_stat)
      call memocc(i_stat, product(shape(vins))*kind(vins), 'VINS', 'allocate_locals_tmat_newsolver')
      vins = 0.0d0
      allocate (aux(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(aux))*kind(aux), 'AUX', 'allocate_locals_tmat_newsolver')
      aux = czero
      allocate (ipiv(lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(ipiv))*kind(ipiv), 'IPIV', 'allocate_locals_tmat_newsolver')
      ipiv = 0
      allocate (tmat0(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(tmat0))*kind(tmat0), 'TMAT0', 'allocate_locals_tmat_newsolver')
      tmat0 = czero
      allocate (tmatll(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(tmatll))*kind(tmatll), 'TMATLL', 'allocate_locals_tmat_newsolver')
      tmatll = czero
      allocate (alpha0(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(alpha0))*kind(alpha0), 'ALPHA0', 'allocate_locals_tmat_newsolver')
      alpha0 = czero
      allocate (dtmatll(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(dtmatll))*kind(dtmatll), 'DTMATLL', 'allocate_locals_tmat_newsolver')
      dtmatll = czero
      allocate (alphall(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(alphall))*kind(alphall), 'ALPHALL', 'allocate_locals_tmat_newsolver')
      alphall = czero
      allocate (dalphall(lmmaxso,lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(dalphall))*kind(dalphall), 'DALPHALL', 'allocate_locals_tmat_newsolver')
      dalphall = czero
      allocate (jlk_index(nsra*lmmaxso), stat=i_stat)
      call memocc(i_stat, product(shape(jlk_index))*kind(jlk_index), 'JLK_INDEX', 'allocate_locals_tmat_newsolver')
      jlk_index = 0

    else ! allocmode/=1: deallocate arrays

      if (nsra==2) then
        deallocate (vnspll, stat=i_stat)
        call memocc(i_stat, -product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'allocate_locals_tmat_newsolver')
      else
        deallocate (vnspll, stat=i_stat)
        call memocc(i_stat, -product(shape(vnspll))*kind(vnspll), 'VNSPLL', 'allocate_locals_tmat_newsolver')
      end if
      deallocate (vnspll0, stat=i_stat)
      call memocc(i_stat, -product(shape(vnspll0))*kind(vnspll0), 'VNSPLL0', 'allocate_locals_tmat_newsolver')
      deallocate (vnspll1, stat=i_stat)
      call memocc(i_stat, -product(shape(vnspll1))*kind(vnspll1), 'VNSPLL1', 'allocate_locals_tmat_newsolver')

      deallocate (hlk, stat=i_stat)
      call memocc(i_stat, -product(shape(hlk))*kind(hlk), 'HLK', 'allocate_locals_tmat_newsolver')
      deallocate (jlk, stat=i_stat)
      call memocc(i_stat, -product(shape(jlk))*kind(jlk), 'JLK', 'allocate_locals_tmat_newsolver')
      deallocate (hlk2, stat=i_stat)
      call memocc(i_stat, -product(shape(hlk2))*kind(hlk2), 'HLK2', 'allocate_locals_tmat_newsolver')
      deallocate (jlk2, stat=i_stat)
      call memocc(i_stat, -product(shape(jlk2))*kind(jlk2), 'JLK2', 'allocate_locals_tmat_newsolver')

      deallocate (tmatsph, stat=i_stat)
      call memocc(i_stat, -product(shape(tmatsph))*kind(tmatsph), 'TMATSPH', 'allocate_locals_tmat_newsolver')

      deallocate (rll, stat=i_stat)
      call memocc(i_stat, -product(shape(rll))*kind(rll), 'RLL', 'allocate_locals_tmat_newsolver')
      deallocate (sll, stat=i_stat)
      call memocc(i_stat, -product(shape(sll))*kind(sll), 'SLL', 'allocate_locals_tmat_newsolver')

      if (opt('XCPL    ') .or. (t_wavefunctions%save_rllleft .or. t_wavefunctions%save_sllleft .or. test('rhoqtest'))) then
        deallocate (rllleft, stat=i_stat)
        call memocc(i_stat, -product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'allocate_locals_tmat_newsolver')
        deallocate (sllleft, stat=i_stat)
        call memocc(i_stat, -product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'allocate_locals_tmat_newsolver')
      else
        deallocate (rllleft, stat=i_stat)
        call memocc(i_stat, -product(shape(rllleft))*kind(rllleft), 'RLLLEFT', 'allocate_locals_tmat_newsolver')
        deallocate (sllleft, stat=i_stat)
        call memocc(i_stat, -product(shape(sllleft))*kind(sllleft), 'SLLLEFT', 'allocate_locals_tmat_newsolver')
      end if                       ! ( opt('XCPL    ') .or. ... )

      deallocate (vins, stat=i_stat)
      call memocc(i_stat, -product(shape(vins))*kind(vins), 'VINS', 'allocate_locals_tmat_newsolver')
      deallocate (aux, stat=i_stat)
      call memocc(i_stat, -product(shape(aux))*kind(aux), 'AUX', 'allocate_locals_tmat_newsolver')
      deallocate (ipiv, stat=i_stat)
      call memocc(i_stat, -product(shape(ipiv))*kind(ipiv), 'IPIV', 'allocate_locals_tmat_newsolver')
      deallocate (tmat0, stat=i_stat)
      call memocc(i_stat, -product(shape(tmat0))*kind(tmat0), 'TMAT0', 'allocate_locals_tmat_newsolver')
      deallocate (tmatll, stat=i_stat)
      call memocc(i_stat, -product(shape(tmatll))*kind(tmatll), 'TMATLL', 'allocate_locals_tmat_newsolver')
      deallocate (alpha0, stat=i_stat)
      call memocc(i_stat, -product(shape(alpha0))*kind(alpha0), 'ALPHA0', 'allocate_locals_tmat_newsolver')
      deallocate (dtmatll, stat=i_stat)
      call memocc(i_stat, -product(shape(dtmatll))*kind(dtmatll), 'DTMATLL', 'allocate_locals_tmat_newsolver')
      deallocate (alphall, stat=i_stat)
      call memocc(i_stat, -product(shape(alphall))*kind(alphall), 'ALPHALL', 'allocate_locals_tmat_newsolver')
      deallocate (dalphall, stat=i_stat)
      call memocc(i_stat, -product(shape(dalphall))*kind(dalphall), 'DALPHALL', 'allocate_locals_tmat_newsolver')
      deallocate (jlk_index, stat=i_stat)
      call memocc(i_stat, -product(shape(jlk_index))*kind(jlk_index), 'JLK_INDEX', 'allocate_locals_tmat_newsolver')

    end if ! allocmode ==1 or /=1

  end subroutine allocate_locals_tmat_newsolver

end module mod_tmatnewsolver
