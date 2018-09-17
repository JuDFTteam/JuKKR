module mod_kloopz1

contains

  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: KLOOPZ1_QDOS
  !> @note
  !> - Jonathan Chico Apr. 2018: Removed inc.p dependencies and rewrote to Fortran90
  ! -------------------------------------------------------------------------------
  subroutine kloopz1_qdos(eryd, gmatll, ins, alat, ie, igf, nshell, naez, nofks, volbz, bzkp, volcub, cls, nacls, naclsmax, ncls, rr, rbasis, ezoa, atom, rcls, icc, ginp, ideci, &
    lefttinvll, righttinvll, vacflag, nlbasis, nrbasis, factl, natomimp, nsymat, dsymll, ratom, rrot, nsh1, nsh2, ijtabsym, ijtabsh, icheck, invmod, refpot, trefll, tsst, msst, &
    cfctor, cfctorinv, crel, rc, rrel, srrel, irrel, nrrel, drotq, symunitary, kmrot, natyp, ncpa, icpa, itcpamax, cpatol, noq, iqat, itoq, conc, iprint, icpaflag, ispin, nspin, &
    tqdos, iqdosrun, &             ! qdos ruess
    dtrefll, dtmatll, dginp, lly_grtr, tracet, lly) ! LLY Lloyd

    use :: mod_types, only: t_inc
    use :: mod_mympi, only: myrank, master
    use :: global_variables
    use :: constants
    use :: mod_profiling
    use :: mod_datatypes, only: dp
    use :: mod_rotgll
    use :: mod_mssinit
    use :: mod_kkrmat01
    use :: mod_gijdmat
    use :: mod_cpamillsx
    use :: mod_cmatstr
    use :: mod_symetrmat
    use :: mod_rotate
    use :: mod_projtau
    use :: mod_cinit

    implicit none

    ! .. Parameters
    integer, parameter :: linmax = 1
    real (kind=dp), parameter :: tolmssq = 1.0d-6
    ! .. Input variables
    integer, intent (in) :: ie
    integer, intent (in) :: lly    !! LLY <> 0 => use Lloyd formula
    integer, intent (in) :: igf    !! Do not print or print (0/1) the KKRFLEX_* files
    integer, intent (in) :: ins    !! 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (in) :: icc    !! Enables the calculation of off-diagonal elements of the GF.(0=SCF/DOS; 1=cluster; -1=custom)
    integer, intent (in) :: naez   !! Number of atoms in unit cell
    integer, intent (in) :: ncpa   !! NCPA = 0/1 CPA flag
    integer, intent (in) :: ncls   !! Number of reference clusters
    integer, intent (in) :: kmrot  !! 0: no rotation of the magnetisation; 1: individual rotation of the magnetisation for every site
    integer, intent (in) :: natyp  !! Number of kinds of atoms in unit cell
    integer, intent (in) :: nofks
    integer, intent (in) :: ideci
    integer, intent (in) :: ispin
    integer, intent (in) :: nspin  !! Counter for spin directions
    integer, intent (in) :: nsymat
    integer, intent (in) :: invmod !! Inversion scheme
    integer, intent (in) :: iprint
    integer, intent (in) :: nlbasis !! Number of basis layers of left host (repeated units)
    integer, intent (in) :: nrbasis !! Number of basis layers of right host (repeated units)
    integer, intent (in) :: natomimp !! Size of the cluster for impurity-calculation output of GF should be 1, if you don't do such a calculation
    integer, intent (in) :: naclsmax
    integer, intent (in) :: iqdosrun !! qdos ruess: counts qdos run
    real (kind=dp), intent (in) :: alat !! Lattice constant in a.u.
    real (kind=dp), intent (in) :: volbz
    real (kind=dp), intent (in) :: cpatol !! Convergency tolerance for CPA-cycle
    complex (kind=dp), intent (in) :: eryd
    complex (kind=dp), intent (in) :: cfctor
    complex (kind=dp), intent (in) :: cfctorinv
    ! .. Input arrays
    integer, dimension (nembd2), intent (in) :: cls !! Cluster around atomic sites
    integer, dimension (naez), intent (in) :: noq !! Number of diff. atom types located
    integer, dimension (nsheld), intent (in) :: nsh1 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (nsheld), intent (in) :: nsh2 !! Corresponding index of the sites I/J in  (NSH1/2) in the unit cell in a shell
    integer, dimension (naez), intent (in) :: icpa !! ICPA = 0/1 site-dependent CPA flag
    integer, dimension (natyp), intent (in) :: iqat !! The site on which an atom is located on a given site
    integer, dimension (nclsd), intent (in) :: nacls !! Number of atoms in cluster
    integer, dimension (0:nsheld), intent (in) :: nshell !! Index of atoms/pairs per shell (ij-pairs); nshell(0) = number of shells
    integer, dimension (nembd2), intent (in) :: refpot !! Ref. pot. card  at position
    integer, dimension (nofgij), intent (in) :: ijtabsh !! Linear pointer, assigns pair (i,j) to a shell in the array GS(*,*,*,NSHELD)
    integer, dimension (nofgij), intent (in) :: ijtabsym !! Linear pointer, assigns pair (i,j) to the rotation bringing GS into Gij
    integer, dimension (naclsd, nembd2), intent (in) :: atom !! Atom at site in cluster
    integer, dimension (naclsd, nembd2), intent (in) :: ezoa !! EZ of atom at site in cluster
    integer, dimension (natyp, naez), intent (in) :: itoq
    integer, dimension (2, lmmaxd), intent (in) :: nrrel
    integer, dimension (naez/nprincd, naez/nprincd), intent (in) :: icheck
    integer, dimension (2, 2, lmmaxd), intent (in) :: irrel
    real (kind=dp), dimension (natyp), intent (in) :: conc !! Concentration of a given atom
    real (kind=dp), dimension (kpoibz), intent (in) :: volcub
    real (kind=dp), dimension (3, 0:nrd), intent (in) :: rr !! Set of real space vectors (in a.u.)
    real (kind=dp), dimension (3, kpoibz), intent (in) :: bzkp
    real (kind=dp), dimension (3, nsheld), intent (in) :: ratom
    real (kind=dp), dimension (3, nembd2), intent (in) :: rbasis !! Position of atoms in the unit cell in units of bravais vectors
    real (kind=dp), dimension (3, naclsd, nclsd), intent (in) :: rcls !! Real space position of atom in cluster
    real (kind=dp), dimension (48, 3, nsheld), intent (in) :: rrot
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rc !! NREL REAL spher. harm. > CMPLX. spher. harm. NREL CMPLX. spher. harm. > REAL spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: crel !! Non-relat. CMPLX. spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. CMPLX. spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: rrel !! Non-relat. REAL spher. harm. > (kappa,mue) (kappa,mue)  > non-relat. REAL spher. harm.
    complex (kind=dp), dimension (lmmaxd, lmmaxd), intent (in) :: factl
    complex (kind=dp), dimension (lmmaxd, lmmaxd, natyp), intent (in) :: tsst
    complex (kind=dp), dimension (lmmaxd, lmmaxd, natyp), intent (in) :: msst
    complex (kind=dp), dimension (2, 2, lmmaxd), intent (in) :: srrel
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: tqdos ! qdos : Read-in inverse t-matrix
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: drotq !! Rotation matrices to change between LOCAL/GLOBAL frame of reference for magnetisation <> Oz or noncollinearity
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nrefd), intent (in) :: trefll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsymaxd), intent (in) :: dsymll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nrefd), intent (in) :: dtrefll !! LLY Lloyd dtref/dE
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez), intent (in) :: dtmatll ! LLY  dt/dE (should be av.-tmatrix in CPA)
    complex (kind=dp), dimension (lmgf0d*naclsmax, lmgf0d, ncls), intent (in) :: ginp !! Cluster GF (ref syst.)
    complex (kind=dp), dimension (lmgf0d*naclsmax, lmgf0d, ncls), intent (in) :: dginp !! LLY Lloyd Energy derivative of GINP
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspin), intent (in) :: lefttinvll
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nembd1, nspin), intent (in) :: righttinvll
    logical, dimension (2), intent (in) :: vacflag
    logical, dimension (nsymaxd), intent (in) :: symunitary !! unitary/antiunitary symmetry flag
    ! .. Output variables
    integer, intent (out) :: icpaflag
    ! .. In/Out variables
    integer, intent (inout) :: itcpamax !! Max. number of CPA iterations
    complex (kind=dp), intent (inout) :: tracet !! \f$Tr\left[ (t-tref)^{-1} \frac{d(t-tref)}{dE} \right]\f$
    complex (kind=dp), intent (inout) :: lly_grtr !! Trace Eq.5.38 PhD Thiess (k-integrated)! LLY Lloyd
    complex (kind=dp), dimension (lmmaxd, lmmaxd, nsheld), intent (inout) :: gmatll !! GMATLL = diagonal elements of the G matrix (system)
    ! .. Local Scalars
    integer :: i_stat, i_all
    integer :: ih, lm1, lm2, ns, nsdia, icall, irec
    integer :: iq, jq, it, i, j, iqtau, icpastart, itcpa, iu, nsmax
    real (kind=dp) :: cpaerrl, cpaerr, cpacorr, cpachng
    complex (kind=dp) :: ez, cnsymat, tauvbz
    logical :: ldia
    character (len=4) :: str4
    character (len=10) :: str10
    ! .. Local Arrays ..
    integer, dimension (naez) :: nkmq
    integer, dimension (naez) :: nlinq
    integer, dimension (nsymaxd) :: isumg
    integer, dimension (linmax) :: ikm1lin
    integer, dimension (linmax) :: ikm2lin
    integer, dimension (nsymaxd, naez) :: isumq
    integer, dimension (nsymaxd, natyp) :: isumt
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: gll
    complex (kind=dp), dimension (linmax, natyp) :: tautlin
    ! .. Effective (site-dependent) Delta_t^(-1) matrix ..
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: xc
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: w1
    complex (kind=dp), dimension (lmmaxd, lmmaxd) :: w2
    complex (kind=dp), dimension (lmmaxd, lmmaxd, naez) :: mssq
    ! .. Local allocatable arrays
    complex (kind=dp), dimension (:, :, :), allocatable :: dmssq
    complex (kind=dp), dimension (:, :, :), allocatable :: taudelq
    complex (kind=dp), dimension (:, :, :), allocatable :: taudelt
    complex (kind=dp), dimension (:, :, :, :), allocatable :: gs

#ifdef CPP_MPI
    integer :: irec0
#endif
    ! .. External Functions
    logical :: test, opt
    external :: test, opt
    ! .. Data statement
    data icall/0/
    ! .. Save statement
    save :: icall, isumg, cnsymat

    ! ----------------------------------------------------------------------------
    if (test('flow    ')) then
      write (1337, *) '>>> KLOOPZ1: invert delta_t and do Fourier transformation'
    end if
    icall = icall + 1

    ! Reinitialise the ICALL counter for second run of kloopz1                    ! qdos ruess
    if ((opt('qdos    ') .and. (iqdosrun==1))) icall = 1 ! qdos ruess

    if (opt('FERMIOUT') .and. myrank==master) then ! fswrt
      write (6801, '(A)') 'energy(ie):' ! fswrt
      write (6801, '(2ES25.16)') eryd ! fswrt
    end if                         ! fswrt
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The arrays ISUM are used in the symmetrisation routine SYMETRMAT
    ! Symmetrising single-site : same matrix for each symmetry
    ! Symmetrising G matrix    : pick G(ISYM) for symmetry ISYM
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if ( ICALL.EQ.1 ) then
    cnsymat = cone/dble(nsymat)
    do iu = 1, nsymaxd
      do it = 1, natyp
        isumt(iu, it) = it
      end do
      do iq = 1, naez
        isumq(iu, iq) = iq
      end do
      isumg(iu) = iu
    end do
    ! end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! TSST in the LOCAL frame is used to set up
    ! MSST = (TSST-TREF)^(-1) in the LOCAL frame to be used in <CPAMILLSX> and < PROJTAU > below

    ! MSSQ = the inverse of the effective (on-site) Delta_t matrix in the GLOBAL frame;
    ! the Average T-matrix Approximation (ATA) (ICPASTART=1) is used
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    icpastart = 1
    call mssinit(ncpa, icpastart, tsst, msst, mssq, trefll, drotq, refpot, iqat, itoq, noq, conc, kmrot, natyp, naez) ! nref was taken out of calling list 1.2.2012

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! VIRTUAL ATOMS:
    ! Be careful! in case of  OPT('VIRATOMS')==1 MSSQ is the Tmatrix
    ! not the inverse T-matrix!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ----------------------------------------------------------------------------
    ! Output now:
    ! MSSQ is the (Delta_t)^(-1) matrix in the GLOBAL frame
    ! and refers to a site IQ occupied by NOQ(IQ) atoms having
    ! the occupancies CONC(ITOQ(1..NOQ(IQ))

    ! TSST is the t-matrix in the LOCAL frame
    ! for each of the NATYP atoms
    ! MSST is the Delta_t^(-1) matrix in the LOCAL frame
    ! for each of the NATYP atoms
    ! ----------------------------------------------------------------------------

    ! Write gref to the TBkkr-container-file
    if (opt('FERMIOUT') .and. myrank==master) then ! fswrt
      write (6801, '(A)') 'GINP(ie):' ! fswrt
      do i = 1, ncls               ! fswrt
        do lm2 = 1, lmgf0d         ! fswrt
          do lm1 = 1, naclsmax*lmgf0d ! fswrt
            write (6801, '(2ES25.16)') ginp(lm1, lm2, i) ! fswrt
          end do                   ! fswrt
        end do                     ! fswrt
      end do                       ! fswrt
    end if                         ! fswrt

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BEGIN CPA - LOOP  (if required)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! ----------------------------------------------------------------------------
    ! ikm1lin,ikm2lin,nlinq --> dummy settings for <PROJTAU>
    ! ----------------------------------------------------------------------------
    ikm1lin(1) = 1
    ikm2lin(1) = 1
    do iq = 1, naez
      nkmq(iq) = lmmaxd
      nlinq(iq) = 1
    end do

    cpaerrl = 1.0d+6
    itcpa = 0
    ez = eryd

    if (ncpa/=0) then
      allocate (dmssq(lmmaxd,lmmaxd,naez), stat=i_stat)
      call memocc(i_stat, product(shape(dmssq))*kind(dmssq), 'DMSSQ', 'kloopz1')
      dmssq(:, :, :) = czero
      call cinit(lmmaxd*lmmaxd*naez, dmssq)
    end if
    ! ----------------------------------------------------------------------------
    !! GLL2K >  incorporated now
    ! ----------------------------------------------------------------------------
    tauvbz = 1.d0/volbz
    ! ----------------------------------------------------------------------------
    ! Convert inverted delta_t-matrices to p.u.
    ! ----------------------------------------------------------------------------
    do i = 1, naez

      if (.not. opt('VIRATOMS')) then
        call zscal(lmmaxd*lmmaxd, cfctor, mssq(1,1,i), 1)
      else
        call zscal(lmmaxd*lmmaxd, cfctorinv, mssq(1,1,i), 1)
      end if                       ! ( .not. OPT('VIRATOMS') )

      if (.not. opt('NEWSOSOL')) then
        if (kmrot==0) then
          do lm2 = 1, lmmaxd
            do lm1 = 1, lm2
              mssq(lm1, lm2, i) = 0.5d0*(mssq(lm1,lm2,i)+mssq(lm2,lm1,i))
              mssq(lm2, lm1, i) = mssq(lm1, lm2, i)
            end do
          end do
        end if
      end if
    end do

    do i = 1, natyp
      call zscal(lmmaxd*lmmaxd, cfctor, msst(1,1,i), 1)
    end do
    ! ----------------------------------------------------------------------------
    ! Symmetrise the delta_t^(-1) matrices
    ! ----------------------------------------------------------------------------
    if (.not. opt('NEWSOSOL')) then
      if ((krel==1) .or. (ins/=0)) then
        do iq = 1, naez
          call symetrmat(nsymat, cnsymat, dsymll, symunitary, mssq, isumq(1,iq), mssq(1,1,iq), lmmaxd, nsymaxd)
          if (kmrot==0) then
            do i = 1, noq(iq)
              it = itoq(i, iq)
              call symetrmat(nsymat, cnsymat, dsymll, symunitary, msst, isumt(1,it), msst(1,1,it), lmmaxd, nsymaxd)
            end do
          end if
        end do
      end if
    end if

    ez = ez*cfctor*cfctor
    nsdia = max(naez, natyp)

    allocate (taudelq(lmmaxd,lmmaxd,nshell(0)), stat=i_stat)
    call memocc(i_stat, product(shape(taudelq))*kind(taudelq), 'TAUDELQ', 'kloopz1')
    taudelq(:, :, :) = czero
    ! ----------------------------------------------------------------------------
    ! Jonathan Chico: This possibly can be removed by a do while
100 continue
    ! ----------------------------------------------------------------------------
    itcpa = itcpa + 1

    allocate (gs(lmmaxd,lmmaxd,nsymaxd,nshell(0)), stat=i_stat)
    call memocc(i_stat, product(shape(gs))*kind(gs), 'GS', 'kloopz1')
    gs(:, :, :, :) = czero

    ! copy read-in cpa t-matrix but only after fort.37 was created in first run   !qdos ruess
    if (opt('readcpa ') .or. (opt('qdos    ') .and. (iqdosrun==1))) then ! qdos ruess
      mssq(:, :, :) = tqdos(:, :, :) ! lmmaxd,lmmaxd,naez                            !qdos ruess
    end if

    call kkrmat01(bzkp, nofks, gs, volcub, mssq, rrot, nshell(0), nsdia, alat, nsymat, naez, cls, nacls, naclsmax, rr, ezoa, atom, nsh1, nsh2, ginp, rbasis, rcls, &
      lefttinvll(1,1,1,ispin), righttinvll(1,1,1,ispin), vacflag, nlbasis, nrbasis, factl, icheck, invmod, ideci, srrel, irrel, nrrel, dtrefll, dtmatll, dginp, refpot, lly_grtr, &
      tracet, cfctor, lly)         ! LLY

    nsmax = nshell(0)
    ! ----------------------------------------------------------------------------
    ! NS=1,NSMAX
    ! ----------------------------------------------------------------------------
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! For qdos calculation, do not symmetrize the GF matrix, but call routine
    ! symetrmat anyway because of factor tauvbz:
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ns = 1, nsmax
      ! -------------------------------------------------------------------------
      ! symmetrise GS, get GLL as sum over all symmetry wedges of GS
      ! (isumg(i) = i, see above)
      ! -------------------------------------------------------------------------
      call symetrmat(nsymat, tauvbz, dsymll, symunitary, gs(1,1,1,ns), isumg, gll, lmmaxd, nsymaxd)

      if (ns<=natyp) then
        iqtau = iqat(ns)
      else
        iqtau = ns
      end if

      taudelq(1:lmmaxd, 1:lmmaxd, iqtau) = -gll(1:lmmaxd, 1:lmmaxd)
    end do
    ! ----------------------------------------------------------------------------
    ! NS=1,NSMAX
    ! ----------------------------------------------------------------------------
    i_all = -product(shape(gs))*kind(gs)
    deallocate (gs, stat=i_stat)
    call memocc(i_stat, i_all, 'GS', 'kloopz1')
    ! ----------------------------------------------------------------------------
    if (ncpa>0) then
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! do one CPA iteration, the output is a new MSSQ = (Delta_t)^(-1)
      ! in the GLOBAL frame
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      icpaflag = 0
      if (opt('readcpa ') .or. (opt('qdos    ') .and. iqdosrun>0)) then
        ! copy read-in cpa t-matrix in second run
        call cpamillsx(itcpa, cpaerr, cpacorr, cpachng, iprint, icpa, naez, nkmq, noq, itoq, conc, mssq, msst, taudelq, dmssq, kmrot, drotq, natyp, naez, lmmaxd)
        mssq(:, :, :) = tqdos(:, :, :) ! lmmaxd,lmmaxd,naez   ! qdos
        itcpamax = 0
        cpaerr = 0.d0
      else
        call cpamillsx(itcpa, cpaerr, cpacorr, cpachng, iprint, icpa, naez, nkmq, noq, itoq, conc, mssq, msst, taudelq, dmssq, kmrot, drotq, natyp, naez, lmmaxd)
        ! ----------------------------------------------------------------------
        ! Symmetrise m-CPA
        ! ----------------------------------------------------------------------
        if (.not. opt('NEWSOSOL')) then
          do iq = 1, naez
            if (icpa(iq)/=0) then
              call symetrmat(nsymat, cnsymat, dsymll, symunitary, mssq, isumq(1,iq), mssq(1,1,iq), lmmaxd, nsymaxd)
            end if
          end do
        end if

        if (iprint>=1 .and. (t_inc%i_write>0)) then
          write (1337, 140) cpaerr, cpacorr, cpachng
        end if

        if (cpaerr<=cpatol) then
          if (iprint>0 .and. (t_inc%i_write>0)) then
            write (1337, 110) itcpa, cpaerr, cpacorr, cpachng
          end if
        else if (itcpa>itcpamax) then
          if (t_inc%i_write>0) then
            write (1337, 120) itcpa, cpaerr, cpacorr, cpachng
          end if
          icpaflag = 1
        else if (cpaerr>20*cpaerrl) then
          if (t_inc%i_write>0) write (1337, 130) itcpa
          if (t_inc%i_write>0) write (1337, 140) cpaerr, cpacorr, cpachng
          icpaflag = 2
        else
          ! -------------------------------------------------------------------
          ! Go to the next CPA iteration if not converged
          ! -------------------------------------------------------------------
          cpaerrl = cpaerr
          go to 100
        end if                     ! ( CPAERR.LE.CPATOL )
      end if                       ! (OPT('readcpa ').OR.OPT('qdos    '))
      i_all = -product(shape(dmssq))*kind(dmssq)
      deallocate (dmssq, stat=i_stat)
      call memocc(i_stat, i_all, 'DMSSQ', 'kloopz1')

    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END CPA - LOOP
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! ----------------------------------------------------------------------------
    ! The inverse of the Delta_t(CPA)-matrix is calculated, now write
    ! this on a file in case of decimation output.
    ! attention: new format Dec. 2004
    ! ----------------------------------------------------------------------------

    ! In first qdos run write out t matrix which is read in for the calculation for every k point
    if (opt('deci-out') .or. (iqdosrun==0)) then ! qdos ruess
      do ih = 1, naez
#ifdef CPP_MPI
        irec0 = lmmaxd**2*(ih-1) + lmmaxd**2*naez*(ie-1) + lmmaxd**2*t_inc%ielast*naez*(ispin-1)
        if (opt('deci-out')) write (37, 150) ie, eryd, ih
#else
        write (37, 150) ie, eryd, ih
#endif
        do lm1 = 1, lmmaxd
          do lm2 = 1, lmmaxd
            if (lm1==lm2) then
#ifdef CPP_MPI
              if (.not. opt('deci-out')) then
                irec = irec0 + lm2 + lmmaxd*(lm1-1)
                write (37, rec=irec) mssq(lm1, lm2, ih)*cfctorinv
              else
                write (37, 160) lm1, lm2, mssq(lm1, lm2, ih)*cfctorinv
              end if
#else
              write (37, 160) lm1, lm2, mssq(lm1, lm2, ih)*cfctorinv
#endif
            else
#ifdef CPP_MPI
              irec = irec0 + lm2 + lmmaxd*(lm1-1)
              if (abs(mssq(lm1,lm2,ih)/mssq(lm1,lm1,ih))>tolmssq) then
                if (.not. opt('deci-out')) then
                  write (37, rec=irec) mssq(lm1, lm2, ih)*cfctorinv
                else
                  write (37, 160) lm1, lm2, mssq(lm1, lm2, ih)*cfctorinv
                end if
              else
                if (.not. opt('deci-out')) then
                  write (37, rec=irec)(0.d0, 0.d0)
                end if
              end if
#else
              if (abs(mssq(lm1,lm2,ih)/mssq(lm1,lm1,ih))>tolmssq) then
                write (37, 160) lm1, lm2, mssq(lm1, lm2, ih)*cfctorinv
              end if
#endif
            end if
          end do
        end do
      end do
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! calculate the component-projected site-diagonal
    ! TAU-matrices TAUDELT(IT).
    ! - there are NSMAX = NAEZ/NATYP (NCPA=0/1) site-diagonal
    ! elements TAUDELQ(1..NSMAX) in the GLOBAL (crystal) frame of
    ! reference
    ! - they are NSMAX site-diagonal elements projected on atomic types
    ! in the array TAUDELT(1..NSMAX) in the LOCAL frame of reference
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate (taudelt(lmmaxd,lmmaxd,nshell(0)), stat=i_stat)
    call memocc(i_stat, product(shape(taudelt))*kind(taudelt), 'TAUDELT', 'kloopz1')
    taudelt(:, :, :) = czero
    if (ncpa>0) then
      nsmax = natyp
      call projtau(icpaflag, cpachng, kmrot, .false., .false., 9, czero, natyp, naez, nkmq, msst, mssq, nlinq, iqat, conc, taudelq, taudelt, tautlin, ikm1lin, ikm2lin, drotq, &
        nsheld, naez, lmmaxd, linmax)
    else
      nsmax = naez
      do ns = 1, nsmax
        call zcopy(lmmaxd*lmmaxd, taudelq(1,1,ns), 1, taudelt(1,1,ns), 1)
        if (kmrot/=0) then
          do j = 1, lmmaxd
            call zcopy(lmmaxd, taudelt(1,j,ns), 1, w1(1,j), 1)
          end do
          call rotate(w1, 'G->L', taudelt(1,1,ns), lmmaxd, drotq(1,1,ns), lmmaxd)
        end if
      end do
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! off-diagonal elements (if present) are in the same array
    ! TAUDELQ in the range (NSMAX+1,..,NSHELL(0)).
    ! The G_ij's are left UNPROJECTED and in the GLOBAL frame of
    ! reference (remember this for further use!), i.e. in case of CPA,
    ! G_ij(CPA) is stored. The component-projected G-elements can
    ! be obtained through the projection matrices (see below)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do ns = nsmax + 1, nshell(0)
      call zcopy(lmmaxd*lmmaxd, taudelq(1,1,ns), 1, taudelt(1,1,ns), 1)
    end do
    ! ----------------------------------------------------------------------------
    ! Switch from TAU to G multiplying with MSST on the site-diagonal,
    ! system-matrix (1..NSMAX) and with MSSQ on the rest of the elements
    ! ----------------------------------------------------------------------------

    if (test('Gmat    ') .and. (t_inc%i_write>0)) then
      write (1337, '(/,4X,70("-"),/,4X,A,I4)') 'system G_ii matrix for i = 1,', nsmax
    end if
    do ns = 1, nshell(0)

      ldia = (abs(ratom(1,ns)**2+ratom(2,ns)**2+ratom(3,ns)**2)<1d-6)
      ! -------------------------------------------------------------------------
      ! GLL = -TAU
      ! -------------------------------------------------------------------------
      call zcopy(lmmaxd*lmmaxd, taudelt(1,1,ns), 1, gll, 1)
      call zscal(lmmaxd*lmmaxd, -cone, gll, 1)
      ! -------------------------------------------------------------------------
      ! XC = (Delta_t(I))^(-1) * GLL
      ! -------------------------------------------------------------------------
      if (ns<=nsmax) then
        ! ----------------------------------------------------------------------
        ! deal with the site-diagonal GFUN of the system - use MSST(I)
        ! on both sides of TAU multiplication
        ! ----------------------------------------------------------------------
        do j = 1, lmmaxd
          call zcopy(lmmaxd, msst(1,j,ns), 1, w1(1,j), 1)
        end do
        call zcopy(lmmaxd*lmmaxd, w1, 1, w2, 1)
      else
        ! ----------------------------------------------------------------------
        ! deal with the Gij of the cluster - use MSSQ(IQ) and MSSQ(JQ)
        ! ----------------------------------------------------------------------
        iq = nsh1(ns)
        jq = nsh2(ns)
        do j = 1, lmmaxd
          call zcopy(lmmaxd, mssq(1,j,iq), 1, w1(1,j), 1)
        end do
        do j = 1, lmmaxd
          call zcopy(lmmaxd, mssq(1,j,jq), 1, w2(1,j), 1)
        end do
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! NS.LT.NSMAX
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (.not. opt('VIRATOMS')) then
        if (.not. test('testgmat')) then
          call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, cone, w1, lmmaxd, gll, lmmaxd, czero, xc, lmmaxd)

          if (ldia) then
            ! ----------------------------------------------------------------
            ! GLL = - (Delta_t(I))^(-1) - (Delta_t(I))^(-1) * GLL * (Delta_t(I))^(-1)
            ! ----------------------------------------------------------------
            call zcopy(lmmaxd*lmmaxd, w2, 1, gll, 1)
            call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, -cone, xc, lmmaxd, w2, lmmaxd, -cone, gll, lmmaxd)
          else
            ! ----------------------------------------------------------------
            ! GLL =  - (Delta_t(I))^(-1)  * GLL * (Delta_t(J))^(-1)
            ! ----------------------------------------------------------------
            call zgemm('N', 'N', lmmaxd, lmmaxd, lmmaxd, -cone, xc, lmmaxd, w2, lmmaxd, czero, gll, lmmaxd)
          end if
        end if                     ! ( .not. TEST('testgmat') ) THEN
      end if                       ! ( .not. OPT('VIRATOMS') ) THEN
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! LDIA
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! -------------------------------------------------------------------------
      ! GMATLL = GLL/RFCTOR
      ! -------------------------------------------------------------------------
      gmatll(1:lmmaxd, 1:lmmaxd, ns) = gll(1:lmmaxd, 1:lmmaxd)*cfctorinv

      if ((ns<=nsmax) .and. (test('Gmat    '))) then
        write (str4, '(I4)') ns
        str10 = '   i =' // str4(1:4)
        call cmatstr(str10, 10, gmatll(1,1,ns), lmmaxd, lmmaxd, 2*krel+1, 2*krel+1, 0, 1.0e-8_dp, 6)
      end if
    end do
    ! ----------------------------------------------------------------------------
    if (test('Gmat    ') .and. (t_inc%i_write>0)) write (1337, '(/,4X,70("-"))')
    ! ----------------------------------------------------------------------------
    ! it calculates the rest of the G n n' matrix from the
    ! knowledge of the representative pairs (shells) using the
    ! real space symmetries (added 23.2.2000)
    ! ----------------------------------------------------------------------------
    if (icc>0) then
      irec = 1 + ie + t_inc%ielast*(ispin-1) ! added for mpi run
      call rotgll(gmatll, natomimp, ijtabsym, ijtabsh, dsymll, symunitary, igf, rc, crel, rrel, krel, lmmaxd, irec)
    end if
    ! IF ( OPT('VIRATOMS') ) THEN
    ! write(*,*) 'VIRTUAL ATOM OPTION : stop calculation '
    ! stop
    ! END IF
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! In the case of NCPA.NE.0 and NSHELL(0).GT.NATYP the projection
    ! matrices DMAT and DTIL which are used to get
    ! ij            ij    _
    ! G     =  D  * G    * D
    ! ab       a    CPA    b
    ! - with a/b the atom of type a/b sitting on site i/j - are calculated
    ! and stored for later use.  the allocated work space for
    ! TSST (DMAT) and MSST (DTIL) is used.
    ! for an atom having occupancy 1, DMAT/DTIL = unit matrix
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((ncpa/=0) .and. (nshell(0)>nsmax)) then
      do it = 1, natyp
        iq = iqat(it)
        ih = refpot(iq)

        if (kmrot/=0) then
          call rotate(tsst(1,1,it), 'L->G', w1, lmmaxd, drotq(1,1,iq), lmmaxd)
        else
          call zcopy(lmmaxd*lmmaxd, tsst(1,1,it), 1, w1, 1)
        end if

        do j = 1, lmmaxd
          call zaxpy(lmmaxd, -cone, trefll(1,j,ih), 1, w1(1,j), 1)
        end do

        call gijdmat(taudelq(1,1,iq), w1, mssq(1,1,iq), tsst(1,1,it), msst(1,1,it), cfctorinv, iprint, ie, it, krel, lmmaxd)
      end do
    end if
    ! ----------------------------------------------------------------------------
    if (test('flow    ') .and. (t_inc%i_write>0)) write (1337, *) '<<< KLOOPZ1'

    i_all = -product(shape(taudelq))*kind(taudelq)
    deallocate (taudelq, stat=i_stat)
    call memocc(i_stat, i_all, 'TAUDELQ', 'kloopz1')
    i_all = -product(shape(taudelt))*kind(taudelt)
    deallocate (taudelt, stat=i_stat)
    call memocc(i_stat, i_all, 'TAUDELT', 'kloopz1')

    return

110 format (' CPA converged after', i3, ' iterations   ', 'ERR:', f9.6, ' CORR:', f9.6, ' CHNG:', f9.6)
120 format (10x, 10('!'), ' CPA-cycle  NOT  converged after ', i3, ' iterations ', 10('!'), /, 14x, 'ERROR ', f12.8, ' CORRECTION ', f15.8, ' CHANGE ', f15.8)
130 format (' CPA: ERROR increased by more than ', '20*TOL for ITCPA=', i4, ' >>>> iteration stopped ', 5x, 10('!'))
140 format (' CPA: ERROR ', f12.8, '    CORRECTION ', f15.8, '    CHANGE ', f15.8)

150 format (/, 80('*'), /, 'ENERGY ', i5, 2d16.8, ' SITE ', i3)
160 format (2i5, 1p, 2d22.14)

  end subroutine kloopz1_qdos

end module mod_kloopz1
