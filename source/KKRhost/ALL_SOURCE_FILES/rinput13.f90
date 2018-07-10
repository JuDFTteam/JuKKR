module rinput

  implicit none

contains
  ! -------------------------------------------------------------------------------
  ! SUBROUTINE: RINPUT13
  ! > @brief Routine to read the information from the input file
  ! > @author Bernd Zimmermann
  ! > @note Jonathan Chico: Added calls for the different parameters options,
  ! > that are being changed from the inc.p to now the reader. Of this way if
  ! there
  ! > is a problem the user only has to add a line or change a value of the
  ! inc.p
  ! > not recompile the code
  ! -------------------------------------------------------------------------------
  subroutine rinput13(kte, igf, kxc, lly, icc, ins, kws, ipe, ipf, ipfe, &
    icst, imix, lpot, naez, nemb, nref, ncls, npol, lmax, kcor, kefg, khyp, &
    kpre, kvmad, lmmax, lmpot, ncheb, nleft, ifile, kvrel, nspin, natyp, &
    nineq, npnt1, npnt2, npnt3, kfrozn, ishift, n1semi, n2semi, n3semi, &
    nsteps, insref, kshape, itdbry, nright, kforce, ivshift, khfield, nlbasis, &
    nrbasis, intervx, intervy, intervz, npan_eq, npan_log, npolsemi, tk, fcm, &
    emin, emax, rmax, gmax, alat, r_log, rcutz, rcutxy, eshift, qbound, &
    hfield, mixing, abasis, bbasis, cbasis, vconst, tksemi, tolrdif, emusemi, &
    ebotsemi, fsemicore, lambda_xc, deltae, lrhosym, linipol, lcartesian, imt, &
    cls, lmxc, irns, irws, ntcell, refpot, inipol, ixipol, hostimp, kfg, vbc, &
    zperleft, zperight, bravais, rmt, zat, rws, mtfac, rmtref, rmtnew, &
    rmtrefat, fpradius, tleft, tright, rbasis, socscale, cscl, socscl, solver, &
    i12, i13, i19, i25, i40, txc, drotq, ncpa, itcpamax, cpatol, noq, iqat, &
    icpa, kaoez, conc, kmrot, qmtet, qmphi, kreadldau, lopt, ueff, jeff, &
    erefldau)

    use :: profiling
    use :: constants
    use :: mod_wunfiles, only: t_params
    use :: memoryhandling
    use :: global_variables
    use :: mod_types, only: t_inc
    use :: mod_save_wavefun, only: t_wavefunctions
    use :: mod_version_info
    use :: mod_datatypes
    use :: godfrin, only: t_godfrin ! GODFRIN Flaviano

    implicit none
    ! ..
    ! > @note VP : there should be some crosscheck of competing options
    ! >            e.g., XCPL and CONDUCT cannot be done simultaneously
    ! >            neither SOC1 and SOC2 manipulation etc.
    ! ..
    ! .. External Subroutines ..
    external :: rcstop
    ! ..
    ! .. Intrinsic Functions ..
    intrinsic :: min
    ! ..
    ! .. Scalar Arguments ..
    integer, intent (inout) :: kte ! < Calculation of the total energy On/Off
                                   ! (1/0)
    integer, intent (inout) :: igf ! < Do not print or print (0/1) the
                                   ! KKRFLEX_* files
    integer, intent (inout) :: kxc ! < Type of xc-potential 0=vBH 1=MJW 2=VWN
                                   ! 3=PW91
    integer, intent (inout) :: lly ! < LLY <> 0 : apply Lloyds formula
    integer, intent (inout) :: icc ! < Enables the calculation of off-diagonal
                                   ! elements of the GF.(0=SCF/DOS; 1=cluster;
                                   ! -1=custom)
    integer, intent (inout) :: ins ! < 0 (MT), 1(ASA), 2(Full Potential)
    integer, intent (inout) :: kws ! < 0 (MT), 1(ASA)
    integer, intent (inout) :: ipe ! < Not real used, IPFE should be 0
    integer, intent (inout) :: ipf ! < Not real used, IPFE should be 0
    integer, intent (inout) :: ipfe ! < Not real used, IPFE should be 0
    integer, intent (inout) :: icst ! < Number of Born approximation
    integer, intent (inout) :: imix ! < Type of mixing scheme used
                                    ! (0=straight, 4=Broyden 2nd, 5=Anderson)
    integer, intent (inout) :: lpot ! < Maximum l component in potential
                                    ! expansion
    integer, intent (inout) :: naez ! < Number of atoms in unit cell
    integer, intent (inout) :: nemb ! < Number of 'embedding' positions
    integer, intent (inout) :: nref ! < Number of diff. ref. potentials
    integer, intent (inout) :: ncls ! < Number of reference clusters
    integer, intent (inout) :: npol ! < Number of Matsubara Pols (EMESHT)
    integer, intent (inout) :: lmax ! < Maximum l component in wave function
                                    ! expansion
    integer, intent (inout) :: kcor
    integer, intent (inout) :: kefg
    integer, intent (inout) :: khyp
    integer, intent (inout) :: kpre
    integer, intent (inout) :: kvmad
    integer, intent (inout) :: lmmax
    integer, intent (inout) :: lmpot
    integer, intent (inout) :: ncheb ! < Number of Chebychev pannels for the
                                     ! new solver
    integer, intent (inout) :: nleft ! < Number of repeated basis for left
                                     ! host to get converged  electrostatic
                                     ! potentials
    integer, intent (inout) :: ifile ! < Unit specifier for potential card
    integer, intent (inout) :: kvrel ! < 0,1 : non / scalar relat. calculation
    integer, intent (inout) :: nspin ! < Counter for spin directions
    integer, intent (inout) :: natyp ! < Number of kinds of atoms in unit cell
    integer, intent (inout) :: nineq ! < Number of ineq. positions in unit
                                     ! cell
    integer, intent (inout) :: npnt1 ! < number of E points (EMESHT) for the
                                     ! contour integration
    integer, intent (inout) :: npnt2 ! < number of E points (EMESHT) for the
                                     ! contour integration
    integer, intent (inout) :: npnt3 ! < number of E points (EMESHT) for the
                                     ! contour integration
    integer, intent (inout) :: kfrozn
    integer, intent (inout) :: ishift
    integer, intent (inout) :: n1semi ! < Number of energy points for the
                                      ! semicore contour
    integer, intent (inout) :: n2semi ! < Number of energy points for the
                                      ! semicore contour
    integer, intent (inout) :: n3semi ! < Number of energy points for the
                                      ! semicore contour
    integer, intent (inout) :: nsteps ! < number of iterations
    integer, intent (inout) :: insref ! < INS for reference pot. (usual 0)
    integer, intent (inout) :: kshape ! < Exact treatment of WS cell
    integer, intent (inout) :: itdbry ! < Number of SCF steps to remember for
                                      ! the Broyden mixing
    integer, intent (inout) :: nright ! < Number of repeated basis for right
                                      ! host to get converged  electrostatic
                                      ! potentials
    integer, intent (inout) :: kforce ! < Calculation of the forces
    integer, intent (inout) :: ivshift
    integer, intent (inout) :: khfield ! < 0,1: no / yes external magnetic
                                       ! field
    integer, intent (inout) :: nlbasis ! < Number of basis layers of left host
                                       ! (repeated units)
    integer, intent (inout) :: nrbasis ! < Number of basis layers of right
                                       ! host (repeated units)
    integer, intent (inout) :: intervx ! < Number of intervals in x-direction
                                       ! for k-net in IB of the BZ
    integer, intent (inout) :: intervy ! < Number of intervals in y-direction
                                       ! for k-net in IB of the BZ
    integer, intent (inout) :: intervz ! < Number of intervals in z-direction
                                       ! for k-net in IB of the BZ
    integer, intent (inout) :: npan_eq ! < Number of intervals from [R_LOG] to
                                       ! muffin-tin radius Used in conjunction
                                       ! with runopt NEWSOSOL
    integer, intent (inout) :: npan_log ! < Number of intervals from nucleus
                                        ! to [R_LOG] Used in conjunction with
                                        ! runopt NEWSOSOL
    integer, intent (inout) :: npolsemi ! < Number of poles for the semicore
                                        ! contour
    real (kind=dp), intent (inout) :: tk ! < Temperature
    real (kind=dp), intent (inout) :: fcm ! < Factor for increased linear
                                          ! mixing of magnetic part of
                                          ! potential compared to non-magnetic
                                          ! part.
    real (kind=dp), intent (inout) :: emin ! < Lower value (in Ryd) for the
                                           ! energy contour
    real (kind=dp), intent (inout) :: emax ! < Maximum value (in Ryd) for the
                                           ! DOS calculation Controls also
                                           ! [NPT2] in some cases
    real (kind=dp), intent (inout) :: rmax ! < Ewald summation cutoff
                                           ! parameter for real space
                                           ! summation
    real (kind=dp), intent (inout) :: gmax ! < Ewald summation cutoff
                                           ! parameter for reciprocal space
                                           ! summation
    real (kind=dp), intent (inout) :: alat ! < Lattice constant (in a.u.)
    real (kind=dp), intent (inout) :: r_log ! < Radius up to which log-rule is
                                            ! used for interval width. Used in
                                            ! conjunction with runopt NEWSOSOL
    real (kind=dp), intent (inout) :: rcutz ! < Parameter for the screening
                                            ! cluster along the z-direction
    real (kind=dp), intent (inout) :: rcutxy ! < Parameter for the screening
                                             ! cluster along the x-y plane
    real (kind=dp), intent (inout) :: eshift
    real (kind=dp), intent (inout) :: qbound ! < Convergence parameter for the
                                             ! potential
    real (kind=dp), intent (inout) :: hfield ! < External magnetic field, for
                                             ! initial potential shift in spin
                                             ! polarised case
    real (kind=dp), intent (inout) :: mixing ! < Magnitude of the mixing
                                             ! parameter
    real (kind=dp), intent (inout) :: abasis ! < Scaling factors for rbasis
    real (kind=dp), intent (inout) :: bbasis ! < Scaling factors for rbasis
    real (kind=dp), intent (inout) :: cbasis ! < Scaling factors for rbasis
    real (kind=dp), intent (inout) :: vconst ! < Potential shift in the first
                                             ! iteration
    real (kind=dp), intent (inout) :: tksemi ! < Temperature for semi-core
                                             ! contour
    real (kind=dp), intent (inout) :: tolrdif ! < For distance between
                                              ! scattering-centers smaller
                                              ! than [<TOLRDIF>], free GF is
                                              ! set to zero. Units are Bohr
                                              ! radii.
    real (kind=dp), intent (inout) :: emusemi ! < Top of semicore contour in
                                              ! Ryd.
    real (kind=dp), intent (inout) :: ebotsemi ! < Bottom of semicore contour
                                               ! in Ryd
    real (kind=dp), intent (inout) :: fsemicore ! < Initial normalization
                                                ! factor for semicore states
                                                ! (approx. 1.)
    real (kind=dp), intent (inout) :: lambda_xc ! < Scale magnetic moment (0 <
                                                ! Lambda_XC < 1,0=zero moment,
                                                ! 1= full moment)
    complex (kind=dp), intent (inout) :: deltae ! < LLY Energy difference for
                                                ! numerical derivative
    logical, intent (inout) :: lrhosym
    logical, intent (inout) :: linipol ! < True: Initial spin polarization;
                                       ! false: no initial spin polarization
    logical, intent (inout) :: lcartesian ! < True: Basis in cartesian coords;
                                          ! false: in internal coords
    ! .. Array Arguments ..
    integer, dimension (:), allocatable, intent (out) :: imt ! < R point at MT
                                                             ! radius
    integer, dimension (:), allocatable, intent (out) :: cls ! < Cluster
                                                             ! around atomic
                                                             ! sites
    integer, dimension (:), allocatable, intent (out) :: lmxc
    integer, dimension (:), allocatable, intent (out) :: irns ! < Position of
                                                              ! atoms in the
                                                              ! unit cell in
                                                              ! units of
                                                              ! bravais
                                                              ! vectors
    integer, dimension (:), allocatable, intent (out) :: irws ! < R point at
                                                              ! WS radius
    integer, dimension (:), allocatable, intent (out) :: ntcell ! < Index for
                                                                ! WS cell
    integer, dimension (:), allocatable, intent (out) :: refpot ! < Ref. pot.
                                                                ! card  at
                                                                ! position
    integer, dimension (:), allocatable, intent (out) :: inipol ! < Initial
                                                                ! spin
                                                                ! polarisation
    integer, dimension (:), allocatable, intent (out) :: ixipol ! < Constraint
                                                                ! of spin pol.
    integer, dimension (:), allocatable, intent (out) :: hostimp
    integer, dimension (:, :), allocatable, intent (out) :: kfg
    real (kind=dp), dimension (2), intent (inout) :: vbc ! < Potential
                                                         ! constants
    real (kind=dp), dimension (3), intent (inout) :: zperleft ! < Vector to
                                                              ! define how to
                                                              ! repeat the
                                                              ! basis of the
                                                              ! left host
    real (kind=dp), dimension (3), intent (inout) :: zperight ! < Vector to
                                                              ! define how to
                                                              ! repeat the
                                                              ! basis of the
                                                              ! right host
    real (kind=dp), dimension (3, 3), intent (inout) :: bravais ! < Bravais
                                                                ! lattice
                                                                ! vectors
    real (kind=dp), dimension (:), allocatable, intent (out) :: rmt ! <
                                                                    ! Muffin-tin
                                                                    ! radius
                                                                    ! of true
                                                                    ! system
    real (kind=dp), dimension (:), allocatable, intent (out) :: zat ! <
                                                                    ! Nuclear
                                                                    ! charge
    real (kind=dp), dimension (:), allocatable, intent (out) :: rws ! < Wigner
                                                                    ! Seitz
                                                                    ! radius
    real (kind=dp), dimension (:), allocatable, intent (out) :: mtfac ! <
                                                                      ! Scaling
                                                                      ! factor
                                                                      ! for
                                                                      ! radius
                                                                      ! MT
    real (kind=dp), dimension (:), allocatable, intent (out) :: rmtref ! <
                                                                       ! Muffin-tin
                                                                       ! radius
                                                                       ! of
                                                                       ! reference
                                                                       ! system
    real (kind=dp), dimension (:), allocatable, intent (out) :: rmtnew ! <
                                                                       ! Adapted
                                                                       ! muffin-tin
                                                                       ! radius
    real (kind=dp), dimension (:), allocatable, intent (out) :: rmtrefat
    real (kind=dp), dimension (:), allocatable, intent (out) :: fpradius ! < R
                                                                         ! point
                                                                         ! at
                                                                         ! which
                                                                         ! full-potential
                                                                         ! treatment
                                                                         ! starts
    real (kind=dp), dimension (:, :), allocatable, intent (out) :: tleft ! <
                                                                         ! Vectors
                                                                         ! of
                                                                         ! the
                                                                         ! basis
                                                                         ! for
                                                                         ! the
                                                                         ! left
                                                                         ! host
    real (kind=dp), dimension (:, :), allocatable, intent (out) :: tright ! <
                                                                          ! vectors
                                                                          ! of
                                                                          ! the
                                                                          ! basis
                                                                          ! for
                                                                          ! the
                                                                          ! right
                                                                          ! host
    real (kind=dp), dimension (:, :), allocatable, intent (out) :: rbasis ! <
                                                                          ! Position
                                                                          ! of
                                                                          ! atoms
                                                                          ! in
                                                                          ! the
                                                                          ! unit
                                                                          ! cell
                                                                          ! in
                                                                          ! units
                                                                          ! of
                                                                          ! bravais
                                                                          ! vectors
    ! variables for spin-orbit/speed of light scaling
    real (kind=dp), dimension (:), allocatable, intent (out) :: socscale ! <
                                                                         ! Spin-orbit
                                                                         ! scaling
    real (kind=dp), dimension (:, :), allocatable, intent (out) :: cscl ! <
                                                                        ! Speed
                                                                        ! of
                                                                        ! light
                                                                        ! scaling
    real (kind=dp), dimension (:, :), allocatable, intent (out) :: socscl
    character (len=10), intent (inout) :: solver                           ! < Type of solver


    character (len=40), intent (inout) :: i12                           ! < File identifiers


    character (len=40), intent (inout) :: i13                           ! < Potential file name


    character (len=40), intent (inout) :: i19                           ! < Shape function file name


    character (len=40), intent (inout) :: i25                           ! < Scoef file name


    character (len=40), intent (inout) :: i40                           ! < File identifiers


    character (len=124), dimension (6), intent (inout) :: txc
    complex (kind=dp), dimension (:, :, :), allocatable, intent (out) :: drotq
    ! < Rotation matrices to change between LOCAL/GLOBAL frame of reference
    ! for magnetisation <> Oz or noncollinearity
    ! ----------------------------------------------------------------------------
    ! > @note CPA variables. Routine has been modified to look for
    ! >     the token ATOMINFOC and only afterwards, if not found, for the
    ! >     old token ATOMINFO. The only necessary extra information
    ! >     required is the site IQAT(IATOM) on which the atom IATOM
    ! >     is located and the occupancy (concentration) CONC(IATOM).
    ! >     The rest of CPA variables are deduced from these two.
    ! >     The tolerance for the CPA-cycle and the number of CPA iterations
    ! >     can be modified adding the token <CPAINFO> in the input file.
    ! ----------------------------------------------------------------------------
    integer, intent (inout) :: ncpa ! < ncpa = 0/1 CPA flag
    integer, intent (inout) :: itcpamax ! < max. number of CPA iterations
    real (kind=dp), intent (inout) :: cpatol ! < convergency tolerance for
                                             ! CPA-cycle
    integer, dimension (:), allocatable, intent (out) :: noq ! < number of
                                                             ! diff. atom
                                                             ! types located
    integer, dimension (:), allocatable, intent (out) :: iqat ! < the site on
                                                              ! which an atom
                                                              ! is located on
                                                              ! a given site
    integer, dimension (:), allocatable, intent (out) :: icpa ! < icpa = 0/1
                                                              ! site-dependent
                                                              ! CPA flag
    integer, dimension (:, :), allocatable, intent (out) :: kaoez ! < atom
                                                                  ! types
                                                                  ! located at
                                                                  ! a given
                                                                  ! site
    real (kind=dp), dimension (:), allocatable, intent (out) :: conc ! <
                                                                     ! concentration
                                                                     ! of a
                                                                     ! given
                                                                     ! atom

    ! ----------------------------------------------------------------------------
    ! > @note Variables storing the magnetization direction information.
    ! >     QMTET/QMPHI(NAEZ) give the angles to which the magnetic moment
    ! >     on a given site is rotated against the z-axis. Default values
    ! >     0.0 and 0.0, i.e., magnetic moment parallel to the z-axis.
    ! >     The angles are read in after the token RBASISANG is found
    ! >     (sought in input file prior to RBASIS token)
    ! >
    ! >   *  KMROT                                                           *
    ! >   *  0: no rotation of the magnetisation                             *
    ! >   *  1: individual rotation of the magnetisation for every site      *
    ! >   ( see also the routine < FINDGROUP > and ff)
    ! ----------------------------------------------------------------------------
    integer, intent (inout) :: kmrot ! < 0: no rotation of the magnetisation;
                                     ! 1: individual rotation of the
                                     ! magnetisation for every site
    real (kind=dp), dimension (:), allocatable, intent (out) :: qmtet ! < \f$
                                                                      ! \theta\f$
                                                                      ! angle
                                                                      ! of the
                                                                      ! agnetization
                                                                      ! with
                                                                      ! respect
                                                                      ! to the
                                                                      ! z-axis
    real (kind=dp), dimension (:), allocatable, intent (out) :: qmphi ! < \f$
                                                                      ! \phi\f$
                                                                      ! angle
                                                                      ! of the
                                                                      ! agnetization
                                                                      ! with
                                                                      ! respect
                                                                      ! to the
                                                                      ! z-axis
    ! ---------------------------------------------------------------------------
    ! LDA+U
    integer, intent (inout) :: kreadldau ! < LDA+U arrays available
    integer, dimension (:), allocatable, intent (inout) :: lopt ! < angular
                                                                ! momentum
                                                                ! QNUM for the
                                                                ! atoms on
                                                                ! which LDA+U
                                                                ! should be
                                                                ! applied (-1
                                                                ! to switch it
                                                                ! OFF)
    real (kind=dp), dimension (:), allocatable, intent (out) :: ueff ! < input
                                                                     ! U
                                                                     ! parameter
                                                                     ! for
                                                                     ! each
                                                                     ! atom
    real (kind=dp), dimension (:), allocatable, intent (out) :: jeff ! < input
                                                                     ! J
                                                                     ! parameter
                                                                     ! for
                                                                     ! each
                                                                     ! atom
    real (kind=dp), dimension (:), allocatable, intent (out) :: erefldau ! <
                                                                         ! the
                                                                         ! energies
                                                                         ! of
                                                                         ! the
                                                                         ! projector's
                                                                         ! wave
                                                                         ! functions
                                                                         ! (REAL)
    ! LDA+U
    ! ---------------------------------------------------------------------------

    ! ----------------------------------------------------------------------------
    ! Local variables
    ! ----------------------------------------------------------------------------
    ! for OPERATOR option
    logical :: lexist, operator_imp
    ! IVSHIFT test option
    logical :: test, opt
    external :: test, opt
    ! ..
    ! .. Local Scalars ..
    real (kind=dp), parameter :: eps=10D-13
    integer :: ndim                ! < Dimension for the Bravais lattice for
                                   ! slab or bulk (2/3)
    integer :: nasoc
    integer :: i, il, j, ier, ier2, i1, ii, ir, idosemicore, i_stat, i_all
    real (kind=dp) :: soscale, ctlscale
    real (kind=dp) :: brymix, strmix, tx, ty, tz
    character (len=43) :: tshape
    character (len=256) :: uio                             ! NCOLIO=256


    logical :: lnew                ! < Logical variable for old/new treatment
                                   ! of left and right host
    logical :: mansoc
    logical :: manctl
    logical :: latominfo           ! < Logical variable for old/new treatment
                                   ! of the ATOMINFO
    ! .. Local CPA variables
    integer :: io, ia, iq, iprint
    real (kind=dp) :: sum1
    character (len=3), dimension (0:1) :: cpaflag

    ! .. Local Arrays ..
    integer, dimension (:), allocatable :: isp
    integer, dimension (:), allocatable :: imansoc
    real (kind=dp), dimension (10) :: dvec
    character (len=4), dimension (3) :: tspin
    character (len=8), dimension (3) :: tkws
    character (len=2), dimension (-2:-1) :: socii
    character (len=43), dimension (0:3) :: tins
    character (len=43), dimension (0:3) :: tkcor
    character (len=43), dimension (0:2) :: tvrel
    ! ..
    ! .. Data statements ..
    data tspin/'non-', '    ', '    '/
    data tshape/' exact cell treatment (shape correction)  '/
    data tvrel/' non relativistic calculation              ', &
      ' s.r.a. calculation                        ', &
      ' fully relativistic calculation            '/
    data tkcor/' frozen core approximation                 ', &
      ' core relaxation s.r.a.                    ', &
      ' core relaxation nonsra                    ', &
      ' core relaxation                           '/
    data tins/' spherical averaged input potential        ', &
      ' non spherical input potential for cluster ', &
      ' non spherical input potential for cluster ', &
      ' non spherical input potential             '/
    data tkws/' full mt', '   ws   ', ' full ws'/

    data cpaflag/' NO', 'YES'/
    data socii/'xy', 'zz'/
    ! ..

    ! ------------ array set up and definition of input parameter -----------

    ! concatenate name & serial number
    txc(1) = ' Morruzi,Janak,Williams  #serial: ' // serialnr
    txc(2) = ' von Barth,Hedin         #serial: ' // serialnr
    txc(3) = ' Vosko,Wilk,Nusair       #serial: ' // serialnr
    txc(4) = ' GGA PW91                #serial: ' // serialnr
    txc(5) = ' GGA PBE                 #serial: ' // serialnr
    txc(6) = ' GGA PBEsol              #serial: ' // serialnr

    iprint = 0

    open (111, file='inputcard_generated.txt') ! Write out found or assumed
                                               ! values
    call version_print_header(111)

    nemb = 0

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read RUNNING options
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ioinput('RUNOPT          ', uio, 1, 7, ier)
    if (ier/=0) then
      write (111, *) 'RUNOPT not found'
    else
      read (unit=uio, fmt=130)(t_params%optc(i), i=1, 8)
      write (111, fmt='(A6)') 'RUNOPT'
      write (111, fmt=130)(t_params%optc(i), i=1, 8)
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read TEST options
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ioinput('TESTOPT         ', uio, 1, 7, ier)
    if (ier/=0) then
      write (111, *) 'TESTOPT not found'
    else
      read (unit=uio, fmt=130)(t_params%testc(i), i=1, 8)
      call ioinput('TESTOPT         ', uio, 2, 7, ier)
      read (unit=uio, fmt=130)(t_params%testc(8+i), i=1, 8)
      write (111, fmt='(A7)') 'TESTOPT'
      write (111, fmt=130)(t_params%testc(i), i=1, 8)
      write (111, fmt=130)(t_params%testc(8+i), i=1, 8)
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin lattice structure definition
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ioinput('ALATBASIS       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) alat
      write (111, *) 'ALATBASIS=', alat
    else
      write (111, *) 'ALATBASIS not found in inputcard'
      write (*, *) 'rinput13: ALATBASIS not found in inputcard'
      stop 'rinput13: ALATBASIS not found in inputcard'
    end if

    ! Set 2-d or 3-d geometry
    linterface = .false.
    call ioinput('INTERFACE       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) linterface
      write (111, *) 'INTERFACE=', linterface
    else
      write (111, *) 'Default INTERFACE= ', linterface
    end if

    ndim = 3
    if (linterface) ndim = 2
    if (.not. linterface .and. .not. opt('SUPRCELL')) then
      write (1337, *) &
        '3D-calculation, adding run-option "full inv" for full inversion.'
      call addopt('full inv')
    end if

    if (opt('WRTGREEN')) then
      write (1337, *) 'WRTGREEN option found'
      write (1337, *) 'adding run-opt "full inv" for full inversion.'
      write (1337, *) 'adding run-opt "fix mesh"'
      call addopt('full inv')
      call addopt('fix mesh')
    end if

    write (111, *) 'Bravais vectors in units of ALAT'
    bravais(1:3, 1:3) = 0d0
    do i = 1, ndim
      call ioinput('BRAVAIS         ', uio, i, 7, ier)
      if (ier/=0) stop 'RINPUT: BRAVAIS NOT FOUND'
      read (unit=uio, fmt=*)(bravais(j,i), j=1, ndim)
    end do
    write (111, fmt='(A7)') 'BRAVAIS'
    do i = 1, ndim
      write (111, *)(bravais(j,i), j=1, ndim)
    end do

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read the number of atoms in the unit cell
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ioinput('NAEZ            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) naez
      write (111, *) 'NAEZ=', naez
    else
      write (111, *) 'NAEZ not found'
      stop 'NAEZ not found in <RINPUT13>'
    end if
    ! if (NAEZ.GT.NAEZD) then
    ! write(6,*) ' set NAEZD to at least ',NAEZ
    ! stop ' in < RINPUT13 > '
    ! end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Read the atom types, if no CPA NATYP=NAEZ
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    natyp = naez
    call ioinput('NATYP           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) natyp
      write (111, *) 'NATYP=', natyp
    else
      write (111, *) 'Default NATYP= ', naez
    end if
    ! if (NATYP.GT.NATYPD) then
    ! write(6,*) 'RINPUT13: NATYP > NATYPD',NATYP,NATYPD
    ! stop ' IN < RINPUT13 > '
    ! end if
    if (natyp<naez) then
      write (6, *) 'RINPUT13: NATYP < NAEZ ', natyp, naez
      stop ' IN < RINPUT13 > '
    end if

    allocate (isp(natyp), stat=i_stat)
    call memocc(i_stat, product(shape(isp))*kind(isp), 'ISP', 'rinput13')
    isp = 0

    lcartesian = .false.
    ier = 0
    call ioinput('CARTESIAN       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) lcartesian
      write (111, *) 'CARTESIAN= ', lcartesian
    else
      write (111, *) 'Default CARTESIAN= ', lcartesian
    end if

    ! Jonathan Chico: This call needs to be done before the rest as one needs
    ! to
    ! find out the value of NEMB to be able to allocate several arrays
    if (linterface) then
      write (1337, 770)

      nright = 10
      call ioinput('NRIGHTHO        ', uio, 1, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) nright
        write (111, *) 'NRIGHTHO=', nright
      else
        write (111, *) 'Default NRIGHTHO=', nright
      end if

      nleft = 10
      call ioinput('NLEFTHOS        ', uio, 1, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) nleft
        write (111, *) 'NLEFTHOS=', nleft
      else
        write (111, *) 'Default NLEFTHOS=', nleft
      end if

      call ioinput('<NLBASIS>       ', uio, 1, 7, ier)
      if (ier/=0) then
        write (1337, *) 'rinput13: <NLBASIS> not found in inputcard'
        ier = 0
        call ioinput('NLBASIS         ', uio, 1, 7, ier)
        if (ier/=0) then
          write (*, *) 'rinput13: NLBASIS also not found in inputcard'
          stop 'rinput13: NLBASIS not found in inputcard'
        end if
      end if
      if (ier==0) then
        read (unit=uio, fmt=*) nlbasis
        write (111, *) '<NLBASIS>=', nlbasis
      end if

      call ioinput('<NRBASIS>       ', uio, 1, 7, ier)
      if (ier/=0) then
        write (1337, *) 'rinput13: <NRBASIS> not found in inputcard'
        ier = 0
        call ioinput('NRBASIS         ', uio, 1, 7, ier)
        if (ier/=0) then
          write (*, *) 'rinput13: NRBASIS also not found in inputcard'
          stop 'rinput13: NRBASIS not found in inputcard'
        end if
      end if
      if (ier==0) then
        read (unit=uio, fmt=*) nrbasis
        write (111, *) '<NRBASIS>=', nrbasis
      end if

      nemb = nlbasis + nrbasis
      write (1337, *) 'Number of embedded atoms NEMB=NLBASIS + NRBASIS=', nemb
      ! if(NEMB.GT.NEMBD) then
      ! write(6,*) 'Please, increase the parameter nembd (',nembd,') in inc.p
      ! to',nemb
      ! stop 'ERROR in NEMBD.'
      ! endif

      ier = 0
      ! Check if the keywords exist for old/new treatment of left and right
      ! host
      call ioinput('LEFTBASIS       ', uio, 1, 7, ier)
      if (ier==0) then
        lnew = .false.
      else
        lnew = .true.
        ier = 0
        call ioinput('<RBLEFT>        ', uio, 1, 7, ier)
      end if
      if (ier/=0) then
        write (*, *) 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
        stop 'rinput13: LEFTBASIS or <RBLEFT> not found in inputcard'
      end if
      ier = 0
      call ioinput('RIGHBASIS       ', uio, 1, 7, ier)
      if (ier==0) then
        lnew = .false.
      else
        lnew = .true.
        ier = 0
        call ioinput('<RBRIGHT>       ', uio, 1, 7, ier)
      end if
      if (ier/=0) then
        write (*, *) 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
        stop 'rinput13: RIGHBASIS or <RBRIGHT> not found in inputcard'
      end if
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate the unit cell arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_cell(1, naez, nemb, natyp, cls, imt, irws, irns, ntcell, &
      refpot, kfg, kaoez, rmt, zat, rws, mtfac, rmtref, rmtrefat, rmtnew, &
      rbasis, lmxc)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation of the unit cell arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate the right and left hosts for slab calculation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_semi_inf_host(1, nemb, tleft, tright)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation of the right and left hosts for slab calculation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Basis atoms
    write (111, fmt='(A16)') '<RBASIS>        '
    do i = 1, naez
      call ioinput('<RBASIS>        ', uio, i, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*)(rbasis(j,i), j=1, 3)
        write (111, fmt='(3E24.12)')(rbasis(j,i), j=1, 3)
      else
        ier = 0
        call ioinput('RBASIS          ', uio, i, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*)(rbasis(j,i), j=1, 3)
          write (111, fmt='(3E24.12)')(rbasis(j,i), j=1, 3)
        else
          write (*, *) &
            'RINPUT13: Keyword <RBASIS> or RBASIS not found. Stopping.'
          stop 'RINPUT13: RBASIS'
        end if
      end if
    end do                         ! I=1,NAEZ
    call idreals(rbasis(1,1), 3*naez, iprint)

    dvec(1:3) = 1.d0
    call ioinput('BASISCALE       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*)(dvec(i), i=1, 3)
      write (111, fmt='(A10,3E12.4)') 'BASISCALE=', dvec(1:3)
    else
      write (111, fmt='(A18,3E12.4)') 'Default BASISCALE=', dvec(1:3)
    end if

    call idreals(dvec(1), 3, iprint)
    abasis = dvec(1)
    bbasis = dvec(2)
    cbasis = dvec(3)

    write (1337, 220) abasis, bbasis, cbasis
    write (1337, 360)
    write (1337, 180) alat

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin read left- and right-host information in 2d-case.
    ! Set up the embeding positions
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (linterface) then
      ! -------------------------------------------------------------------------
      ! > @note In leftbasis and rightbasis, kaoez is used only in decimation
      ! case.
      ! > Then it indicates the correspondence of the atom-coordinate given
      ! > by leftbasis and rightbasis to the left- and right-host t-matrix
      ! read in
      ! > by decimaread. For the slab case, kaoez is not used in the embedded
      ! positions.
      ! -------------------------------------------------------------------------
      if (lnew) then

        write (111, fmt='(A82)') '<RBLEFT>                     &
          &                                 <RMTREFL>   <KAOEZL>'
        do i = 1, nlbasis
          call ioinput('<RBLEFT>        ', uio, i, 7, ier)
          read (unit=uio, fmt=*)(tleft(i1,i), i1=1, 3)
          kaoez(1, naez+i) = i     ! Default
          call ioinput('<KAOEZL>        ', uio, i, 7, ier)
          if (ier==0) read (unit=uio, fmt=*) kaoez(1, naez+i)
          call ioinput('<RMTREFL>       ', uio, i, 7, ier)
          if (ier==0) read (unit=uio, fmt=*) rmtrefat(naez+i)
          write (111, fmt='(3E20.12,3X,F9.6,3X,I5)')(tleft(i1,i), i1=1, 3), &
            rmtrefat(naez+i), kaoez(1, naez+i)
        end do
        write (111, fmt='(A82)') '<RBRIGHT>                    &
          &                                 <RMTREFR>   <KAOEZL>'
        do i = 1, nrbasis
          call ioinput('<RBRIGHT>       ', uio, i, 7, ier)
          read (unit=uio, fmt=*)(tright(i1,i), i1=1, 3)
          kaoez(1, naez+nlbasis+i) = i ! Default
          call ioinput('<KAOEZR>        ', uio, i, 7, ier)
          if (ier==0) read (unit=uio, fmt=*) kaoez(1, naez+nlbasis+i)
          call ioinput('<RMTREFR>       ', uio, i, 7, ier)
          if (ier==0) read (unit=uio, fmt=*) rmtrefat(naez+nlbasis+i)
          write (111, fmt='(3E20.12,3X,F9.6,3X,I5)')(tright(i1,i), i1=1, 3), &
            rmtrefat(naez+nlbasis+i), kaoez(1, naez+nlbasis+i)
        end do

      else                         ! (LNEW) now old-style input

        do i = 1, nlbasis
          call ioinput('LEFTBASIS       ', uio, i, 7, ier)
          read (unit=uio, fmt=*)(tleft(i1,i), i1=1, 3), ii, ir
          kaoez(1, naez+i) = ii    ! changed 1.11.99
          refpot(naez+i) = ir
        end do
        do i = 1, nrbasis
          call ioinput('RIGHBASIS       ', uio, i, 7, ier)
          read (unit=uio, fmt=*)(tright(i1,i), i1=1, 3), ii, ir
          kaoez(1, naez+nlbasis+i) = ii ! changed 1.11.99
          refpot(naez+nlbasis+i) = ir
        end do
      end if

      call idreals(tleft, 3*(nemb+1), iprint)
      call idreals(tright, 3*(nemb+1), iprint)


      ! Put The additional atoms in the "embeding" positions

      do i = 1, nlbasis
        rbasis(1:3, naez+i) = tleft(1:3, i)
      end do
      do i = 1, nrbasis
        rbasis(1:3, naez+nlbasis+i) = tright(1:3, i)
      end do
      ! -------------------------------------------------------------------------
      ! In RBASIS we have first the basis atoms or the interface
      ! atoms then the left host then the right host the host
      ! goes in the NEMB positions

      ! IN CASE OF CPA the host is treated as an effective
      ! CPA medium, that is, there is only one kind of atom
      ! occupying a crystallographic site.

      ! -------------------------------------------------------------------------
      call ioinput('ZPERIODL        ', uio, 1, 7, ier)
      if (ier/=0) then
        write (*, *) 'rimput13: ZPERIODL not found in inputcard'
        stop 'rimput13: ZPERIODL not found in inputcard'
      else
        read (unit=uio, fmt=*)(zperleft(i1), i1=1, 3)
        write (111, fmt='(A9,3E20.12)') 'ZPERIODL=', (zperleft(i1), i1=1, 3)
      end if
      call idreals(zperleft(1), 3, iprint)

      call ioinput('ZPERIODR        ', uio, 1, 7, ier)
      if (ier/=0) then
        write (*, *) 'rimput13: ZPERIODR not found in inputcard'
        stop 'rimput13: ZPERIODR not found in inputcard'
      else
        read (unit=uio, fmt=*)(zperight(i1), i1=1, 3)
        write (111, fmt='(A9,3E20.12)') 'ZPERIODR=', (zperight(i1), i1=1, 3)
      end if
      call idreals(zperight(1), 3, iprint)

      write (1337, 790) nleft, nlbasis
      write (1337, 800) nright, nrbasis
      write (1337, 810)(zperleft(i1), i1=1, 3)
      write (1337, 820)(zperight(i1), i1=1, 3)
      write (1337, 830)
      write (1337, 840)
      do i = nleft, 1, -1
        do i1 = nlbasis, 1, -1
          tx = tleft(1, i1) + (i-1)*zperleft(1)
          ty = tleft(2, i1) + (i-1)*zperleft(2)
          tz = tleft(3, i1) + (i-1)*zperleft(3)
          write (1337, 780)(i-1)*nlbasis + i1, tx, ty, tz, kaoez(1, i1)
        end do
      end do
      write (1337, 850)
      do i = 1, naez
        write (1337, 780) i, (rbasis(i1,i), i1=1, 3)
      end do
      write (1337, 860)
      do i = 1, nright
        do i1 = 1, nrbasis
          tx = tright(1, i1) + (i-1)*zperight(1)
          ty = tright(2, i1) + (i-1)*zperight(2)
          tz = tright(3, i1) + (i-1)*zperight(3)
          write (1337, 780)(i-1)*nrbasis + i1, tx, ty, tz, kaoez(1, i1)
        end do
      end do

    end if                         ! LINTERFACE
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End read left- and right-host information in 2d-case.
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Although NSPIN is fixed to 1 in REL mode,
    ! NSPIN should be used as 1 or 2 at this stage
    ! to indicate a non- or spin-polarised potential
    ! that has to be read in. NSPIN is set to 1 before
    ! being passed to the subsequent programs.
    ! < TESTDIM > has been accordingly modified
    call ioinput('NSPIN           ', uio, 1, 7, ier)
    if (ier/=0) then
      write (111, *) 'NSPIN not found'
      stop 'NSPIN not found'
    else
      read (unit=uio, fmt=*) nspin
      write (111, *) 'NSPIN=', nspin
    end if

    write (1337, 150) nspin
    write (1337, 350)

    ! Atomic number
    call ioinput('<ZATOM>         ', uio, 1, 7, ier)
    if (ier==0) then
      write (111, '(A10)') '<ZATOM>   '
      do i = 1, natyp
        call ioinput('<ZATOM>         ', uio, i, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) zat(i)
          write (111, fmt='(F6.3)') zat(i)
        end if
      end do
    else
      write (111, *) 'zatom will be read in from pot-file'
    end if

    ! Angular momentum cutoff
    call ioinput('LMAX            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) lmax
      write (111, *) 'LMAX=', lmax
    else
      stop 'LMAX not found'
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation of CPA arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_cpa(1, naez, natyp, noq, icpa, iqat, hostimp, conc)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation of CPA arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i = 1, naez
      icpa(i) = 0
      noq(i) = 1
    end do
    ncpa = 0

    do i = 1, naez
      kaoez(1, i) = i              ! default
      iqat(i) = i                  ! Basis-Site of atom I
    end do
    if (natyp==naez) conc(1:natyp) = 1.d0

    ! CPA calculation, read concentrations
    if (natyp>naez) then

      ncpa = 1
      noq(1:naez) = 0              ! re-initialize

      ier = 0
      ier2 = 0
      call ioinput('<SITE>          ', uio, 1, 7, ier)
      call ioinput('<CPA-CONC>      ', uio, 1, 7, ier2)
      if (ier/=0 .or. ier2/=0) then
        write (1337, *) &
          '<SITE> or <CPA-CONC> not found, will search for ATOMINFOC'
      else

        write (111, fmt='(A18)') '<SITE>  <CPA-CONC>'
        do i = 1, natyp
          call ioinput('<SITE>          ', uio, i, 7, ier)
          read (unit=uio, fmt=*) iqat(i)
          call ioinput('<CPA-CONC>      ', uio, i, 7, ier)
          read (unit=uio, fmt=*) conc(i)
          write (111, fmt='(I5,4X,E16.8)') iqat(i), conc(i)
        end do

        do i = 1, natyp
          iq = iqat(i)
          noq(iq) = noq(iq) + 1
          if (noq(iq)>1) icpa(iq) = 1
          kaoez(noq(iq), iq) = i
        end do

        do iq = 1, naez
          sum1 = 0d0
          if (noq(iq)<1) then
            write (6, *) 'RINPUT13: CPA: SITE', iq, 'HAS NO ASSIGNED ATOM'
            stop 'RINPUT13: CPA'
          end if
          do io = 1, noq(iq)
            sum1 = sum1 + conc(kaoez(io,iq))
          end do
          if (abs(sum1-1.d0)>1d-6) then
            write (6, *) ' SITE ', iq, ' CONCENTRATION <> 1.0 !'
            write (6, *) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
            stop ' IN <RINPUT99>'
          end if
        end do

      end if
    end if                         ! (NATYP.GT.NAEZ)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End atom type information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin relativistic treatment information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    kcor = 2
    ! call IoInput('KCOR      ',UIO,1,7,IER)
    ! read (UNIT=UIO,FMT=*) kcor

    kvrel = 1                      ! 0=Schroedinger / 1=SRA / 2=Dirac
    call ioinput('KVREL           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kvrel
      write (111, *) 'KVREL= ', kvrel
    else
      write (111, *) 'Default KVREL= ', kvrel
    end if

    if (opt('NEWSOSOL')) KORBIT = 1
    call ioinput('KORBIT          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) korbit
      write (111, *) 'KORBIT= ', korbit
    else
      write (111, *) 'Default KORBIT= ', korbit
    end if

    ! ----------------------------------------------------------------------------
    ! Start of the reading of variables that used to be in the inc.p
    ! ----------------------------------------------------------------------------
    ! > @note JC: Read the IRM value from the inputcard. This in principle can
    ! be determined from
    ! > the potential file, hence maybe it is best to do it that way instead
    call ioinput('IRMD            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) irmd
      write (111, *) 'IRMD= ', irmd
    else
      write (111, *) 'Default IRMD= ', irmd
    end if

    call ioinput('IRNSD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) irnsd
      write (111, *) 'IRNSD= ', irnsd
    else
      write (111, *) 'Default IRNSD= ', irnsd
    end if

    call ioinput('NSHELD          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) nsheld
      write (111, *) 'NSHELD= ', nsheld
    else
      write (111, *) 'Default NSHELD= ', nsheld
    end if

    call ioinput('KNOSPH          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) knosph
      write (111, *) 'KNOSPH= ', knosph
    else
      write (111, *) 'Default KNOSPH= ', knosph
    end if

    call ioinput('IEMXD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) iemxd
      write (111, *) 'IEMXD= ', iemxd
    else
      write (111, *) 'Default IEMXD= ', iemxd
    end if

    call ioinput('NRMESH          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) nrd
      write (111, *) 'NRMESH= ', nrd
    else
      write (111, *) 'Default NRMESH= ', nrd
    end if

    kvrel = 1                      ! 0=Schroedinger / 1=SRA / 2=Dirac
    call ioinput('KVREL           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kvrel
      write (111, *) 'KVREL= ', kvrel
    else
      write (111, *) 'Default KVREL= ', kvrel
    end if
    ! store KVREL to be used later on
    t_inc%kvrel = kvrel

    call ioinput('KPOIBZ          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kpoibz
      write (111, *) 'KPOIBZ= ', kpoibz
    else
      write (111, *) 'Default KPOIBZ= ', kpoibz
    end if

    call ioinput('NMAXD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) nmaxd
      write (111, *) 'NMAXD= ', nmaxd
    else
      write (111, *) 'Default NMAXD= ', nmaxd
    end if

    call ioinput('ISHLD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ishld
      write (111, *) 'ISHLD= ', ishld
    else
      write (111, *) 'Default ISHLD= ', ishld
    end if

    call ioinput('KNOCO           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) knoco
      write (111, *) 'KNOCO= ', knoco
    else
      write (111, *) 'Default KNOCO= ', knoco
    end if

    call ioinput('NTREFD          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ntrefd
      write (111, *) 'NTREFD= ', ntrefd
    else
      write (111, *) 'Default NTREFD= ', ntrefd
    end if

    call ioinput('NATOMIMPD       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) natomimpd
      write (111, *) 'NATOMIMPD= ', natomimpd
    else
      write (111, *) 'Default NATOMIMPD= ', natomimpd
    end if

    call ioinput('NPRINCD         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) nprincd
      write (111, *) 'NPRINCD= ', nprincd
    else
      write (111, *) 'Default NPRINCD= ', nprincd
    end if

    call ioinput('IPAND           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ipand
      write (111, *) 'IPAND= ', ipand
    else
      write (111, *) 'Default IPAND= ', ipand
    end if

    call ioinput('NFUND           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) nfund
      write (111, *) 'NFUND= ', nfund
    else
      write (111, *) 'Default NFUND= ', nfund
    end if

    call ioinput('IRID            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) irid
      write (111, *) 'IRID= ', irid
    else
      write (111, *) 'Default IRID= ', irid
    end if

    call ioinput('NGSHD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ngshd
      write (111, *) 'NGHSD= ', ngshd
    else
      write (111, *) 'Default NGSHD= ', ngshd
    end if

    call ioinput('WLENGTH         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) wlength
      write (111, *) 'WLENGTH= ', wlength
    else
      write (111, *) 'Default WLENGTH= ', wlength
    end if

    call ioinput('NACLSD          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) naclsd
      write (111, *) 'NACLSD= ', naclsd
    else
      write (111, *) 'Default NACLSD= ', naclsd
    end if

    call ioinput('NTOTD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ntotd
      write (111, *) 'NTOTD= ', ntotd
    else
      write (111, *) 'Default NTOTD= ', ntotd
    end if
    ! ----------------------------------------------------------------------------
    ! End of variables that used to be in the inc.
    ! ----------------------------------------------------------------------------
    ! ----------------------------------------------------------------------------
    ! Calculate derived parameters
    ! ----------------------------------------------------------------------------
    lm2d = (2*lmax+1)**2
    nclsd = naez + nemb
    mmaxd = 2*lmax + 1
    ncleb = (lmax*2+1)**2*(lmax+1)**2
    nspind = krel + (1-krel)*2     ! (KSP+1) where KSP is always 1
    npotd = (2*(krel+korbit)+(1-(krel+korbit))*nspind)*natyp
    lmmaxd = (krel+korbit+1)*(lmax+1)**2
    lmgf0d = (lmax+1)**2
    lassld = 4*lmax
    nembd1 = nemb + 1
    irmind = irmd - irnsd
    nofgij = natomimpd**2 + 1
    ntperd = natyp - ntrefd
    nspindd = nspind - korbit
    nsatypd = (natyp-1)*knosph + 1
    nspotd = (2*krel+(1-krel)*nspind)*nsatypd
    if (krel/=0 .or. korbit/=0 .or. knoco/=0) then
      lnc = .true.
    else
      lnc = .false.
    end if
    ! ----------------------------------------------------------------------------
    ! End of calculation of the derived parameters
    ! ----------------------------------------------------------------------------

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocation of SOC arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_soc(1, krel, natyp, lmax, socscale, cscl, socscl)
    allocate (imansoc(natyp), stat=i_stat)
    call memocc(i_stat, product(shape(imansoc))*kind(imansoc), 'IMANSOC', &
      'rinput13')
    imansoc = 0
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation of SOC arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (opt('NEWSOSOL')) then      ! Spin-orbit
      if (opt('NEWSOSOL') .and. (nspin/=2)) stop &
        ' set NSPIN = 2 for SOC solver in inputcard'
      npan_log = 30
      npan_eq = 30
      ncheb = 10
      r_log = 0.1d0
      call ioinput('NPAN_LOG        ', uio, 1, 7, ier)
      if (ier==0) read (unit=uio, fmt=*) npan_log
      call ioinput('NPAN_EQ         ', uio, 1, 7, ier)
      if (ier==0) read (unit=uio, fmt=*) npan_eq
      call ioinput('NCHEB           ', uio, 1, 7, ier)
      if (ier==0) read (unit=uio, fmt=*) ncheb
      call ioinput('R_LOG           ', uio, 1, 7, ier)
      if (ier==0) read (unit=uio, fmt=*) r_log
      write (111, *) 'NPAN_LOG= ', npan_log
      write (111, *) 'NPAN_EQ= ', npan_eq
      write (111, *) 'NCHEB= ', ncheb
      write (111, *) 'R_LOG= ', r_log
    end if

    call ioinput('<SOCSCL>        ', uio, 1, 7, ier)
    if (ier==0) then
      write (111, '(A10)') '<SOCSCL>  '
      do i = 1, natyp
        call ioinput('<SOCSCL>        ', uio, i, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) socscale(i)
          write (111, fmt='(F6.3)') socscale(i)
        end if
      end do
      ! read (UNIT=UIO,FMT=*) (SOCSCALE(I1),I1=1,NATYP)
      ! !Bernd - old way
      ! write(111,FMT='(A10,50E10.2)') '<SOCSCL>= ',(SOCSCALE(I1),I1=1,NATYP)
      ! !Bernd - old way
    else
      write (111, fmt='(A18,50E10.2)') 'Default <SOCSCL>= ', &
        (socscale(i1), i1=1, natyp)
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End relativistic treatment information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin cell control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call ioinput('<FPRADIUS>      ', uio, 1, 7, ier)
    if (ier==0) then
      write (111, '(A10)') '<FPRADIUS>'
      do i = 1, natyp
        call ioinput('<FPRADIUS>      ', uio, i, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) fpradius(i)
        end if
        write (111, fmt='(F6.3)') fpradius(i)
      end do
    else
      write (111, *) 'fpradius will be read in from pot-file'
    end if


    ins = 1
    call ioinput('INS             ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ins
      write (111, *) 'INS= ', ins
    else
      write (111, *) 'Default INS= ', ins
    end if

    kshape = 2
    if (ins==0) kshape = 0
    call ioinput('KSHAPE          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kshape
      write (111, *) 'KSHAPE= ', kshape
    else
      write (111, *) 'Default KSHAPE= ', kshape
    end if

    if ((krel==1) .and. (kshape/=0)) then
      write (1337, *) ' WARNING : KSHAPE set to ZERO for REL case'
      write (111, *) ' WARNING : kshape set to ZERO for REL case'
      kshape = 0
    end if

    ! Read cell information
    write (1337, *) 'Cell information <SHAPE>:'
    write (111, fmt='(A16)') '<SHAPE>         '
    do i = 1, natyp
      ntcell(i) = iqat(i)          ! Default: Different shape function per
                                   ! atom
      call ioinput('<SHAPE>         ', uio, i, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) ntcell(i)
        write (111, fmt='(I6)') ntcell(i)
      end if
    end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End cell control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin exchange correlation treatment information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    kxc = 2                        ! 0=vBH 1=MJW 2=VWN 3=PW91
    call ioinput('KEXCOR          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kxc
      write (111, *) 'KEXCOR= ', kxc
    else
      write (111, *) 'Default KEXCOR= ', kxc
    end if

    ! Scale magnetic moment (0 < Lambda_XC < 1,  0=zero moment, 1= full
    ! moment)
    lambda_xc = 1.d0
    call ioinput('LAMBDA_XC       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) lambda_xc
      write (111, *) 'LAMBDA_XC= ', lambda_xc
    else
      write (111, *) 'Default LAMBDA_XC= ', lambda_xc
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDA+U treatment
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate the LDA+U arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_ldau(1, natyp, lopt, ueff, jeff, erefldau)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of LDA+U array allocation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (opt('LDA+U   ')) then

      ! Check for LDA+U consistency -- if INS=0 suppress it
      if ((ins==0)) then
        write (1337, *)
        write (1337, *) ' WARNING: LDA+U should be used only in NON-SPHERICAL' &
          , ' case (INS=1) '
        write (1337, *) ' Running option LDA+U will be ignored'
        write (1337, *)
        do i = 1, 32
          if (t_params%optc(i)(1:8)=='LDA+U   ') t_params%optc(i) = '        '
        end do
      end if

      ! -> get number of atoms for lda+u:

      ier = 0
      call ioinput('NAT_LDAU        ', uio, 1, 7, ier)
      if (ier/=0) then
        nasoc = natyp
      else
        read (unit=uio, fmt=*) nasoc
        if (nasoc>natyp) stop ' main0: NAT_LDAU > NATYP'
      end if

      ! -> read in UEFF,JEFF,LOPT,EREFLDAU for the desired atoms

      il = 0
      do i = 1, nasoc
        ier = 0
        call ioinput('LDAU_PARA       ', uio, i, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) i1, lopt(i1), ueff(i1), jeff(i1), &
            erefldau(i1)
          il = il + 1
        end if
      end do
      if (il/=nasoc) then
        write (6, *) ' ERROR: LDA+U invoked for ', nasoc, ' atoms'
        write (6, *) &
          '        Some (all) parameters are missing in the input-file'
        stop
      end if
      kreadldau = 0
      ier = 0
      call ioinput('KREADLDAU       ', uio, 1, 7, ier)
      if (ier==0) read (unit=uio, fmt=*) kreadldau

    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End exchange correlation treatment information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin external field control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    khfield = 0
    hfield = 0.d0
    call ioinput('HFIELD          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) hfield
      if (abs(hfield)>eps) then
        khfield = 1
        write (*, *) 'WARNING: HFIELD>0.0 found, set KHFIELD to 1'
        write (1337, *) 'WARNING: HFIELD>0.0 found, set KHFIELD to 1'
      end if
      write (111, *) 'HFIELD= ', hfield
    else
      write (111, *) 'Default HFIELD= ', hfield
    end if

    vconst = 0.d0
    call ioinput('VCONST          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) vconst
      write (111, *) 'VCONST= ', vconst
    else
      write (111, *) 'Default VCONST= ', vconst
    end if


    if (test('atptshft')) then
      write (1337, *) 'read IN IVSHIFT'
      call ioinput('IVSHIFT         ', uio, 1, 7, ier)
      read (unit=uio, fmt=*) ivshift
    end if

    ! Initial polarization
    linipol = .false.
    call ioinput('LINIPOL         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) linipol
      write (111, *) 'LINIPOL= ', linipol
    else
      write (111, *) 'Default: LINIPOL= ', linipol
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Allocate magnetization arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call allocate_magnetization(1, naez, natyp, lmmaxd, inipol, ixipol, qmtet, &
      qmphi, drotq)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End of allocation of magnetization arrays
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (linipol) then
      inipol(1:natyp) = 1
      call ioinput('XINIPOL         ', uio, 1, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*)(inipol(i), i=1, natyp)
        write (111, fmt='(A10,80I2)') 'XINIPOL=  ', (inipol(i), i=1, natyp)
      else
        write (111, fmt='(A18,80I2)') 'Default XINIPOL=  ', &
          (inipol(i), i=1, natyp)
      end if
    end if


    write (1337, 230)(inipol(i), i=1, natyp)
    write (1337, 340)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End external field control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin Green function calculation control (diag./non-diag)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    igf = 0
    call ioinput('IGREENFUN       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) igf
      write (111, *) 'IGREENFUN= ', igf
    else
      write (111, *) 'Default IGREENFUN= ', igf
    end if

    if (opt('OPERATOR')) then
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
    if (opt('KKRFLEX ') .or. opt('WRTGREEN') .or. opt('GREENIMP') .or. &
      operator_imp) then
      write (1337, *) &
        'Setting IGREENFUN=1 for KKRFLEX/WRTGREEN/GREENIMP/OPERATOR options'
      igf = 1
    end if

    icc = 0
    call ioinput('ICC             ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) icc
      write (111, *) 'ICC= ', icc
    else
      write (111, *) 'Default ICC= ', icc
    end if
    if (opt('KKRFLEX ') .or. opt('WRTGREEN') .or. opt('GREENIMP') .or. &
      operator_imp) then
      write (1337, *) &
        'Setting ICC=1 for KKRFLEX/WRTGREEN/GREENIMP/OPERATOR  options'
      icc = 1
    end if
    if ((opt('XCPL    ')) .or. (opt('CONDUCT '))) icc = -1

    if (icc/=0 .and. igf==0) igf = 1
    if (icc==0 .and. igf/=0) icc = -1
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Green function calculation control (diag./non-diag)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin accuracy parameters
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Brilloun zone mesh
    intervx = 10
    intervy = 10
    intervz = 10
    call ioinput('BZDIVIDE        ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) intervx, intervy, intervz
      write (111, fmt='(A9,3I5)') 'BZDIVIDE=', intervx, intervy, intervz
    else
      write (111, fmt='(A17,3I5)') 'Default BZDIVIDE=', intervx, intervy, &
        intervz
    end if
    write (1337, 350)
    write (1337, 190) intervx, intervy, intervz
    write (1337, 330)

    if (opt('GREENIMP')) then
      write (*, *) &
        'WARNING! Found option GREENIMP: resetting BZDIVIDE to 1,1,1'
      write (1337, *) &
        'WARNING! Found option GREENIMP: resetting BZDIVIDE to 1,1,1'
      intervx = 1
      intervy = 1
      intervz = 1
    end if

    ! Energy contour
    npol = 7
    ! if (OPT('dos     ').OR.OPT('DOS     ')) NPOL = 0
    call ioinput('NPOL            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) npol
      write (111, *) 'NPOL=', npol
    else
      write (111, *) 'Default NPOL=', npol
    end if

    call ioinput('EMIN            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) emin
      write (111, *) 'EMIN= ', emin
    else if (npol==0) then
      emin = -1.d0
      write (111, *) 'Default for DOS: EMIN= ', emin
    else
      write (1337, *) 'Error in rinput13: EMIN not found'
      write (111, *) 'Error in rinput13: EMIN not found'
      stop 'Error in rinput13: EMIN not found'
    end if

    emax = 1.d0
    call ioinput('EMAX            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) emax
      write (111, *) ' EMAX=', emax
    else
      write (111, *) 'Default  EMAX=', emax
    end if

    tk = 800.d0
    call ioinput('TEMPR           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) tk
      write (111, *) 'TEMPR=', tk
    else
      write (111, *) 'Default TEMPR=', tk
    end if

    npnt1 = 3
    if (npol==0) npnt1 = 0
    call ioinput('NPT1            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) npnt1
      write (111, *) ' NPT1=', npnt1
    else
      write (111, *) 'Default  NPT1=', npnt1
    end if

    npnt2 = nint((emax-emin)*20.d0) ! 20 pts/Ryd
    if (npol==0) npnt2 = nint((emax-emin)*100.d0) ! For dos, 100 pts/Ryd
    call ioinput('NPT2            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) npnt2
      write (111, *) ' NPT2=', npnt2
    else
      write (111, *) 'Default  NPT2=', npnt2
    end if

    npnt3 = 3
    if (npol==0) npnt3 = 0
    call ioinput('NPT3            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) npnt3
      write (111, *) ' NPT3=', npnt3
    else
      write (111, *) 'Default  NPT3=', npnt3
    end if

    ! -> semicore
    ! initialise variables
    idosemicore = 0
    ebotsemi = emin
    emusemi = ebotsemi
    npolsemi = 0
    n1semi = 0
    n2semi = 0
    n3semi = 0
    fsemicore = 1.d0

    ier = 0
    if (opt('SEMICORE')) then
      call ioinput('EBOTSEMI        ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) ebotsemi
      call ioinput('EMUSEMI         ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) emusemi

      ! -> EMUSEMI < EBOT
      if (emusemi>=emin) go to 100
      call ioinput('TKSEMI          ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) tksemi

      call ioinput('NPOLSEMI        ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) npolsemi
      call ioinput('N1SEMI          ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) n1semi
      call ioinput('N2SEMI          ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) n2semi
      call ioinput('N3SEMI          ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) n3semi
      call ioinput('FSEMICORE       ', uio, 1, 7, ier)
      if (ier/=0) go to 100
      read (unit=uio, fmt=*) fsemicore
      idosemicore = 1
100   continue
      if (idosemicore==0) then
        write (1337, *)
        write (1337, *) ' WARNING: SEMICORE used', &
          ' with incomplete/incorrect contour description'
        write (1337, *) ' Running option SEMICORE will be ignored'
        write (111, *)
        write (111, *) ' WARNING: SEMICORE used', &
          ' with incomplete/incorrect contour description'
        write (111, *) ' Running option SEMICORE will be ignored'
        do i = 1, 32
          if (t_params%optc(i)(1:8)=='SEMICORE') t_params%optc(i) = '        '
        end do
      end if
    end if

    ! CPA convergence parameters
    cpatol = 1d-4
    itcpamax = 20
    call ioinput('CPAINFO         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) cpatol, itcpamax
    else
      write (111, *) 'Default cpainfo:'
    end if
    write (111, fmt='(A7)') 'CPAINFO'
    write (111, fmt='(E12.4,I5)') cpatol, itcpamax

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin screening cluster information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rcutz = 11.d0/alat             ! Default 11 Bohr radii
    call ioinput('RCLUSTZ         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) rcutz
      write (111, *) 'RCLUSTZ=', rcutz
    else
      write (111, *) 'Default RCLUSTZ=', rcutz
    end if

    rcutxy = rcutz
    call ioinput('RCLUSTXY        ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) rcutxy
      write (111, *) 'RCLUSTXY=', rcutxy
    else
      write (111, *) 'Default RCLUSTXY=', rcutxy
    end if

    write (1337, *) 'Parameters used for the cluster calculation'
    if (abs(rcutz-rcutxy)<1.d-4) then
      write (1337, *) 'Clusters inside spheres with radius R = ', rcutz
    else
      write (1337, *) 'Clusters inside cylinders with '
      write (1337, *) 'Rz = ', rcutz, ' Rxy = ', rcutxy
    end if
    write (1337, 350)
    write (1337, 210)              ! rbasis
    write (1337, 320)
    do i = 1, naez
      write (1337, 260) i, (rbasis(j,i), j=1, 3), qmtet(i), qmphi(i), icpa(i), &
        noq(i), (kaoez(j,i), j=1, noq(i))
    end do

    do i = 1, naez
      call ioinput('<RMTREF>        ', uio, i, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) rmtrefat(i)
      end if
    end do
    if (ier==0) then
      write (111, fmt='(A18)') '        <RMTREF>  '
    else
      write (111, fmt='(A18)') 'Default <RMTREF>  '
    end if
    do i = 1, naez
      write (111, fmt='(9X,F9.6)') rmtrefat(i)
    end do

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End screening cluster information
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Number of Born iterations
    icst = 2
    call ioinput('ICST            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) icst
      write (111, *) 'ICST=', icst
    else
      write (111, *) 'Default ICST=', icst
    end if

    ! Usage of Lloyd's formula
    lly = 0                        ! LLY Default=0 : do not apply Lloyds
                                   ! formula
    if (opt('LLOYD   ') .or. opt('Lloyd   ') .or. opt('lloyd   ')) lly = 1
    call ioinput('<LLOYD>         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) lly
      write (111, *) '<LLOYD>=', lly
    else
      write (111, *) 'Default <LLOYD>=', lly
    end if
    if (lly/=0) write (1337, *) 'Applying Lloyds formula, LLY=', lly

    deltae = (1.d-5, 0.d0)         ! Difference for numer. derivative in
                                   ! Lloyds formula
    call ioinput('<DELTAE>        ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) deltae
      write (111, *) '<DELTAE>=', deltae
    else
      write (111, *) 'Default <DELTAE>=', deltae
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End accuracy parameters
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin old-type of ATOMINFO
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    latominfo = .false.
    ! Initialize all clusters to 1
    cls(1:naez+nemb) = 1
    write (1337, *) 'ATOMINFOC or ATOMINFO:'
    do i = 1, natyp
      call ioinput('ATOMINFOC       ', uio, i+1, 7, ier)
      ia = 1
      if (ier==0) then
        latominfo = .true.
        read (unit=uio, fmt=*) zat(i), lmxc(i), (kfg(j,i), j=1, 4), j, ier, &
          ntcell(i), mtfac(i), irns(i), rmtref(ier), iqat(i), conc(i)
        iq = iqat(i)
        refpot(iq) = ier
        rmtrefat(i) = rmtref(ier)
        cls(iq) = j
        noq(iq) = noq(iq) + 1
        if (noq(iq)>1) then
          icpa(iq) = 1
          ncpa = 1
        end if
        kaoez(noq(iq), iq) = i
      else
        ier = 0
        call ioinput('ATOMINFO        ', uio, i+1, 7, ier)
        if (ier==0) then
          latominfo = .true.
          read (unit=uio, fmt=*) zat(i), lmxc(i), (kfg(j,i), j=1, 4), j, &
            refpot(i), ntcell(i), mtfac(i), irns(i), rmtref(refpot(i))
          iqat(i) = i
          rmtrefat(i) = rmtref(refpot(i))
          cls(i) = j
          conc(i) = 1.d0
          noq(i) = 1
          kaoez(1, i) = i
        end if
      end if
    end do

    ! If old-style ATOMINFO is present, and if a 2-dim calculation is
    ! performed,
    ! and also if the RMTREF of the "outside region" is not read in explicitly
    ! (LNEW is false) then assign the RMTREF of the outside region according
    ! to
    ! the already-read-in REFPOT under LEFTBASIS  and RIGHBASIS.
    if (latominfo .and. linterface .and. .not. lnew) then
      do i = naez + 1, naez + nemb
        rmtrefat(i) = rmtref(refpot(i))
      end do
    end if


    ! Determine total number of clusters
    ncls = 0
    do i = 1, natyp
      ncls = max(ncls, cls(iqat(i)))
    end do

    ! Determine total number of different reference potentials
    nref = 0
    do i = 1, naez + nemb
      nref = max(nref, refpot(i))
    end do

    ! in line 1792  this is done: NINEQ = NAEZ, so here NINEQ is still
    ! undefinded
    ! so we move this writeout back

    ! write(6,2016) NCLS,NREF,NINEQ
    ! write(6,2110)
    ! write(6,2103)

    do iq = 1, naez
      sum1 = 0d0
      if (noq(iq)<1) then
        write (6, *) 'RINPUT13: CPA: SITE', iq, 'HAS NO ASSIGNED ATOM'
        stop 'RINPUT13: CPA'
      end if
      do io = 1, noq(iq)
        sum1 = sum1 + conc(kaoez(io,iq))
      end do
      if (abs(sum1-1.d0)>1d-6) then
        write (6, *) ' SITE ', iq, ' CONCENTRATION <> 1.0 !'
        write (6, *) ' CHECK YOUR <ATOMINFO-CPA> INPUT '
        stop ' IN <RINPUT99>'
      end if
    end do
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End old-type of ATOMINFO
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Write out atominfo
    write (1337, 270) natyp
    write (1337, 350)
    write (1337, 140)(zat(i), lmxc(i), (kfg(j,i),j=1,4), cls(iqat(i)), refpot( &
      iqat(i)), ntcell(i), mtfac(i), irns(i), iqat(i), conc(i), i=1, natyp)
    write (1337, 370)
    write (1337, 350)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin SCF convergence control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    nsteps = 1
    call ioinput('NSTEPS          ', uio, 1, 7, ier)
    if (ier/=0) then
      write (111, *) 'Default NSTEPS=', nsteps
    else
      read (unit=uio, fmt=*) nsteps
    end if
    if (npol==0) then
      nsteps = 1
      write (1337, *) 'NPOL=0, setting NSTEPS to 1'
    end if
    if (igf/=0) then
      nsteps = 1
      write (1337, *) 'IGF.NE.0, setting NSTEPS to 1'
    end if
    if (icc/=0) then
      nsteps = 1
      write (1337, *) 'ICC.NE.0, setting NSTEPS to 1'
    end if
    if (opt('XCPL    ')) then
      nsteps = 1
      write (1337, *) 'RUNOPT XCPL used, setting NSTEPS to 1'
    end if
    if (opt('KKRFLEX ')) then
      nsteps = 1
      write (1337, *) 'RUNOPT KKRFLEX used, setting NSTEPS to 1'
    end if

    write (1337, 160) nsteps
    write (1337, 350)

    imix = 0
    call ioinput('IMIX            ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) imix
      write (111, *) 'IMIX= ', imix
    else
      write (111, *) 'Default IMIX= ', imix
    end if
    if (npol==0) then
      write (1337, *) 'NPOL=0, setting IMIX= 0'
      imix = 0
    end if

    strmix = 0.01d0
    call ioinput('STRMIX          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) strmix
      write (111, *) 'STRMIX= ', strmix
    else
      write (111, *) 'Default STRMIX= ', strmix
    end if
    if (npol==0) then
      write (1337, *) 'NPOL=0, setting STRMIX= 0.'
      strmix = 0
    end if

    itdbry = 40
    call ioinput('ITDBRY          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) itdbry
      write (111, *) 'ITDBRY= ', itdbry
    else
      write (111, *) 'Default ITDBRY= ', itdbry
    end if

    fcm = 20.d0
    call ioinput('FCM             ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) fcm
      write (111, *) 'FCM= ', fcm
    else
      write (111, *) 'Default FCM= ', fcm
    end if

    qbound = 1.d-7
    call ioinput('QBOUND          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) qbound
      write (111, *) 'QBOUND= ', qbound
    else
      write (111, *) 'Default QBOUND= ', qbound
    end if

    brymix = 0.01d0
    call ioinput('BRYMIX          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) brymix
      write (111, *) 'BRYMIX= ', brymix
    else
      write (111, *) 'Default BRYMIX= ', brymix
    end if

    call ioinput('RMAX            ', uio, 1, 7, ier)
    if (ier/=0) stop 'rinput13: RMAX not in the inputcard'
    read (unit=uio, fmt=*) rmax
    write (111, *) 'RMAX= ', rmax

    call ioinput('GMAX            ', uio, 1, 7, ier)
    if (ier/=0) stop 'rinput13: GMAX not in the inputcard'
    read (unit=uio, fmt=*) gmax
    write (111, *) 'GMAX= ', gmax
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End SCF convergence control
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Begin file name definitions
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    il = 1
    call ioinput('FILES           ', uio, il, 7, ier)
    if (ier==0) then
      call ioinput('FILES           ', uio, il, 7, ier)
      read (unit=uio, fmt='(A40)') i12
      call ioinput('FILES           ', uio, il+1, 7, ier)
      read (unit=uio, fmt='(A40)') i13
      call ioinput('FILES           ', uio, il+2, 7, ier)
      read (unit=uio, fmt='(A40)') i40
      call ioinput('FILES           ', uio, il+3, 7, ier)
      read (unit=uio, fmt='(A40)') i19
      call ioinput('FILES           ', uio, il+4, 7, ier)
      read (unit=uio, fmt='(A40)') i25
    else
      i13 = 'potential                               ' ! 40 chars
      i19 = 'shapefun                                ' ! 40 chars
      i25 = 'scoef                                   ' ! 40 chars
      i12 = '                                        ' ! 40 chars (not used)
      i40 = '                                        ' ! 40 chars (not used)
    end if

    write (1337, *) 'I12="', i12, '"'
    write (1337, *) 'I13="', i13, '"'
    write (1337, *) 'I40="', i40, '"'
    write (1337, *) 'I19="', i19, '"'
    write (1337, *) 'I25="', i25, '"'

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End file name definitions
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ifile = 13
    call ioinput('<IFILE>         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ifile
      write (111, *) '<IFILE>= ', ifile
    else
      write (111, *) 'Default <IFILE>= ', ifile
    end if

    ipe = 1                        ! Used to print out in calrmt
    ishift = 0
    call ioinput('ISHIFT          ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) ishift
      write (111, *) 'ISHIFT= ', ishift
    else
      write (111, *) 'Default ISHIFT= ', ishift
    end if
    if (opt('rigid-ef') .or. opt('DECIMATE')) then
      ishift = 2
      write (1337, *) ' Rigid Fermi Energy, ISHIFT is set to ', ishift
      write (111, *) ' Rigid Fermi Energy, ishift is set to ', ishift
    end if
    if (test('no-neutr') .or. opt('no-neutr')) then
      ishift = 1
      write (1337, *) 'No charge neutrality required, ISHIFT is set to', &
        ishift
      write (111, *) 'No charge neutrality required, ISHIFT is set to', ishift
    end if

    eshift = 0.d0
    insref = 0
    kws = 2
    khyp = 0

    tolrdif = 0.5d0                ! Set free GF to zero for r<tolrdif
                                   ! (a.u.)(vir. atoms)
    call ioinput('<TOLRDIF>       ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) tolrdif
      write (111, *) '<TOLRDIF>=', tolrdif
    else
      write (111, *) 'Default <TOLRDIF>=', tolrdif
    end if

    ! -------------------------------------------------
    kte = 1
    ! call IoInput('KTE       ',UIO,1,7,IER)
    ! read (UNIT=UIO,FMT=*) kte

    kpre = 1
    ! call IoInput('KPRE      ',UIO,1,7,IER)
    ! read (UNIT=UIO,FMT=*) kpre

    kefg = 0
    ! call IoInput('KEFG      ',UIO,1,7,IER)
    ! read (UNIT=UIO,FMT=*) kefg

    kvmad = 0
    call ioinput('KVMAD           ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) kvmad
      write (111, *) 'KVMAD= ', kvmad
    else
      write (111, *) 'Default KVMAD= ', kvmad
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Determination of properties at Fermi level
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (opt('GF-EF   ')) then
      igf = 1
      if (npol>0) npol = 0
      if (npol<0) then
        npnt1 = 0
        npnt3 = 0
      end if
      npnt2 = 1
    end if

    if (opt('DOS-EF  ')) then
      npol = 0
      npnt2 = 1
    end if
    ! ----------------------------------------------------------------------
    ! ---------------------------------------------------------------------

    kforce = 0
    if (ins>0) then
      call ioinput('KFORCE          ', uio, 1, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) kforce
        write (111, *) 'KFORCE= ', kforce
      else
        write (111, *) 'Default KFORCE= ', kforce
      end if
    end if

    kfrozn = kcor
    if (kcor==0) kcor = 2

    ! ------------------------------------------------------------------------
    write (1337, 560) lmax
    write (1337, 680)
    write (1337, 570) emin, emax, tk
    write (1337, 690)
    write (1337, 580) npol, npnt1, npnt2, npnt3
    write (1337, 710)
    write (1337, 700)
    write (1337, 590) ifile, ipe, ishift, eshift
    write (1337, 720)
    write (1337, 600) kshape, irmd, ins, icst, insref
    write (1337, 760)
    write (1337, 610) kcor, kvrel, kws, khyp, khfield, kxc
    write (1337, 730)
    write (1337, 670) kte, kpre, kefg, kvmad
    write (1337, 760)
    write (1337, 630) imix, igf, icc
    write (1337, 710)
    write (1337, 640) itdbry
    write (1337, 740)
    write (1337, 650) strmix, fcm, qbound
    write (1337, 690)
    write (1337, 660) brymix
    write (1337, 750)
    write (1337, 620) hfield, vconst
    ! ------------------------------------------------------------------------

    ! ------------------------------------------------------------------------
    ipf = 1337
    ipfe = ipf + 3

    if (opt('SEARCHEF')) then
      imix = 0
      mixing = 0.0d0
      strmix = mixing
      itdbry = 1
      qbound = 1.0d-10
      write (1337, '(1X,A)') 'Option SEARCHEF used overriding INPUT for'
      write (1337, '(1X,A)') 'IMIX,MIX,QBOUND,ITDBRY: 0, 0.0, 1E-10, 1'
      write (1337, *)
    end if

    if (imix>2) then
      fcm = 1.0d0
      mixing = brymix
    else
      mixing = strmix
    end if

    if (imix>=6) write (1337, fmt=490)(imix-5), itdbry - 1

    write (1337, fmt=450) mixing, qbound
    ! --------------------------------------------------------
    write (1337, fmt=460) cpaflag(ncpa)
    if (ncpa/=0) write (1337, 470) itcpamax, cpatol
    ! --------------------------------------------------------

    lmmax = (lmax+1)**2
    lpot = 2*lmax
    lmpot = (lpot+1)*(lpot+1)
    lmxspd = (2*lpot+1)**2

    write (1337, fmt=400) lmax, lmax, natyp, natyp, irmd, irmd, nspin, nspind

    if (ins>0) then
      write (1337, fmt=510)
      write (1337, fmt=520)
      do i = 1, natyp
        write (1337, fmt=530) i, irns(i), irnsd
        if (irns(i)>irnsd) call rcstop('19      ')
      end do
    end if

    write (1337, fmt=510)

    if (khfield==1) write (1337, fmt=410) hfield
    if (kvrel<=1) then
      write (1337, fmt=420) tspin(nspin)
    else
      write (1337, fmt=420) tspin(nspin+1)
    end if
    write (1337, fmt=550) tvrel(kvrel)
    write (1337, fmt=550) tkcor(kfrozn)
    if (kshape==0) then
      write (1337, fmt=430) tkws(kws+1)
    else
      write (1337, fmt=550) tshape
    end if

    write (1337, fmt=480) txc(kxc+1)
    if (ins>0) write (1337, fmt=540) tins(ins), icst
    write (1337, fmt=440)

    vbc(1) = vconst
    vbc(2) = vbc(1)

    lrhosym = .false.
    call ioinput('LRHOSYM         ', uio, 1, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) lrhosym
      write (111, *) 'LRHOSYM= ', lrhosym
    else
      write (111, *) 'Default LRHOSYM= ', lrhosym
    end if

    if ((ncpa/=0) .and. lrhosym) then
      write (1337, *) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
      write (1337, *) '        YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
      write (111, *) ' WARNING : CHARGE SYMMETRISATION NOT ALLOWED FOR CPA '
      write (111, *) '    YOUR SETTING IN INPUT FILE IS OVERRIDDEN'
      lrhosym = .false.
    end if

    if (lrhosym) then
      call ioinput('IXIPOL          ', uio, 1, 7, ier)
      read (unit=uio, fmt=*)(ixipol(i), i=1, natyp)
      write (1337, 240)(ixipol(i), i=1, natyp)
      write (1337, 340)
      do i = 1, natyp
        if (ixipol(i)/=0 .and. abs(ixipol(abs(ixipol(i))))/=i) then
          write (6, *) 'Error in IXIPOL at atom ', i, '.'
          stop 'IXIPOL'
        end if
      end do
    else
      do i = 1, natyp
        ixipol(i) = 0
      end do
      write (1337, 240)(ixipol(i), i=1, natyp)
      write (1337, 340)
    end if
    write (1337, 250) naez, nemb
    write (1337, 380)

    nineq = naez
    write (1337, 200) ncls, nref, nineq
    write (1337, 380)
    write (1337, 340)

    ! ----------------------------------------------------------------------------
    kmrot = 0

    do i = 1, naez
      ! -------------------------------------------------------------------------
      ! Atoms equivalent by inversional symmetry
      ! -------------------------------------------------------------------------
      qmtet(i) = 0d0
      qmphi(i) = 0d0
      ier = 0
      call ioinput('RBASISANG       ', uio, i, 7, ier)

      if (ier==0) then
        read (unit=uio, fmt=*)(rbasis(j,i), j=1, 3), qmtet(i), qmphi(i)
        if (abs(qmtet(i))>1d-6) kmrot = 1
        if (abs(qmphi(i))>1d-6) kmrot = 1
      end if
    end do                         ! I=1,NAEZ
    call idreals(rbasis(1,1), 3*naez, iprint)
    ! ----------------------------------------------------------------------------

    if (nemb>0) write (1337, *)
    write (1337, 280)((rbasis(j,i),j=1,3), i, refpot(i), i=naez+1, naez+nemb)

    ! ------------------------------------------------------------------------
    if (.not. opt('VIRATOMS')) then
      do i = 1, naez
        do io = 1, noq(i)
          if (kaoez(io,i)<1) stop 'Error in KAOEZ'
        end do
      end do
    end if
    ! ------------------------------------------------------------------------
    write (1337, 390)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check for DECIMATE consistency
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (opt('DECIMATE')) then
      if (mod(nprincd,nlbasis)/=0) then
        write (6, *) ' Decimation cannot continue '
        write (6, *) 'NPRINCD=', nprincd, ' NLBASIS=', nlbasis
        stop
      end if
      if (mod(nprincd,nrbasis)/=0) then
        write (6, *) ' Decimation cannot continue '
        write (6, *) 'NPRINCD=', nprincd, ' NRBASIS=', nrbasis
        stop
      end if
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check for ITERMDIR consistency -- if KMROT=0 suppress it
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if ((opt('ITERMDIR')) .and. (kmrot==0)) then
      write (1337, *)
      write (1337, *) ' WARNING: ITERMDIR running option used with collinear/' &
        , 'parallel Oz starting'
      write (1337, *) '          system (KMROT = 0 ). Please check token', &
        ' RBASISANG in your input'
      write (1337, *) ' Running option ITERMDIR will be ignored'
      write (1337, *)
      do i = 1, 32
        if (t_params%optc(i)(1:8)=='ITERMDIR') t_params%optc(i) = '        '
      end do
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check for XCPL consistency
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    manctl = (kmrot==0) .and. (krel==0) .and. (nspin>1)
    if ((opt('XCPL    ')) .and. (.not. manctl)) then
      write (1337, *)
      write (1337, *) ' WARNING: XCPL running option requires collinear ', &
        'magnetic systems'
      write (1337, *) ' in a NON/SCALAR/SCALAR+SOC relativistic mode (KREL=0)'
      write (1337, *) ' Running option XCPL will be ignored'
      write (1337, *)
      do i = 1, 32
        if (t_params%optc(i)(1:8)=='XCPL    ') t_params%optc(i) = '        '
      end do
    end if

    write (1337, 110)(t_params%optc(i), i=1, 8)
110 format (79('-'), /, ' EXECUTION OPTIONS:', /, 1x, a8, 7('//',a8), /, &
      79('-'))
    write (1337, 120)(t_params%testc(i), i=1, 16)
120 format (79('-'), /, ' TEST OPTIONS:', /, 2(1x,a8,7('//', &
      a8),/), /, 79('-'))
130 format (8a8)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Initialise SOLVER, SOC and CTL parameters in REL case
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cscl(:, :) = cvlight
    mansoc = .false.
    manctl = .false.

    if (krel==1) then
      solver = 'BS        '
      call ioinput('SOLVER          ', uio, 0, 7, ier)
      if (ier==0) then
        read (unit=uio, fmt=*) solver
        if (solver(1:2)=='BS') then
          solver = 'BS        '
        else
          if (solver/='ABM-OP    ') solver = 'ABM-OP    '
        end if
      end if

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! SOC-MAN
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! For Dirac-ASA
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (opt('SOC     ')) then
        call ioinput('SOSCALE         ', uio, 0, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) soscale
          if (soscale>-2.5d0) then
            if (soscale>=0.0d0) then ! SOC-I
              solver = 'ABM-SOC   '
              mansoc = .true.
            else                   ! SOC-II
              solver = 'ABM-SOC-II'
              mansoc = .true.
              do i = 1, natyp
                socscl(1:lmax+1, i) = soscale
              end do
              write (1337, 960) socii(nint(soscale))
            end if
          else
            write (1337, 870) '< SOC >'
            write (1337, 890)
          end if
        else
          write (1337, 880) '< SOC >'
          write (1337, 890)
        end if

        if (mansoc .and. (soscale>=0d0)) then
          imansoc(1:natyp) = 1
          ! -------------------------------------------------------------------
          ! Now look for a possible include/exclude list (SOCLIST= +/- NASOC)
          ! if SOCLIST is not found, ALL the atoms will have SOC modified with
          ! SOSCALE (+NASOC=only NASOC atoms, -NASOC=all but these NASOC
          ! atoms)
          ! Note that this is allowed only for SOC-I manipulation
          ! -------------------------------------------------------------------
          call ioinput('SOCLIST         ', uio, 0, 7, ier)
          if (ier==0) then
            read (unit=uio, fmt=*) nasoc, (isp(i), i=1, abs(nasoc))

            if (nasoc/=0) then
              if (nasoc<0) then    ! exclude this atoms
                do i = 1, -nasoc
                  imansoc(isp(i)) = 0
                end do
              else
                imansoc(1:natyp) = 0
                do i = 1, nasoc
                  imansoc(isp(i)) = 1
                end do
              end if
            end if
          end if

          write (1337, 310)
          do i = 1, natyp
            if (imansoc(i)==1) then
              socscl(1:lmax+1, i) = soscale
            end if
          end do
          write (1337, 900, advance='no')
          if (nasoc==0) write (1337, 910)
          if (nasoc>0) then
            write (1337, 920)
            write (1337, 940)(isp(i), i=1, nasoc)
          end if
          if (nasoc<0) then
            write (1337, 930)
            write (1337, 940)(isp(i), i=1, abs(nasoc))
          end if
          write (1337, 950) soscale
          write (1337, 310)
        end if
      end if
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! SOC-MAN
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      write (1337, '('' SOLVER used for the DIRAC equation : '',2X,A)') solver
      write (1337, 310)

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CTL-MAN
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (opt('CSCALE  ')) then
        call ioinput('CTLSCALE        ', uio, 0, 7, ier)
        if (ier==0) then
          read (unit=uio, fmt=*) ctlscale
          if (ctlscale>=1d-12) then
            manctl = .true.
          else
            write (1337, 870) '< CSCALE >'
            write (1337, 970)
          end if
        else
          write (1337, 880) '< CSCALE >'
          write (1337, 970)
        end if

        if (manctl) then
          cscl(:, :) = cscl(:, :)/dsqrt(ctlscale)
          write (1337, 980, advance='no')
          write (1337, 910)
          write (1337, 950) 1.d0/dsqrt(ctlscale)
        end if
        write (1337, 310)
      end if

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! CTL-MAN
      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! LDA+U
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (opt('qdos    ')) then
      allocate (t_params%qdos_atomselect(natyp), stat=i_stat) ! INTEGER
      call memocc(i_stat, product(shape(t_params%qdos_atomselect))*kind( &
        t_params%qdos_atomselect), 't_params%qdos_atomselect', 'rinput13')

      t_params%qdos_atomselect(:) = 1
      ! for now this is not used. Later this should be used to speed up the
      ! qdos calculations if not all atoms are supposed to be calculated Then
      ! if fullinv was not chosen then tmatrix is only needed for the
      ! principle layer of the atom of interest and the calculation of G(k)
      ! can be done only on that subblock.
      ! CALL IoInput('qdosatoms       ',UIO,1,7,IER)
      ! IF (IER.EQ.0) THEN
      ! READ (UNIT=UIO,FMT=*) (t_params%qdos_atomselect(I),I=1,NATYP)
      ! WRITE(111,FMT='(A10,80I2)') 'qdosatoms=  ',
      ! (t_params%qdos_atomselect(I),I=1,NATYP)
      ! ELSE
      ! WRITE(111,FMT='(A18,80I2)') 'Default qdosatoms=  ',
      ! (t_params%qdos_atomselect(I),I=1,NATYP)
      ! ENDIF

      ! WRITE (1337,'(A)') 'atom selective writeout for qdos:'
      ! WRITE (1337,'(A,1000I5)') 'qdosatoms=',
      ! (t_params%qdos_atomselect(I),I=1,NATYP)

      if (.not. test('MPIatom ')) then
        ! enforce MPIenerg since this is usually faster for qdos option
        call addtest('MPIenerg')
      end if

    end if

    ! =============================================================         !
    ! fswrt
    ! check and correct some settings automatically for FERMIOUT writeout   !
    ! fswrt
    if (opt('FERMIOUT') .or. opt('OPERATOR')) then ! fswrt
      if (nsteps/=1) then          ! fswrt
        write (6, 170)             ! fswrt
        nsteps = 1                 ! fswrt
      end if                       ! fswrt
      if (.not. test('STOP1B  ')) call addtest('STOP1B  ') ! fswrt
      if (.not. test('STOP1B  ')) stop 'addtest failed for STOP1B' ! fswrt
    end if                         ! fswrt
    ! =============================================================         !
    ! fswrt

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! WF_SAVE
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call ioinput('MEMWFSAVE       ', uio, 0, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) t_wavefunctions%maxmem_number
      write (1337, *) '< MEMWFSAVE >', t_wavefunctions%maxmem_number
      write (111, *) 'MEMWFSAVE=', t_wavefunctions%maxmem_number
    else
      t_wavefunctions%maxmem_number = 0
      write (1337, *) '< MEMWFSAVE >, use default:', &
        t_wavefunctions%maxmem_number
      write (111, *) 'Default MEMWFSAVE= ', t_wavefunctions%maxmem_number
    end if
    call ioinput('UNITMEMWFSAVE   ', uio, 0, 7, ier)
    if (ier==0) then
      read (unit=uio, fmt=*) t_wavefunctions%maxmem_units
      write (1337, *) '< UNITMEMWFSAVE >', t_wavefunctions%maxmem_units, &
        ' (max memory= UNITMEMWFSAVE*1024**MEMWFSAVE)'
      write (111, *) 'UNITMEMWFSAVE=', t_wavefunctions%maxmem_units
    else
      t_wavefunctions%maxmem_units = 2
      write (1337, *) '< UNITMEMWFSAVE >, use default:', &
        t_wavefunctions%maxmem_units, &
        '(MB) (max memory= MEMWFSAVE*1024**UNITMEMWFSAVE)'
      write (111, *) 'Default UNITMEMWFSAVE= ', t_wavefunctions%maxmem_units, &
        '(MB)'
    end if

    ! default flags: save only rll from main1a>tmatnewsolver since left
    ! solutions can be calculated always in main1c>rhovalnew and sll is not
    ! used
    t_wavefunctions%save_rll = .true.
    t_wavefunctions%save_sll = .false.
    t_wavefunctions%save_rllleft = .false.
    t_wavefunctions%save_sllleft = .false.

    if (opt('OPERATOR')) then
      write (1337, *) 'Found option "OPERATOR"'
      write (1337, *) 'Overwrite MEMWFSAVE input with big numbers'
      t_wavefunctions%maxmem_number = 5
      t_wavefunctions%maxmem_units = 3
    end if
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! END WF_SAVE
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Begin Godfrin inversion scheme control                       ! GODFRIN Flaviano
    !==========================================================
    IF(OPT('godfrin ')) THEN
       WRITE(111 ,*) 'Godfrin inversion scheme parameters'
       WRITE(1337,*) 'Godfrin inversion scheme parameters'
   
       t_godfrin%na = NAEZD
       CALL IOINPUT('GODFRIN         ',UIO,2,7,IER)
       IF (IER.NE.0) STOP 'RINPUT: GODFRIN not found!'
       READ (UNIT=UIO,FMT=*) t_godfrin%nb, t_godfrin%ldiag, t_godfrin%lper, t_godfrin%lpardiso
   
       CALL IOINPUT('GODFRIN         ',UIO,4,7,IER)
       ALLOCATE(t_godfrin%bdims(t_godfrin%nb))
       READ (UNIT=UIO,FMT=*) t_godfrin%bdims(:)
   
      !Inconsistency check
      IF( t_godfrin%na /= sum(t_godfrin%bdims) ) stop 'godfrin: na /= sum(bdims)'
#ifdef __INTEL_COMPILER
      ! can only use pardiso solver with intel mkl at the moment, probably only
      ! a linking issue that should be solved in the future
      if( t_godfrin%lpardiso) stop 'No pardiso library available. Try the intel compiler or fix the linking issues'
#endif
     
      WRITE(111 ,FMT='(A100)') 'na, nb, ldiag, lper, lpardiso; then bdims(1:nb)'
      WRITE(1337,FMT='(A100)') 'na, nb, ldiag, lper, lpardiso; then bdims(1:nb)'
      WRITE(111 ,*) t_godfrin%na, t_godfrin%nb, t_godfrin%ldiag, t_godfrin%lper, t_godfrin%lpardiso
      WRITE(1337,*) t_godfrin%na, t_godfrin%nb, t_godfrin%ldiag, t_godfrin%lper, t_godfrin%lpardiso
      WRITE(111 ,FMT='(50(I0," "))') t_godfrin%bdims(:)
      WRITE(1337,FMT='(50(I0," "))') t_godfrin%bdims(:)
     
      !multiply blocks by angular momentum dimension
      t_godfrin%na    = t_godfrin%na*lmmaxd
      t_godfrin%bdims = t_godfrin%bdims*lmmaxd
     
      if(ICC/=0 .and. t_godfrin%ldiag) then
        t_godfrin%ldiag = .false.
        write(111 ,fmt='(A100)') 'rinput13: Warning! ICC/=0. Setting ldiag = T'
        write(1337,fmt='(A100)') 'rinput13: Warning! ICC/=0. Setting ldiag = T'
      end if 
      
    END IF
    ! End Godfrin inversion scheme control                         ! GODFRIN Flaviano
    !==========================================================

    write (1337, 310)
    write (1337, 300) kmrot
    write (1337, 380)
    write (1337, *) ' >>>>>>>>> RINPUT13 EXITS NOW <<<<<<<<<< '

    close (111)                    ! Close file inputcard_generated.txt

    if (allocated(imansoc)) then
      i_all = -product(shape(imansoc))*kind(imansoc)
      deallocate (imansoc, stat=i_stat)
      call memocc(i_stat, i_all, 'IMANSOC', 'rinput13')
    end if

    return
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! INPUT END
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
140 format ((f4.0,i4,4x,4i1,3i4,f8.4,i4,i5,1x,f8.5))
    ! ------------------------------------------------------------------------
150 format (' NSPIN ', /, i4)
160 format (' NSTEPS', /, i4)
170 format (' WARINING: Setting NSTEPS to 1 for runoption FERMOUT')
180 format ('          ALAT = ', f15.8)
190 format ('   INTERVX   INTERVY   INTERVZ', /, 3i10)
200 format ('    NCLS    NREF   NINEQ', /, 3i8)
210 format (' RBASIS', /, 'SITE                BASIS VECTORS                 ' &
      , 'THETA   PHI CPA OCC KAOEZ')
220 format ('         ABASIS         BBASIS         CBASIS', /, 3f15.8)
230 format (' INIPOL', /, (10i4))
240 format (' IXIPOL', /, (10i4))
250 format ('    NAEZ    NEMB  ', /, 2i8)
260 format ((i4,3f15.8,2f6.1,2(1x,i3),4i3))
270 format (' NATYP ', /, i4, /, &
      '   Z lmx     KFG cls pot ntc  MTFAC irns SITE  CONC')
280 format ((3f15.8,2i6))
290 format (' NTCELLR', /, (10i4))
300 format (' KMROT', /, 4i8)
    ! ------------------------------------------------------------------------
310 format (79('-'))
320 format (3('-'), '+', 3(14('-'),'+'), 30('-'))
330 format (3(9('-'),'+'), 49('-'))
340 format (10(3('-'),'+'), 39('-'))
350 format (3('-'), '+', 75('-'))
360 format (3(14('-'),'+'), 34('-'))
370 format (2(3('-'),'+'), 7('-'), '+', 3(3('-'),'+'), 7('-'), '+', 3('-'), &
      '+', 39('-'))
380 format (3(7('-'),'+'), 55('-'))
390 format (7(7('-'),'+'), 23('-'))
400 format (/, 33x, 'check of dimension-data consistency', /, 33x, 35('-'), /, &
      40x, 'lmax   : (', i6, ',', i6, ')', /, 40x, 'natyp  : (', i6, ',', i6, &
      ')', /, 40x, 'irm    : (', i6, ',', i6, ')', /, 40x, 'nspin  : (', i6, &
      ',', i6, ')', /)
410 format (1x, 10('*'), ' external magnetic field applied hfield=', f8.5)
420 format (20x, a4, 'spin polarized calculation')
430 format (1x, 20x, ' calculation with', a8, '-potential')
440 format (1x, 79('*'))
450 format (' mixing factor used           :', f15.6, /, &
      ' convergence quality required :', 1p, d15.2)
460 format (' make use of CPA algorithm    :', 1x, a14)
470 format ('         max. iterations      :', i15, /, &
      '         req. CPA convergency :', 1p, d15.2)
480 format (1x, 20x, a24, 'exchange-correlation potential')
490 format (/, 20x, 'broyden"s method # :', i3, &
      ' is used up to iteration-      ', /, 20x, 'depth :', i3, &
      '  then jacobian is fixed and potential      ', /, 20x, &
      'is updated using that jacobian')
500 format (13x, ' in case of calculating non - spherical wavefcts ', &
      'the parameter lmaxd has to be set equal lmax ')
510 format (/)
520 format (20x, 'full potential calculation ', &
      '- cut off of non spherical potential', /, ' >', /)
530 format (31x, 'representive atom no.', i3, ' irns :', i5, ' irnsd :', i5)
540 format (21x, a43, /, 21x, ' using', i3, '-th. born approximation ')
550 format (21x, a43)
560 format (' lmax', /, i4)
570 format ('          EMIN        EMAX        TK', /, 3f12.6)
580 format ('   NPOL  NPNT1  NPNT2  NPNT3', /, 4i7)
590 format ('  IFILE    IPE ISHIFT ESHIFT', /, 3i7, f12.6)
600 format (' KSHAPE    IRM    INS   ICST INSREF', /, 5i7)
610 format ('   KCOR  KVREL    KWS   KHYP KHFIELD   KXC', /, 6i7)
620 format (' external magnetic hfield     :', f15.4, /, &
      ' VCONST                       :', f15.6)
630 format ('   IMIX    IGF    ICC', /, 3i7)
640 format (' ITDBRY', /, i7)
650 format ('      STRMIX        FCM       QBOUND', /, 3f12.6)
660 format ('      BRYMIX', /, f12.6)
670 format ('    KTE   KPRE   KEFG  KVMAD ', /, 5i7)
680 format (3('-'), '+', 75('-'))
690 format (3(11('-'),'+'), 43('-'))
700 format (3(6('-'),'+'), 58('-'))
710 format (4(6('-'),'+'), 51('-'))
720 format (3(6('-'),'+'), 11('-'), '+', 46('-'))
730 format (6(6('-'),'+'), 37('-'))
740 format (6('-'), '+', 72('-'))
750 format (11('-'), '+', 67('-'))
760 format (5(6('-'),'+'), 44('-'))
770 format ('*** SLAB - INTERFACE CALCULATION ***', /)
780 format (i5, 3f14.8, i5)
790 format ('Number of LEFT  Host Layers : ', i5, ' with ', i5, ' basis')
800 format ('Number of RIGHT Host Layers : ', i5, ' with ', i5, ' basis')
810 format ('Left  side periodicity : ', 3f10.5)
820 format ('Right side periodicity : ', 3f10.5)
830 format ('    Geommetry used : ', /, &
      ' ATOM       TX          TY          TZ ')
840 format ('--------------- Left  Host -------------- ')
850 format ('---------------   S L A B  -------------- ')
860 format ('--------------- Right Host -------------- ')
870 format (/, 1x, 'WARNING: Option ', a, ' used with an INVALID ', &
      'scaling parameter.')
880 format (/, 1x, 'WARNING: Option ', a, ' found but NO value given for the', &
      ' scaling parameter.')
890 format (15x, '++++++++++   SOC option will be IGNORED   ++++++++++', /, &
      1x, 'Please use SOCSCALE= XXX (real>-2.5) in the inputcard', &
      ' to make your option valid ', /)
900 format (1x, 'The SOC will be SCALED')
910 format (' for ALL the atoms in the unit cell.')
920 format (' for the FOLLOWING atoms in the unit cell :')
930 format (' for all the atoms in the unit cell EXCLUDING :')
940 format (1x, 6(2x,i3))
950 format (1x, 'Scaling factor = ', 1p, d9.2)
960 format (1x, 'The SOC is manipulated', ' -- part of the SOC kept: ', a)
970 format (15x, '+++++++++  CSCALE option will be IGNORED  ++++++++++', /, &
      1x, 'Please use CTLSCALE= X (real>=1D-12) in the inputcard', &
      ' to make your option valid ', /)
980 format (1x, 'The CLIGHT will be SCALED')

  end subroutine rinput13
  ! ---------------------------------------------------------------------
  ! ---------------------------------------------------------------------

  subroutine addopt(string)
    use :: mod_wunfiles, only: t_params

    implicit none

    integer, parameter :: noptd = 32
    character (len=8) :: string
    integer :: ii
    logical, external :: opt

    if (.not. opt('        ')) then
      write (*, *) 'Error in ADDOPT for ', string, &
        ' : No free slots in array OPTC.'
      stop 'Error in ADDOPT: No free slots in array OPTC.'
    end if

    if (.not. opt(string)) then
      ii = 1
      do while (ii<=noptd)
        if (t_params%optc(ii)=='        ') then
          t_params%optc(ii) = string
          ii = noptd + 1
        end if
        ii = ii + 1
      end do
    end if

  end subroutine addopt

  subroutine addtest(string)
    use :: mod_types, only: t_inc
    use :: mod_wunfiles, only: t_params
    implicit none
    integer, parameter :: ntstd = 64
    character (len=8) :: string
    integer :: ii
    logical, external :: test

    if (t_inc%i_write>0) write (1337, *) 'in ADDTEST: adding option ', string

    if (.not. test('        ')) then
      write (*, *) 'Error in ADDTEST for ', string, &
        ' : No free slots in array TESTC.'
      stop 'Error in ADDTEST: No free slots in array TESTC.'
    end if

    if (.not. test(string)) then
      ii = 1
      do while (ii<=ntstd)
        if (t_params%testc(ii)=='        ') then
          t_params%testc(ii) = string
          ii = ntstd + 1
        end if
        ii = ii + 1
      end do
    end if

  end subroutine addtest

end module rinput
