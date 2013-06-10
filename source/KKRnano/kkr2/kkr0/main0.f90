  program MAIN0

! Explanation of most variables follows below

!     ALAT                     : lattice constant (in a.u.)
!     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
!     E1,E2,                   : energies needed in EMESHT
!     HFIELD                   : external magnetic field, for
!                              : initial potential shift in
!                              : spin polarised case
!     TK,                      : temperature
!     VCONST,                  : potential shift
!     BRAVAIS(3,3),            : bravais lattice vectors
!     RECBV(3,3),              : reciprocal basis vectors
!     RMTREF(NREFD),           : muffin-tin radius of reference system
!     RBASIS(3,NAEZD),         : position of atoms in the unit cell
!                              : in units of bravais vectors
!     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
!     RR(3,0:NRD)              : set of real space vectors (in a.u.)
!     VBC(2),                  : potential constants
!     WG(LASSLD),              : integr. weights for Legendre polynomials
!     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
!     ZAT(NAEZD)               : nuclear charge
!     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction
!                              : for k-net in IB of the BZ
!     ICST,                    : number of Born approximation
!     IEND,                    : number of nonzero gaunt coeffizients
!     IFILE,                   : unit specifier for potential card
!     IPE,IPF,IPFE,            : not real used, IPFE should be 0
!     KHFELD,                  : 0,1: no / yes external magnetic field
!     KVREL,                   : 0,1 : non / scalar relat. calculation
!     LMAX,                    : maximum l component in
!                              : wave function expansion
!     LPOT,                    : maximum l component in
!                              : potential expansion
!     NAEZ,                    : number of atoms in unit cell
!     NCLS,                    : number of reference clusters
!     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
!     NPOL,                    : number of Matsubara Pols (EMESHT)
!     NR,                      : number of real space vectors rr
!     NREF,                    : number of diff. ref. potentials
!     NSPIN,                   : counter for spin directions
!     IGUESS                   : 0,1 : no / yes (sc) initial guess, set
!                              : IGUESSD to 1 in inc.p if needed
!     BCP                      : 0,1 : no / yes bc-preconditioning, set
!                              : BCPD to 1 in inc.p if needed
!     QBOUND                   : exit condition for self-consistent
!                              : iteration
!     QMRBOUND                 : exit condition for QMR iterations
!     SCFSTEPS                 : number of scf iterations
!     INIPOL(NAEZD),           : initial spin polarisation

!     ATOM(NACLSD,NAEZD),      : atom at site in cluster
!     CLS(NAEZD),              : cluster around atom
!     NACLS(NCLSD),            : number of atoms in cluster
!     EZOA(NACLSD,NAEZD),      : EZ of atom at site in cluster
!     RMT(NAEZD)               : muffin-tin radius of true system
!     RMTNEW(NAEZD)            : adapted muffin-tin radius
!     RWS(NAEZD)               : Wigner Seitz radius
!     ICST                     : the regular non spherical wavefunctions, the
!                              : alpha matrix and the t-matrix in the ICST-th. born approx
!     IMT(NAEZD),              : r point at MT radius
!     IPAN(NAEZD),             : number of panels in non-MT-region
!     IRC(NAEZD),              : r point for potential cutting
!     IRCUT(0:IPAND,NAEZD),    : r points of panel borders
!     IRMIN(NAEZD),            : max r for spherical treatment
!     IRNS(NAEZD)              : number r points for non spher. treatm.
!     IRWS(NAEZD),             : r point at WS radius
!     LMSP(NAEZD,LMXSPD)       : 0,1 : non/-vanishing lm=(l,m) component
!                              : of non-spherical potential
!     LLMSP(NAEZD,NFUND)       : lm=(l,m) of 'nfund'th nonvanishing
!                              : component of non-spherical pot.
!     NTCELL(NAEZD),           : index for WS cell

!     A(NAEZD),B(NAEZD)        : contants for exponential r mesh
!     R(IRMD,NAEZD)            : radial mesh ( in units a Bohr)
!     DRDI(IRMD,NAEZD)         : derivative dr/di
!     THETAS(IRID,NFUND,NCELLD): shape function
!                              :         ( 0 outer space
!                              : THETA = (
!                              :         ( 1 inside WS cell
!                              : in spherical harmonics expansion
!     LCORE(20,NPOTD)          : angular momentum of core states
!     NCORE(NPOTD)             : number of core states
! ----------------------------------------------------------------------

    use Config_Reader
    use InputParams_mod
    use DimParams_mod
    use Main2Arrays_mod

    implicit none

!   double precision KIND constant
    integer, parameter :: DP = kind(1.0D0)

!     .. Parameters ..

    integer :: NSYMAXD
    integer :: MAXMSHD
!
!     Maximal number of Brillouin zone symmetries, 48 is largest
!     possible number
    parameter (NSYMAXD=48)
!     Maximal number of k-meshes used
    parameter (MAXMSHD=8)

!     .. Energy Mesh ..
    double precision :: EFERMI

    integer :: IELAST

    complex(kind=DP), dimension(:), allocatable :: EZ
    complex(kind=DP), dimension(:), allocatable :: WEZ

    double precision :: VOLUME0
    double precision :: RECBV(3,3)

    integer, dimension(:),   allocatable :: NTCELL
    double precision, dimension(:), allocatable :: RMTref

!     .. auxillary variables, not passed to kkr2
    double precision :: PI

    integer :: IE

    integer :: ierror

    complex(kind=DP), dimension(:), allocatable :: DEZ ! needed for EMESHT

    integer ::   IEMXD
    integer, parameter :: KREL = 0

    integer :: EKMD

    type (InputParams)    :: params
    type (DimParams)      :: dims
    type (Main2Arrays)    :: arrays

! ------------ end of declarations ---------------------------------

    call createDimParamsFromConf(dims)

    ierror = getInputParamsValues("input.conf", params)
    if (ierror /= 0) stop

    if (params%NPOL /= 0) then
      IEMXD = params%NPOL + params%NPNT1 + params%NPNT2 + params%NPNT3
    else
      ! DOS-calculation
      IEMXD = params%NPNT2
    end if

    dims%IEMXD = IEMXD

    ! important: determine IEMXD before creating arrays
    call createMain2Arrays(arrays, dims)

!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
    allocate(NTCELL(dims%NAEZ))
    allocate(RMTref(dims%naez))
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

     call RINPUTNEW99(arrays%RBASIS, NTCELL, &
                      arrays%NAEZ,arrays%ZAT, &
                      arrays%RMTref)

!     in case of a LDA+U calculation - read file 'ldauinfo'
!     and write 'wldau.unf', if it does not exist already
!    if (LDAU) then
!      call ldauinfo_read(LMAXD, NSPIND, ZAT, NAEZD)
!    end if

!===================================================================

    ! read starting potential and shapefunctions
    call STARTB1_wrapper(params%alat, &
                 dims%LPOT,dims%NSPIND,NTCELL, &
                 EFERMI, arrays%ZAT, arrays%RMTref, &
                 dims%IPAND, dims%IRID, dims%NFUND, dims%IRMD, dims%NCELLD, &
                 dims%NAEZ, dims%IRNSD)

! ----------------------------------------------------------------------
! update Fermi energy, adjust energy window according to running options

    IELAST = IEMXD

    ! Energy mesh
    allocate(EZ(IEMXD), stat=ierror)
    allocate(WEZ(IEMXD), stat=ierror)

    !   auxillary
    allocate(DEZ(IEMXD), stat=ierror)

! for non-DOS calculation upper energy bound corresponds to Fermi energy
! BAD: params%Emax is changed
    if ( params%NPOL /= 0 ) params%Emax = EFERMI

! --> set up energy contour
    PI = 4.0D0*ATAN(1.0D0)

    call EMESHT(EZ,DEZ,IELAST,params%Emin,params%Emax,EFERMI,params%tempr, &
    params%NPOL,params%NPNT1,params%NPNT2,params%NPNT3,IEMXD)
    do IE = 1,IELAST
      WEZ(IE) = -2.D0/PI*DEZ(IE)
    end do

! ================================================ deal with the lattice

    arrays%BRAVAIS(:,1) = params%bravais_a
    arrays%BRAVAIS(:,2) = params%bravais_b
    arrays%BRAVAIS(:,3) = params%bravais_c

    ! only for informative purposes - prints info about lattice
    call LATTIX99(params%ALAT,arrays%BRAVAIS,RECBV,VOLUME0, .true.)


    call SCALEVEC(arrays%RBASIS,params%basisscale(1), &
                  params%basisscale(2),params%basisscale(3), &
                  arrays%NAEZ,arrays%BRAVAIS,params%CARTESIAN)

! ======================================================================
!     setting up kpoints
! ======================================================================

    call BZKINT0(arrays%NAEZ, &
                 arrays%RBASIS,arrays%BRAVAIS,RECBV, &
                 arrays%NSYMAT,arrays%ISYMINDEX, &
                 arrays%DSYMLL, &
                 params%bzdivide(1),params%bzdivide(2),params%bzdivide(3), &
                 IELAST,EZ,arrays%KMESH,arrays%MAXMESH,MAXMSHD, &
                 arrays%LMAXD, IEMXD, KREL, arrays%KPOIBZ, EKMD)

    ! after return from bzkint0, EKMD contains the right value
    dims%EKMD = EKMD

    call readKpointsFile(arrays%BZKP, arrays%MAXMESH, arrays%NOFKS, &
                         arrays%VOLBZ, arrays%VOLCUB)

!    if (BCPD == 1 .and. NATBLD*XDIM*YDIM*ZDIM /= NAEZD) then
!      write(*,*) "ERROR: When BCPD==1 then NATBLD*XDIM*YDIM*ZDIM has to be equal to NAEZD."
!      stop
!    endif

!     Conversion of RMAX and GMAX to units of ALAT
    params%RMAX = params%RMAX*params%ALAT
    params%GMAX = params%GMAX/params%ALAT

    call TESTDIMLAT(params%ALAT,arrays%BRAVAIS,RECBV,params%RMAX,params%GMAX, &
                    dims%NMAXD, dims%ISHLD)

    call writeDimParams(dims)
    ierror = writeInputParamsToFile('input.unf', params)

    call writeMain2Arrays(arrays, 'arrays.unf')

    ! write energy mesh
    open (67,FILE='energy_mesh',FORM='unformatted')
    write (67) IELAST,EZ,WEZ,params%Emin,params%Emax
    write (67) params%NPOL,params%tempr,params%NPNT1,params%NPNT2,params%NPNT3
    write (67) EFERMI
    close (67)

    write(*,*)
    write(*,*) "WARNING: Assumes identical reference clusters!"

! ======================================================================

! deallocations

    ! Energy mesh
    deallocate(EZ, stat=ierror)
    deallocate(WEZ, stat=ierror)

    !   auxillary
    deallocate(DEZ, stat=ierror)
    deallocate(NTCELL, stat=ierror)
    deallocate(RMTref, stat=ierror)

! -------------- Helper routine -----------------------------------------------

    CONTAINS

      !------------------------------------------------------------------------
      ! Read k-mesh file
      subroutine readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)
        implicit none
        double precision :: BZKP(:,:,:)
        integer :: MAXMESH
        integer :: NOFKS(:)
        double precision :: VOLBZ(:)
        double precision :: VOLCUB(:,:)

        ! -----------------------------
        integer :: I
        integer :: ID
        integer :: L

        open (52,file='kpoints',form='formatted')
        rewind (52)

        do L = 1,MAXMESH
          read (52,fmt='(I8,f15.10)') NOFKS(L),VOLBZ(L)
          read (52,fmt=*) (BZKP(ID,1,L),ID=1,3),VOLCUB(1,L)
          do I=2,NOFKS(L)
            read (52,fmt=*) (BZKP(ID,I,L),ID=1,3),VOLCUB(I,L)
          end do
        end do

        close (52)
      end subroutine

end program

