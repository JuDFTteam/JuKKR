#include "macros.h"

! TODO: if energy_mesh.0 is not generated, check if settings in
! input.conf match those in energy_mesh.0
! PROBLEM: Fermi-Energy specified in 'potential' file - E-mesh depends on
! E-Fermi -> k-mesh depends on E-mesh => k-mesh depends on EFERMI
! NOTE: k-mesh ALWAYS DEPENDS ON THE FERMI ENRGY FROM 'potential' NOT ON THE ACTUAL ONE!!!!!! - BUG?
!  program MAIN0
module kkr0_mod
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: main0
  
  contains
  
  subroutine main0()

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

    use InputParams_mod, only: InputParams, getInputParamsValues, writeInputParamsToFile
    use DimParams_mod, only: DimParams, createDimParamsFromFile, writeDimParams
    use Main2Arrays_mod, only: Main2Arrays, createMain2Arrays, writeMain2Arrays
    use Warnings_mod, only: get_number_of_warnings, show_warning_lines

    integer, parameter :: NSYMAXD=48 ! Maximal number of Brillouin zone symmetries, 48 is largest possible number
    integer, parameter :: MAXMSHD=8 ! Maximal number of k-meshes used

!     .. Energy Mesh ..
    double precision :: EFERMI

    integer :: IELAST
    integer :: iesemicore

    double complex, allocatable :: EZ(:)
    double complex, allocatable :: WEZ(:)

    double precision :: VOLUME0
    double precision :: RECBV(3,3)

    integer,   allocatable :: NTCELL(:)
    double precision, allocatable :: radius_muffin_tin(:)

!     .. auxillary variables, not passed to kkr2
    double precision :: PI
    integer :: IE
    integer :: ierror
    double complex, allocatable :: DEZ(:) ! needed for EMESHT
    integer :: IEMXD
    integer, parameter :: KREL = 0
    integer :: EKMD
    logical :: startpot_exists
    type (InputParams)    :: params
    type (DimParams)      :: dims
    type (Main2Arrays)    :: arrays

! ------------ end of declarations ---------------------------------

    call createDimParamsFromFile(dims, "global.conf")

    ierror = getInputParamsValues("input.conf", params)
    if (ierror /= 0) stop

    ! Calculate number of energy points

    ! semicore contour is used
    if (params%use_semicore == 1) then

      write(*,*) "WARNING: Using semicore contour because use_semicore=1 in input.conf. This feature is still subject to beta testing"

      if (params%NPOL /= 0) then
        IEMXD = params%NPOL + params%NPNT1 + params%NPNT2 + params%NPNT3 + params%n1semi + params%n2semi + params%n3semi
      else ! DOS-calculation
        IEMXD = params%NPNT2 + params%n2semi
      endif

    ! semicore contour is not used
    else

      if (params%NPOL /= 0) then
        IEMXD = params%NPOL + params%NPNT1 + params%NPNT2 + params%NPNT3
      else ! DOS-calculation
        IEMXD = params%NPNT2
      endif

    endif


    dims%IEMXD = IEMXD

    ! important: determine IEMXD before creating arrays
    call createMain2Arrays(arrays, dims)

!-----------------------------------------------------------------------------
! Array allocations BEGIN
!-----------------------------------------------------------------------------
    allocate(NTCELL(dims%NAEZ))
    allocate(radius_muffin_tin(dims%naez))
!-----------------------------------------------------------------------------
! Array allocations END
!-----------------------------------------------------------------------------

    call RINPUTNEW99(arrays%RBASIS, arrays%NAEZ)

!     in case of a LDA+U calculation - read file 'ldauinfo'
!     and write 'wldau.unf', if it does not exist already
    if (params%LDAU) then
      call ldauinfo_read(dims%LMAXD, dims%NSPIND, arrays%ZAT, dims%NAEZ)
    endif

!===================================================================

    ! read starting potential and shapefunctions
    startpot_exists = .false.
    inquire(file='potential', exist=startpot_exists)
    ! if energy_mesh.0 file is missing, also regenerate start files
    if (startpot_exists) then

      call STARTB1_wrapper_new(params%alat, dims%NSPIND, NTCELL, &
                               EFERMI, arrays%ZAT, radius_muffin_tin, &
                               dims%NAEZ)

    else
      ! no formatted potential provided
      write(*,*) "WARNING: file 'potential' not found... skipping start potential generation."
      warn(6, "file 'potential' not found... skipping start potential generation.")
      write(*,*) "Trying to read initial, approximate EFermi from EFERMI file..."
      open (67, file='EFERMI', form='formatted')
      read (67, *) EFERMI
      close(67)
    endif

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
    if (params%NPOL /= 0) params%Emax = EFERMI

! --> set up energy contour
! --> set up semicore energy contour if use_semicore == 1
    iesemicore = 0
    if (params%use_semicore == 1) then
! EPATHTB calls EMESHT both for the semicore contour and the valence contour
      call EPATHTB(EZ,DEZ,EFERMI,IELAST,iesemicore,params%use_semicore, &
                  params%emin,params%emax,params%tempr,params%npol,params%npnt1,params%npnt2,params%npnt3, &
                  params%ebotsemi,params%emusemi,params%tempr,params%npol,params%n1semi,params%n2semi,params%n3semi, &
                  IEMXD)
    else
! Call EMESTH for valence contour only (can be included in EPATHTB when semicore contour feature is stable)
      call EMESHT(EZ,DEZ,IELAST,params%Emin,params%Emax,EFERMI,params%tempr, &
                  params%NPOL,params%NPNT1,params%NPNT2,params%NPNT3,IEMXD)
    endif

    PI = 4.0D0*ATAN(1.0D0)

    do IE = 1, IELAST
      WEZ(IE) = -2.D0/PI*DEZ(IE)
      IF (IE <= IESEMICORE) WEZ(IE) = WEZ(IE)*params%FSEMICORE
    enddo ! ie


! ================================================ deal with the lattice

    arrays%BRAVAIS(:,1) = params%bravais_a
    arrays%BRAVAIS(:,2) = params%bravais_b
    arrays%BRAVAIS(:,3) = params%bravais_c

    ! only for informative purposes - prints info about lattice
    call LATTIX99(params%ALAT,arrays%BRAVAIS,RECBV,VOLUME0, .true.)


    call SCALEVEC(arrays%RBASIS, arrays%NAEZ,arrays%BRAVAIS,params%CARTESIAN)

! ======================================================================
!     setting up kpoints
! ======================================================================

    call BZKINT0(arrays%NAEZ, &
                 arrays%RBASIS,arrays%BRAVAIS,RECBV, &
                 arrays%NSYMAT,arrays%ISYMINDEX, &
                 arrays%DSYMLL, &
                 params%bzdivide(1),params%bzdivide(2),params%bzdivide(3), &
                 IELAST,EZ,arrays%KMESH,arrays%MAXMESH,MAXMSHD, &
                 dims%LMAXD, IEMXD, KREL, arrays%KPOIBZ, EKMD)

    ! after return from bzkint0, EKMD contains the right value
    dims%EKMD = EKMD

    ! bzkint0 wrote a file 'kpoints': read this file and use it as k-mesh
    call readKpointsFile(arrays%BZKP, arrays%MAXMESH, arrays%NOFKS, arrays%VOLBZ, arrays%VOLCUB)

!     Conversion of RMAX and GMAX to atomic units
    params%RMAX = params%RMAX*params%ALAT
    params%GMAX = params%GMAX/params%ALAT

    call TESTDIMLAT(params%ALAT,arrays%BRAVAIS,RECBV,params%RMAX,params%GMAX, &
                    dims%NMAXD, dims%ISHLD)

    call writeDimParams(dims, 'inp0.unf')
    ierror = writeInputParamsToFile('input.unf', params)

    call writeMain2Arrays(arrays, 'arrays.unf')

    ! write start energy mesh
    if (params%use_semicore == 1) then
      open  (67,FILE='energy_mesh.0',FORM='unformatted')
      write (67) IELAST,EZ,WEZ,params%Emin,params%Emax
      write (67) params%NPOL,params%tempr,params%NPNT1,params%NPNT2,params%NPNT3
      write (67) EFERMI
      write (67) IESEMICORE,params%FSEMICORE,params%EBOTSEMI
      write (67) params%EMUSEMI
      write (67) params%N1SEMI,params%N2SEMI,params%N3SEMI
      close (67)
    else
      open  (67,FILE='energy_mesh.0',FORM='unformatted')
      write (67) IELAST,EZ,WEZ,params%Emin,params%Emax
      write (67) params%NPOL,params%tempr,params%NPNT1,params%NPNT2,params%NPNT3
      write (67) EFERMI
      close (67)
    endif

! ======================================================================

! deallocations

    ! Energy mesh
    deallocate(EZ, stat=ierror)
    deallocate(WEZ, stat=ierror)

    !   auxillary
    deallocate(DEZ, stat=ierror)
    deallocate(NTCELL, stat=ierror)
    deallocate(radius_muffin_tin, stat=ierror)

    
! -------------- Helper routine -----------------------------------------------

    if (get_number_of_warnings() > 0) ierror = show_warning_lines(unit=6)
    
  endsubroutine ! main0

  !------------------------------------------------------------------------
  ! Read k-mesh file
  subroutine readKpointsFile(BZKP, MAXMESH, NOFKS, VOLBZ, VOLCUB)
    double precision, intent(out) :: BZKP(:,:,:)
    integer, intent(in) :: MAXMESH
    integer, intent(out) :: NOFKS(:)
    double precision, intent(out) :: VOLBZ(:)
    double precision, intent(out) :: VOLCUB(:,:)

    ! -----------------------------
    integer, parameter :: fu = 52 ! file unit
    integer :: I
    integer :: ID
    integer :: L
    logical :: new_kpoints

    new_kpoints = .false.
    inquire(file='new.kpoints',exist=new_kpoints)

    if (.not. new_kpoints) then
      open (fu,file='kpoints',form='formatted')
    else
      ! if file new.kpoints exists - use those kpoints
      write(*,*) "WARNING: rejecting file kpoints - using file new.kpoints instead."
      warn(6, "rejecting file kpoints - using file new.kpoints instead.")
      open (fu,file='new.kpoints',form='formatted')
    endif

    rewind (fu)

    do L = 1, MAXMESH
      read (fu,fmt='(I8,f15.10)') NOFKS(L),VOLBZ(L)
      do I = 1, NOFKS(L)
        read (fu,fmt=*) BZKP(1:3,I,L),VOLCUB(I,L)
      enddo ! i
    enddo ! l

    close (fu)
  endsubroutine ! readKpointsFile

endmodule ! kkr0_mod

