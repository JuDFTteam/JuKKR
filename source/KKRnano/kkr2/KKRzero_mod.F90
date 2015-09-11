
! TODO: if energy_mesh.0 is not generated, check if settings in
! input.conf match those in energy_mesh.0
! PROBLEM: Fermi-Energy specified in 'potential' file - E-mesh depends on
! E-Fermi -> k-mesh depends on E-mesh => k-mesh depends on EFERMI
! NOTE: k-mesh ALWAYS DEPENDS ON THE FERMI ENRGY FROM 'potential' NOT ON THE ACTUAL ONE!!!!!! - BUG?
!  program MAIN0
module KKRzero_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: main0
  
  contains
  
  subroutine main0(checkmode)

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
!     R(IRMD,NAEZD)            : radial mesh (in units a Bohr)
!     DRDI(IRMD,NAEZD)         : derivative dr/di
!     THETAS(IRID,NFUND,NCELLD): shape function
!                              :         (0 outer space
!                              : THETA = (
!                              :         (1 inside WS cell
!                              : in spherical harmonics expansion
!     LCORE(20,NPOTD)          : angular momentum of core states
!     NCORE(NPOTD)             : number of core states
! ----------------------------------------------------------------------

    use InputParams_mod, only: InputParams, getInputParamsValues, writeInputParamsToFile
    use DimParams_mod, only: DimParams, createDimParamsFromFile, writeDimParams
    use Main2Arrays_mod, only: Main2Arrays, createMain2Arrays, writeMain2Arrays
    use BrillouinZone_mod, only: bzkint0
    use Warnings_mod, only: get_number_of_warnings, show_warning_lines
    use Lattice_mod, only: lattix99
    use EnergyMeshHelpers_mod, only: emesht, epathtb
    use Startb1_mod, only: startb1_wrapper_new

    integer, intent(in) :: checkmode ! 0: usual kkr0, >0: checks only
    
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

    integer, allocatable :: NTCELL(:)
    double precision, allocatable :: radius_muffin_tin(:)

!     .. auxillary variables, not passed to kkr2
    double precision :: PI
    integer :: IE
    integer :: ierror
    double complex, allocatable :: DEZ(:) ! needed for EMESHT
    integer :: IEMXD
    integer, parameter :: KREL = 0
    logical :: startpot_exists
    type(InputParams)    :: params
    type(DimParams)      :: dims
    type(Main2Arrays)    :: arrays

! ------------ endof declarations ---------------------------------

    call createDimParamsFromFile(dims, "global.conf")

    ierror = getInputParamsValues("input.conf", params)
    if (ierror /= 0) die_here("failed to read ''input.conf''!")

    ! Calculate number of energy points

    ! semicore contour is used
    if (params%use_semicore == 1) then

      write(*,*) "WARNING: Using semicore contour because use_semicore=1 in input.conf. This feature is still subject to beta testing"
      warn(6, "using semicore contour because use_semicore=1 in input.conf. This feature is still subject to beta testing")
      
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

!     in case of a LDA+U calculation - read file 'ldauinfo' and write 'wldau.unf', if it does not exist already
    if (params%LDAU) then
      call ldauinfo_read(dims%LMAXD, dims%NSPIND, arrays%ZAT, dims%NAEZ)
    endif

!===================================================================

    ! read starting potential and shapefunctions
    startpot_exists = .false.
    inquire(file='potential', exist=startpot_exists)
    ! if energy_mesh.0 file is missing, also regenerate start files
    if (startpot_exists) then

      call STARTB1_wrapper_new(params%alat, dims%NSPIND, NTCELL, EFERMI, arrays%ZAT, radius_muffin_tin, dims%NAEZ, nowrite=(checkmode /= 0))

    else
      ! no formatted potential provided
      write(*,*) "WARNING: file 'potential' not found... skipping start potential generation."
      warn(6, "file 'potential' not found... skipping start potential generation.")
      write(*,*) "Trying to read initial, approximate EFermi from EFERMI file..."
      open (67, file='EFERMI', form='formatted', action='read', status='old')
      read (67, *) EFERMI
      close(67)
    endif

! ----------------------------------------------------------------------
! update Fermi energy, adjust energy window according to running options

    IELAST = IEMXD
    allocate(EZ(IEMXD), WEZ(IEMXD), DEZ(IEMXD), stat=ierror) ! Energy mesh, DEZ is aux.

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
    call lattix99(params%ALAT, arrays%BRAVAIS, RECBV, VOLUME0, .true.)


    call SCALEVEC(arrays%RBASIS, arrays%NAEZ, arrays%BRAVAIS, params%CARTESIAN)

! ======================================================================
!     setting up kpoints
! ======================================================================

    call BZKINT0(arrays%NAEZ, arrays%RBASIS, arrays%BRAVAIS,RECBV, arrays%NSYMAT, arrays%ISYMINDEX, &
                 arrays%DSYMLL, params%bzdivide, IELAST, EZ, arrays%KMESH, arrays%MAXMESH, MAXMSHD, &
                 dims%LMAXD, IEMXD, KREL, arrays%KPOIBZ, dims%EKMD, nowrite=(checkmode /= 0)) ! after return from bzkint0, EKMD contains the right value

    ! bzkint0 wrote a file 'kpoints': read this file and use it as k-mesh
    call readKpointsFile(arrays%MAXMESH, arrays%NOFKS, arrays%BZKP, arrays%VOLCUB, arrays%VOLBZ)
    
!   Conversion of RMAX and GMAX to atomic units
    params%RMAX = params%RMAX*params%ALAT
    params%GMAX = params%GMAX/params%ALAT

    call TESTDIMLAT(params%ALAT, arrays%BRAVAIS, RECBV, params%RMAX, params%GMAX, dims%NMAXD, dims%ISHLD)

    if (checkmode == 0) then 
      ! write binary files that are needed in the main program
    
      call writeDimParams(dims, 'inp0.unf')
      ierror = writeInputParamsToFile(params, 'input.unf')
      call writeMain2Arrays(arrays, 'arrays.unf')

      ! write start energy mesh
        open (67, file='energy_mesh.0', form='unformatted', action='write')
        write(67) IELAST,EZ,WEZ,params%Emin,params%Emax
        write(67) params%NPOL,params%tempr,params%NPNT1,params%NPNT2,params%NPNT3
        write(67) EFERMI
      if (params%use_semicore == 1) then
        write(67) IESEMICORE,params%FSEMICORE,params%EBOTSEMI
        write(67) params%EMUSEMI
        write(67) params%N1SEMI,params%N2SEMI,params%N3SEMI
      endif ! semicore
        close(67)
        
     else  ! checkmode == 0
       write(*,'(A)') "CheckMode: binary files 'inp0.unf', 'input.unf' and arrays.unf' are not created!" ! do we need a warning here?
     endif ! checkmode == 0

! ======================================================================

! deallocations

    deallocate(EZ, WEZ, stat=ierror) ! Energy mesh
    deallocate(DEZ, NTCELL, radius_muffin_tin, stat=ierror) !   auxillary
    
! -------------- Helper routine -----------------------------------------------

    if (get_number_of_warnings() > 0) &
      ierror = show_warning_lines(unit=6)
    
  endsubroutine ! main0

  !------------------------------------------------------------------------
  ! Read k-mesh file
  subroutine readKpointsFile(maxmesh, nofks, bzkp, volcub, volbz)
    integer, intent(in) :: maxmesh
    integer, intent(out) :: nofks(:)
    double precision, intent(out) :: bzkp(:,:,:), volcub(:,:), volbz(:)

    ! -----------------------------
    integer, parameter :: fu = 52 ! file unit
    integer :: i, l, ios

    open (fu, file='new.kpoints', form='formatted', status='old', action='read', iostat=ios)
    if (ios /= 0) then
      ! default to the old kpoint file name
      open (fu, file='kpoints', form='formatted', status='old', action='read')
    else
      ! if file new.kpoints exists - use those kpoints
      write(*,*) "WARNING: rejecting file kpoints - using file new.kpoints instead."
      warn(6, "rejecting file kpoints - using file new.kpoints instead.")
    endif

    rewind (fu)
    do l = 1, maxmesh
      read (fu,fmt='(i8,f15.10)') nofks(l), volbz(l)
      do i = 1, nofks(l)
        read (fu,fmt=*) bzkp(1:3,i,l), volcub(i,l)
      enddo ! i
    enddo ! l
    close (fu)
    
  endsubroutine ! readKpointsFile

  
!>    Reads atominfo and rbasis files.
  subroutine rinputnew99(rbasis, naez)
    double precision, intent(out) :: rbasis(3,*)
    integer, intent(in) :: naez
    
    character(len=9), parameter :: version = "Sept 2015"
    integer :: i

!------------ array set up and definition of input parameter -----------

! ================================================================================
! |                                                                              |
! | KKRnano                                                                      |
! | Massively Parallel Screened Korringa-Kohn-Rostoker Electronic Structure Code |
! | for Bulk                                                                     |
! |                                                                              |
! | established Juelich 2008                          Version : Jun 2013         |
! |                                                                              |
! ================================================================================
      
    write(6, fmt="(/80(1h=))") 
    write(6, fmt="('|',78x,'|')") 
    write(6, fmt="('| KKRnano',70x,'|')") 
    write(6, fmt="('| Massively Parallel Screened Korringa-Kohn-Rostoker Electronic Structure Code |')") 
    write(6, fmt="('| for Bulk',69x,'|')") 
    write(6, fmt="('|',78x,'|')") 
    write(6, fmt="('| established Juelich 2008',26x,'Version : ',a9,8x,'|')") version
    write(6, fmt="('|',78x,'|')") 
    write(6, fmt="(80(1h=))") 
      
    open(77, file='rbasis', form='formatted', action='read', status='old') ! open input file containing the atomic positions
    do i = 1, naez
      read(unit=77, fmt=*) rbasis(1:3,i) ! reading atomic positions
    enddo ! i
    close(77)

    write(6, fmt="(3(7(1h-),1h+) ,55(1h-))")
    write(6, fmt="(10(3(1h-),1h+) ,39(1h-))")
     
#if 0     
!      the file 'atominfo' is no longer needed in kkrnano 
!      (z and rmt are read from 'potential' and ntcell is read from 'shapefun')

!      open(77, file='atominfo', form='formatted')
!      do i = 1, naez
!        read (unit=77,fmt=*) z(i), lmxcdummy, kfgdummy(1:4), cls, refpot, ntcell(i), mtfacdummy, irns_dummy, temp
!        radius_muffin_tin(i) = temp
!      enddo ! i
!      close(77)

!
!     No cleanup of Format statements due to nostalgic reasons.
!
 2010 FORMAT(' NSPIN '/I4)
 2011 FORMAT(' NSTEPS'/I4)
 2014 FORMAT('          ALAT = ',F15.8)
 2015 FORMAT('   INTERVX   INTERVY   INTERVZ'/3I10)
 2016 FORMAT('    NCLS    NREF   '/,2I8)
 2018 FORMAT(' RBASIS'/,'SITE                BASIS VECTORS                 ')
 2019 FORMAT('         ABASIS         BBASIS         CBASIS'/3F15.8)
 2025 FORMAT((i4,3F15.8))
 2028 FORMAT(' NAEZ ',/,I8)
! ------------------------------------------------------------------------
 2100 FORMAT(79(1H-))
 2101 format(  3(1H-),1H+  , 3(14(1H-),1H+),  30(1H-))
 2102 format(3(9(1H-),1H+) ,49(1H-))
 2103 FORMAT(10(3(1H-),1H+) ,39(1H-))
 2104 format(  3(1H-),1H+  ,75(1H-))
 2107 format(3(14(1H-),1H+),34(1H-))
 2110 format(3(7(1H-),1H+) ,55(1H-))
 9000 FORMAT (I2,3X,4I5)
 9010 FORMAT (1X,I1,1X,4I1)
 9020 FORMAT (/,33x,'check of dimension-data consistency',/,33x,35('-'),/,40x,'lmax   : (',i6,',',i6,')',/,40x,'NAEZ  : (',i6,',',i6,')',/,40x,'irm    : (',i6,',',i6,')',/,40x,'nspin  : (',i6,',',i6,')',/)
 9030 FORMAT (1x,10 ('*'),' external magnetic field applied hfield=',f8.5)
 9040 FORMAT (3f12.7)
 9050 FORMAT (20x,a4,'spin polarized calculation')
 9060 FORMAT (8i4)
 9070 FORMAT (1x,20x,' calculation with',a8,'-potential')
 9080 FORMAT (1x,79 ('*'))
 9090 FORMAT (' mixing factor used           :',f15.6,/,' convergence quality required :',1p,d15.2)
 9100 FORMAT (1x,20x,a24,'exchange-correlation potential')
 9110 FORMAT (/,20x,'broyden"s method # :',i3,' is used up to iteration-      ',/,20x,'depth :',i3,'  then jacobian is fixed and potential      ',/,20x,'is updated using that jacobian')
 9120 FORMAT (13x,' in case of calculating non - spherical wavefcts the parameter lmaxd has to be set equal lmax ')
 9130 FORMAT (/)
 9140 FORMAT (20x,'full potential calculation - cut off of non spherical potential',/,' >',/)
 9150 FORMAT (31x,'representive atom no.',i3,' irns :',i5,' irnsd :',i5)
 9160 FORMAT (21x,a43,/,21x,' using',i3,'-th. born approximation ')
 9170 FORMAT (21x,a43)
 9180 FORMAT (2i5)
 9190 FORMAT (3f12.7,/,4i4)
 9200 FORMAT (3i4,1f12.7)
 9210 FORMAT (' lmax'/,i4)
 9220 FORMAT ('          E1          E2          TK'/,3f12.6)
 9230 FORMAT ('   NPOL  NPNT1  NPNT2  NPNT3'/,4i7)
 9250 FORMAT ('  IFILE    IPE ISHIFT'/,3i7)
 9260 FORMAT (' KSHAPE    IRM    ICST'/,3i7)
 9270 FORMAT ('   KCOR  KVREL    KWS   KHFELD    KXC'/,5i7)
 9280 FORMAT (' external magnetic hfield     :',f15.4/,' VCONST                       :',f15.6)
 9290 FORMAT ('   IMIX   '/,i7)
 9300 FORMAT (' ITDBRY'/,i7)
 9310 FORMAT ('      STRMIX        FCM       QBOUND'/,3f12.6)
 9320 FORMAT ('      BRYMIX'/,f12.6)
 9330 FORMAT ('    KTE   KPRE   KVMAD'/,3i7)
 9301 format(  3(1H-),1H+  ,75(1H-))
 9302 format(3(11(1H-),1H+),43(1H-))
 9303 format(3(6(1H-),1H+) ,58(1H-))
 9304 format(4(6(1H-),1H+) ,51(1H-))
 9305 format(3(6(1H-),1H+),11(1H-),1H+ ,46(1H-))
 9306 format(6(6(1H-),1H+) ,37(1H-))
 9307 format(6(1H-),1H+,72(1H-))
 9308 format(11(1H-),1H+,67(1H-))
 9309 format(5(6(1H-),1H+) ,44(1H-))
#endif
  endsubroutine ! rinputnew99
  
  
  
  subroutine testdimlat(alat, bravais, recbv, rmax, gmax, nmaxd, ishld) ! todo: remove nmaxd and ishld from interface
    use Constants_mod, only: pi
! **********************************************************************
! *  modified version of lattice3d.f                                   *
! *  this one only tests the dimension of arrays!                      *
! *                                            Alexander Thiess 2010   *
! **********************************************************************
! *  generate lattice vectors of direct and reciprocal space from      *
! *  basic translation vectors br                                      *
! *                                                                    *
! *  alat            : lattice constant                                *
! *  br(i,j)         : i=x,y,z j= 1,2,3 bravais vectors                *
! *                    *** in a.u. ****                                *
! *  rmax            : maximum radius in real space        (input)     *
! *  gmax            : maximum radius in reciprocal space  (input)     *
! *  ngmax           : Number of reciprocal lattice vectors            *
! *  gn(3,nmaxd)     : x,y,z   of reciprocal lattice vectors           *
! *  nrmax           : Number of real lattice vectors                  *
! *  rm(3,nmaxd)     : x,y,z  of real space vectors                    *
! *  nshlg           : shells in reciprocal space                      *
! *  nshlr           : shells in real space                            *
! *  nsg,nsr         : integer arrays, number of atoms in each shell   *
! *                                                                    *
! *  Dimension of arrays GN,RM changed from (4,*) to (3,*), the 4th    *
! *  one it is used only locally (GNR/RMR)       v.popescu May 2004    *
! *                                                                    *
! **********************************************************************
    integer, intent(in) :: nmaxd, ishld
    double precision, intent(in) :: alat, rmax, gmax
    double precision, intent(in) :: bravais(3,3), recbv(3,3)
    
    integer :: ngmax, nrmax, nshlr, nshlg, i,k,n,n1,ng,nr,nsh,nshl,numg,numgh,numr,numrh, i1, i2, i3, ist
    double precision :: absg2(3), absr2(3), bg(3,3), br(3,3)
    double precision :: absgm, absrm, ar2, da, db, rv(3), av(3), vmin, vmin2, rmax2, gv(3), gmax2, ag2
    double precision, allocatable :: cv(:,:), cd(:) ! cv(1:3,nmaxd), cd(nmaxd)
    
    br = bravais*alat ! --> basic trans. vectors and basis vectors
    bg = recbv*(2.d0*pi/alat) ! --> generate primitive vectors bg of reciprocal space

    ! --> estimate no. of lattice vectors
    do i = 1, 3
      absr2(i) = br(1,i)**2 + br(2,i)**2 + br(3,i)**2
      absg2(i) = bg(1,i)**2 + bg(2,i)**2 + bg(3,i)**2
    enddo ! i

    absrm = 2.0d0*pi/sqrt(maxval(absr2(1:3)))
    absgm = 2.0d0*pi/sqrt(maxval(absg2(1:3)))
    numr = 2*ceiling(rmax/absgm) + 1
    numg = 2*ceiling(gmax/absrm) + 1
    numrh = numr/2 + 1
    numgh = numg/2 + 1

    rmax2 = rmax**2
    gmax2 = gmax**2
    
!     allocate(cv(0:3,numr**3), stat=ist)
    allocate(cd(numr**3), stat=ist)
    if (ist /= 0) stop 'testdimlat out of memory (real space)'
    
!   generate lattice vectors of real space
    nr = 0
    do i1 = 1, numr
      av(1) = dble(i1 - numrh)
      do i2 = 1, numr
        av(2) = dble(i2 - numrh)
        do i3 = 1, numr
          av(3) = dble(i3 - numrh)
          rv(1:3) = av(1)*br(1:3,1) + av(2)*br(1:3,2) + av(3)*br(1:3,3)
          ar2 = rv(1)**2 + rv(2)**2 + rv(3)**2
          if (ar2 <= rmax2) then
            nr = nr + 1
!             cv(0,nr) = sqrt(ar2) ! also store the radius
!             cv(1:3,nr) = rv(1:3)
            cd(nr) = ar2
          endif ! ar <= rmax
        enddo ! i3
      enddo ! i2
    enddo ! i1
    nrmax = nr

    ! --> sort vectors in order of increasing absolute value
    da = 1.d-6
    nsh = 0
    nshl = -1
    do k = 1, nr 
!     vmin = rmax + 1.d0
      vmin2 = rmax**2 + 1.d0
      
      ! find vmin2 = minval(cd(1:nr))
      ! find n1    = minloc(cd(1:nr))
      do n = 1, nr
!       if (cv(0,n) - vmin < 0.d0) then
        if (cd(n) < vmin2) then
!         vmin = cv(0,n)
          vmin2 = cd(n)
          n1 = n
        endif
      enddo ! n
      
      nshl = nshl + 1
      db = sqrt(vmin2)

      if (db > da + 1.d-6) then
        nsh = nsh + 1
        nshl = 0
        da = db
      endif
      
!       cv(0,n1) = rmax + 1.d0
      cd(n1) = rmax**2 + 1.d0
    enddo ! k
    
    nsh = nsh + 1
    nshl = nshl + 1

!   nsr(nsh) = nshl
    nshlr = nsh
    if (nshlr <= 1) die_here("cut-off radius rmax too small")

    deallocate(cd, stat=ist)    
    
    allocate(cv(0:3,numg**3), stat=ist)
    if (ist /= 0) stop 'testdimlat out of memory (reciprocal space)'
    
!   generate lattice vectors of real space
    ng = 0
    do i1 = 1, numg
      av(1) = dble(i1 - numgh)
      do i2 = 1, numg
        av(2) = dble(i2 - numgh)
        do i3 = 1, numg
          av(3) = dble(i3 - numgh)
          gv(1:3) = av(1)*bg(1:3,1) + av(2)*bg(1:3,2) + av(3)*bg(1:3,3)
          ag2 = gv(1)**2 + gv(2)**2 + gv(3)**2
          if (ag2 <= gmax2) then
            ng = ng + 1
            cv(0,ng) = sqrt(ag2)
            cv(1:3,ng) = gv(1:3)
          endif
        enddo ! i3
      enddo ! i2
    enddo ! i1
    ngmax = ng

    ! --> sort vectors in order of increasing abs. value
    da = 1.d-6
    nsh = 0
    nshl = -1
    do k = 1, ng
      vmin = gmax + 1.d0
      do n = 1, ng
        if (cv(0,n) < vmin) then
          vmin = cv(0,n)
          n1 = n
        endif
      enddo ! n

      nshl = nshl + 1
!     gn(1:3,k) = cv(1:3,n1)
!     gnr(k) = cv(0,n1)
      db = vmin
      if (db > da + 1.d-7) then
        nsh = nsh + 1

!       nsg(nsh) = nshl
        nshl = 0
        da = db
      endif

      cv(0,n1) = gmax + 1.d0
    enddo ! k

    nsh = nsh + 1
    nshl = nshl + 1

!   nsg(nsh) = nshl
    nshlg = nsh
    if (nshlg <= 1) die_here("cut-off radius gmax too small")

    deallocate(cv, stat=ist)    
    
    write(6,'(79(1h=),/)')
    write(6, fmt="(10x,'R max =',F9.5,' (a.u.)',/,10X,'G max =',f9.5,' (1/a.u.)',/)") rmax, gmax
    write(6,'(79(1h=),/,15x,a)') 'checking lattice for Ewald-summ ........... OK'
    write(6,'(79(1h=),/)')
      
  endsubroutine ! testdimlat
  
  
  
  subroutine scalevec(rbasis, naez, bravais, lcartesian)
    integer, intent(in) :: naez
    double precision, intent(inout) :: rbasis(3,*)
    double precision, intent(in) :: bravais(3,3)
    logical, intent(in) :: lcartesian
    
    integer :: i
    double precision :: rb(3)

    write(6,'(79(1h=))')
    write(6,'(23x,a)') 'SCALEVEC: scale site coordinates'
    write(6,'(23x,a)') '          bring all to CARTESIAN system'
    write(6,'(79(1h=))')
    write(6,*)

! ---> normalization of atomic positions in the unit cell
!
!      if lcartesian is true cartesian coordinates are used
!      else the basis atoms are in units of the lattice vectors
!
    write(6,'(5x,a,$)') 'Position of atoms in the unit cell READ IN as:'
    if (lcartesian) then
      write(6,'(a)') ' CARTESIAN coordinates'
      write(6,'(42x,a,/)') '---> No transformation required'
      write(6, fmt="(12x,51(1h-),/,16x,'    Positions of (ALL) generated sites')")
      write(6, fmt="(16x,'   in CARTESIAN coordinates (ALAT units)')")
      write(6, fmt="(12x,51(1h-),/,15x,a,/,12x,51(1h-))") 'IQ       x           y           z       IT'
      do i = 1, naez
        write(6, fmt="(13x,i5,3f12.6,10i3)") i, rbasis(1:3,i), i
      enddo ! i
      write(6,'(12x,51(1h-),/)')
    else ! rescale lattice
      write(6,'(a)') ' LATTICE VECTORS units'
      write(6,*)
      write(6,'(12x,49(1h-))') 
      write(6,'(13x,a)') 'Input positions transformed to CARTESIAN system'
      write(6,'(12x,49(1h-),/,13x,a,/,12x,49(1h-))') 'IQ        x             y             z        IT'
      do i = 1, naez
        rb(1:3) = rbasis(1:3,i)
        rbasis(1:3,i) = rb(1)*bravais(1:3,1) + rb(2)*bravais(1:3,2) + rb(3)*bravais(1:3,3)
        write(6, fmt="(12x,i3,3f14.8,10i3)") i, rbasis(1:3,i), i
      enddo ! i
      write(6,'(12x,49(1h-),/)')
    endif
    ! todo: last section only differ in terms of the output format --> make same
    
    !**********************************************************************
    ! from now on after < scalevec > rbasis are the basis vectors
    ! in units of au/alat in (xyz) reference
    !**********************************************************************
  endsubroutine scalevec
  
  
endmodule ! kkr0_mod

