!--------------------------------------------------------------------------------
! Copyright (c) 2018 Forschungszentrum Juelich GmbH, Juelich, Germany
! This file is part of KKRnano and available as free software under the conditions
! of the MIT license as expressed in the LICENSE file in more detail.
!--------------------------------------------------------------------------------
module KKRzero_mod
!-------------------------------------------------------------------------------
!> Summary: Preparation of binary data files
!> Author: Marcel Bornemann, Paul F Baumeister, Elias Rabel, Alexander Thiess, Rudolf Zeller
!> Category: KKRnano, initialization, potential, input-output, k-points, geometry
!>
!> ToDo: if energy_mesh.0 is not generated, check if settings in
!>       input.conf match those in energy_mesh.0
!> PROBLEM: Fermi-Energy specified in 'potential' file - E-mesh depends on
!>          E-Fermi -> k-mesh depends on E-mesh => k-mesh depends on EFERMI
!> NOTE: k-mesh ALWAYS DEPENDS ON THE FERMI ENRGY FROM 'potential' NOT ON THE ACTUAL ONE!!!!!! - BUG?
!> program MAIN0
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: main0
  
  contains
  
  subroutine main0(checkmode, voronano)

! Explanation of most variables follows below

!     ALAT                     : lattice constant (in a.u.)
!     ABASIS,BBASIS,CBASIS,    : scaling factors for rbasis
!     E1,E2,                   : energies needed in EMESHT
!     HFIELD                   : external magnetic field, for initial potential shift in spin polarised case
!     TK,                      : temperature
!     VCONST,                  : potential shift
!     BRAVAIS(3,3),            : bravais lattice vectors
!     RECBV(3,3),              : reciprocal basis vectors
!     RMTREF(NREFD),           : muffin-tin radius of reference system
!     RBASIS(3,NAEZD),         : position of atoms in the unit cell in units of bravais vectors
!     RCLS(3,NACLSD,NCLSD),    : real space position of atom in cluster
!     RR(3,0:NRD)              : set of real space vectors (in a.u.)
!     VBC(2),                  : potential constants
!     WG(LASSLD),              : integr. weights for Legendre polynomials
!     YRG(LASSLD,0:LASSLD)     : spherical harmonics (GAUNT2)
!     ZAT(NAEZD)               : nuclear charge
!     INTERVX,INTERVY,INTERVZ, : number of intervals in x,y,z-direction for k-net in IB of the BZ
!     ICST,                    : number of Born approximation
!     IEND,                    : number of nonzero gaunt coeffizients
!     IFILE,                   : unit specifier for potential card
!     IPE,IPF,IPFE,            : not real used, IPFE should be 0
!     KHFELD,                  : 0,1: no / yes external magnetic field
!     KVREL,                   : 0,1 : non / scalar relat. calculation
!     LMAX,                    : maximum l component in wave function expansion
!     LPOT,                    : maximum l component in potential expansion
!     NAEZ,                    : number of atoms in unit cell
!     NCLS,                    : number of reference clusters
!     NPNT1,NPNT2,NPNT3,       : number of E points (EMESHT)
!     NPOL,                    : number of Matsubara Pols (EMESHT)
!     NR,                      : number of real space vectors rr
!     NREF,                    : number of diff. ref. potentials
!     NSPIN,                   : counter for spin directions
!     IGUESS                   : 0,1 : no / yes (sc) initial guess, set IGUESSD to 1 in inc.p if needed
!     BCP                      : 0,1 : no / yes bc-preconditioning, set BCPD to 1 in inc.p if needed
!     QBOUND                   : exit condition for self-consistency iterations
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
!     ICST                     : the regular non spherical wavefunctions, the alpha matrix and the t-matrix in the ICST-th. born approx
!     IMT(NAEZD),              : r point at MT radius
!     IPAN(NAEZD),             : number of panels in non-MT-region
!     IRC(NAEZD),              : r point for potential cutting
!     IRCUT(0:IPAND,NAEZD),    : r points of panel borders
!     IRMIN(NAEZD),            : max r for spherical treatment
!     IRNS(NAEZD)              : number r points for non spher. treatm.
!     IRWS(NAEZD),             : r point at WS radius
!     LMSP(NAEZD,LMXSPD)       : 0,1 : non/-vanishing lm=(l,m) component of non-spherical potential
!     LLMSP(NAEZD,NFUND)       : lm=(l,m) of 'nfund'th nonvanishing component of non-spherical pot.
!     NTCELL(NAEZD),           : index for WS cell

!     A(NAEZD),B(NAEZD)        : contants for exponential r mesh
!     R(IRMD,NAEZD)            : radial mesh (in units a Bohr)
!     DRDI(IRMD,NAEZD)         : derivative dr/di
!     THETAS(IRID,NFUND,NCELLD): shape function: THETA = (0 outer space and 1 inside WS cell) in spherical harmonics expansion
!     LCORE(20,NPOTD)          : angular momentum of core states
!     NCORE(NPOTD)             : number of core states
! ----------------------------------------------------------------------

    use InputParams_mod, only: InputParams, getValues, store
    use DimParams_mod, only: DimParams, parse, store
    use Main2Arrays_mod, only: Main2Arrays, create, store
    use BrillouinZone_mod, only: bzkint0, readKpointsFile
    use BrillouinZoneMesh_mod, only: BrillouinZoneMesh, create, load, store, destroy
    use Warnings_mod, only: show_warning_lines
    use Lattice_mod, only: lattix99
    use Constants_mod, only: pi
    use MadelungCalculator_mod, only: testdimlat
    use EnergyMesh_mod, only: getEnergyMeshSize, EnergyMesh, create, init, update, store, destroy
    use Startb1_mod, only: startb1_wrapper_new
    
    integer, intent(in) :: checkmode ! 0: usual kkr0, >0: checks only, no writing of any files
    integer, intent(in) :: voronano  ! 0: usual kkr0,  1: returns before reading potential and shapefunctions
    
    integer, parameter  :: krel=0

    double precision    :: efermi, recbv(3,3), volume0
    integer             :: ist, ierror
    logical             :: startpot_exists
    
    type(DimParams)     :: dims
    type(InputParams)   :: params
    type(Main2Arrays)   :: arrays
    type(EnergyMesh)    :: emesh
    type(BrillouinZoneMesh) :: kmeshes(8)

    call parse(dims, "global.conf", altfile="input.conf")

#ifndef USE_OLD_MESH
    ! This is related to the issue kkr/jukkr#53 in the IFF GitLab
    if (dims%korbit /= 0) then
      warn(6,"For NOCO calculations (korbit = 1) usage of makefile option 'TYPE=voronoi_mesh' is advised") 
    endif
#endif

    ist = getValues("input.conf", params)
    if (ist /= 0) die_here('failed to read "input.conf"!')

    ! global.conf <-> input.conf consistency checks
    if (dims%KPOIBZ < params%bzdivide(1)*params%bzdivide(2)*params%bzdivide(3)) then
      dims%KPOIBZ = params%bzdivide(1)*params%bzdivide(2)*params%bzdivide(3)
      warn(6,'Kpoint allocation insufficient. KPOIBZ is increased to ' + params%bzdivide(1)*params%bzdivide(2)*params%bzdivide(3)) 
    endif

    if (params%soc) then
      write(*,'(A,F4.2)') 'Spin-orbit coupling is scaled with socscale=', params%socscale
      write(*,*) 'test=', params%socscale
    endif

    dims%iemxd = getEnergyMeshSize(params%npol, [params%npnt1, params%npnt2, params%npnt3], params%npntsemi)
    call create(emesh, dims%iemxd)

    call create(arrays, dims%lmmaxd, dims%lmmaxd_noco, dims%naez, dims%kpoibz, dims%maxmshd)

    call rinputnew99(arrays%rbasis, arrays%zat, dims%naez) ! will modify naez if naez == 0 (auto mode)
    arrays%naez = dims%naez ! store corrected number of all atoms also in arrays%
    
!   in case of a LDA+U calculation - read file 'ldauinfo' and write 'wldau.unf', if it does not exist already
    if (params%LDAU) call ldauinfo_read(dims%lmaxd, dims%nspind, arrays%zat, dims%naez)

!   NOCO calculations require NSPIND==2
    if (dims%korbit == 1) then
       if(dims%nspind .NE. 2) die_here('NSPIND=2 in global.conf is mandatory for SOC calculations')
    else
       if(dims%korbit .NE. 0) die_here('When not using NOCO: KORBIT in global.conf should be zero')    
    endif

    if (voronano == 1) then
      ! Conversion of rmax and gmax to atomic units
      params%rmax = params%rmax*params%alat
      params%gmax = params%gmax/params%alat
      call store(dims, 'bin.dims')
      call store(params, 'bin.input')
      arrays%bravais(:,1) = params%bravais_a(1:3)
      arrays%bravais(:,2) = params%bravais_b(1:3)
      arrays%bravais(:,3) = params%bravais_c(1:3)
      call store(arrays, 'bin.arrays')
      write(*,*) 'voronano == 1: Starting potential and shapefunctions are not read in by kkr0'
      return
    endif

!===================================================================

    ! read starting potential and shapefunctions
    startpot_exists = .false.; inquire(file='potential', exist=startpot_exists)
    ! if energy_mesh.0 file is missing, also regenerate start files
    if (startpot_exists) then

!      write(*,*) 'Entering startb1_wrapper_new'
      call startb1_wrapper_new(params%alat, dims%nspind, efermi, arrays%zat, dims%naez, nowrite=(checkmode /= 0))
!      write(*,*) 'Leaving startb1_wrapper_new'

    else
      ! no formatted potential provided
      write(*,*) "WARNING: file 'potential' not found... skipping start potential generation."
      warn(6, "file 'potential' not found... skipping start potential generation.")
      write(*,*) "Trying to read initial, approximate Fermi energy from EFERMI file..."
      open(67, file='EFERMI', form='formatted', action='read', status='old')
      read(67,*) efermi
      close(67)
    endif

! ----------------------------------------------------------------------
! update fermi energy, adjust energy window according to running options

! for non-dos calculation upper energy bound corresponds to fermi energy
! bad: params%emax is changed
    if (params%npol /= 0) params%emax = efermi

    call init(emesh, efermi, params%emin, params%emax, params%tempr, params%npol, [params%npnt1, params%npnt2, params%npnt3], &
                    params%ebotsemi, params%emusemi, params%npntsemi, params%fsemicore)
    call update(emesh)

! ================================================ deal with the lattice

    arrays%bravais(:,1) = params%bravais_a(1:3)
    arrays%bravais(:,2) = params%bravais_b(1:3)
    arrays%bravais(:,3) = params%bravais_c(1:3)

    ! only for informative purposes - prints info about lattice
    call lattix99(params%alat, arrays%bravais, recbv, volume0, .true.)


    call scalevec(arrays%rbasis, arrays%naez, arrays%bravais, params%cartesian)

! ======================================================================
!     setting up kpoints
! ======================================================================

    call bzkint0(arrays%naez, arrays%rbasis, arrays%bravais, recbv, arrays%nsymat, arrays%isymindex, &
                 arrays%dsymll, params%bzdivide, emesh%ielast, emesh%iesemicore, emesh%ez, dims%iemxd, emesh%kmesh, arrays%maxmesh, &
                 dims%lmaxd, dims%lmmaxd_noco, krel+dims%korbit, dims%ekmd, params%fullbz, dims%korbit, nowrite=(checkmode /= 0), kpms=kmeshes) ! after return from bzkint0, ekmd contains the right value
    
!   Conversion of rmax and gmax to atomic units
    params%rmax = params%rmax*params%alat
    params%gmax = params%gmax/params%alat

    call testdimlat(params%alat, arrays%bravais, recbv, params%rmax, params%gmax) ! recommends mimimum values for NMAXD and ISHLD 
    
    if (checkmode == 0) then 
      ! write binary files that are needed in the main program

      ! bzkint0 wrote a file 'kpoints': read this file and use it as k-mesh
      call readKpointsFile(arrays%maxmesh, arrays%nofks, arrays%bzkp, arrays%volcub, arrays%volbz)
      
      call store(dims, 'bin.dims')
      call store(params, 'bin.input')
      call store(arrays, 'bin.arrays')

      call store(emesh, filename='bin.energy_mesh.0')
        
    else  ! checkmode == 0
      write(*,'(A)') "CheckMode: binary files 'bin.dims', 'bin.input' and 'bin.arrays' are not created!" ! do we need a warning here?
    endif ! checkmode == 0

    ist = show_warning_lines(unit=6)
    
    call destroy(emesh)
    
  endsubroutine ! main0

  
!>    Reads atominfo and rbasis files.
  subroutine rinputnew99(rbasis, zat, naez)
    use PositionReader_mod, only: getAtomData
    include 'mpif.h'

    double precision, allocatable, intent(inout) :: zat(:) ! (naez)
    double precision, allocatable, intent(inout) :: rbasis(:,:) ! (3,naez)
    integer, intent(inout) :: naez !< if naez < 1 (auto) it will be modified to the number of atoms on exit
    
    double precision, allocatable  :: pos(:,:)
    character(len=9), parameter :: version = "May  2018"
    integer :: ist, naez_xyz

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

! read data from .xyz-file 'rbasis.xyz'
    ist = getAtomData('rbasis.xyz', naez_xyz, pos, MPI_COMM_WORLD)
! check if number of atoms from 'global.conf' equals number of atoms from .xyz-file 
    if (naez < 1) then ! automatic mode: naez is determined here
      naez = naez_xyz
      write(6, fmt="(9(a,i0))") "found ",naez," atoms in rbasis.xyz"
      ! reshape the arrays
      deallocate(zat, rbasis, stat=ist)
      allocate(zat(naez), rbasis(3,naez), stat=ist) ; assert(ist == 0)
    else
      if (naez_xyz /= naez) die_here('number of atoms in global.conf ('-naez-') differ from that in rbasis.xyz ('-naez_xyz-')!')
    endif
! assign nuclear charge and rbasis
    zat(1:naez)        = pos( 0 ,1:naez)
    rbasis(1:3,1:naez) = pos(1:3,1:naez)
    !ToDo:: check whether read-in was successful    

    write(6, fmt="(3(7(1h-),1h+) ,55(1h-))")
    write(6, fmt="(10(3(1h-),1h+) ,39(1h-))")
     
  endsubroutine ! rinputnew99
  
  subroutine scalevec(rbasis, naez, bravais, lcartesian)
    double precision, intent(inout) :: rbasis(3,*)
    integer, intent(in) :: naez
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
        write(6, fmt="(13x,i5,3f12.6,'  ',i0)") i, rbasis(1:3,i), i
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
        write(6, fmt="(12x,i3,3f14.8,'  ',i0)") i, rbasis(1:3,i), i
      enddo ! i
      write(6,'(12x,49(1h-),/)')
    endif
    ! todo: last section only differ in terms of the output format --> make same
    
    !**********************************************************************
    ! from now on after < scalevec > rbasis are the basis vectors
    ! in units of au/alat in (xyz) reference
    !**********************************************************************
  endsubroutine ! scalevec
  
  
endmodule ! kkr0_mod
