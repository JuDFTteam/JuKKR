module PotentialConverter_mod
!-------------------------------------------------------------------------------
!> Summary: Converts unformatted potential file to the KKRhost formatted potential file format
!> Author: Elias Rabel, Paul F Baumeister
!> Category: KKRnano, input-output, potential
!> 
!> NOTE: VBC entries are just dummy values!!! RMT, RMTNEW also?
!-------------------------------------------------------------------------------
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  implicit none
  private
  public :: kkrvform

  contains
  
  subroutine kkrvform()
    use DimParams_mod, only: DimParams, load, destroy
    use BasisAtom_mod, only: BasisAtom, load, associateBasisAtomMesh, destroy
    use RadialMeshData_mod, only: RadialMeshData, load, represent, destroy
    use PotentialData_mod, only: represent

    type(DimParams)      :: dims
    type(RadialMeshData) :: mesh
    type(BasisAtom)      :: atomdata
    integer :: iatom
    double precision :: alat = 0.0d0
    double precision :: efermi = 9999.0d0
    double precision :: vbc(2) = 0.0d0
    integer :: kxc = 2
    character(len=:), allocatable :: str
    
    efermi = getFermiEnergy()
    call getStuffFromInputCard(alat, kxc)

    call load(dims, 'bin.dims') ! read dim. parameters from 'bin.dims'

    do iatom = 1, dims%naez

      write(*,*) "Writing potential ", iatom

      call load(atomdata, "bin.atoms", "bin.vpotnew", iatom)
      CHECKASSERT( atomdata%atom_index == iatom )

      call load(mesh, "bin.meshes", iatom)

      call associateBasisAtomMesh(atomdata, mesh)

      ! show data on stdout
      call represent(mesh, str)
      write(*, '(A)') str
      call represent(atomdata%potential, str)
      write(*, '(A)') str

      call writeFormattedPotential(Efermi, ALAT, VBC, KXC, atomdata)

      call destroy(atomdata)
      call destroy(mesh)
    enddo ! iatom

    call destroy(dims)

  endsubroutine ! kkrvform


  !------------------------------------------------------------------------------
  !> Wraps writeFormattedPotentialImpl
  subroutine writeFormattedPotential(Efermi, ALAT, VBC, KXC, atomdata)
    use BasisAtom_mod, only: BasisAtom
    
    double precision, intent(in) :: Efermi
    double precision, intent(in) :: ALAT
    double precision, intent(in) :: VBC(2)
    integer, intent(in) :: KXC
    type(BasisAtom), intent(in) :: atomdata

    integer :: nspind, irnsd

    nspind = atomdata%nspin

#define mesh atomdata%mesh_ptr
    CHECKASSERT( associated(atomdata%mesh_ptr) )

    irnsd = atomdata%potential%irmd - atomdata%potential%irmind

    CHECKASSERT( atomdata%potential%irmd == mesh%irmd )

    call writeFormattedPotentialImpl(Efermi,VBC,NSPIND, &
              KXC,atomdata%potential%LPOT,mesh%A,mesh%B,mesh%IRC, &
              atomdata%potential%VINS,atomdata%potential%VISP,mesh%DRDI,mesh%IRNS,mesh%R,mesh%RWS,mesh%RMT,ALAT, &
              atomdata%core%ECORE,atomdata%core%LCORE(:,1:NSPIND),atomdata%core%NCORE(1:NSPIND),atomdata%Z_nuclear,atomdata%core%ITITLE(:,1:NSPIND), &
              atomdata%atom_index, mesh%irmd, irnsd)
#undef mesh
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> get Fermi energy from file 'energy_mesh'
  double precision function getFermiEnergy() result(efermi)
    double precision :: e1, e2, tk
    integer :: ielast, npnt1, npnt2, npnt3, npol
    double complex, allocatable :: wez(:), ez(:)

    open(67, file='bin.energy_mesh', form='unformatted', action='read', status='old')
    read(67) ielast
    close(67)
    allocate(wez(ielast), ez(ielast))
    open(67, file='bin.energy_mesh', form='unformatted', action='read', status='old')
    read(67) ielast,ez,wez,e1,e2
    read(67) npol,tk,npnt1,npnt2,npnt3
    read(67) efermi
    close(67)
  endfunction ! getFermiEnergy

  !---------------------------------------------------------------------------
  subroutine getStuffFromInputCard(alat, kxc)
    double precision, intent(out) :: alat
    integer, intent(out) :: kxc

!   character(len=80) :: uio
!   integer :: ier

    !CALL IoInput('KEXCOR    ',UIO,1,7,IER)
    !READ (UNIT=UIO,FMT=*) kxc
    !CALL IoInput('ALATBASIS ',UIO,1,7,IER)
    !READ (UNIT=UIO,FMT=*) ALAT

    write(*,*) "WARNING: COULD NOT GET ALAT and KXC - TODO"
    warn(6, "ToBeImplemented: could not get ALAT and KXC")
    alat = 0.d0
    kxc = 2

  endsubroutine ! getStuffFromInputCard
  
  
  !===========================================================================================
  !> Output formatted potential files for each atom.
  subroutine writeFormattedPotentialImpl(fermiEnergy, vbc, nspin, kxc, lpot, a, b, irc, vins, visp, drdi, irns, r, rws, rmt, alat, ecore, lcore, ncore, zat, ititle, rank, irmd, irnsd)
    integer, intent(in) :: irmd, irnsd, rank, nspin, kxc, lpot, irns, irc
    integer, intent(in) :: lcore(20,nspin), ncore(nspin), ititle(20,nspin)
    double precision, intent(in) :: fermiEnergy, alat, vbc(2), a, b, rws, rmt, zat
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(lpot+1)**2,2), visp(irmd,2), drdi(irmd), ecore(20,2), r(irmd)

    character(len=16) :: fname
        
    write(unit=fname, fmt='(a,i5.5)') "vpot.",rank ! generate file name

    open(11, file=fname, form='formatted', action='write')

    call rites(11, nspin, zat, alat, rmt, rmt, rws, ititle, r, drdi, visp, a, b, kxc, irns, lpot, vins, irc, fermiEnergy, vbc, ecore, lcore, ncore, irmd, irnsd)

    close(11)

  endsubroutine ! write
  
  subroutine rites(ifile, nspin, z, alat, rmt, rmtnew, rws, ititle, r, drdi, visp, a, b, kxc, irns, lpot, vins, irc, efermi, vbc, ecore, lcore, ncore, irmd, irnsd)
! ************************************************************************
!      formatted output of potential.
!      this subroutine stores in 'ifile' the necessary results
!      (potentials e.t.c.) to start self-consistency iterations
!
!       if the sm of absolute values of an lm component of vins (non
!       spher. potential) is less than qbound this
!       component will not be stored .
!
!-----------------------------------------------------------------------
    integer, intent(in) :: ifile !< file unit
    integer, intent(in) :: irmd
    integer, intent(in) :: irnsd
    integer, intent(in) :: kxc, lpot, nspin
    integer, intent(in) :: irc
    integer, intent(in) :: irns
    integer, intent(in) :: ititle(20,nspin)
    integer, intent(in) :: lcore(20,nspin)
    integer, intent(in) :: ncore(nspin)
    double precision, intent(in) :: alat
    double precision, intent(in) :: a, b
    double precision, intent(in) :: drdi(irmd), r(irmd)
    double precision, intent(in) :: ecore(20,2)
    double precision, intent(in) :: efermi
    double precision, intent(in) :: rmt, rmtnew, rws
    double precision, intent(in) :: vbc(2)
    double precision, intent(in) :: vins(irmd-irnsd:irmd,(lpot+1)**2,2)
    double precision, intent(in) :: visp(irmd,2)
    double precision, intent(in) :: z

    double precision :: rv, sm
    integer :: ic, ir, irmin, is, lm, lmnr, lmpot, nr
    integer, parameter :: isave = 1, inew = 1
    double precision, parameter :: qbound = 1.d-10
    character(len=24), parameter :: txc(0:3) = [' Morruzi,Janak,Williams', ' von Barth,Hedin       ', ' Vosko,Wilk,Nusair     ', ' GGA PW91              ']
    character(len=*), parameter :: F9000 = "(7a4,6x,'  exc:',a24,3x,a10)", &
      F9010 = "(3f12.8)", F9020 = "(f10.5,/,f10.5,2f15.10)", F9030 = "(i3,/,2d15.8,/,2i2)", &
      F9040 = "(i5,1p,d20.11)", F9060 = "(10i5)", F9070 = "(1p,4d20.13)"
    
    lmpot = (lpot+1)**2

    do is = 1, nspin
      nr = irc
      irmin = nr - irns
      
      write(ifile, fmt=F9000) ititle(1:7,is), txc(kxc)
      write(ifile, fmt=F9010) rmt, alat, rmtnew
      write(ifile, fmt=F9020) z, rws, efermi, vbc(is)
      write(ifile, fmt=F9030) nr, a, b, ncore(is), inew

      if (ncore(is) >= 1) then
        write(ifile, fmt=F9040) (lcore(ic,is), ecore(ic,is), ic=1,ncore(is))
      endif

!--->   store the full potential , but the non spherical contribution
!       only from irns1 up to irws1 ;
!       remember that the lm = 1 contribution is multiplied
!       by a factor 1/sqrt(4 pi)
!
      write(ifile, fmt=F9060) nr, irns, lmpot, isave
      write(ifile, fmt=F9070) visp(1:nr,is)
      if (lpot > 0) then
        lmnr = 1
        do lm = 2, lmpot
          sm = 0.d0
          do ir = irmin, nr
            rv = vins(ir,lm,is)*r(ir)
            sm = sm + rv*rv*drdi(ir)
          enddo ! ir

          if (sqrt(sm) > qbound) then
            lmnr = lmnr + 1
            write(ifile, fmt=F9060) lm
            write(ifile, fmt=F9070) vins(irmin:nr,lm,is)
          endif
        enddo ! lm

!--->     write a one to mark the end
        if (lmnr < lmpot) write(ifile, fmt=F9060) isave
      endif
    enddo ! is

  endsubroutine ! rites
  
endmodule ! PotentialConverter_mod  
