!> Module that contains most of the information and data of a basis atom/site.
!> @author Elias Rabel
!>
!> Example of usage:
!> \verbatim
!>
!>   call createBasisAtom(atomdata, I1, lpot, nspind, irmind, irmd)
!>   call openBasisAtomDAFile(atomdata, 37, "atoms")
!>   call readBasisAtomDA(atomdata, 37, I1)
!>   call closeBasisAtomDAFile(37)
!>
!>   call createShapefunData(cell, irid, (2*LPOT+1)**2, nfund)
!>   call openShapefunDataDAFile(cell, 37 , "cells")
!>   call readShapefunDataDA(cell, 37, atomdata%cell_index)
!>   call closeShapefunDataDAFile(37)
!>
!>   call associateBasisAtomCell(atomdata, cell)
!>
!>   call createRadialMeshData(mesh, irmd, ipand)
!>   call openRadialMeshDataDAFile(mesh, 37 , "meshes")
!>   call readRadialMeshDataDA(mesh, 37, I1)
!>   call closeRadialMeshDataDAFile(37)
!>
!>   call associateBasisAtomMesh(atomdata, mesh)
!>
!> \endverbatim

! move core to potential??

! cell_index, cluster_index --> probably unnecessary, only for additional safety?

! TODO: output of Atom:
! historical reasons: write to two files: vpotnew and atoms

! potential and atom information should be written separately
  ! ==> separation in changing and nonchanging part

! -============= rather not ===============
! Potential not as member of datastructure??? and core as member of potential?
! new potential file format including all core info  (also nonchanging?) ?
! ==========================================
! reference cluster: own file: ref_clusters - or move creation to kkr2

! build in checks for compatibility of shapes with potential etc...

#ifdef TASKLOCAL_FILES
#define FILEWRITE(X,Y) write(X)
#define FILEREAD(X,Y) read(X)
#else
#define FILEWRITE(X,Y) write(X,Y)
#define FILEREAD(X,Y) read(X,Y)
#endif

#define CHECKASSERT(X) if (.not. (X)) then; write(*,*) "ERROR: Check " // #X // " failed. ", __FILE__, __LINE__; STOP; endif

module BasisAtom_mod
#include "macros.h"
  use Exceptions_mod, only: die, launch_warning, operator(-), operator(+)
  use ShapefunData_mod, only: ShapefunData, create, destroy
  use PotentialData_mod, only: PotentialData, create, destroy
  use AtomicCoreData_mod, only: AtomicCoreData, create, destroy
  use RadialMeshData_mod, only: RadialMeshData, create, destroy
  implicit none

  public :: BasisAtom, create, destroy, load
  public :: getMinReclenBasisAtomPotential
  public :: associateBasisAtomCell, associateBasisAtomMesh, resetPotentials, readBasisAtomPotentialHeader
  
  ! direct access file routines
  public :: openBasisAtomDAFile, writeBasisAtomDA, readBasisAtomDA, closeBasisAtomDAFile
  public :: openBasisAtomPotentialDAFile, writeBasisAtomPotentialDA, readBasisAtomPotentialDA, closeBasisAtomPotentialDAFile  
  public :: openBasisAtomPotentialIndexDAFile, writeBasisAtomPotentialIndexDA, readBasisAtomPotentialIndexDA, closeBasisAtomPotentialIndexDAFile

  type BasisAtom
    integer :: atom_index !< position in atominfo/rbasis file
    integer :: cell_index
    double precision :: Z_nuclear
    integer :: nspin
    double precision :: rMTref !< radius of repulsive reference potential
    double precision :: radius_muffin_tin !< user-specified muffin-tin radius
    type(PotentialData) :: potential
    type(AtomicCoreData) :: core
    type(ShapefunData), pointer :: cell_ptr => null()
    type(RadialMeshData), pointer :: mesh_ptr => null()
  endtype

  interface create
    module procedure createBasisAtom
  endinterface
  
  interface load
    module procedure createBasisAtomFromFile
  endinterface
  
  interface destroy
    module procedure destroyBasisAtom
  endinterface
  
  integer, parameter, private :: MAGIC_NUMBER = 385306

  contains

  !----------------------------------------------------------------------------
  subroutine createBasisAtom(atom, atom_index, lpot, nspin, irmind, irmd)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: lpot, nspin, irmind, irmd  !< number of mesh points
    integer, intent(in) :: atom_index

    atom%atom_index = atom_index
    atom%nspin = nspin
    atom%cell_index = -1
    atom%Z_nuclear = 1.d9
    atom%rMTref = 1.d9
    atom%radius_muffin_tin = 1.d9

    call create(atom%potential, lpot, nspin, irmind, irmd) ! createPotentialData

    ! does irmd have to be the same?
    call create(atom%core, irmd) ! createAtomicCoreData

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its cell data.
  !>
  !> This is done that way because different atoms can share the same
  !> cell data.
  subroutine associateBasisAtomCell(atom, cell)
    type(BasisAtom), intent(inout) :: atom
    type(ShapefunData), target, intent(in) :: cell

    atom%cell_ptr => cell

    if (atom%cell_index /= cell%cell_index) die_here("Mismatch in cell indices for atom"+atom%atom_index)

  endsubroutine ! associate

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its radial mesh.
  !>
  !> That way one does not have to keep track of which mesh belongs
  !> to which atom - as soon the relation has been established.
  subroutine associateBasisAtomMesh(atom, mesh)
    type(BasisAtom), intent(inout) :: atom
    type(RadialMeshData), target, intent(in) :: mesh

    atom%mesh_ptr => mesh

    ! check if mesh fits potential dimensions
    CHECKASSERT(mesh%irmd == atom%potential%irmd)
    CHECKASSERT(mesh%irmin == atom%potential%irmind)
    CHECKASSERT(mesh%irns == atom%potential%irnsd)

  endsubroutine ! associate

  !----------------------------------------------------------------------------
  subroutine destroyBasisAtom(atom)
    type(BasisAtom), intent(inout) :: atom

    nullify(atom%cell_ptr)

    call destroy(atom%potential)
    call destroy(atom%core)
  endsubroutine ! destroy


!=========================== I / O ============================================

  !----------------------------------------------------------------------------
  !> creates basis atom AND potential from record 'atom_id' from files
  !> <filename>, <filenamepot>, <filenamepot>.idx
  !> Note: One has to still deal with the mesh!

  !> atom_id = atom_index
  subroutine createBasisAtomFromFile(atom, filename, filenamepot, atom_id)
    type(BasisAtom), intent(inout) :: atom
    character(len=*), intent(in) :: filename, filenamepot
    integer, intent(in) :: atom_id
    
    integer :: irmd, lpot, nspin, irmind, max_reclen
    integer, parameter :: fu=42

    ! index file has extension .idx
    call openBasisAtomPotentialIndexDAFile(atom, fu, filenamepot-".idx", action='read')
    call readBasisAtomPotentialIndexDA(atom, fu, atom_id, lpot, nspin, irmind, irmd, max_reclen)
    call closeBasisAtomPotentialDAFile(fu)

    call createBasisAtom(atom, atom_id, lpot, nspin, irmind, irmd)

    call openBasisAtomDAFile(atom, fu, filename, action='read')
    call readBasisAtomDA(atom, fu, atom_id)
    call closeBasisAtomDAFile(fu)

    call openBasisAtomPotentialDAFile(atom, fu, filenamepot, max_reclen, action='read')
    call readBasisAtomPotentialDA(atom, fu, atom_id)
    call closeBasisAtomPotentialDAFile(fu)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Write basis atom data to direct access file 'fu' at record 'atom_id'.
  subroutine writeBasisAtomDA(atom, fu, atom_id)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fu, atom_id

    write(fu, rec=atom_id) MAGIC_NUMBER, &
                           atom%atom_index, &
                           atom%cell_index, &
                           atom%nspin, &
                           atom%Z_nuclear, &
                           atom%radius_muffin_tin, &
                           atom%rMTref, &
                           atom%core%NCORE, &
                           atom%core%LCORE, &
                           atom%core%ITITLE, &
                           MAGIC_NUMBER
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read basis atom data from direct access file 'fu' at record 'atom_id'.
  subroutine readBasisAtomDA(atom, fu, atom_id)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fu, atom_id

    integer :: magic(2)

    read(fu, rec=atom_id)  magic(1), &
                           atom%atom_index, &
                           atom%cell_index, &
                           atom%nspin, &
                           atom%Z_nuclear, &
                           atom%radius_muffin_tin, &
                           atom%rMTref, &
                           atom%core%NCORE, &
                           atom%core%LCORE, &
                           atom%core%ITITLE, &
                           magic(2)

    if (any(magic /= MAGIC_NUMBER)) die_here("Invalid basis atom data read.")
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens BasisAtom direct access file.
  subroutine openBasisAtomDAFile(atom, fu, filename, action)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fu
    character(len=*), intent(in) :: filename, action

    integer :: rlen

    inquire(iolength=rlen) MAGIC_NUMBER, &
                           atom%atom_index, &
                           atom%cell_index, &
                           atom%nspin, &
                           atom%Z_nuclear, &
                           atom%radius_muffin_tin, &
                           atom%rMTref, &
                           atom%core%NCORE, &
                           atom%core%LCORE, &
                           atom%core%ITITLE, &
                           MAGIC_NUMBER

    open(fu, access='direct', file=filename, recl=rlen, form='unformatted', action=action)
  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes BasisAtom direct access file.
  subroutine closeBasisAtomDAFile(fu)
    integer, intent(in) :: fu

    close(fu)
  endsubroutine ! close


  !----------------------------------------------------------------------------
  !> Write to unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'atom_id'
  !> File has to be opened and closed by user. No checks
  subroutine writeBasisAtomPotentialDA(atom, fu, atom_id)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fu, atom_id

#ifdef TASKLOCAL_FILES
    character(len=16) :: filename
    write(unit=filename, fmt='(A,I7.7)') "bin.pot.",atom_id
    open(fu, file=filename, form='unformatted', action='write')

    call writeBasisAtomPotentialIndexDA(atom, fu, atom_id, 0)
#endif

    FILEWRITE(fu, rec=atom_id) atom%potential%VINS, &
                               atom%potential%VISP, &
                               atom%core%ECORE
#ifdef TASKLOCAL_FILES
    close(fu)
#endif
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Opens unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> WARNING: If number of mesh points is different for each atom, specify
  !> the (otherwise optional) parameter max_reclen !!!
  subroutine openBasisAtomPotentialDAFile(atom, fu, filename, max_reclen, action)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fu
    character(len=*), intent(in) :: filename, action
    integer, intent(in), optional :: max_reclen

#ifndef TASKLOCAL_FILES
    integer :: rlen
    
    if (present(max_reclen)) then
      rlen = max_reclen
    else
      rlen = getMinReclenBasisAtomPotential(atom)
    endif

    open(fu, access='direct', file=filename, recl=rlen, form='unformatted', action=action)
#endif
  endsubroutine ! open

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store potential of this atom
  !>
  !> Note: for a file containing potentials of several atoms, the maximum
  !> of all their record lengths has to be determined
  integer function getMinReclenBasisAtomPotential(atom) result(rlen)
    type(BasisAtom), intent(in) :: atom

    inquire (iolength=rlen) atom%potential%VINS, &
                              atom%potential%VISP, &
                              atom%core%ECORE
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Reads from unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'atom_id'
  !> File has to be opened and closed by user. No checks
  subroutine readBasisAtomPotentialDA(atom, fu, atom_id)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fu, atom_id

#ifdef TASKLOCAL_FILES
    character(len=16) :: filename
    integer :: lpot, nspin, irmind, irmd, max_reclen

    write(unit=filename, fmt='(A,I7.7)') 'bin.pot.',atom_id
    open(fu, file=filename, form='unformatted', action="read", status="old")

    !  skip header at beginning of file
    call readBasisAtomPotentialHeader(atom, fu, atom_id, lpot, nspin, irmind, irmd, max_reclen)
#endif

    FILEREAD(fu, rec=atom_id)  atom%potential%VINS, &
                               atom%potential%VISP, &
                               atom%core%ECORE
#ifdef TASKLOCAL_FILES
    close(fu)
#endif
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Closes unformatted direct-access (DA) potential-file.
  !>
  !> The format is compatible to the 'vpotnew' file
  subroutine closeBasisAtomPotentialDAFile(fu)
    integer, intent(in) :: fu

#ifndef TASKLOCAL_FILES
    close(fu)
#endif
  endsubroutine ! close

  !===========  Index file ======================================================

  !----------------------------------------------------------------------------
  !> Write potential dimension data to direct access file 'fu' at record 'atom_id'
  subroutine writeBasisAtomPotentialIndexDA(atom, fu, atom_id, max_reclen)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fu, atom_id, max_reclen

    FILEWRITE(fu, rec=atom_id) atom%potential%lpot, &
                               atom%potential%nspin, &
                               atom%potential%irmind, &
                               atom%potential%irmd, &
                               max_reclen, &
                               MAGIC_NUMBER + atom_id

  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read potential dimension data from direct access file 'fu' at record 'atom_id'
  !>
  !> Returns dimensions irmd and ipand
  subroutine readBasisAtomPotentialIndexDA(atom, fu, atom_id, lpot, nspin, irmind, irmd, max_reclen)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fu, atom_id
    integer, intent(out) :: lpot, nspin, irmind, irmd, max_reclen
    
#ifdef TASKLOCAL_FILES
    character(len=16) :: filename

    write(unit=filename, fmt='(A,I7.7)') "bin.pot.",atom_id
    open(fu, file=filename, form='unformatted', action='read', status="old")
#endif

    call readBasisAtomPotentialHeader(atom, fu, atom_id, lpot, nspin, irmind, irmd, max_reclen)

#ifdef TASKLOCAL_FILES
    close(fu)
#endif

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens BasisAtomPotential index file.
  subroutine openBasisAtomPotentialIndexDAFile(atom, fu, filename, action)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fu
    character(len=*), intent(in) :: filename, action

#ifndef TASKLOCAL_FILES
    integer :: rlen, max_reclen
    
    max_reclen = 0
    inquire(iolength=rlen) atom%potential%lpot, &
                           atom%potential%nspin, &
                           atom%potential%irmind, &
                           atom%potential%irmd, &
                           max_reclen, &
                           MAGIC_NUMBER

    open(fu, access='direct', file=filename, recl=rlen, form='unformatted', action=action)
#endif
  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes BasisAtomPotential index file.
  subroutine closeBasisAtomPotentialIndexDAFile(fu)
    integer, intent(in) :: fu

#ifndef TASKLOCAL_FILES
    close(fu)
#endif
  endsubroutine ! close

  !----------------------------------------------------------------------------
  !> Copy output potential to input potential for new iteration.
  subroutine resetPotentials(atom)
    type(BasisAtom), intent(inout) :: atom

#define  mesh atom%mesh_ptr
    call resetPotentialsImpl(mesh%IRC, mesh%IRMD, mesh%IRMIN, &
                         atom%potential%IRMIND, atom%potential%LMPOT, &
                         atom%potential%NSPIN, atom%potential%VINS, &
                         atom%potential%VISP, atom%potential%VONS)
#undef mesh
  endsubroutine ! reset

!================ Private helper functions ====================================

  subroutine readBasisAtomPotentialHeader(atom, fu, atom_id, lpot, nspin, irmind, irmd, max_reclen)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fu, atom_id
    integer, intent(out) :: lpot, nspin, irmind, irmd, max_reclen

    integer :: magic, checkmagic

    checkmagic = MAGIC_NUMBER + atom_id

    FILEREAD(fu, rec=atom_id)  lpot, nspin, irmind, irmd, max_reclen, magic

    if (magic /= checkmagic) die_here("Invalid mesh index data read.")

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Copy output potential to input potential for new iteration.
  subroutine resetPotentialsImpl(irc1, irmd, irmin1, irmind, lmpotd, nspin, vins, visp, vons)
    integer, intent(in) :: irc1, irmin1, irmd, irmind, lmpotd, nspin
    double precision, intent(out) :: vins(irmind:irmd,lmpotd,2)
    double precision, intent(out) :: visp(irmd,2)
    double precision, intent(inout) :: vons(irmd,lmpotd,2)

    integer :: ispin

    if (irmin1 < irmind) die_here("resetpotentials:"+irmin1+"= irmin1 < irmind ="+irmind) ! debug
    if (irc1 > irmd) die_here("resetpotentials:"+irc1+"= irc1 > irmd ="+irmd) ! debug

    vins(:,:,:) = 0.d0 ! initialise vins
    visp(:,:) = 0.d0 ! initialise visp

    ! copy output potential to input potential for new iteration
    ! note: nonspherical part for indices < irmin1 is thrown away!
    do ispin = 1, nspin
      visp(1:irc1,ispin) = vons(1:irc1,1,ispin)
      if (lmpotd > 1) vins(irmin1:irc1,2:lmpotd,ispin) = vons(irmin1:irc1,2:lmpotd,ispin)
    enddo ! ispin
    
  endsubroutine ! resetPotentialsImpl

endmodule ! BasisAtom_mod
