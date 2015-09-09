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
!>   call createCellData(cell, irid, (2*LPOT+1)**2, nfund)
!>   call openCellDataDAFile(cell, 37 , "cells")
!>   call readCellDataDA(cell, 37, getCellIndex(atomdata))
!>   call closeCellDataDAFile(37)
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
  use CellData_mod, only: CellData
  use PotentialData_mod, only: PotentialData
  use AtomicCoreData_mod, only: AtomicCoreData
  use RadialMeshData_mod, only: RadialMeshData
  implicit none
  
  public :: BasisAtom, create, destroy
  public :: getMinReclenBasisAtomPotential, getCellIndex
  public :: createBasisAtom, destroyBasisAtom, createBasisAtomFromFile ! deprecated
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
    double precision :: RMTref !< radius of repulsive reference potential
    double precision :: radius_muffin_tin !< user-specified muffin-tin radius
    type(PotentialData) :: potential
    type(AtomicCoreData) :: core
    type(CellData), pointer :: cell_ptr => null()
    type(RadialMeshData), pointer :: mesh_ptr => null()
  endtype

  interface create
    module procedure createBasisAtom, createBasisAtomFromFile
  endinterface
  
  interface destroy
    module procedure destroyBasisAtom
  endinterface
  
  integer, parameter, private :: MAGIC_NUMBER = 385306
  
  contains

  !----------------------------------------------------------------------------
  subroutine createBasisAtom(atom, atom_index, lpot, nspin, irmind, irmd)
    use PotentialData_mod, only: createPotentialData
    use AtomicCoreData_mod, only: createAtomicCoreData

    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: lpot, nspin, irmind, irmd  !< number of mesh points
    integer, intent(in) :: atom_index

    atom%atom_index = atom_index
    atom%nspin = nspin
    atom%cell_index = -1
    atom%Z_nuclear = 1.d9
    atom%RMTref = 1.d9
    atom%radius_muffin_tin = 1.d9

    call createPotentialData(atom%potential, lpot, nspin, irmind, irmd)

    ! does irmd have to be the same?
    call createAtomicCoreData(atom%core, irmd)

  endsubroutine ! create

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its cell data.
  !>
  !> This is done that way because different atoms can share the same
  !> cell data.
  subroutine associateBasisAtomCell(atom, cell)
    use CellData_mod, only: CellData
    type(BasisAtom), intent(inout) :: atom
    type(CellData), target, intent(in) :: cell
    
    type(CellData), pointer :: cell_ptr

    cell_ptr => cell
    atom%cell_ptr => cell_ptr

    if (atom%cell_index /= cell_ptr%cell_index) &
      die_here("Mismatch in cell indices for atom"+atom%atom_index)

  endsubroutine ! associate

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its radial mesh.
  !>
  !> That way one does not have to keep track of which mesh belongs
  !> to which atom - as soon the relation has been established.
  subroutine associateBasisAtomMesh(atom, mesh)
    type(BasisAtom), intent(inout) :: atom
    type(RadialMeshData), target, intent(in) :: mesh

    type(RadialMeshData), pointer :: mesh_ptr

    mesh_ptr => mesh
    atom%mesh_ptr => mesh_ptr

    ! check if mesh fits potential dimensions
    CHECKASSERT(mesh%irmd == atom%potential%irmd)
    CHECKASSERT(mesh%irmin == atom%potential%irmind)
    CHECKASSERT(mesh%irns == atom%potential%irnsd)

  endsubroutine ! associate

  !----------------------------------------------------------------------------
  subroutine destroyBasisAtom(atom)
    use PotentialData_mod, only: destroyPotentialData
    use AtomicCoreData_mod, only: destroyAtomicCoreData

    type(BasisAtom), intent(inout) :: atom

    nullify(atom%cell_ptr)

    call destroyPotentialData(atom%potential)
    call destroyAtomicCoreData(atom%core)
  endsubroutine ! destroy

  !----------------------------------------------------------------------------
  !> Returns index of corresponding cell.
  pure integer function getCellIndex(atom)
    type(BasisAtom), intent(in) :: atom

    getCellIndex = atom%cell_index
  endfunction getCellIndex

!=========================== I / O ============================================

  !----------------------------------------------------------------------------
  !> creates basis atom AND potential from record 'recnr' from files
  !> <filename>, <filenamepot>, <filenamepot>.idx
  !> Note: One has to still deal with the mesh!

  !> recnr = atom_index
  subroutine createBasisAtomFromFile(atom, filename, filenamepot, recnr)
    type(BasisAtom), intent(inout) :: atom
    character(len=*), intent(in) :: filename, filenamepot
    integer, intent(in) :: recnr
    
    integer :: irmd, lpot, nspin, irmind, max_reclen
    integer, parameter :: FILEUNIT = 42

    ! index file has extension .idx
    call openBasisAtomPotentialIndexDAFile(atom, FILEUNIT, filenamepot-".idx")
    call readBasisAtomPotentialIndexDA(atom, FILEUNIT, recnr, lpot, nspin, irmind, irmd, max_reclen)
    call closeBasisAtomPotentialDAFile(FILEUNIT)

    call createBasisAtom(atom, recnr, lpot, nspin, irmind, irmd)

    call openBasisAtomDAFile(atom, FILEUNIT, filename)
    call readBasisAtomDA(atom, FILEUNIT, recnr)
    call closeBasisAtomDAFile(FILEUNIT)

    call openBasisAtomPotentialDAFile(atom, FILEUNIT, filenamepot, max_reclen)
    call readBasisAtomPotentialDA(atom, FILEUNIT, recnr)
    call closeBasisAtomPotentialDAFile(FILEUNIT)

  endsubroutine

  !----------------------------------------------------------------------------
  !> Write basis atom data to direct access file 'fileunit' at record 'recnr'.
  subroutine writeBasisAtomDA(atom, fileunit, recnr)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit, recnr

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%radius_muffin_tin, &
                                atom%RMTref, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read basis atom data from direct access file 'fileunit' at record 'recnr'.
  subroutine readBasisAtomDA(atom, fileunit, recnr)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit, recnr

    integer :: magic, magic2

    read  (fileunit, rec=recnr) magic, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%radius_muffin_tin, &
                                atom%RMTref, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                magic2

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) die_here("Invalid basis atom data read.")
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens BasisAtom direct access file.
  subroutine openBasisAtomDAFile(atom, fileunit, filename)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    integer :: reclen

    inquire (iolength=reclen) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%radius_muffin_tin, &
                                atom%RMTref, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes BasisAtom direct access file.
  subroutine closeBasisAtomDAFile(fileunit)
    integer, intent(in) :: fileunit

    close(fileunit)
    
  endsubroutine


  !----------------------------------------------------------------------------
  !> Write to unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine writeBasisAtomPotentialDA(atom, fileunit, recnr)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit, recnr

#ifdef TASKLOCAL_FILES
    character(len=16) :: num
    write(unit=num, fmt='(A,I7.7)') "pot.",recnr
    open(fileunit, file=num, form='unformatted', action='write')

    call writeBasisAtomPotentialIndexDA(atom, FILEUNIT, recnr, 0)
#endif

    FILEWRITE (fileunit, rec=recnr) atom%potential%VINS, &
                                    atom%potential%VISP, &
                                    atom%core%ECORE

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Opens unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> WARNING: If number of mesh points is different for each atom, specify
  !> the (otherwise optional) parameter max_reclen !!!
  subroutine openBasisAtomPotentialDAFile(atom, fileunit, filename, max_reclen)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: max_reclen

#ifndef TASKLOCAL_FILES
    integer :: reclen
    
    if (present(max_reclen)) then
      reclen = max_reclen
    else
      reclen = getMinReclenBasisAtomPotential(atom)
    endif

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
  endsubroutine ! open

  !---------------------------------------------------------------------------
  !> Return MINIMUM record length needed to store potential of this atom
  !>
  !> Note: for a file containing potentials of several atoms, the maximum
  !> of all their record lengths has to be determined
  integer function getMinReclenBasisAtomPotential(atom) result(reclen)
    type(BasisAtom), intent(in) :: atom

    inquire (iolength=reclen) atom%potential%VINS, &
                              atom%potential%VISP, &
                              atom%core%ECORE
  endfunction ! get

  !----------------------------------------------------------------------------
  !> Reads from unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine readBasisAtomPotentialDA(atom, fileunit, recnr)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit, recnr

#ifdef TASKLOCAL_FILES
    character(len=16) :: num
    integer :: lpot, nspin, irmind, irmd, max_reclen

    write(unit=num, fmt='(A,I7.7)') 'pot.',recnr
    open(fileunit, file=num, form='unformatted', action="read", status="old")

    !  skip header at beginning of file
    call readBasisAtomPotentialHeader(atom, fileunit, recnr, lpot, nspin, irmind, irmd, max_reclen)
#endif

    FILEREAD (fileunit, rec=recnr) atom%potential%VINS, &
                                   atom%potential%VISP, &
                                   atom%core%ECORE

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Closes unformatted direct-access (DA) potential-file.
  !>
  !> The format is compatible to the 'vpotnew' file
  subroutine closeBasisAtomPotentialDAFile(fileunit)
    integer, intent(in) :: fileunit

#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! close

  !===========  Index file ======================================================

  !----------------------------------------------------------------------------
  !> Write potential dimension data to direct access file 'fileunit' at record 'recnr'
  subroutine writeBasisAtomPotentialIndexDA(atom, fileunit, recnr, max_reclen)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit, recnr, max_reclen

    FILEWRITE (fileunit, rec=recnr) atom%potential%lpot, &
                                    atom%potential%nspin, &
                                    atom%potential%irmind, &
                                    atom%potential%irmd, &
                                    max_reclen, &
                                    MAGIC_NUMBER + recnr

  endsubroutine ! write

  !----------------------------------------------------------------------------
  !> Read potential dimension data from direct access file 'fileunit' at record 'recnr'
  !>
  !> Returns dimensions irmd and ipand
  subroutine readBasisAtomPotentialIndexDA(atom, fileunit, recnr, lpot, nspin, irmind, irmd, max_reclen)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit, recnr
    integer, intent(out) :: lpot, nspin, irmind, irmd, max_reclen
    
#ifdef TASKLOCAL_FILES
    character(len=16) :: num

    write(unit=num, fmt='(A,I7.7)') "pot.",recnr
    open(fileunit, file=num, form='unformatted', action='read', status="old")
#endif

    call readBasisAtomPotentialHeader(atom, fileunit, recnr, lpot, nspin, irmind, irmd, max_reclen)

#ifdef TASKLOCAL_FILES
    close(fileunit)
#endif

  endsubroutine ! read

  !----------------------------------------------------------------------------
  !> Opens BasisAtomPotential index file.
  subroutine openBasisAtomPotentialIndexDAFile(atom, fileunit, filename)
    type(BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

#ifndef TASKLOCAL_FILES
    integer :: reclen, max_reclen
    
    max_reclen = 0
    inquire (iolength=reclen) atom%potential%lpot, &
                              atom%potential%nspin, &
                              atom%potential%irmind, &
                              atom%potential%irmd, &
                              max_reclen, &
                              MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')
#endif
  endsubroutine ! open

  !----------------------------------------------------------------------------
  !> Closes BasisAtomPotential index file.
  subroutine closeBasisAtomPotentialIndexDAFile(fileunit)
    integer, intent(in) :: fileunit

#ifndef TASKLOCAL_FILES
    close(fileunit)
#endif
  endsubroutine ! close

  !----------------------------------------------------------------------------
  !> Copy output potential to input potential for new iteration.
  subroutine resetPotentials(atom)
    type(BasisAtom), intent(inout) :: atom

    type(RadialMeshData), pointer :: mesh

    mesh => atom%mesh_ptr
    call resetPotentialsImpl(mesh%IRC, mesh%IRMD, mesh%IRMIN, &
                         atom%potential%IRMIND, atom%potential%LMPOT, &
                         atom%potential%NSPIN, atom%potential%VINS, &
                         atom%potential%VISP, atom%potential%VONS)

  endsubroutine ! reset

!================ Private helper functions ====================================

  subroutine readBasisAtomPotentialHeader(atom, fileunit, recnr, lpot, nspin, irmind, irmd, max_reclen)
    type(BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit, recnr
    integer, intent(out) :: lpot, nspin, irmind, irmd, max_reclen

    integer :: magic, checkmagic

    checkmagic = MAGIC_NUMBER + recnr

    FILEREAD  (fileunit, rec=recnr) lpot, nspin, irmind, irmd, max_reclen, magic

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

      call dcopy(irc1,vons(1,1,ispin),1,visp(1,ispin),1)

      if (lmpotd > 1) vins(irmin1:irc1,2:lmpotd,ispin) = vons(irmin1:irc1,2:lmpotd,ispin)
    enddo
  endsubroutine ! resetPotentialsImpl

endmodule BasisAtom_mod
