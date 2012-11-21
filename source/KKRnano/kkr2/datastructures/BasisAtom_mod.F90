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


! eliminates: ZAT, VREF, REFPOT?
! TODO: establish BasisAtomGroup???
! references to t-matrix and G_ref ???

! build in checks for compatibility of shapes with potential etc...

module BasisAtom_mod
  use CellData_mod
  use RefClusterData_mod
  use PotentialData_mod
  use AtomicCoreData_mod
  use RadialMeshData_mod
  implicit none

  type BasisAtom

    integer :: atom_index !todo (position in atominfo/rbasis file)
    integer :: cell_index
    ! TODO: do the same with cluster
    integer :: cluster_index
    double precision :: Z_nuclear !todo
    integer :: nspin
    !double precision :: Vref      !todo  ! or move to reference cluster? yes
    !double precision :: RMTref    ! move to ref cluster
    ! double precision, dimension(3) :: position !todo
    ! position, etc...
    ! lmax ???

    type (PotentialData) :: potential
    type (AtomicCoreData) :: core

    type (CellData), pointer :: cell_ptr => null()
    type (RefClusterData), pointer :: cluster_ptr => null()
    type (RadialMeshData), pointer :: mesh_ptr => null()

  end type

CONTAINS

  !----------------------------------------------------------------------------
  subroutine createBasisAtom(atom, atom_index, lpot, nspin, irmind, irmd)
    use PotentialData_mod
    use AtomicCoreData_mod
    implicit none

    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: lpot
    integer, intent(in) :: nspin
    integer, intent(in) :: irmind
    integer, intent(in) :: irmd  !< number of mesh points
    integer, intent(in) :: atom_index

    atom%atom_index = atom_index
    atom%nspin = nspin
    atom%cell_index = -1
    atom%cluster_index = -1
    atom%Z_nuclear = 1.d9

    call createPotentialData(atom%potential, lpot, nspin, irmind, irmd)

    ! does irmd have to be the same?
    call createAtomicCoreData(atom%core, irmd)

  end subroutine

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its cell data.
  !>
  !> This is done that way because different atoms can share the same
  !> cell data.
  subroutine associateBasisAtomCell(atom, cell)
    use CellData_mod
    implicit none
    type (BasisAtom), intent(inout) :: atom
    type (CellData), target, intent(in) :: cell

    type (CellData), pointer :: cell_ptr

    cell_ptr => cell
    atom%cell_ptr => cell_ptr

    if (atom%cell_index /= cell_ptr%cell_index) then
      write(*,*) "ERROR: Mismatch in cell indices for atom ", atom%atom_index
      stop
    endif

  end subroutine

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its radial mesh.
  !>
  !> That way one does not have to keep track of which mesh belongs
  !> to which atom - as soon the relation has been established.
  subroutine associateBasisAtomMesh(atom, mesh)
    use RadialMeshData_mod
    implicit none
    type (BasisAtom), intent(inout) :: atom
    type (RadialMeshData), target, intent(in) :: mesh

    type (RadialMeshData), pointer :: mesh_ptr

    mesh_ptr => mesh
    atom%mesh_ptr => mesh_ptr

  end subroutine

  !----------------------------------------------------------------------------
  !> Associates a basis atom with its reference cluster.
  !>
  !> This is done that way because different atoms can share the same
  !> reference cluster data.
!  subroutine associateBasisAtomRefCluster(atom, cluster_ptr)
!    use RefClusterData_mod
!    implicit none
!    type (BasisAtom), intent(inout) :: atom
!    type (RefClusterData), pointer :: cluster_ptr
!
!    atom%cluster_ptr => cluster_ptr
!    atom%cluster_index = cluster_ptr%cluster_index  ! TODO
!
!  end subroutine

  !----------------------------------------------------------------------------
  subroutine destroyBasisAtom(atom)
    use PotentialData_mod
    use AtomicCoreData_mod
    implicit none

    type (BasisAtom), intent(inout) :: atom

    nullify(atom%cell_ptr)
    nullify(atom%cluster_ptr)

    call destroyPotentialData(atom%potential)
    call destroyAtomicCoreData(atom%core)
  end subroutine

  !----------------------------------------------------------------------------
  !> Returns index of corresponding cell.
  pure integer function getCellIndex(atom)
     implicit none

    type (BasisAtom), intent(in) :: atom

    getCellIndex = atom%cell_index
  end function

  !----------------------------------------------------------------------------
  !> Write to unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine writeBasisAtomPotentialDA(atom, fileunit, recnr)
    implicit none
    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    write (fileunit, rec=recnr) atom%potential%VINS, &
                                atom%potential%VISP, &
                                atom%core%ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Opens unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  subroutine openBasisAtomPotentialDAFile(atom, fileunit, filename)
    implicit none
    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    integer :: reclen

    inquire (iolength = reclen)   atom%potential%VINS, &
                                  atom%potential%VISP, &
                                  atom%core%ECORE

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Reads from unformatted direct-access (DA) file.
  !>
  !> The format is compatible to the 'vpotnew' file
  !> A record number has to be specified - argument 'recnr'
  !> File has to be opened and closed by user. No checks
  subroutine readBasisAtomPotentialDA(atom, fileunit, recnr)
    implicit none
    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    read (fileunit, rec=recnr) atom%potential%VINS, &
                                atom%potential%VISP, &
                                atom%core%ECORE

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes unformatted direct-access (DA) potential-file.
  !>
  !> The format is compatible to the 'vpotnew' file
  subroutine closeBasisAtomPotentialDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine

  !----------------------------------------------------------------------------
  !> Write basis atom data to direct access file 'fileunit' at record 'recnr'.
  subroutine writeBasisAtomDA(atom, fileunit, recnr)

    implicit none
    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = 385306

    write (fileunit, rec=recnr) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER
  end subroutine

  !----------------------------------------------------------------------------
  !> Read basis atom data from direct access file 'fileunit' at record 'recnr'.
  subroutine readBasisAtomDA(atom, fileunit, recnr)
    implicit none

    type (BasisAtom), intent(inout) :: atom
    integer, intent(in) :: fileunit
    integer, intent(in) :: recnr

    integer, parameter :: MAGIC_NUMBER = 385306
    integer :: magic, magic2

    read  (fileunit, rec=recnr) magic, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                magic2

    if (magic /= MAGIC_NUMBER .or. magic2 /= MAGIC_NUMBER) then
      write (*,*) "ERROR: Invalid basis atom data read. ", __FILE__, __LINE__
      STOP
    end if
  end subroutine

  !----------------------------------------------------------------------------
  !> Opens BasisAtom direct access file.
  subroutine openBasisAtomDAFile(atom, fileunit, filename)
    implicit none

    type (BasisAtom), intent(in) :: atom
    integer, intent(in) :: fileunit
    character(len=*), intent(in) :: filename

    !------
    integer :: reclen

    integer, parameter :: MAGIC_NUMBER = 385306

    inquire (iolength = reclen) MAGIC_NUMBER, &
                                atom%atom_index, &
                                atom%cell_index, &
                                atom%cluster_index, &
                                atom%nspin, &
                                atom%Z_nuclear, &
                                atom%core%NCORE, &
                                atom%core%LCORE, &
                                atom%core%ITITLE, &
                                MAGIC_NUMBER

    open(fileunit, access='direct', file=filename, recl=reclen, form='unformatted')

  end subroutine

  !----------------------------------------------------------------------------
  !> Closes BasisAtom direct access file.
  subroutine closeBasisAtomDAFile(fileunit)
    implicit none
    integer, intent(in) :: fileunit

    close(fileunit)

  end subroutine


end module BasisAtom_mod
